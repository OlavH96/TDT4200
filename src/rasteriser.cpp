#include "rasteriser.hpp"
#include "utilities/lodepng.h"
#include <vector>
#include <chrono>
#include "Timer.cpp"
#include <math.h>
#include <algorithm>

// --- Overview ---

// I'm going to assume most of you who are reading through this file have never worked with computer graphics before.
// Let me therefore start by explaining the general idea of this file a little bit, because you'll need to at least _slightly_
// understand what is going on in this program if you are to optimise it effectively.

// We start off with two lists; one with vertices, and one with indices.
// Vertices are three-dimensional coordinates (even though in this file they have the type float4, where the last coordinate is always 1. The reason why is a long story).
// Indices represent the indices of vertices in the vertex buffer that should be connected together,
// much like you do in "connect the dots". Except here every three indices define a triangle.
// So for instance, index 0, 1, and 2 reference the three vertices of the first triangle. Indices 3, 4, and 5 the next. And so on.

// Rendering an image is done in two steps.
// First, all vertices are multiplied with a so-called "transformation matrix" in a section / program called the "vertex shader".
// How it works is not really relevant now. All you should know is that its purpose is to move objects around on the scene, and apply a "perspective" effect
// that causes objects to appear smaller the further they move away from the camera.
// The entire process only requires that all vertices are multiplied by a 4x4 matrix, using regular matrix multiplication (4x4 matrix with a 4x1 vertex)
// The second main stage involves rasterising the transformed triangles
// A "depth buffer" is used to make sure that triangles that should appear in front of others
// A "fragment shader" is run on each rendered pixel, which determines the colour the pixel should receive.
// The rendered image is stored in the so-called "framebuffer".


/**
 * Executes the vertex shader, transforms vertices and normals of the mesh object
 * @param mesh                    Mesh object with all vertices and normals
 * @param transformedVertexBuffer returned transformed vertices
 * @param transformedNormalBuffer returned transformed normals
 */
void runVertexShader( Mesh &mesh,
					  std::vector<float4> &transformedVertexBuffer,
					  std::vector<float4> &transformedNormalBuffer )
{
	// The & in front of the variable names cause the function to modify variables from the function
	// calling this one, rather than making a copy of them.

	// The matrices defined below are the ones used to transform the vertices and normals.

	// This projection matrix assumes a 16:9 aspect ratio, and an field of view (FOV) of 90 degrees.
	mat4x4 projectionMatrix(
		0.347270,   0, 			0, 		0,
		0,	  		0.617370, 	0,		0,
		0,	  		0,			-1, 	-0.2f,
		0,	  		0,			-1,		0);

	mat4x4 viewMatrix(
		0.5, 	0, 		0, 		5,
		0, 		-0.5, 	0, 		30,
		0, 		0, 		0.5, 	-55,
		0, 		0, 		0, 		1);

	mat4x4 normalMatrix(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	mat4x4 MVP = projectionMatrix * viewMatrix;

	for(unsigned int i = 0; i < transformedVertexBuffer.size(); i++) {
		float4 transformed = MVP * mesh.vertices[i];
		transformed.x /= transformed.w;
		transformed.y /= transformed.w;
		transformed.z /= transformed.w;
		transformedVertexBuffer.at(i) = transformed;
	}

	for(unsigned int j = 0; j < transformedNormalBuffer.size(); j++) {
		float4 normal;
		normal.x = mesh.normals[j].x;
		normal.y = mesh.normals[j].y;
		normal.z = mesh.normals[j].z;
		normal.w = 1;
		transformedNormalBuffer.at(j) = normalMatrix * normal;
	}
}

/**
 * Executes the fragment shader to calculate the colour of the pixel from its
 * normal
 * @param normal triangle pixel normal
 * @return		 colour of the pixel in RGBA
 */
std::vector<unsigned char> runFragmentShader( float3 const normal )
{
	std::vector<unsigned char> pixelColour(4);
	const float3 lightDirection(0.0f, 0.0f, 1.0f);

	// Computing the dot product between the surface normal and a light
	// direction gives a diffuse-like reflection. It looks more than
	// good enough for a few static images.
	float colour = normal.x * lightDirection.x +
		normal.y * lightDirection.y +
		normal.z * lightDirection.z;

	// We first scale the colour value from a range between 0 and 1,
	// to between 0 and 255.
	// Since single bytes are only able to go between 0 and 255,
	// we subsequently clamp the colour to lie within that range.
	unsigned char colourByte = (unsigned char) std::min(255.0f,
		std::max(colour * 255.0f, 0.0f));

	// And this writes the pixel to the pixelColor vector. The first three
	// channels are red, green, and blue. The fourth represents transparency.
	pixelColour.at(0) = colourByte;
	pixelColour.at(1) = colourByte;
	pixelColour.at(2) = colourByte;
	pixelColour.at(3) = 255;
	// This colour vector is supposed to go into the frame buffer
	return pixelColour;
}

/**
 * interpolates given normals using the barycentric weights
 *
 * Using float3 to return all three normals at once
 *
 * @param  n0 normal of triangle vertex
 * @param  n1 normal of triangle vertex
 * @param  n2 normal of triangle vertex
 * @param  w0 barycentric weight
 * @param  w1 barycentric weight
 * @param  w2 barycentric weight
 * @return    interpolated normal
 */
inline float3 interpolateNormals( float4 const n0,
						   float4 const n1,
					 	   float4 const n2,
					 	   float const w0,
						   float const w1,
						   float const w2 )
{
	float3 res;
	res.x = w0 * n0.x + w1 * n1.x + w2 * n2.x;
	res.y = w0 * n0.y + w1 * n1.y + w2 * n2.y;
	res.z = w0 * n0.z + w1 * n1.z + w2 * n2.z;
	return res;
}

/**
 * converts a vertex from clipping space to screen pixel coordinates
 * @param  vertex vertex
 * @param  width  screen width
 * @param  height screen height
 * @return        vertex in screen pixel coordinates
 */
inline float4 convertClippingSpace( float4 const vertex,
							 unsigned int const width,
							 unsigned int const height )
{
	float4 res;
	res.x = (vertex.x + 0.5f) * float(width);
	res.y = (vertex.y + 0.5f) * float(height);
	res.z = vertex.z;
	res.w = vertex.w;
	return res;
}

/**
 * Computing barycentric weights. This is both for determining whether the point
 * lies within the triangle, as well as interpolating coordinates needed for
 * rendering. If you need an explanation, you can find one here:
 * https://codeplea.com/triangular-interpolation
 * Using float3 to return all weights at once (x -> w0, y -> w1, z -> w2)
 *
 * @param  v0 triangle vertex
 * @param  v1 triangle vertex
 * @param  v2 triangle vertex
 * @param  x  screen pixel x-coordinate
 * @param  y  screen pixel y-coordinate
 * @return    barycentric weights of the pixel in relation to the triangle vertices
 */
// float3 getTriangleBarycentricWeights( float4 const v0,
// 									  float4 const v1,
// 									  float4 const v2,
// 									  unsigned int const x,
// 									  unsigned int const y )
// {
// 	float3 res;
// 	auto var0 = v1.y - v2.y;
// 	auto var4 = v0.x - v2.x;
// 	auto var5 = v2.x - v1.x;
// 	auto divisor = (var0 * var4) + (var5 * (v0.y - v2.y));
// 	auto var1 = x - v2.x;
// 	auto var2 = y - v2.y;
//
// 	res.x = ((var0 * var1) + (var5 * var2)) / divisor;
// 	res.y = (((v2.y - v0.y) * var1) + (var4 * var2)) / divisor;
// 	res.z = 1 - res.x - res.y;
// 	return res;
// }
float4 precalculateTriangleBarycentricWeights(float4 const &v0,
									  float4 const &v1,
									  float4 const &v2){

auto var0 = v1.y - v2.y;
auto var4 = v0.x - v2.x;
auto var5 = v2.x - v1.x;
auto divisor = (var0 * var4) + (var5 * (v0.y - v2.y));

float4 res;
res.x = var0;
res.y = var4;
res.z = var5;
res.w = 1/divisor;
return res;

}
float3 continueCalcuateTriangleBarycentricWeights(const float4 &precalc, const float4 &v0,
	const float4 &v2, unsigned int const &x, unsigned int const &y ){

	auto var1 = x - v2.x;
	auto var2 = y - v2.y;

	float3 res;
	res.x = ((precalc.x * var1) + (precalc.z * var2)) * precalc.w;
	res.y = (((v2.y - v0.y) * var1) + (precalc.y * var2)) * precalc.w;
	res.z = 1 - res.x - res.y;
	return res;

}

/**
 * Calculates the triangle's pixel depth from the barycentric weights
 * @param  v0 triangle vertex
 * @param  v1 triangle vertex
 * @param  v2 triangle vertex
 * @param  w0 barycentric weight
 * @param  w1 barycentric weight
 * @param  w2 barycentric weight
 * @return    pixel depth (z)
 */
inline float getTrianglePixelDepth( float4 const &v0,
							 float4 const &v1,
							 float4 const &v2,
							 float const &w0,
							 float const &w1,
							 float const &w2 )
{
	return w0 * v0.z + w1 * v1.z + w2 * v2.z;
}

/**
 * The main procedure which rasterises all triangles on the framebuffer
 * @param mesh                    Mesh object
 * @param transformedVertexBuffer transformed vertices from the mesh obj
 * @param transformedNormalBuffer transformed normals from the mesh obj
 * @param frameBuffer             frame buffer for the rendered image
 * @param depthBuffer             depth buffer for every pixel on the image
 * @param width                   width of the image
 * @param height                  height of the image
 */
void rasteriseTriangles( Mesh &mesh,
                         std::vector<float4> &transformedVertexBuffer,
                         std::vector<float4> &transformedNormalBuffer,
                         std::vector<unsigned char> &frameBuffer,
                         std::vector<float> depthBuffer,
                         const unsigned int &width,
                         const unsigned int &height )
{
	Timer timer;
	timer.start();

	Timer pixelTimer;
	Timer innerForLoopTimer;
	Timer outerForLoopTimer;

	// We rasterise one triangle at a time
	unsigned int triangleCount = mesh.indexCount / 3;
	for(unsigned int triangleIndex = 0; triangleIndex < triangleCount; triangleIndex++) {
		// '\r' returns to the beginning of the current line
		//std::cout << "Rasterizing triangle "<<triangleIndex << "/" << triangleCount<<'\r';
		// As vertices are commonly reused within a model, rendering libraries use an
		// index buffer which specifies the indices of the vertices in the vertex buffer
		// which together make up the specific triangle.
		const int base = 3 * triangleIndex;
		unsigned int index0 = mesh.indices[base + 0];
		unsigned int index1 = mesh.indices[base + 1];
		unsigned int index2 = mesh.indices[base + 2];

		// We look up those triangles here
		float4 vertex0 = float4(transformedVertexBuffer.at(index0));
		float4 vertex1 = float4(transformedVertexBuffer.at(index1));
		float4 vertex2 = float4(transformedVertexBuffer.at(index2));

		// These triangles are still in so-called "clipping space". We first convert them
		// to screen pixel coordinates
		vertex0 = convertClippingSpace(vertex0, width, height);
		vertex1 = convertClippingSpace(vertex1, width, height);
		vertex2 = convertClippingSpace(vertex2, width, height);

		float4 precalc = precalculateTriangleBarycentricWeights(vertex0, vertex1, vertex2);
		if (precalc.w >= 0.0) {
			continue;
		}

		// Read the normals belonging to each vertex
		float4 normal0 = float4(transformedNormalBuffer.at(index0));
		float4 normal1 = float4(transformedNormalBuffer.at(index1));
		float4 normal2 = float4(transformedNormalBuffer.at(index2));

		const int xcoordinates[] = {ceil(vertex0.x), ceil(vertex1.x), ceil(vertex2.x)};
		const int ycoordinates[] = {ceil(vertex0.y), ceil(vertex1.y), ceil(vertex2.y)};

		const unsigned int xmax = *std::max_element(xcoordinates, xcoordinates+3);
		const unsigned int ymax = *std::max_element(ycoordinates, ycoordinates+3);

		const unsigned int xmin = *std::min_element(xcoordinates, xcoordinates+3);
		const unsigned int ymin = *std::min_element(ycoordinates, ycoordinates+3);

		if (!pixelTimer.isRunning()) {
			pixelTimer.start();
		}
		// We iterate over each pixel on the screen
		for(unsigned int x = xmin; x < xmax; x++) {
			if (!outerForLoopTimer.isRunning()) {
				outerForLoopTimer.start();
			}
			for(unsigned int y = ymin; y < ymax; y++) {

				if (!innerForLoopTimer.isRunning()) {
					innerForLoopTimer.start();
				}

				// Calculating the barycentric weights of the pixel in relation to the triangle
				auto tbw = continueCalcuateTriangleBarycentricWeights(precalc, vertex0, vertex2, x, y); //getTriangleBarycentricWeights(vertex0, vertex1, vertex2, x, y);
				float weight0 = tbw.x;
				float weight1 = tbw.y;
				float weight2 = tbw.z;

				// Skip pixels outside frame asap
				// The weights have the nice property that if only one is negative, the pixel lies outside the triangle
				if(!(weight0 >= 0 && weight1 >= 0 && weight2 >= 0)) {
					continue;
				}
				// Now we can determine the depth of our pixel
				// Z-clipping discards pixels too close or too far from the camera
				float pixelDepth = getTrianglePixelDepth(vertex0, vertex1, vertex2, weight0, weight1, weight2);
				// Skip pixels too close or too far away asap
				if(!(pixelDepth >= -1 && pixelDepth <= 1)) {
						continue;
				}
				//Have we drawn a pixel above the current?
					// This pixel is going into the frame buffer,
					// save its depth to skip all next pixels underneath it
				int position = y * width + x;
				// Prevents the position going out of bounds for large meshes
				if (position > depthBuffer.size()) {
					continue;
				}
				if(!(pixelDepth < depthBuffer.at(position))){
					continue;
				}


				// But since a pixel can lie anywhere between the vertices, we compute an approximated normal
				// at the pixel location by interpolating the ones from the vertices.
				float3 interpolatedNormal = interpolateNormals(normal0, normal1, normal2, weight0, weight1, weight2);

				// This process can slightly change the length, so we normalise it here to make sure the lighting calculations
				// appear correct.
				float normalLength = std::sqrt( interpolatedNormal.x * interpolatedNormal.x +
					interpolatedNormal.y * interpolatedNormal.y +
					interpolatedNormal.z * interpolatedNormal.z );

				// Turn 3 divisions into 1 divisions and 3 multiplications 30 -> 19
				float normalLengthInverse = 1/normalLength;

				interpolatedNormal.x *= normalLengthInverse;
				interpolatedNormal.y *= normalLengthInverse;
				interpolatedNormal.z *= normalLengthInverse;

				// And we can now execute the fragment shader to compute this pixel's colour.
				std::vector<unsigned char> pixelColour = runFragmentShader(interpolatedNormal);

				depthBuffer.at(position) = pixelDepth;
				// Copy the calculated pixel colour into the frame buffer - RGBA

				//Coordinate of the current pixel in the framebuffer, remember RGBA color code
				unsigned int pixelBaseCoordinate = (position) << 2;

				frameBuffer.at(pixelBaseCoordinate + 0) = pixelColour.at(0);
				frameBuffer.at(pixelBaseCoordinate + 1) = pixelColour.at(1);
				frameBuffer.at(pixelBaseCoordinate + 2) = pixelColour.at(2);
				frameBuffer.at(pixelBaseCoordinate + 3) = pixelColour.at(3);

				innerForLoopTimer.iterate();
			}
			outerForLoopTimer.iterate();
		}
		pixelTimer.iterate();

	}
	std::cout << std::endl;
	innerForLoopTimer.stop();
	timer.stop();
	pixelTimer.stop();
	outerForLoopTimer.stop();

	auto pixelTime = pixelTimer.getAverageTime();
	auto pixelTotalTime = pixelTimer.getTotalTime();
	std::cout << "PixelIteration took " << pixelTime << "ns on average, "<<pixelTotalTime <<"ns in total" << std::endl;

	auto outerTime = outerForLoopTimer.getAverageTime();
	auto outerTotalTime = outerForLoopTimer.getTotalTime();
	std::cout << "outerForLoopTimer took " << outerTime << "ns on average, "<<outerTotalTime << "ns in total" << std::endl;


	auto innerTime = innerForLoopTimer.getAverageTime();
	auto innerTotalTime = innerForLoopTimer.getTotalTime();
	std::cout << "innerForLoopTimer took " << innerTime << "ns on average, "<<innerTotalTime << "ns in total" << std::endl;

	auto duration = timer.getTimeElapsedMilliseconds();
	std::cout << "rasterizeTriangles used " << duration << "ns";
	// finish the progress output with a new line
	std::cout << std::endl;
}

/**
 * Procedure to kick of the rasterisation process
 * @param mesh            Mesh object
 * @param outputImageFile path of the output image
 * @param width           width of the output image
 * @param height          height of the output image
 */
void rasterise(Mesh mesh, std::string outputImageFile, unsigned int width, unsigned int height) {
	// We first need to allocate some buffers.
	// Timer rasteriseTimer;
	// rasteriseTimer.start();
	// The framebuffer contains the image being rendered.
	std::vector<unsigned char> frameBuffer;
	const int number_of_pixels = width * height;
	frameBuffer.resize(number_of_pixels * 4, 0);

	// The depth buffer is used to make sure that objects closer to the camera occlude/obscure objects that are behind it
	std::vector<float> depthBuffer;
	depthBuffer.resize(number_of_pixels, 1);

	// And these two buffers store vertices and normals processed by the vertex shader.
	std::vector<float4> transformedVertexBuffer;
	transformedVertexBuffer.resize(mesh.vertexCount);

	std::vector<float4> transformedNormalBuffer;
	transformedNormalBuffer.resize(mesh.vertexCount);

	// Initializing the framebuffer with RGBA (0,0,0,255), black, no transparency
	for (unsigned int x = 0; x < width; x++) {
		for(unsigned int y = 0; y < height; y++) {
			const int size = (x + y * width) << 2;
			frameBuffer.at(size + 0) = 0;
			frameBuffer.at(size + 1) = 0;
			frameBuffer.at(size + 2) = 0;
			frameBuffer.at(size + 3) = 255;
		}
	}

	std::cout << "Running the vertex shader... ";
	runVertexShader(mesh, transformedVertexBuffer, transformedNormalBuffer);

	std::cout << "complete!" << std::endl;

	rasteriseTriangles(mesh, transformedVertexBuffer, transformedNormalBuffer, frameBuffer, depthBuffer, width, height);

	std::cout << "Finished rendering!" << std::endl;

	// rasteriseTimer.stop();

	std::cout << "Writing image to '" << outputImageFile << "'..." << std::endl;

	unsigned error = lodepng::encode(outputImageFile, frameBuffer, width, height);

	if(error)
	{
		std::cout << "An error occurred while writing the image file: " << error << ": " << lodepng_error_text(error) << std::endl;
	}
	// std::cout << "rasteriseTimer took " << rasteriseTimer.getTimeElapsedMilliseconds() << "ns" << std::endl;
}
