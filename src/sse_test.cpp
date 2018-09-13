#include "sse_test.hpp"
#include <stdlib.h>
#include <time.h>

inline float randFloat() {
    return (float) (rand()) / (float) (RAND_MAX);
}

inline sse_float4 randFloat4() {
    sse_float4 res;
    res.elements[0] = randFloat();
    res.elements[1] = randFloat();
    res.elements[2] = randFloat();
    res.elements[3] = randFloat();
    return res;
}

void sse_test(const Mesh* mesh) {
    //Not allowed to change:
    unsigned int const loadFactor = 1000;

    // Mesh vertices container
    std::vector<sse_float4> vertices;
    vertices.resize(mesh->vertexCount);

    // Containers of random floats
    std::vector<sse_float4> rand1;
    rand1.resize(mesh->vertexCount);
    std::vector<sse_float4> rand2;
    rand2.resize(mesh->vertexCount);
    std::vector<sse_float4> rand3;
    rand3.resize(mesh->vertexCount);
    std::vector<sse_float4> rand4;
    rand4.resize(mesh->vertexCount);

    srand(time(NULL));

    std::cout << "SSE_TEST: Initializing vectors... " << std::flush;
    for (unsigned int i=0; i < vertices.size(); i++) {
        float4 tempVertices = mesh->vertices[i];
        vertices[i].elements[0] = tempVertices.x;
        vertices[i].elements[1] = tempVertices.y;
        vertices[i].elements[2] = tempVertices.z;
        vertices[i].elements[3] = tempVertices.w;

        rand1.push_back(randFloat4());
        rand2.push_back(randFloat4());
        rand3.push_back(randFloat4());
        rand4.push_back(randFloat4());    
    }
    std::cout << "finished!"  << std::endl;
    for (unsigned int loadIterator = 0; loadIterator < loadFactor; loadIterator++) {
        std::cout << "SSE_TEST: " << (loadIterator+1) << "/" << loadFactor << " Crunching numbers on " << vertices.size() << " vertices... " << "\r" << std::flush;
        for (unsigned int i = 0; i < vertices.size(); i++) {
            vertices[i].vector = (vertices[i].vector + rand1[i].vector - rand2[i].vector) * rand3[i].vector;

            if (rand4[i].vector[0] != 0 && rand4[i].elements[1] != 0
                && rand4[i].elements[2] != 0 && rand4[i].elements[3] != 0) {
                vertices[i].vector = vertices[i].vector / rand4[i].vector;
            }
        }
    }
    std::cout << std::endl;
}
