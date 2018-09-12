#pragma once

#include <string>
#include "utilities/OBJLoader.hpp"

void rasterise(Mesh &mesh, const std::string &outputImageFile, const unsigned int &width, const unsigned int &height);
