#pragma once

#include <vector>
#include <stddef.h>

std::vector<std::vector<unsigned>> determineBestViews(const std::vector<std::vector<double>>& scores, size_t topK, double viewChangeCost, size_t* changes = NULL, double* totalCost = NULL);
