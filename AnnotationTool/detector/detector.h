#pragma once
#define _DEBUG_OUT
#include "settings.h"
#include "filters.h"
#include "../util/util.h"
#include "../util/timer.h"

#include <opencv2/imgproc/imgproc.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>

const double GABOR_CANNY_LOWER = 0.2;
const double GABOR_CANNY_UPPER = 0.3;
const int GABOR_NUMBER_OF_DIRECTION = 18;

bool detector(std::string original_image_path, std::string mask_path, std::string output_yml_name, std::string output_ori_png_name);