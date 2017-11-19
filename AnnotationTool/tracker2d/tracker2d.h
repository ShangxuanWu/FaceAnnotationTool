#pragma once
#include "../util/util.h"
#include "../util/cutil.h"

#include <boost/program_options.hpp>

#include <map>
#include <fstream>

#include <stdio.h>

bool tracker2d(std::string input_yml, std::string ori_png, std::string output_fn);