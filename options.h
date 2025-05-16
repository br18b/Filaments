#ifndef __OPTIONS__
#define __OPTIONS__

#include <filesystem>
#include <iostream>
#include <string>
#include <set>
#include <map>

#include "argparser.h"

void show_options();
bool check_pars(int argc, char *argv[], std::string& inFilaments, std::string& inDensity, std::string& inHx, std::string& inHy, std::string& inHz, std::string& mainName, bool& magnetized, std::string& message);

#endif