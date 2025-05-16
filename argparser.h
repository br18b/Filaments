#ifndef __ARGPARSER__
#define __ARGPARSER__

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

std::map<std::string, std::string> load_default_pars();
std::vector<std::string> parse(int argc, char *argv[], std::vector<std::string> options);
void sanitize_filename(std::string full_name, std::string& main, std::string& postfix, std::string& extension);
void handle_pars(int argc, char *argv[], std::string mainName, std::string& outDir, std::string& outName, double& dL, double& dR, double& line_rms_threshold, double& line_length_threshold, double& radial_frac, double& radial_fit_threshold, int& resolution, bool& saveFailed);

#endif