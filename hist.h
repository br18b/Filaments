#ifndef __HIST__
#define __HIST__

#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <functional>

#include <filesystem>
#include <algorithm>
#include <regex>

namespace fs = std::filesystem;

std::string pad(int number, int places, char fill);
std::string pad(int number);

bool is_int(const std::string& s);
bool is_float(const std::string& string);
bool is_path_valid(const std::string& path);

std::vector<int> get_frames(const std::string& path);

bool validate_input(int argc, char *argv[], std::string& sim_path, std::string& path_out, std::vector<int>& frames);

std::pair<double, double> find_bounds(const std::vector<double>& data);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbin, std::pair<double, double> bounds);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbin);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins, std::pair<double, double> bounds);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins);

std::vector<std::vector<double>> make_joint_histogram(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& x_bins, const std::vector<double>& y_bins);
std::vector<std::vector<double>> make_joint_histogram(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& x_bins, const std::vector<double>& y_bins, std::function<double(double, double)> weight_fun);

std::vector<std::vector<std::vector<double>>> make_joint_histogram(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& data_z, const std::vector<double>& x_bins, const std::vector<double>& y_bins, const std::vector<double>& z_bins);
std::vector<std::vector<std::vector<double>>> make_joint_histogram(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& data_z, const std::vector<double>& x_bins, const std::vector<double>& y_bins, const std::vector<double>& z_bins, std::function<double(double, double, double)> weight_fun);

void make_histograms(const std::string& path, const std::vector<int>& frames, const std::string& path_out);

#endif