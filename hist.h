#ifndef __HIST__
#define __HIST__

#include <vector>
#include <set>
#include <iostream>
#include <functional>

std::pair<double, double> find_bounds(const std::vector<double>& data);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbin, std::pair<double, double> bounds);
std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins, std::pair<double, double> bounds);


#endif