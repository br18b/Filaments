#include "hist.h"

std::pair<double, double> find_bounds(const std::vector<double>& data) {
    double xmin = 1e10;
    double xmax = -1e10;
    for (const double& d: data) {
        if (xmin > d) xmin = d;
        if (xmax < d) xmax = d;
    }
    return std::make_pair(xmin, xmax);
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbins, std::pair<double, double> bounds) {
    std::pair<std::vector<double>, std::vector<double>> res;
    if (nbins > 1) {
        int n = data.size();
        int npt = n / nbins;
        int remainder = n - npt * nbins;
        std::multiset<double> sorted;
        double xmin = bounds.first;
        double xmax = bounds.second;
        for (const double& d: data) {
            if ((d >= xmin) && (d <= xmax)) sorted.insert(d);
        }
        std::vector<double> bin_dividers, weights;
        bin_dividers.push_back(xmin);
        int count = 0; int bin_index = 0;
        std::multiset<double>::iterator it;
        for (it = sorted.begin(); it != sorted.end(); ++it) {
            if (count < npt) {
                count++;
            }
            else {
                double val_old = *std::prev(it);
                double val = *it;
                bin_dividers.push_back(0.5 * (val_old + val)); weights.push_back(npt);
                bin_index++;
                count = 1;
            }
        }
        if (remainder < 0.5 * npt) {
            bin_dividers.back() = xmax; weights.back() += remainder;
        }
        else {
            bin_dividers.push_back(xmax); weights.push_back(remainder);
        }
        for (auto& w: weights) w /= n;
        res.first = bin_dividers; res.second = weights;
        return res;
    }
    else {
        std::cout << "Number of bins must be > 1!" << std::endl;
        return res;
    }
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins, std::pair<double, double> bounds) {
    std::pair<std::vector<double>, std::vector<double>> res;
    if (nbins > 1) {
        int n = data.size();
        double total_weight = 0;
        for (const double& d : data) {
            total_weight+= weight_fun(d);
        }
        double bin_weight = total_weight / nbins;
        std::multiset<double> sorted;
        double xmin = bounds.first;
        double xmax = bounds.second;
        for (const double& d: data) {
            if ((d >= xmin) && (d <= xmax)) sorted.insert(d);
        }
        std::vector<double> bin_dividers, weights;
        bin_dividers.push_back(xmin);
        double accu_weight = 0; int bin_index = 0;
        std::multiset<double>::iterator it;
        for (it = sorted.begin(); it != sorted.end(); ++it) {
            if (accu_weight < bin_weight) {
                accu_weight += weight_fun(*it);
            }
            else {
                double val_old = *std::prev(it);
                double val = *it;
                bin_dividers.push_back(0.5 * (val_old + val)); weights.push_back(accu_weight);
                bin_index++;
                accu_weight = weight_fun(val);
            }
        }
        if (accu_weight < 0.5 * bin_weight) {
            bin_dividers.back() = xmax; weights.back() += accu_weight;
        }
        else {
            bin_dividers.push_back(xmax); weights.push_back(accu_weight);
        }
        for (auto& w: weights) w /= total_weight;
        res.first = bin_dividers; res.second = weights;
        return res;
    }
    else {
        std::cout << "Number of bins must be > 1!" << std::endl;
        return res;
    }
}