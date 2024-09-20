#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <cstdio>
#include <math.h>
#include <map>

#include "string_pad.h"
#include "files.h"
#include "vec3D.h"
#include "hist.h"
#include "paths_frames.h"

static std::function<double(double)> weight_angle3D = [](double x) {
    if (x != 0) return 1 / sin(std::abs(x));
    else return 0.0;
};

void make_histograms(std::string path, std::string stats_filename, std::string prefix, int frame) {
    std::vector<vec3D> normal, havg;
    std::vector<double> length, radius, aspect, thetaLH, thetaLx, thetaHx;
    std::string frame_str = string_pad_left(std::to_string(frame), "0", 4);
    std::ifstream fin; fin.open(path + stats_filename);
    std::string line;
    while(getline(fin, line)) {
        double a, b, c, d;
        std::istringstream ss(line);
        vec3D P1, P2, h;
        ss >> P1.x >> P1.y >> P1.z >> P2.x >> P2.y >> P2.z >> h.x >> h.y >> h.z >> a >> b >> c >> d;
        vec3D n = P2 - P1; normalize(n);
        normal.push_back(n); havg.push_back(h); vec3D h_norm = h; normalize(h_norm);
        length.push_back(a); radius.push_back(d); aspect.push_back(d / a);
        thetaLH.push_back(acos(n * h_norm));
        thetaLx.push_back(acos(n.x));
        thetaHx.push_back(acos(h_norm.x));
    }
    fin.close(); int nbins = 10;
    std::pair<double, double> bounds_halfangle = std::make_pair(0, M_PI / 2);
    std::pair<double, double> bounds_angle = std::make_pair(0, M_PI);
    std::pair<double, double> bounds_aspect = std::make_pair(0, 1);

    auto hist_thetaLH_weighted = make_histogram(thetaLH, weight_angle3D, nbins, bounds_halfangle);
    auto hist_thetaLx_weighted = make_histogram(thetaLx, weight_angle3D, nbins, bounds_angle);
    auto hist_thetaHx_weighted = make_histogram(thetaHx, weight_angle3D, nbins, bounds_angle);
    save_histogram(hist_thetaLH_weighted, path + "data" + frame_str + "_thetaLH_weighted.txt");
    save_histogram(hist_thetaLx_weighted, path + "data" + frame_str + "_thetaLx_weighted.txt");
    save_histogram(hist_thetaHx_weighted, path + "data" + frame_str + "_thetaHx_weighted.txt");

    auto hist_thetaLH_raw = make_histogram(thetaLH, nbins, bounds_halfangle);
    auto hist_thetaLx_raw = make_histogram(thetaLx, nbins, bounds_angle);
    auto hist_thetaHx_raw = make_histogram(thetaHx, nbins, bounds_angle);
    save_histogram(hist_thetaLH_raw, path + "data" + frame_str + "_thetaLH.txt");
    save_histogram(hist_thetaLx_raw, path + "data" + frame_str + "_thetaLx.txt");
    save_histogram(hist_thetaHx_raw, path + "data" + frame_str + "_thetaHx.txt");

    auto hist_aspect = make_histogram(aspect, nbins, bounds_aspect);
    save_histogram(hist_aspect, path + "data" + frame_str + "_aspect.txt");
}

void make_histograms(std::string prefix, int frame, std::string pers_threshold) {
    std::string frame_str = string_pad_left(std::to_string(frame), "0", 4);
    std::string path = P71_path + prefix + "/";
    std::string stats_filename = "data" + frame_str + "_filament_c" + pers_threshold + "_stats.txt";
    make_histograms(path, stats_filename, prefix, frame);
}

int main(int argc, char *argv[]) {
    std::string sim = "1_half";
    std::string pers_threshold = "0.5";
    if (argc > 1) sim = argv[1];
    int frame = last_frames[sim];
    if (argc > 2) frame = std::stoi(argv[2]);
    if (argc > 3) pers_threshold = argv[3];
    make_histograms(sim, frame, pers_threshold);
}