#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
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

static std::function<double(double,double)> weight_angles3D = [](double x, double y) {
    if ((x != 0) && (y != 0)) return 1 / (sin(std::abs(x)) * sin(std::abs(y)));
    else return 0.0;
};

static double clamp(double x, double min, double max) {
    return std::min(max, std::max(min, x));
}

double fourierLH(const std::vector<double>& angle, int n) {
    double res = 0;
    for (auto& theta: angle) {
        res += sin(n * theta);
    }
    res /= angle.size();
    return res;
}

void fourierLH(const std::vector<double>& angle, const std::string& path_out, int nmax) {
    std::ofstream file_out; file_out.open(path_out + "fourier_odd_thetaLH.txt");
    for (int n = 1; n <= nmax; n+=2) {
        double coeff = fourierLH(angle, n);
        file_out << n << " " << coeff << std::endl;
    }
    file_out.close();
}

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> points(n);
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        points[i] = a + i * step;
    }
    return points;
}

std::vector<double> linspace_bin_friendly(double A, double B, int n) {
    std::vector<double> points(n);
    double a = (2 * A * n - 3 * A - B) / (2 * (n - 2));
    double b = (2 * B * n - A - 3 * B) / (2 * (n - 2));
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        points[i] = a + i * step;
    }
    return points;
}

long double mode_zero(long double u, int i) {
    if ((i + 1) % 3 == 0) {
        int n = 2 * (i + 1) / 3;
        return sqrt(2) * sin(n * M_PI * u);
    }
    else if ((i + 1) % 3 == 1) {
        int n = 2 * i / 3 + 1;
        return sqrt(2) * cos(n * M_PI * u);
    }
    else {
        int n = 2 * (i - 1) / 3 + 2;
        return sqrt(2) * cos(n * M_PI * u);
    }
}

long double mode_full(long double u, int i) {
    if (i == 0) {
        return 1;
    }
    else if (i % 2 == 1) {
        int n = (i - 1) / 2 + 1;
        return cos(n * M_PI * u / 2);
    }
    else {
        int n = i / 2;
        return sin(n * M_PI * u / 2);
    }
}

void extract_fourier(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, int N, std::vector<std::vector<double>>& F1D, std::vector<double>& F3D) {
    F1D.resize(3, std::vector<double>(N, 0.0));
    for (int comp = 0; comp < 3; comp++) {
        F1D[comp].resize(N);
        for (int i = 0; i < N; i++) {
            F1D[comp][i] = 0;
        }
    }
    F3D.resize(N * N * N, 0.0);
    for (int t = 0; t < u.size(); t++) {
        std::vector<double> mu(N), mv(N), mw(N);
        for (int i = 0; i < N; i++) {
            mu[i] = cos(i * M_PI * u[t]);
            mv[i] = cos(i * M_PI * v[t]);
            mw[i] = sin((2 * i + 1) * w[t] / 2);
        }
        for (int i = 0; i < N; i++) {
            F1D[0][i] += mu[i];
            F1D[1][i] += mv[i];
            F1D[2][i] += mw[i];
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    F3D[N * N * i + N * j + k] += mu[i] * mv[j] * mw[k];
                }
            }
        }
    }
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < N; i++) {
            F1D[comp][i] /= u.size();
        }
    }
    for (int i = 0; i < N * N * N; i++) {
        F3D[i] /= u.size();
    }
    
}

void make_histograms_uvw(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& vbar, const std::vector<double>& w, int nbins, const std::string& path_out, std::string axis) {
    auto hist_u = make_histogram(u, nbins, {0,1});
    auto hist_v = make_histogram(v, nbins, {0,1});
    auto hist_v2 = make_histogram(vbar, nbins, {0,1});
    auto hist_w = make_histogram(w, nbins, {0,1});
    save_histogram(hist_u, path_out + "hist_u" + axis + ".txt");
    save_histogram(hist_v, path_out + "hist_v" + axis + ".txt");
    save_histogram(hist_v2, path_out + "hist_vbar" + axis + ".txt");
    save_histogram(hist_w, path_out + "hist_w" + axis + ".txt");
    
    auto bins = linspace_bin_friendly(0, 1, nbins);
    
    auto hist_uv = make_joint_histogram(u, v, bins, bins);
    auto hist_uw = make_joint_histogram(u, w, bins, bins);
    auto hist_vw = make_joint_histogram(v, w, bins, bins);
    
    auto hist_uv2 = make_joint_histogram(u, vbar, bins, bins);
    auto hist_vw2 = make_joint_histogram(vbar, w, bins, bins);
    
    save_joint_histogram(bins, bins, hist_uv, path_out + "hist_uv" + axis + ".txt");
    save_joint_histogram(bins, bins, hist_uw, path_out + "hist_uw" + axis + ".txt");
    save_joint_histogram(bins, bins, hist_vw, path_out + "hist_vw" + axis + ".txt");
    
    save_joint_histogram(bins, bins, hist_uv2, path_out + "hist_uvbar" + axis + ".txt");
    save_joint_histogram(bins, bins, hist_vw2, path_out + "hist_vbarw" + axis + ".txt");
    
    bins = linspace_bin_friendly(0, 1, 7);
    auto hist3D = make_joint_histogram(u, v, w, bins, bins, bins);
    auto hist3D2 = make_joint_histogram(u, vbar, w, bins, bins, bins);
    save_joint_histogram_sparse(bins, bins, bins, hist3D, path_out + "hist_uvw" + axis + ".txt");
    save_joint_histogram_sparse(bins, bins, bins, hist3D2, path_out + "hist_uvbarw" + axis + ".txt");
}

void make_histograms(const std::vector<double>& aspect, const std::vector<double>& len, const std::vector<double>& len_proj, const std::vector<double>& thetaLH, const std::vector<double>& thetaL, const std::vector<double>& thetaH, const std::vector<double>& psiLH, std::pair<double, double> bounds_halfangle, std::pair<double, double> bounds_angle, std::pair<double, double> bounds_fullangle, std::pair<double, double> bounds_aspect, int nbins, const std::string& path_out, std::string axis) {
    auto hist_thetaLH_weighted = make_histogram(thetaLH, weight_angle3D, nbins, bounds_halfangle);
    auto hist_thetaL_weighted = make_histogram(thetaL, weight_angle3D, nbins, bounds_angle);
    auto hist_thetaH_weighted = make_histogram(thetaH, weight_angle3D, nbins, bounds_angle);
    save_histogram(hist_thetaLH_weighted, path_out + "hist_thetaLH_weighted.txt");
    save_histogram(hist_thetaL_weighted, path_out + "hist_thetaL" + axis + "_weighted.txt");
    save_histogram(hist_thetaH_weighted, path_out + "hist_thetaH" + axis + "_weighted.txt");

    auto hist_thetaLH_raw = make_histogram(thetaLH, nbins, bounds_halfangle);
    auto hist_thetaL_raw = make_histogram(thetaL, nbins, bounds_angle);
    auto hist_thetaH_raw = make_histogram(thetaH, nbins, bounds_angle);
    auto hist_psiLH_raw = make_histogram(psiLH, nbins, bounds_angle);
    save_histogram(hist_thetaLH_raw, path_out + "hist_thetaLH.txt");
    save_histogram(hist_thetaL_raw, path_out + "hist_thetaL" + axis + ".txt");
    save_histogram(hist_psiLH_raw, path_out + "hist_psiLH" + axis + ".txt");

    auto hist_aspect = make_histogram(aspect, nbins, bounds_aspect);
    save_histogram(hist_aspect, path_out + "hist_aspect.txt");

    auto hist_len = make_histogram(len, nbins);
    save_histogram(hist_len, path_out + "hist_length.txt");

    auto hist_len_proj = make_histogram(len_proj, nbins);
    save_histogram(hist_len_proj, path_out + "hist_length_proj.txt");

    // now 2D histograms
    int nbins_smaller = 10;
    hist_thetaL_raw = make_histogram(thetaL, nbins_smaller, bounds_angle);
    auto binsL = hist_thetaL_raw.first;
    //auto binsL = linspace(0, M_PI, 10);
    auto binsL_weighted = hist_thetaL_weighted.first;
    
    hist_thetaH_raw = make_histogram(thetaH, nbins_smaller, bounds_angle);
    auto binsH = hist_thetaH_raw.first;
    //auto binsH = linspace(0, M_PI, 10);
    auto binsH_weighted = hist_thetaH_weighted.first;
    
    hist_thetaLH_raw = make_histogram(thetaLH, nbins_smaller, bounds_halfangle);
    auto binsLH = hist_thetaLH_raw.first;
    //auto binsLH = linspace(0, M_PI / 2, 10);
    auto binsLH_weighted = hist_thetaLH_weighted.first;
    
    hist_psiLH_raw = make_histogram(psiLH, nbins_smaller, bounds_angle);
    auto binsPsi = hist_psiLH_raw.first;
    
    auto hist_thetaL_thetaH_raw = make_joint_histogram(thetaL, thetaH, binsL, binsH);
    auto hist_thetaL_thetaH_weighted = make_joint_histogram(thetaL, thetaH, binsL_weighted, binsH_weighted, weight_angles3D);
    
    auto hist_thetaL_thetaLH_raw = make_joint_histogram(thetaL, thetaLH, binsL, binsLH);
    auto hist_thetaL_thetaLH_weighted = make_joint_histogram(thetaL, thetaLH, binsL_weighted, binsLH_weighted, weight_angles3D);
    
    auto hist_thetaH_thetaLH_raw = make_joint_histogram(thetaH, thetaLH, binsH, binsLH);
    auto hist_thetaH_thetaLH_weighted = make_joint_histogram(thetaH, thetaLH, binsH_weighted, binsLH_weighted, weight_angles3D);
    
    auto hist_thetaL_psiLH_raw = make_joint_histogram(thetaL, psiLH, binsL, binsPsi);
    auto hist_thetaH_psiLH_raw = make_joint_histogram(thetaH, psiLH, binsH, binsPsi);
    auto hist_thetaLH_psiLH_raw = make_joint_histogram(thetaLH, psiLH, binsLH, binsPsi);

    save_joint_histogram(binsL, binsH, hist_thetaL_thetaH_raw, path_out + "hist_thetaL" + axis + "_thetaH" + axis + ".txt");
    save_joint_histogram(binsL_weighted, binsH_weighted, hist_thetaL_thetaH_weighted, path_out + "hist_thetaL" + axis + "_thetaH" + axis + "_weighted.txt");

    save_joint_histogram(binsL, binsLH, hist_thetaL_thetaLH_raw, path_out + "hist_thetaL" + axis + "_thetaLH.txt");
    save_joint_histogram(binsL_weighted, binsLH_weighted, hist_thetaL_thetaLH_weighted, path_out + "hist_thetaL" + axis + "_thetaLH_weighted.txt");

    save_joint_histogram(binsH, binsLH, hist_thetaH_thetaLH_raw, path_out + "hist_thetaH" + axis + "_thetaLH.txt");
    save_joint_histogram(binsH_weighted, binsLH_weighted, hist_thetaH_thetaLH_weighted, path_out + "hist_thetaH" + axis + "_thetaLH_weighted.txt");
    

    save_joint_histogram(binsL, binsPsi, hist_thetaL_psiLH_raw, path_out + "hist_thetaL" + axis + "_psiLH" + axis + ".txt");
    save_joint_histogram(binsH, binsPsi, hist_thetaH_psiLH_raw, path_out + "hist_thetaH" + axis + "_psiLH" + axis + ".txt");
    save_joint_histogram(binsLH, binsPsi, hist_thetaLH_psiLH_raw, path_out + "hist_thetaLH_psiLH" + axis + ".txt");
    
    auto hist3D = make_joint_histogram(thetaL, thetaH, thetaLH, binsL, binsH, binsLH);
    save_joint_histogram(binsL, binsH, binsLH, hist3D, path_out + "hist_thetaL" + axis + "_thetaH" + axis + "_thetaLH.txt");

    // histograms involving alpha = thetaL + thetaH and beta = thetaL - thetaH
    std::vector<double> alpha, beta;
    std::pair<double, double> bounds_alpha = std::make_pair(0, 2 * M_PI);
    std::pair<double, double> bounds_beta = std::make_pair(-M_PI, M_PI);
    for (int i = 0; i < thetaL.size(); i++) {
        alpha.push_back(thetaL[i] + thetaH[i]);
        beta.push_back(thetaL[i] - thetaH[i]);
    }

    auto hist_alpha = make_histogram(alpha, nbins, bounds_alpha);
    auto hist_beta = make_histogram(beta, nbins, bounds_beta);

    save_histogram(hist_alpha, path_out + "hist_alpha" + axis + ".txt");
    save_histogram(hist_beta, path_out + "hist_beta" + axis + ".txt");

    auto bins_alpha = linspace(bounds_alpha.first, bounds_alpha.second, 20);
    auto bins_beta = linspace(bounds_beta.first, bounds_beta.second, 20);
    
    auto hist_alpha_beta = make_joint_histogram(alpha, beta, bins_alpha, bins_beta);
    auto hist_alpha_LH = make_joint_histogram(alpha, thetaLH, bins_alpha, binsLH);
    auto hist_beta_LH = make_joint_histogram(beta, thetaLH, bins_beta, binsLH);
    
    save_joint_histogram(bins_alpha, bins_beta, hist_alpha_beta, path_out + "hist_alpha" + axis + "_beta" + axis + ".txt");
    save_joint_histogram(bins_alpha, binsLH, hist_alpha_LH, path_out + "hist_alpha" + axis + "_thetaLH.txt");
    save_joint_histogram(bins_beta, binsLH, hist_beta_LH, path_out + "hist_beta" + axis + "_thetaLH.txt");
}

void make_histograms_dry(const std::vector<double>& aspect, const std::vector<double>& len, const std::vector<double>& thetaLx, std::pair<double, double> bounds_angle, std::pair<double, double> bounds_aspect, int nbins, const std::string& path_out, std::string axis) {
    auto hist_thetaLx_weighted = make_histogram(thetaLx, weight_angle3D, nbins, bounds_angle);
    save_histogram(hist_thetaLx_weighted, path_out + "hist_thetaL" + axis + "_weighted.txt");

    auto hist_thetaLx_raw = make_histogram(thetaLx, nbins, bounds_angle);
    save_histogram(hist_thetaLx_raw, path_out + "hist_thetaL" + axis + ".txt");

    auto hist_aspect = make_histogram(aspect, nbins, bounds_aspect);
    save_histogram(hist_aspect, path_out + "hist_aspect.txt");

    auto hist_len = make_histogram(len, nbins);
    save_histogram(hist_len, path_out + "hist_length.txt");
}

void load_data(const std::string& fil_stats, std::vector<double>& aspect, std::vector<double>& thetaLH, std::vector<double>& thetaLx, std::vector<double>& thetaHx) {
    std::ifstream fin; fin.open(fil_stats);
    std::string line;
    aspect.resize(0);
    thetaLH.resize(0);
    thetaLx.resize(0);
    thetaHx.resize(0);
    while(getline(fin, line)) {
        double a, b, c, d;
        std::istringstream ss(line);
        vec3D P1, P2, h;
        ss >> P1.x >> P1.y >> P1.z >> P2.x >> P2.y >> P2.z >> h.x >> h.y >> h.z >> a >> b >> c >> d;
        vec3D n = P2 - P1; normalize(n);
        //normal.push_back(n); havg.push_back(h);
        //length.push_back(a); radius.push_back(d);
        vec3D h_norm = h; normalize(h_norm);
        aspect.push_back(d / a);
        thetaLH.push_back(acos(n * h_norm));
        thetaLx.push_back(acos(n.x));
        thetaHx.push_back(acos(h_norm.x));
    }
}

void load_data(const std::string& fil_stats, std::vector<double>& aspect, std::vector<double>& len, std::vector<double>& len_proj, std::vector<double>& thetaLH, std::vector<double>& thetaLn, std::vector<double>& thetaHn, std::vector<double>& psiLHn, int axis) {
    std::ifstream fin; fin.open(fil_stats);
    std::string line;
    aspect.resize(0);
    len.resize(0);
    len_proj.resize(0);
    thetaLH.resize(0);
    thetaLn.resize(0);
    thetaHn.resize(0);
    psiLHn.resize(0);
    while(getline(fin, line)) {
        double a, b, c, d;
        std::istringstream ss(line);
        vec3D P1, P2, h;
        ss >> P1.x >> P1.y >> P1.z >> P2.x >> P2.y >> P2.z >> h.x >> h.y >> h.z >> a >> b >> c >> d;
        vec3D n = P2 - P1;
        vec3D La = P2 - P1;
        normalize(n);
        //normal.push_back(n); havg.push_back(h);
        //length.push_back(a); radius.push_back(d);
        vec3D h_norm = h; normalize(h_norm);
        aspect.push_back(d / a);
        len.push_back(a);
        thetaLH.push_back(acos(n * h_norm));
        double sign, psi;
        if (axis == 0) {
            thetaLn.push_back(acos(n.x));
            thetaHn.push_back(acos(h_norm.x));
            n.x = 0; h_norm.x = 0;
            sign = cross(n, h_norm).x;
            La.x = 0;
        }
        else if (axis == 1) {
            thetaLn.push_back(acos(n.y));
            thetaHn.push_back(acos(h_norm.y));
            n.y = 0; h_norm.y = 0;
            sign = cross(n, h_norm).y;
            La.y = 0;
        }
        else {
            thetaLn.push_back(acos(n.z));
            thetaHn.push_back(acos(h_norm.z));
            n.z = 0; h_norm.z = 0;
            sign = cross(n, h_norm).z;
            La.z = 0;
        }
        normalize(n); normalize(h_norm);
        psi = acos(clamp(n * h_norm, -1, 1));
        len_proj.push_back(std::sqrt(mag2(La)));
        if (sign >= 0) {
            psiLHn.push_back(psi);
        }
        else {
             // this is to distinguish the mutual orientation of filament and mag. field projection. Comment out if PDF(psi) is assumed symmetric around psi = pi
            //psi = 2 * M_PI - psi;
            psiLHn.push_back(psi);
        }
    }
}

void load_data_dry(const std::string& fil_stats, std::vector<double>& aspect, std::vector<double>& thetaLx) {
    std::ifstream fin; fin.open(fil_stats);
    std::string line;
    aspect.resize(0);
    thetaLx.resize(0);
    while(getline(fin, line)) {
        double a, b, c, d;
        std::istringstream ss(line);
        vec3D P1, P2;
        ss >> P1.x >> P1.y >> P1.z >> P2.x >> P2.y >> P2.z >> a >> b >> c >> d;
        vec3D n = P2 - P1; normalize(n);
        //normal.push_back(n); havg.push_back(h);
        //length.push_back(a); radius.push_back(d);
        aspect.push_back(d / a);
        thetaLx.push_back(acos(n.x));
    }
}

void load_data_dry(const std::string& fil_stats, std::vector<double>& aspect, std::vector<double>& len, std::vector<double>& thetaLn, int axis) {
    std::ifstream fin; fin.open(fil_stats);
    std::string line;
    aspect.resize(0);
    len.resize(0);
    thetaLn.resize(0);
    while(getline(fin, line)) {
        double a, b, c, d;
        std::istringstream ss(line);
        vec3D P1, P2;
        ss >> P1.x >> P1.y >> P1.z >> P2.x >> P2.y >> P2.z >> a >> b >> c >> d;
        vec3D n = P2 - P1; normalize(n);
        //normal.push_back(n); havg.push_back(h);
        //length.push_back(a); radius.push_back(d);
        aspect.push_back(d / a);
        len.push_back(a);
        if (axis == 0) {
            thetaLn.push_back(acos(n.x));
        }
        else if (axis == 1) {
            thetaLn.push_back(acos(n.y));
        }
        else {
            thetaLn.push_back(acos(n.z));
        }
    }
}

bool is_magnetized(const std::string& fil_path) {
    std::regex maPattern(R"(Ma(-?\d+\.\d+))"); // Matches "Ma" followed by a float
    std::smatch match;

    bool magnetized = false;
    
    if (std::regex_search(fil_path, match, maPattern)) {
        double maValue = std::stod(match[1].str()); // Convert the matched value to double
        magnetized = (maValue != 0.0);
    }
    else {
        std::cerr << "Error: 'Ma' pattern not found in file path!" << std::endl;
    }

    return magnetized;
}

void make_histograms(const std::string& fil_path, const std::string& path_out, const std::vector<int>& frames) {
    std::pair<double, double> bounds_halfangle = std::make_pair(0, M_PI / 2);
    std::pair<double, double> bounds_angle = std::make_pair(0, M_PI);
    std::pair<double, double> bounds_fullangle = std::make_pair(0, 2 * M_PI);
    std::pair<double, double> bounds_aspect = std::make_pair(0, 1);

    int nbins = 10;
    int nbins_ens = 25;

    std::vector<std::string> axis = {"x","y","z"};
    for (int i = 0; i <= 2; i++) {
        double umin = 1e6;
        double vmin = 1e6;
        double wmin = 1e6;
        double umax = -1e6;
        double vmax = -1e6;
        double wmax = -1e6;
        int excluded_count = 0;
        std::vector<double> aspect_ens, len_ens, len_proj_ens, thetaLH_ens, thetaLn_ens, thetaHn_ens, psiLHn_ens, u_ens, v_ens, vbar_ens, w_ens, F3D;
        std::vector<std::vector<double>> F1D;
        for (auto frame : frames) {
            std::vector<double> aspect, len, len_proj, thetaLH, thetaLn, thetaHn, psiLHn, u, v, vbar, w;
            std::string stats_name = fil_path + "DD" + pad(frame) + "/cube" + pad(frame) + "_filament_stats.txt";
            std::string output_dir = path_out + "DD" + pad(frame) + "/";
            load_data(stats_name, aspect, len, len_proj, thetaLH, thetaLn, thetaHn, psiLHn, i);

            aspect_ens.insert(aspect_ens.end(), aspect.begin(), aspect.end());
            len_ens.insert(len_ens.end(), len.begin(), len.end());
            len_proj_ens.insert(len_proj_ens.end(), len_proj.begin(), len_proj.end());
            thetaLH_ens.insert(thetaLH_ens.end(), thetaLH.begin(), thetaLH.end());
            thetaLn_ens.insert(thetaLn_ens.end(), thetaLn.begin(), thetaLn.end());
            thetaHn_ens.insert(thetaHn_ens.end(), thetaHn.begin(), thetaHn.end());
            psiLHn_ens.insert(psiLHn_ens.end(), psiLHn.begin(), psiLHn.end());

            for (int i = 0; i < thetaLn.size(); i++) {
                if ((thetaLH[i] > 0) && (thetaLH[i] < M_PI)) {
                    auto alpha = thetaLn[i] + thetaHn[i] - thetaLH[i];
                    auto beta = thetaLn[i] - thetaHn[i] + thetaLH[i];
                    auto beta2 = thetaHn[i] - thetaLn[i] + thetaLH[i];
                    u.push_back(alpha / (2 * (M_PI - thetaLH[i])));
                    v.push_back(beta / (2 * thetaLH[i]));
                    vbar.push_back(beta2 / (2 * thetaLH[i]));
                    w.push_back(2 * thetaLH[i] / M_PI);
                    umin = std::min(umin, u.back());
                    vmin = std::min(vmin, v.back());
                    wmin = std::min(wmin, w.back());
                    umax = std::max(umax, u.back());
                    vmax = std::max(vmax, v.back());
                    wmax = std::max(wmax, w.back());
                }
                else excluded_count++;
            }

            u_ens.insert(u_ens.end(), u.begin(), u.end());
            v_ens.insert(v_ens.end(), v.begin(), v.end());
            vbar_ens.insert(vbar_ens.end(), vbar.begin(), vbar.end());
            w_ens.insert(w_ens.end(), w.begin(), w.end());

            //std::cout << output_dir << " " << aspect.size() << std::endl;

            make_histograms(aspect, len, len_proj, thetaLH, thetaLn, thetaHn, psiLHn, bounds_halfangle, bounds_angle, bounds_fullangle, bounds_aspect, nbins, output_dir, axis[i]);
            make_histograms_uvw(u, v, vbar, w, nbins, path_out, axis[i]);
        }

        std::ofstream file_bounds(path_out + "uvw_" + axis[i] + "_bounds.txt");
        file_bounds << umin << " " << umax << std::endl;
        file_bounds << vmin << " " << vmax << std::endl;
        file_bounds << wmin << " " << wmax;
        file_bounds.close();

        std::ofstream file_excluded_count(path_out + "uvw_" + axis[i] + "_excluded.txt");
        file_excluded_count << excluded_count;
        file_excluded_count.close();

        //std::cout << aspect_ens.size() << std::endl;

        make_histograms(aspect_ens, len_ens, len_proj_ens, thetaLH_ens, thetaLn_ens, thetaHn_ens, psiLHn_ens, bounds_halfangle, bounds_angle, bounds_fullangle, bounds_aspect, nbins_ens, path_out, axis[i]);
        make_histograms_uvw(u_ens, v_ens, vbar_ens, w_ens, nbins, path_out, axis[i]);

        fourierLH(thetaLH_ens, path_out, 150);
        extract_fourier(u_ens, v_ens, w_ens, 32, F1D, F3D);

        save1D(F1D[0], path_out + "Fourier_u" + axis[i] + ".txt");
        save1D(F1D[1], path_out + "Fourier_v" + axis[i] + ".txt");
        save1D(F1D[2], path_out + "Fourier_w" + axis[i] + ".txt");

        save1D(F3D, path_out + "Fourier_uvw" + axis[i] + ".txt");
    }
}

void make_histograms_dry(const std::string& fil_path, const std::string& path_out, const std::vector<int>& frames) {
    std::pair<double, double> bounds_angle = std::make_pair(0, M_PI);
    std::pair<double, double> bounds_aspect = std::make_pair(0, 1);

    int nbins = 10;
    int nbins_ens = 25;

    std::vector<std::string> axis = {"x","y","z"};
    for (int i = 0; i <= 2; i++) {
        std::vector<double> aspect_ens, len_ens, thetaLn_ens;
        for (auto frame : frames) {
            std::vector<double> aspect, len, thetaLn;
            std::string stats_name = fil_path + "DD" + pad(frame) + "/cube" + pad(frame) + "_filament_stats.txt";
            std::string output_dir = path_out + "DD" + pad(frame) + "/";
            load_data_dry(stats_name, aspect, len, thetaLn, i);

            aspect_ens.insert(aspect_ens.end(), aspect.begin(), aspect.end());
            len_ens.insert(len_ens.end(), len.begin(), len.end());
            thetaLn_ens.insert(thetaLn_ens.end(), thetaLn.begin(), thetaLn.end());

            //std::cout << output_dir << " " << aspect.size() << std::endl;

            make_histograms_dry(aspect, len, thetaLn, bounds_angle, bounds_aspect, nbins, output_dir, axis[i]);
        }

        //std::cout << aspect_ens.size() << std::endl;

        make_histograms_dry(aspect_ens, len_ens, thetaLn_ens, bounds_angle, bounds_aspect, nbins_ens, path_out, axis[i]);
    }
}

int main(int argc, char *argv[]) {
    std::string fil_path, path_out;
    std::vector<int> frames;
    if (validate_input(argc, argv, fil_path, path_out, frames)) {
        std::cout << "Now analyzing " << fil_path << std::endl;
        if (is_magnetized(fil_path)) make_histograms(fil_path, path_out, frames);
        else make_histograms_dry(fil_path, path_out, frames);
        //if (is_magnetized(fil_path)) std::cout << "Is magnetized!";
        //else std::cout << "Is not magnetized!";
    }
}