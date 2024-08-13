#include <iostream>
#include <fstream>
#include <functional>

#include "filament.h"
#include "fitting.h"
#include "files.h"

auto fit_exp_lin = [](double r, const std::vector<double>& param) {
    double a = param[0];
    return std::exp(- (r / a));
    };

auto fit_exp_sq = [](double r, const std::vector<double>& param) {
    double a = param[0];
    return std::exp(-(r / a) * (r / a));
    };

auto fit_exp_pow = [](double r, const std::vector<double>& param) {
    double a = param[0];
    double alpha = param[1];
    return std::exp(-std::pow((r / a), alpha));
    };

auto fit_parabola = [](double r, const std::vector<double>& param) {
    double a = param[0];
    return 1 - r * r / (a * a);
    };

double get_radius(double a, double alpha) {
    return a * std::tgamma(1 + 2 / alpha) / (2 * std::tgamma(1 + 1 / alpha));
}

void extract_filaments() {
    std::string name = "6_half_density_128_c0.44_filament";
    std::ofstream out; out.open("H:\\scratch\\projects\\filaments\\plots\\filaments\\" + name + "_lines.txt");
    double line_rms_threshold = 0.65;
    double line_length_threshold = 0.05;
    for (int i = 1; i <= 6397; i++) {
        std::cout << i << std::endl;
        std::string filename = name + "_" + std::to_string(i);
        //Filament2D fil("H:\\scratch\\projects\\filaments\\plots\\filaments\\" + filename); fil.period = 512;
        Filament3D fil("H:\\scratch\\projects\\filaments\\plots\\filaments\\" + filename); fil.period = 128; fil.threshold = 2;
        fil.connect();
        fil.write("H:\\scratch\\projects\\filaments\\plots\\filaments\\" + filename + "_sorted");
        for (int seg = 0; seg < fil.Nseg; seg++) {
            //std::cout << "Segment: " << seg << " " << fil.segments[seg].index.size() << std::endl;
            find_fit_progressive(fil, seg, out, line_rms_threshold, line_length_threshold);
        }
    }
    out.close();
}

void extract_stats() {
    std::string path = "H:\\scratch\\projects\\filaments\\plots\\filaments\\";
    std::string name = "6_half_density_128";
    std::string name_lines = "6_half_density_128_c0.44_filament_lines";
    int N = 128;
    double*** rho = new double** [N];
    for (int i = 0; i < N; i++) {
        rho[i] = new double* [N];
        for (int j = 0; j < N; j++) {
            rho[i][j] = new double[N];
        }
    }
    load3D(rho, N, path + name + ".txt", std::exp);
    auto lines = load_lines(path + name_lines + ".txt");
    double dL = 0.25;
    auto ca_test = cylinder_radial_avg_norm(vec3D(0,0,0), vec3D(10,20,30), dL, -1, -0.25, rho, N);
    save2D(ca_test, path + "fil_bogus_cyl_rad_avg_norm.txt");
    for (int i = 0; i < lines.size(); i++) {
        auto P1 = lines[i].first; auto P2 = lines[i].second;
        double L = std::sqrt(::dist2(P1, P2));
        auto sp = spine(P1, P2, rho, N, dL);
        // auto ca = cylinder_accumulant_normalized(P1, P2, dL, -1, -0.25, rho, N);
        auto ca = cylinder_radial_avg_norm_constrained(P1, P2, dL, -0.25, 0.75, rho, N);
        save_spine(sp, path + "line_" + std::to_string(i) + "_spine.txt");
        //save2D(ca, path + "fil_" + std::to_string(i) + "_cyl_acc.txt");
        save2D(ca, path + "fil_" + std::to_string(i) + "_cyl_rad_avg_norm.txt");
        std::vector<double> param = { 1, 2 }; auto& a = param[0]; auto& alpha = param[1]; int n = 3;
        double dist = fit_progressive(ca, fit_exp_pow, param, { {0,L}, {0, 1000} }, dist_max_rel_perc, 0.1, n);
        double r;
        if (dist < 0) r = -1;
        else r = get_radius(a, alpha);
        std::cout << i << " " << dist << " " << a << " " << r << " " << r / L << " " << n << std::endl;
        save1D(param, path + "fil_" + std::to_string(i) + "_fitpars_exp_pow.txt");
        save_num(r, L, r / L, path + "fil_" + std::to_string(i) + "_dimensions.txt");
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            delete rho[i][j];
        }
        delete[] rho[i];
    }
    delete[] rho;
}

void extract_all_stats() {
    std::string path = "H:\\scratch\\projects\\filaments\\plots\\filaments\\";
    std::string prefix = "6_half";
    int period = 128; // number of cells in each dimension of the box
    double dL = 0.25; // sampling distance for filaments (quarter of a cell)
    double dR = -0.25; // sampling distance for the radial cylindrical direction (fraction of dL if negative)
    double line_rms_threshold = 0.65; // how well we want to fit a filament with a line - smaller is stricter
    double line_length_threshold = 0.05; // min line length as a fraction of the box size
    double radial_threshold = 0.75; // how far into the radial profile we look
    auto distfun = dist_max_rel_perc; double radial_fit_threshold = 0.1; // in percent - % maximum relative metric used - this can be pretty small
    bool save_failed_filaments = false;
    auto rho = load3D(period, path + prefix + "_density_" + std::to_string(period) + ".txt", std::exp); // density was saved as logrho!
    auto Hx = load3D(period, path + prefix + "_Hx_" + std::to_string(period) + ".fits.txt");
    auto Hy = load3D(period, path + prefix + "_Hy_" + std::to_string(period) + ".fits.txt");
    auto Hz = load3D(period, path + prefix + "_Hz_" + std::to_string(period) + ".fits.txt");
    auto filament_points = load_filaments(period, path + prefix + "_density_" + std::to_string(period) + ".fits_c0.44.up.NDskl.BRK.ASMB.fits.txt");

    std::vector<double> param = { 1, 2 }; auto& a = param[0]; auto& alpha = param[1]; int n = 3; double r;
    std::ofstream out; out.open("H:\\scratch\\projects\\filaments\\plots\\filaments\\6_half_all_stats.txt");
    int fil_index = 1; int total_lines = 0;
    for (const auto& fil_pt : filament_points) {
        std::cout << "Fil. " << fil_index << "/" << filament_points.size() << "...";
        auto fil = Filament3D(fil_pt); fil.period = period;
        fil.connect();
        std::cout << " " << fil.Nseg << " seg."; int linecount = 0;
        for (int seg = 0; seg < fil.Nseg; seg++) {
            auto lines = find_fit_progressive(fil, seg, line_rms_threshold, line_length_threshold); linecount += lines.size();
            for (int i = 0; i < lines.size(); i++) {
                auto P1 = lines[i].first; auto P2 = lines[i].second;
                double L = std::sqrt(::dist2(P1, P2));
                auto filament_profile = cylinder_radial_avg_norm_constrained(P1, P2, dL, dR, radial_threshold, rho, period);
                double dist = fit_progressive(filament_profile, fit_exp_pow, param, { {0,L}, {0, 1000} }, distfun, radial_fit_threshold, n);
                auto Havg = mean_field(P1, P2, Hx, Hy, Hz, period, dL);
                if ((P2 - P1) * Havg < 0) swap(P1, P2);
                if (dist < 0) r = -1;
                else r = get_radius(a, alpha);
                if (save_failed_filaments || r > 0) save_num({ P1.x, P1.y, P1.z, P2.x, P2.y, P2.z, Havg.x, Havg.y, Havg.z, L, a, alpha, r }, out);
            }
        }
        total_lines += linecount;
        if (linecount == 0) std::cout << ", no suitable lines found";
        else if (linecount == 1) std::cout << ", 1 line (" << total_lines << " total lines)";
        else std::cout << ", " << linecount << " lines (" << total_lines << " total lines)";
        std::cout << std::endl;
        fil_index++;
    }
    out.close();
}

int main() {
    extract_all_stats();
}