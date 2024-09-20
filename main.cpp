#include <iostream>
#include <fstream>
#include <functional>
#include <cstdio>

#include "filament.h"
#include "fitting.h"
#include "files.h"
#include "string_pad.h"
#include "vec3D.h"
#include "analyze.h"
#include "paths_frames.h"

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

double get_max_a(double L, double alpha) {
    return 2 * L * std::tgamma(1 + 1 / alpha) / std::tgamma(1 + 2 / alpha);
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
        double L = std::sqrt(::dist2(P1, P2)); double amax = get_max_a(L, 3);
        auto sp = spine(P1, P2, rho, N, dL);
        // auto ca = cylinder_accumulant_normalized(P1, P2, dL, -1, -0.25, rho, N);
        auto ca = cylinder_radial_avg_norm_constrained(P1, P2, dL, -0.25, 0.75, rho, N);
        save_spine(sp, path + "line_" + std::to_string(i) + "_spine.txt");
        //save2D(ca, path + "fil_" + std::to_string(i) + "_cyl_acc.txt");
        save2D(ca, path + "fil_" + std::to_string(i) + "_cyl_rad_avg_norm.txt");
        std::vector<double> param = { 1, 2 }; auto& a = param[0]; auto& alpha = param[1]; int n = 3;
        double dist = fit_progressive(ca, fit_exp_pow, param, { {0, amax}, {0, 1000} }, dist_max_rel_perc, 0.1, n);
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

void extract_all_stats_brano_windows_machine() {
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
                double L = std::sqrt(::dist2(P1, P2)); double amax = get_max_a(L, 3);
                auto filament_profile = cylinder_radial_avg_norm_constrained(P1, P2, dL, dR, radial_threshold, rho, period);
                double dist = fit_progressive(filament_profile, fit_exp_pow, param, { {0, amax}, {0, 1000} }, distfun, radial_fit_threshold, n);
                auto Havg = mean_field(P1, P2, Hx, Hy, Hz, period, dL);
                auto rhoavg = mean_field(P1, P2, rho, period, dL);
                if ((P2 - P1) * Havg < 0) swap(P1, P2);
                if (dist < 0) r = -1;
                else r = get_radius(a, alpha);
                if (save_failed_filaments || r > 0) save_num({ P1.x, P1.y, P1.z, P2.x, P2.y, P2.z, Havg.x, Havg.y, Havg.z, L, a, alpha, r, rhoavg }, out);
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

void P71_stats(std::string path, std::string path2, std::string filament_filename, std::string stats_out_filename, std::string prefix, int frame, std::string pers_threshold) {
    std::string frame_str = string_pad_left(std::to_string(frame), "0", 4);
    int field_size = getSize(path + "data" + frame_str + ".cube.h5");
    int period = std::cbrt(field_size); // number of cells in each dimension of the box
    double dL = 0.25; // sampling distance for filaments (quarter of a cell)
    double dR = -0.25; // sampling distance for the radial cylindrical direction (fraction of dL if negative)
    double line_rms_threshold = 0.65; // how well we want to fit a filament with a line - smaller is stricter
    double line_length_threshold = 0.05; // min line length as a fraction of the box size
    double radial_threshold = 0.75; // how far into the radial profile we look
    auto distfun = dist_max_rel_perc; double radial_fit_threshold = 0.1; // in percent - % maximum relative metric used - this can be pretty small
    bool save_failed_filaments = false;
    
    double* rho = new double[field_size];
    double* Hx = new double[field_size];
    double* Hy = new double[field_size];
    double* Hz = new double[field_size];
    double* filaments = new double[field_size];
    read_from_hdf5(path + "data" + frame_str + ".cube.h5", "Density", rho, field_size);
    read_from_hdf5(path + "data" + frame_str + ".cube.h5", "Bx", Hx, field_size);
    read_from_hdf5(path + "data" + frame_str + ".cube.h5", "By", Hy, field_size);
    read_from_hdf5(path + "data" + frame_str + ".cube.h5", "Bz", Hz, field_size);
    read_fits_array(path2 + filament_filename, filaments);
    auto filament_points = load_filaments(filaments, period);
    
    std::vector<double> param = { 1, 2 }; auto& a = param[0]; auto& alpha = param[1]; int n = 3; double r;
    std::ofstream out; out.open(path2 + stats_out_filename);
    int fil_index = 1; int total_lines = 0;
    for (const auto& fil_pt : filament_points) {
        std::cout << "Fil. " << fil_index << "/" << filament_points.size() << ":";
        auto fil = Filament3D(fil_pt); fil.period = period;
        fil.connect();
        std::cout << " " << fil.Nseg << " seg."; int linecount = 0;
        for (int seg = 0; seg < fil.Nseg; seg++) {
            auto lines = find_fit_progressive(fil, seg, line_rms_threshold, line_length_threshold); linecount += lines.size();
            for (int i = 0; i < lines.size(); i++) {
                auto P1 = lines[i].first; auto P2 = lines[i].second;
                double L = std::sqrt(::dist2(P1, P2)); double amax = get_max_a(L, 3);
                auto filament_profile = cylinder_radial_avg_norm_constrained(P1, P2, dL, dR, radial_threshold, rho, period);
                double dist = fit_progressive(filament_profile, fit_exp_pow, param, { {0, amax}, {0, 1000} }, distfun, radial_fit_threshold, n);
                vec3D Havg = mean_field(P1, P2, Hx, Hy, Hz, period, dL);
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

    std::cout << total_lines << " total lines" << std::endl;

    out.close();
    delete rho; delete Hx; delete Hy; delete Hz; delete filaments;
}

void P71_stats(std::string prefix, int frame, std::string pers_threshold) {
    std::string frame_str = string_pad_left(std::to_string(frame), "0", 4);
    std::string path = downsample_path + prefix + "/";
    std::string path2 = P71_path + prefix + "/";
    std::string filament_filename = "data" + frame_str + ".filaments_c" + pers_threshold + ".fits";
    std::string stats_out_filename = "data" + frame_str + "_filament_c" + pers_threshold + "_stats.txt";
    P71_stats(path, path2, filament_filename, stats_out_filename, prefix, frame, pers_threshold);
}

void test_open_hdf5() {
    std::string filename = "data0045.cube.h5";
    std::string fieldname = "Bx";
    int field_size = getSize(filename); std::cout << "Size (should be 128**3): " << field_size << std::endl;
    double* Bx = new double[field_size];

    read_from_hdf5(filename, fieldname, Bx, field_size);
    std::cout << "Bx[0,0,0]: " << Bx[0] << " " << get_value(Bx, 128, 0, 0, 0) << std::endl;

    delete Bx;
}

void test_saddle() {
    int N = 3;
    double* field = new double[N * N * N];
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                field[index(i,j,k,N)] = 0;
            }
        }
    }
    field[index(N,1,1,1)] = 1;

    std::cout << is_max_x(field, 1, 1, 1, N) << std::endl;
    
    field[index(N,0,1,1)] = 2;
    field[index(N,2,1,1)] = 2;
    
    std::cout << is_min_x(field, 1, 1, 1, N) << std::endl;
    std::cout << is_saddle6(field, 1, 1, 1, N) << std::endl;

    delete field;
}

void skeleton() {
    int N = 128;
    double* pairs = new double[N * N * N];
    double* rho = new double[N * N * N];
    read_fits_array("6_half_density_128.fits", rho);
    read_fits_array("6_half_density_128.fits.up.NDskl.fits", pairs);
    auto pair_array = analyze_pairs3D(pairs, N);

    std::cout << pair_array.size() << std::endl;
    std::cout << pair_array[0].size() << std::endl;
    int ind = 0; int total = 0;
    for (auto& pa: pair_array) {
        std::string line; bool write1 = false; bool write2 = false;
        line = line + std::to_string(ind) + " ";
        int index_max = 0; int index_saddle = 0;
        int index2 = 0;
        for (auto& p: pa) {
            bool ismax = is_maximum6(rho, p, N);
            bool ismin = is_minimum6(rho, p, N);
            bool issaddle = is_saddle6(rho, p, N);
            if (ismax) {
                line = line + "max ";
                write1 = true;
                index_max = index2;

            }
            else if (ismin) line = line + "min ";
            else if (issaddle) {
                line = line + "s ";
                write2 = true;
                index_saddle = index2;
            }
            else line = line + "0 ";
            index2++;
        }
        if (write1 && write2) {
            std::cout << line << " " << rho[index(N, pa[index_max])] - rho[index(N, pa[index_saddle])] << std::endl;
            total++;
        }
        ind++;
    }
    std::cout << total << "/" << pair_array.size() << std::endl;
    delete pairs; delete rho;
}

int main(int argc, char *argv[]) {
    //P71_stats(downsample_path + "6_half/", P71_path + "6_half/test/", "6_half_density_128.fits_c1t0.up.NDskl.BRK.ASMB.fits", "6_half_c1t0_stats.txt", "6_half", 45, "1");
    //P71_stats(downsample_path + "6_half/", P71_path + "6_half/test/", "6_half_density_128.fits_c2t0.up.NDskl.BRK.ASMB.fits", "6_half_c2t0_stats.txt", "6_half", 45, "1");
    std::string sim = "1_half";
    std::string pers_threshold = "0.5";
    if (argc > 1) sim = argv[1];
    int frame = last_frames[sim];
    if (argc > 2) frame = std::stoi(argv[2]);
    if (argc > 3) pers_threshold = argv[3];
    P71_stats(sim, frame, pers_threshold);
}