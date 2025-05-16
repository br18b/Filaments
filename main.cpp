#include <iostream>
#include <fstream>
#include <functional>
#include <cstdio>
#include <map>
#include <filesystem>


#include "filaments_main_functions.h"
#include "argparser.h"
#include "options.h"

int main(int argc, char *argv[]) {
    std::cout << "\033[?7l";
    std::string inFilaments, inDensity, inHx, inHy, inHz, mainName; bool magnetized; std::string err_message;
    if (check_pars(argc, argv, inFilaments, inDensity, inHx, inHy, inHz, mainName, magnetized, err_message)) {
        std::string outDir, outName; double dL, dR, line_rms_threshold, line_length_threshold, radial_frac, radial_fit_threshold; int resolution; bool save_failed_filaments;
        handle_pars(argc, argv, mainName, outDir, outName, dL, dR, line_rms_threshold, line_length_threshold, radial_frac, radial_fit_threshold, resolution, save_failed_filaments);
        if (magnetized)
            P83_stats(inFilaments, inDensity, inHx, inHy, inHz, outDir, outName, dL, dR, line_rms_threshold, line_length_threshold, radial_frac, radial_fit_threshold, resolution, save_failed_filaments);
        else
            P83_stats_dry(inFilaments, inDensity, outDir, outName, dL, dR, line_rms_threshold, line_length_threshold, radial_frac, radial_fit_threshold, resolution, save_failed_filaments);
    }
    else {
        std::cout << err_message << std::endl;
    }
    std::cout << "\033[?7h";
    return 0;
}