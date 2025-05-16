#ifndef __MAIN_FUNS__
#define __MAIN_FUNS__

#include <iostream>
#include <fstream>
#include <functional>
#include <cstdio>
#include <map>
#include <filesystem>

#include "filament.h"
#include "fitting.h"
#include "files.h"
#include "string_pad.h"
#include "vec3D.h"
#include "analyze.h"
#include "paths_frames.h"

double get_radius(double a, double alpha);
double get_max_a(double L, double alpha);

//legacy
void extract_filaments();
void extract_stats();
void extract_all_stats_brano_windows_machine();

//test
void test_saddle();
void skeleton();

//current
void P83_stats(std::string inFilaments, std::string inDensity, std::string inHx, std::string inHy, std::string inHz, std::string outDir, std::string outName, double dL, double dR, double line_rms_threshold, double line_length_threshold, double radial_frac, double radial_fit_threshold, int resolution, bool save_failed_filaments);
void P83_stats_dry(std::string inFilaments, std::string inDensity, std::string outDir, std::string outName, double dL, double dR, double line_rms_threshold, double line_length_threshold, double radial_frac, double radial_fit_threshold, int resolution, bool save_failed_filaments);


#endif