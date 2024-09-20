#pragma once
#ifndef __FITTING__
#define __FITTING__

#include "vec3D.h"
#include "filament.h"
#include "interpolation.h"

#include <vector>
#include <map>
#include <random>
#include <functional>
#include <cmath>

extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_rms;
extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_rms_perc;
extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_abs;
extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max;
extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max_rel;
extern std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max_rel_perc;

double line_dist2(const vec3D& P1, const vec3D& P2, const Point3D& p);
double line_dist2_unconstrained(const vec3D& P1, const vec3D& P2, const Point3D& p);
double multilinear3D_avg(const Filament3D& fil, int segment_index, std::vector<vec3D>& points);
double multilinear3D_max(const Filament3D& fil, int segment_index, std::vector<vec3D>& points);

double multilinear3D_max(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2);
double multilinear3D_max_unconstrained(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2);
double multilinear3D_rms_unconstrained(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2);

std::map<double, int> arclen(const Filament3D& fil, int segment_index);

bool box_condition(const vec3D& p, const double& box_size);

void shuffle_point(vec3D& p, const double& sigma);
void shuffle_point(vec3D& p, const double& sigma, const double& box_size);

double iterate_distance(const Filament3D& fil, int segment_index, std::vector<vec3D>& points, int n_seg);

void extend_line(vec3D& P1, vec3D& P2, const Point3D& p1, const Point3D& p2);
vec3D find_avg(const Filament3D& fil, int start_index, int end_index);
double iterate_distance_rms(const Filament3D& fil, int start_index, int end_index, vec3D& P1, vec3D& P2);

void save_points(std::vector<vec3D>& points, std::string filename);

void initial_multi(const Filament3D& fil, int segment_index, std::map<double, int> arc_length, std::vector<vec3D>& points, int n_seg);
void find_fit_global(const Filament3D& fil, int segment_index, std::string filename);
void find_fit_progressive(const Filament3D& fil, int segment_index, std::ofstream& out, double line_rms_threshold, double line_length_threshold);
std::vector<std::pair<vec3D, vec3D>> find_fit_progressive(const Filament3D& fil, int segment_index, double line_rms_threshold, double line_length_threshold);

void sanitize_point(vec3D& p, int period);
void sanitize_point(vec3D& p, int Nx, int Ny, int Nz);

std::vector<double> spine(const vec3D& P1, const vec3D& P2, const std::vector<std::vector<std::vector<double>>>& field, int period, double dL);
std::vector<double> spine(const vec3D& P1, const vec3D& P2, double*** field, int Nx, int Ny, int Nz, double dL);
std::vector<double> spine(const vec3D& P1, const vec3D& P2, double*** field, int N, double dL);

std::vector<double> spine(const vec3D& P1, const vec3D& P2, double* field, int Nx, int Ny, int Nz, double dL);
std::vector<double> spine(const vec3D& P1, const vec3D& P2, double* field, int N, double dL);

vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, double r, double theta);
vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u);
vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, double r, double theta, int period);
vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, int period);

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double normalizing_factor, double*** field, int period);
double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double*** field, int period);

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double normalizing_factor, double* field, int period);
double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double* field, int period);

std::vector<std::vector<double>> cylinder_accumulant(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period);
std::vector<std::vector<double>> cylinder_accumulant_normalized(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period);
std::vector<std::pair<double, double>> cylinder_accumulant_avg(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period);
std::vector<std::pair<double, double>> cylinder_accumulant_avg_norm(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period);
std::vector<std::pair<double, double>> cylinder_radial_avg_norm(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period);
std::vector<std::pair<double, double>> cylinder_radial_avg_norm_constrained(const vec3D& P1, const vec3D& P2, double dL, double dr, double threshold, double*** field, int period);

std::vector<std::pair<double, double>> cylinder_radial_avg_norm_constrained(const vec3D& P1, const vec3D& P2, double dL, double dr, double threshold, double* field, int period);

void generate_params(const std::vector<double>& param, std::vector<double>& new_param, const std::vector<double>& sigma, const std::vector<std::pair<double, double>>& constraint);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<double>& starting_values, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint);
double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param);

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<double>& starting_values, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N);
double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N);
double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N);
double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N);
double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, double threshold, int& N);
double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, double threshold, int& N);

double mean_field(const vec3D& P1, const vec3D& P2, double*** field, double period, double dL);
vec3D mean_field(const vec3D& P1, const vec3D& P2, double*** Fx, double*** Fy, double*** Fz, double period, double dL);

double mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double*** field, double period, double dL, double dr, double Nphi);
vec3D mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double*** Fx, double*** Fy, double*** Fz, double period, double dL, double dr, double Nphi);

double mean_field(const vec3D& P1, const vec3D& P2, double*** field, double period, double dL);
vec3D mean_field(const vec3D& P1, const vec3D& P2, double* Fx, double* Fy, double* Fz, double period, double dL);

double mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double* field, double period, double dL, double dr, double Nphi);
vec3D mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double* Fx, double* Fy, double* Fz, double period, double dL, double dr, double Nphi);

#endif
