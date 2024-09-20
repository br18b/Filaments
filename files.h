#pragma once
#ifndef __LOAD__
#define __LOAD__

#include "vec3D.h"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdarg>

#include "H5Cpp.h"

#include <fitsio.h>

using namespace H5;

void load1D(double* field, int count, std::string filename);
void load1D(std::vector<double>& field, std::string filename);
void load2D(double** field, int countX, int countY, std::string filename);
void load2D(double** field, int count, std::string filename);
void load3D(double*** field, int countX, int countY, int countZ, std::string filename);
void load3D(double*** field, int countX, int countY, int countZ, std::string filename, double (*transform)(double));
void load3D(double*** field, int count, std::string filename);
void load3D(double*** field, int count, std::string filename, double (*transform)(double));

double* load1D(int count, std::string filename);
double** load2D(int countX, int countY, std::string filename);
double** load2D(int count, std::string filename);
double*** load3D(int countX, int countY, int countZ, std::string filename, double (*transform)(double));
double*** load3D(int count, std::string filename, double (*transform)(double));
double*** load3D(int countX, int countY, int countZ, std::string filename);
double*** load3D(int count, std::string filename);

void load_lines(std::vector<std::pair<vec3D, vec3D>>& lines, std::string filename);
std::vector<std::pair<vec3D, vec3D>> load_lines(std::string filename);

void save_spine(const std::vector<double>& spine, std::string filename);

void save_num(double x, std::string filename);
void save_num(double x, double y, std::string filename);
void save_num(double x, double y, double z, std::string filename);
void save_num(double x1, double x2, double x3, double x4, std::string filename);
void save_num(double x1, double x2, double x3, double x4, double x5, std::string filename);
void save_num(double x1, double x2, double x3, double x4, double x5, double x6, std::string filename);
void save_num(const std::vector<double>& x, std::string filename);
void save_num(const std::vector<double>& x, std::ofstream& file_out);

void save_histogram(const std::pair<std::vector<double>, std::vector<double>>& hist, std::string filename);

void save1D(const std::vector<double>& field, std::string filename);
void save2D(const std::vector<std::vector<double>>& field, std::string filename);
void save2D(const std::vector<std::pair<double, double>>& field, std::string filename);

int getSize(std::string filename, std::string groupname);
int getSize(std::string filename);
int read_from_hdf5(std::string filename, std::string groupname, std::string fieldname, long double res[], int size);
int read_from_hdf5(std::string filename, std::string fieldname, double res[], int size);

bool read_fits_array(const std::string& filename, double data[]);

#endif
