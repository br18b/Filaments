#pragma once
#ifndef __INTERP__
#define __INTERP__

#include <vector>
#include <math.h>

#include "vec3D.h"

double modf(double param, int* intpart);

double get_value(const std::vector<double>& field, double x);
double get_value(const std::vector<std::vector<double>>& field, double x, double y);
double get_value(const std::vector<std::vector<std::vector<double>>>& field, double x, double y, double z);

double get_value(const std::vector<std::vector<std::vector<double>>>& field, vec3D P);

double get_value(double* field, int Nx, double x);
double get_value(double** field, int Nx, int Ny, double x, double y);
double get_value(double*** field, int Nx, int Ny, int Nz, double x, double y, double z);

double get_value(double*** field, int Nx, int Ny, int Nz, vec3D p);

double get_value(double** field, int N, double x, double y);
double get_value(double*** field, int N, double x, double y, double z);
double get_value(double*** field, int N, vec3D p);

vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int Nx, int Ny, int Nz, double x, double y, double z);
vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int N, double x, double y, double z);
vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int N, vec3D p);

#endif