#ifndef __ANALYZE__
#define __ANALYZE__

#include <vector>
#include <iostream>

#include "interpolation.h"
#include "vec3D.h"

std::vector<std::vector<Point3D>> analyze_pairs3D(double* pairs, int period);

inline int moveR(int i, int N) {
    int res = i + 1;
    if (res < N) return res;
    else return 0;
}

inline int moveL(int i, int N) {
    int res = i - 1;
    if (res >= 0) return res;
    else return N - 1;
}

bool is_min_x(double* field, int i, int j, int k, int N);
bool is_max_x(double* field, int i, int j, int k, int N);
bool is_min_y(double* field, int i, int j, int k, int N);
bool is_max_y(double* field, int i, int j, int k, int N);
bool is_min_z(double* field, int i, int j, int k, int N);
bool is_max_z(double* field, int i, int j, int k, int N);
bool is_extr_x(double* field, int i, int j, int k, int N);
bool is_extr_y(double* field, int i, int j, int k, int N);
bool is_extr_z(double* field, int i, int j, int k, int N);
bool is_maximum6(double* field, int i, int j, int k, int N);
bool is_maximum6(double* field, Point3D p, int N);
bool is_minimum6(double* field, int i, int j, int k, int N);
bool is_minimum6(double* field, Point3D p, int N);
bool is_saddle6(double* field, int i, int j, int k, int N);
bool is_saddle6(double* field, Point3D p, int N);

#endif