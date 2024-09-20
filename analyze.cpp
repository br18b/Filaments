#include "analyze.h"

std::vector<std::vector<Point3D>> analyze_pairs3D(double* pairs, int period) {
    std::vector<std::vector<Point3D>> res;
    for (int i = 0; i < period; i++) {
        for (int j = 0; j < period; j++) {
            for (int k = 0; k < period; k++) {
                int pair_index = pairs[index(period, i, j, k)];
                if (pair_index > 0) {
                    if (res.size() < pair_index) res.resize(pair_index);
                    Point3D p(i,j,k);
                    res[pair_index - 1].push_back(p);
                }
            }
        }
    }
    return res;
}

bool is_min_x(double* field, int i, int j, int k, int N) {
    int iL = moveL(i, N); int iR = moveR(i, N);
    double c = field[index(N, i, j, k)];
    
    return (c < field[index(N, iL, j, k)]) && (c < field[index(N, iR, j, k)]);
}

bool is_max_x(double* field, int i, int j, int k, int N) {
    int iL = moveL(i, N); int iR = moveR(i, N);
    double c = field[index(N, i, j, k)];
    
    return (c > field[index(N, iL, j, k)]) && (c > field[index(N, iR, j, k)]);
}

bool is_min_y(double* field, int i, int j, int k, int N) {
    int jL = moveL(j, N); int jR = moveR(j, N);
    double c = field[index(N, i, j, k)];
    
    return (c < field[index(N, i, jL, k)]) && (c < field[index(N, i, jR, k)]);
}

bool is_max_y(double* field, int i, int j, int k, int N) {
    int jL = moveL(j, N); int jR = moveR(j, N);
    double c = field[index(N, i, j, k)];
    
    return (c > field[index(N, i, jL, k)]) && (c > field[index(N, i, jR, k)]);
}

bool is_min_z(double* field, int i, int j, int k, int N) {
    int kL = moveL(k, N); int kR = moveR(k, N);
    double c = field[index(N, i, j, k)];
    
    return (c < field[index(N, i, j, kL)]) && (c < field[index(N, i, j, kR)]);
}

bool is_max_z(double* field, int i, int j, int k, int N) {
    int kL = moveL(k, N); int kR = moveR(k, N);
    double c = field[index(N, i, j, k)];
    
    return (c > field[index(N, i, j, kL)]) && (c > field[index(N, i, j, kR)]);
}

bool is_extr_x(double* field, int i, int j, int k, int N) {
    return is_min_x(field, i, j, k, N) || is_max_x(field, i, j, k, N);
}

bool is_extr_y(double* field, int i, int j, int k, int N) {
    return is_min_y(field, i, j, k, N) || is_max_y(field, i, j, k, N);
}

bool is_extr_z(double* field, int i, int j, int k, int N) {
    return is_min_z(field, i, j, k, N) || is_max_z(field, i, j, k, N);
}

bool is_maximum6(double* field, int i, int j, int k, int N) {
    return is_max_x(field, i, j, k, N) && is_max_y(field, i, j, k, N) && is_max_z(field, i, j, k, N);
}

bool is_maximum6(double* field, Point3D p, int N) {
    return is_maximum6(field, p.x, p.y, p.z, N);
}

bool is_minimum6(double* field, int i, int j, int k, int N) {
    return is_min_x(field, i, j, k, N) && is_min_y(field, i, j, k, N) && is_min_z(field, i, j, k, N);
}

bool is_minimum6(double* field, Point3D p, int N) {
    return is_minimum6(field, p.x, p.y, p.z, N);
}

bool is_saddle6(double* field, int i, int j, int k, int N) { // all directions are extremal but not all max or all min
    return (is_extr_x(field, i, j, k, N) && is_extr_y(field, i, j, k, N) && is_extr_z(field, i, j, k, N)) &&
    !is_maximum6(field, i, j, k, N) && !is_minimum6(field, i, j, k, N);
}

bool is_saddle6(double* field, Point3D p, int N) {
    return is_saddle6(field, p.x, p.y, p.z, N);
}