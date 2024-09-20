#include "interpolation.h"

double modf(double param, int* intpart) {
	double i;
	double fracpart = modf(param, &i);
	*intpart = (int)i;
	return fracpart;
}

double get_value(const std::vector<double>& field, double x) {
	int i;
	int N = field.size();
	double u = modf(x, &i);
	int i2 = i + 1;
	if (i2 == N) i2 = 0;
	return field[i] * (1 - u) + field[i2] * u;
}

double get_value(const std::vector<std::vector<double>>& field, double x, double y) {
	int i, j;
	int Nx = field.size();
	int Ny = field[0].size();
	double u = modf(x, &i);
	double v = modf(y, &j);
	int i2 = i + 1; if (i2 == Nx) i2 = 0;
	int j2 = j + 1; if (j2 == Ny) j2 = 0;
	return field[i][j] * (1 - u) * (1 - v) + field[i2][j] * u * (1 - v) + field[i][j2] * (1 - u) * v + field[i2][j2] * u * v;
}

double get_value(const std::vector<std::vector<std::vector<double>>>& field, double x, double y, double z) {
	int i, j, k;
	int Nx = field.size();
	int Ny = field[0].size();
	int Nz = field[0][0].size();
	double u = modf(x, &i);
	double v = modf(y, &j);
	double w = modf(z, &k);
	int i2 = i + 1; if (i2 == Nx) i2 = 0;
	int j2 = j + 1; if (j2 == Ny) j2 = 0;
	int k2 = k + 1; if (k2 == Nz) k2 = 0;
	return field[i][j][k] * (1 - u) * (1 - v) * (1 - w)
		+ field[i2][j][k] * u * (1 - v) * (1 - w)
		+ field[i][j2][k] * (1 - u) * v * (1 - w)
		+ field[i][j][k2] * (1 - u) * (1 - v) * w
		+ field[i][j2][k2] * (1 - u) * v * w
		+ field[i2][j][k2] * u * (1 - v) * w
		+ field[i2][j2][k] * u * v * (1 - w)
		+ field[i2][j2][k2] * u * v * w;
}

double get_value(const std::vector<std::vector<std::vector<double>>>& field, vec3D p) {
	return get_value(field, p.x, p.y, p.z);
}

double get_value(double* field, int N, double x) {
	int i;
	double u = modf(x, &i);
	int i2 = i + 1;
	if (i2 == N) i2 = 0;
	return field[i] * (1 - u) + field[i2] * u;
}

double get_value(double** field, int Nx, int Ny, double x, double y) {
	int i, j;
	double u = modf(x, &i);
	double v = modf(y, &j);
	int i2 = i + 1; if (i2 == Nx) i2 = 0;
	int j2 = j + 1; if (j2 == Ny) j2 = 0;
	return field[i][j] * (1 - u) * (1 - v) + field[i2][j] * u * (1 - v) + field[i][j2] * (1 - u) * v + field[i2][j2] * u * v;
}

double get_value(double*** field, int Nx, int Ny, int Nz, double x, double y, double z) {
	int i, j, k;
	double u2 = modf(x, &i); double u = 1 - u2; // modf is the fractional part,
	double v2 = modf(y, &j); double v = 1 - v2; // second argument gets filled with floor
	double w2 = modf(z, &k); double w = 1 - w2; // 
	int i2 = i + 1; if (i2 == Nx) i2 = 0;
	int j2 = j + 1; if (j2 == Ny) j2 = 0;
	int k2 = k + 1; if (k2 == Nz) k2 = 0;
	return field[i][j][k] * u * v * w
		+ field[i2][j][k] * u2 * v * w
		+ field[i][j2][k] * u * v2 * w
		+ field[i][j][k2] * u * v * w2
		+ field[i][j2][k2] * u * v2 * w2
		+ field[i2][j][k2] * u2 * v * w2
		+ field[i2][j2][k] * u2 * v2 * w
		+ field[i2][j2][k2] * u2 * v2 * w2;
}

double get_value(double*** field, int Nx, int Ny, int Nz, vec3D p) {
	return get_value(field, Nx, Ny, Nz, p.x, p.y, p.z);
}

double get_value(double** field, int N, double x, double y) {
	return get_value(field, N, N, x, y);
}

double get_value(double*** field, int N, double x, double y, double z) {
	return get_value(field, N, N, N, x, y, z);
}

double get_value(double*** field, int N, vec3D p) {
	return get_value(field, N, N, N, p.x, p.y, p.z);
}

double get_value(double* field, int Nx, int Ny, int Nz, double x, double y, double z) {
	int i, j, k;
	double u2 = modf(x, &i); double u = 1 - u2; // modf is the fractional part,
	double v2 = modf(y, &j); double v = 1 - v2; // second argument gets filled with floor
	double w2 = modf(z, &k); double w = 1 - w2; // 
	int i2 = i + 1; if (i2 == Nx) i2 = 0;
	int j2 = j + 1; if (j2 == Ny) j2 = 0;
	int k2 = k + 1; if (k2 == Nz) k2 = 0;
	return field[index(Ny, Nz, i, j, k)] * u * v * w
		+ field[index(Ny, Nz, i2, j, k)] * u2 * v * w
		+ field[index(Ny, Nz, i, j2, k)] * u * v2 * w
		+ field[index(Ny, Nz, i, j, k2)] * u * v * w2
		+ field[index(Ny, Nz, i, j2, k2)] * u * v2 * w2
		+ field[index(Ny, Nz, i2, j, k2)] * u2 * v * w2
		+ field[index(Ny, Nz, i2, j2, k)] * u2 * v2 * w
		+ field[index(Ny, Nz, i2, j2, k2)] * u2 * v2 * w2;
}

double get_value(double* field, int N, double x, double y, double z) {
	return get_value(field, N, N, N, x, y, z);
}

double get_value(double* field, int Nx, int Ny, int Nz, vec3D p) {
	return get_value(field, Nx, Ny, Nz, p.x, p.y, p.z);
}

double get_value(double* field, int N, vec3D p) {
	return get_value(field, N, p.x, p.y, p.z);
}

vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int Nx, int Ny, int Nz, double x, double y, double z) {
	double fx = get_value(Fx, Nx, Ny, Nz, x, y, z);
	double fy = get_value(Fy, Nx, Ny, Nz, x, y, z);
	double fz = get_value(Fz, Nx, Ny, Nz, x, y, z);
	return vec3D(fx, fy, fz);
}

vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int N, double x, double y, double z) {
	return get_value(Fx, Fy, Fz, N, N, N, x, y, z);
}

vec3D get_value(double*** Fx, double*** Fy, double*** Fz, int N, vec3D p) {
	return get_value(Fx, Fy, Fz, N, p.x, p.y, p.z);
}

vec3D get_value(double* Fx, double* Fy, double* Fz, int Nx, int Ny, int Nz, double x, double y, double z) {
	double fx = get_value(Fx, Nx, Ny, Nz, x, y, z);
	double fy = get_value(Fy, Nx, Ny, Nz, x, y, z);
	double fz = get_value(Fz, Nx, Ny, Nz, x, y, z);
	return vec3D(fx, fy, fz);
}

vec3D get_value(double* Fx, double* Fy, double* Fz, int N, double x, double y, double z) {
	return get_value(Fx, Fy, Fz, N, N, N, x, y, z);
}

vec3D get_value(double* Fx, double* Fy, double* Fz, int N, vec3D p) {
	return get_value(Fx, Fy, Fz, N, p.x, p.y, p.z);
}