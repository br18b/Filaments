#pragma once
#ifndef __VEC__
#define __VEC__

#include <cmath>
#include <algorithm>

struct Point3D {
	int x = 0;
	int y = 0;
	int z = 0;
	int next = -1;
	int prev = -1;
	Point3D operator+(Point3D const& other) const;
	Point3D operator-(Point3D const& other) const;
	Point3D operator-() const;
	int operator*(Point3D const& other) const;

	Point3D() { x = 0; y = 0; z = 0; };
	Point3D(int X, int Y, int Z) { x = X; y = Y; z = Z; };
};

struct vec3D {
	double x = 0;
	double y = 0;
	double z = 0;
	vec3D operator+(vec3D const& other) const;
	vec3D operator-(vec3D const& other) const;
	vec3D operator+(Point3D const& other);
	vec3D operator-(Point3D const& other) const;
	vec3D operator-() const;
	vec3D operator/(double a);
	vec3D& operator*=(double a);
	vec3D& operator/=(double a);
	vec3D& operator+=(vec3D other);
	vec3D& operator-=(vec3D other);
	double operator*(vec3D const& other) const;
	vec3D() { x = 0; y = 0; z = 0; };
	vec3D(Point3D p) { x = p.x; y = p.y; z = p.z; };
	vec3D(double X, double Y, double Z) { x = X; y = Y; z = Z; };
};

vec3D operator*(vec3D v, double other);
vec3D operator*(double other, vec3D v);

vec3D cross(vec3D const& u, vec3D const& v);

double mag2(const vec3D& v);
void normalize(vec3D& v);
double dist2(const vec3D& v1, const vec3D& v2);
double dist2(const Point3D& v1, const vec3D& v2);
double dist2(const vec3D& v1, const Point3D& v2);

int dist2(const Point3D& p1, const Point3D& p2);
int dist2(const Point3D& p1, const Point3D& p2, int period);

void swap(vec3D& v1, vec3D& v2);

#endif