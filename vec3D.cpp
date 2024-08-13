#include "vec3D.h"

Point3D Point3D::operator+(Point3D const& other) const {
	return Point3D(x + other.x, y + other.y, z + other.z);
}

Point3D Point3D::operator-(Point3D const& other) const {
	return Point3D(x - other.x, y - other.y, z - other.z);
}

Point3D Point3D::operator-() const {
	return Point3D(-x, -y, -z);
}

int Point3D::operator*(Point3D const& other) const {
	return x * other.x + y * other.y + z * other.z;
}

vec3D vec3D::operator+(vec3D const& other) const {
	return vec3D(x + other.x, y + other.y, z + other.z);
}

vec3D vec3D::operator-(vec3D const& other) const {
	return vec3D(x - other.x, y - other.y, z - other.z);
}

vec3D vec3D::operator+(Point3D const& other) {
	return vec3D(x + other.x, y + other.y, z + other.z);
}

vec3D vec3D::operator-(Point3D const& other) const {
	return vec3D(x - other.x, y - other.y, z - other.z);
}

vec3D vec3D::operator-() const {
	return vec3D(-x, -y, -z);
}

vec3D vec3D::operator/(double a) {
	return vec3D(x / a, y / a, z / a);
}

vec3D& vec3D::operator*=(double a) {
	this->x = x * a;
	this->y = y * a;
	this->z = z * a;
	return *this;
}

vec3D& vec3D::operator/=(double a) {
	this->x = x / a;
	this->y = y / a;
	this->z = z / a;
	return *this;
}

vec3D& vec3D::operator+=(vec3D other) {
	this->x = x + other.x;
	this->y = y + other.y;
	this->z = z + other.z;
	return *this;
}

vec3D& vec3D::operator-=(vec3D other) {
	this->x = x - other.x;
	this->y = y - other.y;
	this->z = z - other.z;
	return *this;
}

double vec3D::operator*(vec3D const& other) const {
	return x * other.x + y * other.y + z * other.z;
}

vec3D operator*(vec3D v, double other) {
	return vec3D(v.x * other, v.y * other, v.z * other);
}

vec3D operator*(double other, vec3D v) {
	return vec3D(v.x * other, v.y * other, v.z * other);
}

vec3D cross(vec3D const& u, vec3D const& v) {
	vec3D res;
	res.x = u.y * v.z - u.z * v.y;
	res.y = u.z * v.x - u.x * v.z;
	res.z = u.x * v.y - u.y * v.x;
	return res;
}

double mag2(const vec3D& v) {
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

void normalize(vec3D& v) {
	double mag = std::sqrt(mag2(v));
	v = v / mag;
}

double dist2(const vec3D& v1, const vec3D& v2) {
	return (v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z);
}

double dist2(const Point3D& v1, const vec3D& v2) {
	return (v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z);
}

double dist2(const vec3D& v1, const Point3D& v2) {
	return (v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z);
}

int dist2(const Point3D& p1, const Point3D& p2) {
	return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z);
}

int dist2(const Point3D& p1, const Point3D& p2, int period) {
	int dx = std::abs(p1.x - p2.x); int dy = std::abs(p1.y - p2.y); int dz = std::abs(p1.z - p2.z);
	int d1 = dx * dx + dy * dy + dz * dz;
	int d2 = (dx - period) * (dx - period) + dy * dy + dz * dz;
	int d3 = dx * dx + (dy - period) * (dy - period) + dz * dz;
	int d4 = dx * dx + dy * dy + (dz - period) * (dz - period);
	int d5 = (dx - period) * (dx - period) + (dy - period) * (dy - period) + dz * dz;
	int d6 = (dx - period) * (dx - period) + dy * dy + (dz - period) * (dz - period);
	int d7 = dx * dx + (dy - period) * (dy - period) + (dz - period) * (dz - period);
	int d8 = (dx - period) * (dx - period) + (dy - period) * (dy - period) + (dz - period) * (dz - period);
	return std::min({ d1, d2, d3, d4, d5, d6, d7, d8 });
}

void swap(vec3D& v1, vec3D& v2) {
	auto v = v1;
	v1 = v2;
	v2 = v;
}