#ifndef __FIL__
#define __FIL__

#include "vec3D.h"
#include "files.h"
#include <set>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

struct Point2D {
	int x = 0;
	int y = 0;
	int next = -1;
	int prev = -1;
	Point2D operator+(Point2D const& other) const;
	Point2D operator-(Point2D const& other) const;
	Point2D operator-() const;
	int operator*(Point2D const& other) const;
};

int dist2(const Point2D& p1, const Point2D& p2);
int dist2(const Point2D& p1, const Point2D& p2, int period);

class Segment {
public:
	int start = 0;
	int end = 0;
	std::set<int> index;

	Segment() {
		start = 0; end = 0;
		index.clear();
	}
};

class Filament2D {
private:
	int threshold = 2;

	std::vector<std::set<int>> connected_components();
	void connected_components(int index, std::set<int>& to_resolve, std::vector<std::set<int>>& components);

	void flip(std::vector<std::set<int>>& components);
	void flip();

	void find_segments();
	void find_segments(int index, int segment_index, int direction, std::set<int>& resolved);

	//void continuous_segments();
	//void continuous_segments(int index, std::set<int>& traversed);

	void fuse_segments();
	int dist2(int i, int j);

	void sow_ss(int i, int j);
	void sow_ee(int i, int j);
	void sow_es(int i, int j);
	void sow_se(int i, int j);

	int validate_segment(int seg_index);
public:
	std::vector<Point2D> p;
	std::vector<Segment> segments;
	int N = 0;
	int period = 512;
	int start = 0;
	int end = 0;
	int Nseg = 0;

	Filament2D();
	Filament2D(int* x, int* y, int N);
	Filament2D(std::string filename);

	friend std::ostream& operator<<(std::ostream& os, const Filament2D& fil);

	void write(std::string filename);

	void connect();
};

class Filament3D {
private:

	std::vector<std::set<int>> connected_components();
	void connected_components(int index, std::set<int>& to_resolve, std::vector<std::set<int>>& components);

	void flip(std::vector<std::set<int>>& components);
	void flip();

	void find_segments();
	void find_segments(int index, int segment_index, int direction, std::set<int>& resolved);

	//void continuous_segments();
	//void continuous_segments(int index, std::set<int>& traversed);

	void fuse_segments();
	int dist2(int i, int j);

	void sow_ss(int i, int j);
	void sow_ee(int i, int j);
	void sow_es(int i, int j);
	void sow_se(int i, int j);

	int validate_segment(int seg_index);
public:
	std::vector<Point3D> p;
	std::vector<Segment> segments;
	int threshold = 2;
	int N = 0;
	int period = 512;
	int start = 0;
	int end = 0;
	int Nseg = 0;

	Filament3D();
	Filament3D(int* x, int* y, int* z, int N);
	Filament3D(std::string filename);
	Filament3D(const std::vector<Point3D>& pts);

	friend std::ostream& operator<<(std::ostream& os, const Filament3D& fil);

	void write(std::string filename);

	void connect();
};

std::vector<std::vector<Point3D>> load_filaments(int period, std::string filename);

#endif