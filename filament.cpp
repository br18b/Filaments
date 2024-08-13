#include "filament.h"

Point2D Point2D::operator+(Point2D const& other) const {
	Point2D res;
	res.x = x + other.x;
	res.y = y + other.y;
	return res;
}

Point2D Point2D::operator-(Point2D const& other) const {
	Point2D res;
	res.x = x - other.x;
	res.y = y - other.y;
	return res;
}

Point2D Point2D::operator-() const {
	Point2D res;
	res.x = -x;
	res.y = -y;
	return res;
}

int Point2D::operator*(Point2D const& other) const {
	return x * other.x + y * other.y;
}

Filament2D::Filament2D() {
	N = 0;
	p.resize(0);
}

Filament3D::Filament3D() {
	N = 0;
	p.resize(0);
}

Filament2D::Filament2D(int* x, int* y, int n) {
	N = n;
	for (int i = 0; i < N; i++) {
		Point2D pt;
		pt.x = x[i];
		pt.y = y[i];

		p.push_back(pt);
	}
}

Filament3D::Filament3D(int* x, int* y, int* z, int n) {
	N = n;
	for (int i = 0; i < N; i++) {
		Point3D pt;
		pt.x = x[i];
		pt.y = y[i];
		pt.z = z[i];

		p.push_back(pt);
	}
}

Filament3D::Filament3D(const std::vector<Point3D>& pts) {
	N = pts.size();
	p = pts;
}

Filament2D::Filament2D(std::string filename) {
	std::ifstream file_in; file_in.open(filename + ".txt", std::ios_base::in);
	double x_in, y_in;
	while (file_in >> x_in >> y_in) {
		Point2D pt;
		pt.x = x_in;
		pt.y = y_in;
		p.push_back(pt);
	}
	N = p.size();
	file_in.close();
}

Filament3D::Filament3D(std::string filename) {
	std::ifstream file_in; file_in.open(filename + ".txt", std::ios_base::in);
	double x_in, y_in, z_in;
	while (file_in >> x_in >> y_in >> z_in) {
		Point3D pt;
		pt.x = x_in;
		pt.y = y_in;
		pt.z = z_in;
		p.push_back(pt);
	}
	N = p.size();
	file_in.close();
}

std::ostream& operator<<(std::ostream& os, const Filament2D& fil) {
	int index = fil.start;
	while (index != fil.end) {
		os << fil.p[index].x << " " << fil.p[index].y << std::endl;
		index = fil.p[index].next;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const Filament3D& fil) {
	int index = fil.start;
	while (index != fil.end) {
		os << fil.p[index].x << " " << fil.p[index].y << " " << fil.p[index].z << std::endl;
		index = fil.p[index].next;
	}
	return os;
}

void Filament2D::write(std::string filename) {
	int seg_index = 1;
	for (int i = 0; i < Nseg; i++) {
		auto& seg = segments[i];
		if (seg.index.size() >= 3) {
			std::ofstream out; out.open(filename + "_" + std::to_string(seg_index) + ".txt", std::ios_base::out);
			int index = seg.start;
			while (index >= 0) {
				out << p[index].x << " " << p[index].y << std::endl;
				index = p[index].next;
			}
			out.close();
			seg_index++;
		}
	}
}

void Filament3D::write(std::string filename) {
	int seg_index = 1;
	for (int i = 0; i < Nseg; i++) {
		auto& seg = segments[i];
		if (seg.index.size() >= 3) {
			std::ofstream out; out.open(filename + "_" + std::to_string(seg_index) + ".txt", std::ios_base::out);
			int index = seg.start;
			while (index >= 0) {
				out << p[index].x << " " << p[index].y << " " << p[index].z << std::endl;
				index = p[index].next;
			}
			out.close();
			seg_index++;
		}
	}
}

int dist2(const Point2D& p1, const Point2D& p2, int period) {
	int d1 = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
	int d2 = (std::abs(p1.x - p2.x) - period) * (std::abs(p1.x - p2.x) - period) + (p1.y - p2.y) * (p1.y - p2.y);
	int d3 = (p1.x - p2.x) * (p1.x - p2.x) + (std::abs(p1.y - p2.y) - period) * (std::abs(p1.y - p2.y) - period);
	int d4 = (std::abs(p1.x - p2.x) - period) * (std::abs(p1.x - p2.x) - period) + (std::abs(p1.y - p2.y) - period) * (std::abs(p1.y - p2.y) - period);
	return std::min({ d1, d2, d3, d4 });
}

int Filament2D::dist2(int i, int j) {
	return ::dist2(p[i], p[j]);
}

int Filament3D::dist2(int i, int j) {
	return ::dist2(p[i], p[j]);
}

int dist2(const Point2D& p1, const Point2D& p2) {
	return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

void Filament2D::find_segments(int index, int segment_index, int direction, std::set<int>& to_resolve) {
	auto& seg = segments[segment_index];
	bool found = true;
	while (found) {
		seg.index.insert(index);
		found = false;
		to_resolve.erase(index);
		for (auto& i : to_resolve) {
			int d2 = dist2(index, i);
			if (d2 <= 1) {
				if (direction > 0) {
					p[index].next = i;
					p[i].prev = index;
				}
				else {
					p[index].prev = i;
					p[i].next = index;
				}
				index = i;
				found = true;
				break;
			}
		}
	}
	if (direction > 0) seg.end = index;
	else seg.start = index;
}

void Filament3D::find_segments(int index, int segment_index, int direction, std::set<int>& to_resolve) {
	auto& seg = segments[segment_index];
	bool found = true;
	while (found) {
		seg.index.insert(index);
		found = false;
		to_resolve.erase(index);
		for (auto& i : to_resolve) {
			int d2 = dist2(index, i);
			if (d2 <= 1) {
				if (direction > 0) {
					p[index].next = i;
					p[i].prev = index;
				}
				else {
					p[index].prev = i;
					p[i].next = index;
				}
				index = i;
				found = true;
				break;
			}
		}
	}
	if (direction > 0) seg.end = index;
	else seg.start = index;
}

int Filament2D::validate_segment(int seg_index) {
	const auto& seg = segments[seg_index];
	if (seg.index.size() == 0) {
		if ((seg.start < 0) && (seg.end < 0)) return 0;
		else {
			if ((seg.start >= 0) && (seg.end >= 0)) {
				std::cout << "Segment contains no points, but both start and end point to something: " << seg.start << ", " << seg.end << std::endl;
				return -1;
			}
			if (seg.start >= 0) {
				std::cout << "Segment contains no points, but the start points to something: " << seg.start << std::endl;
				return -1;
			}
			else {
				std::cout << "Segment contains no points, but the end points to something: " << seg.end << std::endl;
				return -1;
			}
		}
	}
	else {
		int start = seg.start;
		int end = seg.end;
		if (p[start].prev >= 0) {
			std::cout << "The starting point contains previous: " << p[start].prev << std::endl;
			return -1;
		}
		else if (p[end].next >= 0) {
			std::cout << "The ending point contains next: " << p[end].next << std::endl;
			return -1;
		}
		else {
			std::set<int> indices;
			int index = start;
			while (index >= 0) {
				if (indices.find(index) != indices.end()) {
					std::cout << "There's a loop within the points!" << std::endl;
					return -1;
				}
				else {
					indices.insert(index);
					index = p[index].next;
				}
			}
			if (indices == seg.index) {
				return 0;
			}
			else {
				std::cout << "The looped points do not match the segment list!" << std::endl;
				return -1;
			}
		}
	}
}

int Filament3D::validate_segment(int seg_index) {
	const auto& seg = segments[seg_index];
	if (seg.index.size() == 0) {
		if ((seg.start < 0) && (seg.end < 0)) return 0;
		else {
			if ((seg.start >= 0) && (seg.end >= 0)) {
				std::cout << "Segment contains no points, but both start and end point to something: " << seg.start << ", " << seg.end << std::endl;
				return -1;
			}
			if (seg.start >= 0) {
				std::cout << "Segment contains no points, but the start points to something: " << seg.start << std::endl;
				return -1;
			}
			else {
				std::cout << "Segment contains no points, but the end points to something: " << seg.end << std::endl;
				return -1;
			}
		}
	}
	else {
		int start = seg.start;
		int end = seg.end;
		if (p[start].prev >= 0) {
			std::cout << "The starting point contains previous: " << p[start].prev << std::endl;
			return -1;
		}
		else if (p[end].next >= 0) {
			std::cout << "The ending point contains next: " << p[end].next << std::endl;
			return -1;
		}
		else {
			std::set<int> indices;
			int index = start;
			while (index >= 0) {
				if (indices.find(index) != indices.end()) {
					std::cout << "There's a loop within the points! (forward)" << std::endl;
					return -1;
				}
				else {
					indices.insert(index);
					index = p[index].next;
				}
			}
			if (indices != seg.index) {
				std::cout << "The looped points do not match the segment list! (forward)" << std::endl;
				return -1;
			}
			else {
				index = end;
				indices.clear();
				while (index >= 0) {
					if (indices.find(index) != indices.end()) {
						std::cout << "There's a loop within the points! (backward)" << std::endl;
						return -1;
					}
					else {
						indices.insert(index);
						index = p[index].prev;
					}
				}
				if (indices != seg.index) {
					std::cout << "The looped points do not match the segment list! (backward)" << std::endl;
					return -1;
				}
				else {
					return 0;
				}
			}
		}
	}
}

void Filament2D::find_segments() {
	std::set<int> to_resolve;
	for (int i = 0; i < N; i++) to_resolve.insert(i);
	while (!to_resolve.empty()) {
		int index = *to_resolve.begin();
		Segment seg; seg.start = index; seg.end = index;
		segments.push_back(seg); Nseg++;
		find_segments(index, Nseg - 1, 1, to_resolve);
		to_resolve.insert(index);
		find_segments(index, Nseg - 1, -1, to_resolve);
		validate_segment(Nseg - 1);
	}
}

void Filament3D::find_segments() {
	std::set<int> to_resolve;
	for (int i = 0; i < N; i++) to_resolve.insert(i);
	while (!to_resolve.empty()) {
		int index = *to_resolve.begin();
		Segment seg; seg.start = index; seg.end = index;
		segments.push_back(seg); Nseg++;
		find_segments(index, Nseg - 1, 1, to_resolve);
		to_resolve.insert(index);
		find_segments(index, Nseg - 1, -1, to_resolve);
		validate_segment(Nseg - 1);
	}
}

void Filament2D::connected_components(int index, std::set<int>& to_resolve, std::vector<std::set<int>>& components) {
	if (to_resolve.find(index) != to_resolve.end()) {
		to_resolve.erase(index);
		int component_index = 0;
		for (auto& comp : components) {
			if (comp.find(index) != comp.end()) {
				break;
			}
			else component_index++;
		}
		for (int i = 0; i < N; i++) {
			int d2 = ::dist2(p[index], p[i]);
			if (d2 <= 2) {
				components[component_index].insert(i);
				connected_components(i, to_resolve, components);
			}
		}
	}
}

void Filament3D::connected_components(int index, std::set<int>& to_resolve, std::vector<std::set<int>>& components) {
	if (to_resolve.find(index) != to_resolve.end()) {
		to_resolve.erase(index);
		int component_index = 0;
		for (auto& comp : components) {
			if (comp.find(index) != comp.end()) {
				break;
			}
			else component_index++;
		}
		for (int i = 0; i < N; i++) {
			int d2 = ::dist2(p[index], p[i]);
			if (d2 <= 3) {
				components[component_index].insert(i);
				connected_components(i, to_resolve, components);
			}
		}
	}
}

std::vector<std::set<int>> Filament2D::connected_components() {
	std::set<int> to_resolve;
	for (int i = 0; i < N; i++) {
		to_resolve.insert(i);
	}
	std::vector<std::set<int>> components;
	while (!to_resolve.empty()) {
		int index = *to_resolve.begin();
		std::set<int> connected; connected.insert(index);
		components.push_back(connected);
		connected_components(index, to_resolve, components);
	}
	return components;
}

std::vector<std::set<int>> Filament3D::connected_components() {
	std::set<int> to_resolve;
	for (int i = 0; i < N; i++) {
		to_resolve.insert(i);
	}
	std::vector<std::set<int>> components;
	while (!to_resolve.empty()) {
		int index = *to_resolve.begin();
		std::set<int> connected; connected.insert(index);
		components.push_back(connected);
		connected_components(index, to_resolve, components);
	}
	return components;
}

void Filament2D::flip(std::vector<std::set<int>>& components) {
	/* find the largest component - we'll keep that one within the frame */
	int largest_index = 0; int largest_size = 0;
	for (int i = 0; i < components.size(); i++) {
		if (components[i].size() > largest_size) {
			largest_size = components[i].size();
			largest_index = i;
		}
	}
	auto& largest_component = components[largest_index];

	/* flip all other components till they match up with the largest one */
	for (int i = 0; i < components.size(); i++) {
		if (i != largest_index) {
			auto& component = components[i];
			/* find the sewing point */
			bool found = false;
			for (auto& j : largest_component) {
				if (found) break;
				for (auto& k : component) {
					int d2 = ::dist2(p[j], p[k], period);
					if (d2 <= 2) { /* found the sewing point! */
						int dx = std::abs(p[j].x - p[k].x);
						int dy = std::abs(p[j].y - p[k].y);
						if ((std::abs(dx - period) <= 1) && (std::abs(dy - period) <= 1)) { /* diagonal */
							if (p[k].x <= 1) dx = period;
							else dx = -period;
							if (p[k].y <= 1) dy = period;
							else dy = -period;
							for (auto& l : component) {
								p[l].x += dx;
								p[l].y += dy;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
						else if (std::abs(dx - period) <= 1) { /* left/right */
							if (p[k].x <= 1) dx = period;
							else dx = -period;
							for (auto& l : component) {
								p[l].x += dx;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
						else if (std::abs(dy - period) <= 1) { /* up/down */
							if (p[k].y <= 1) dy = period;
							else dy = -period;
							for (auto& l : component) {
								p[l].y += dy;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
						else {
							std::cout << "Trouble in paradise!";
						}
					}
				}
			}
		}
	}
}

void Filament3D::flip(std::vector<std::set<int>>& components) {
	/* find the largest component - we'll keep that one within the frame */
	int largest_index = 0; int largest_size = 0;
	int flip_threshold = 2;
	for (int i = 0; i < components.size(); i++) {
		if (components[i].size() > largest_size) {
			largest_size = components[i].size();
			largest_index = i;
		}
	}
	auto& largest_component = components[largest_index];

	/* flip all other components till they match up with the largest one */
	for (int i = 0; i < components.size(); i++) {
		if (i != largest_index) {
			auto& component = components[i];
			/* find the sewing point */
			bool found = false;
			for (auto& j : largest_component) {
				if (found) break;
				for (auto& k : component) {
					int d2 = ::dist2(p[j], p[k], period);
					if (d2 <= 3 * flip_threshold * flip_threshold) { /* found the sewing point! */
						int dx = std::abs(p[j].x - p[k].x);
						int dy = std::abs(p[j].y - p[k].y);
						int dz = std::abs(p[j].z - p[k].z);
						if (std::abs(dx - period) <= flip_threshold) { /* up/down */
							if (p[k].x <= 1) dx = period;
							else dx = -period;
							for (auto& l : component) {
								p[l].x += dx;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
						if (std::abs(dy - period) <= flip_threshold) { /* up/down */
							if (p[k].y <= 1) dy = period;
							else dy = -period;
							for (auto& l : component) {
								p[l].y += dy;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
						if (std::abs(dz - period) <= flip_threshold) { /* up/down */
							if (p[k].z <= 1) dz = period;
							else dz = -period;
							for (auto& l : component) {
								p[l].z += dz;
								largest_component.insert(l);
							}
							found = true;
							break;
						}
					}
				}
			}
		}
	}
}

void Filament2D::flip() {
	auto components = connected_components();
	flip(components);
}

void Filament3D::flip() {
	auto components = connected_components();
	flip(components);
}

/*
void Filament2D::continuous_segments(int index, std::set<int>& traversed) {
	if (traversed.find(index) == traversed.end()) {
		traversed.insert(index);
		int next = p[index].next;
		int prev = p[index].prev;
		int segment_index = 0;
		for (auto& seg : segments) {
			if (seg.index.find(index) != seg.index.end()) break;
			else segment_index++;
		}
		auto& seg = segments[segment_index];

		if (next >= 0) {
			int d2 = dist2(index, next);
			if (d2 <= threshold) {
				seg.index.insert(next);
				continuous_segments(next, traversed);
			}
			else {
				seg.end = index;
			}
		}

		if (prev >= 0) {
			int d2 = dist2(index, prev);
			if (d2 <= threshold) {
				seg.index.insert(prev);
				continuous_segments(prev, traversed);
			}
			else {
				seg.start = index;
			}
		}
	}
}

void Filament2D::continuous_segments() {
	std::set<int> traversed;
	while (traversed.size() < N) {
		int index = 0;
		for (int i = 0; i < N; i++) {
			if (traversed.find(i) == traversed.end()) {
				segments.push_back(Segment2D()); segments.back().index.insert(i); Nseg++;
				segments.back().start = i; segments.back().end = i;
				index = i;
				break;
			}
		}
		continuous_segments(index, traversed);
	}
	for (int i = 0; i < Nseg; i++) {
		p[segments[i].start].prev = -1;
		p[segments[i].end].next = -1;
	}
}
*/

void Filament2D::sow_ss(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];
	int size1 = seg1.index.size();
	int size2 = seg1.index.size();

	int start1 = seg1.start;
	int start2 = seg2.start;

	{
		int is_valid = validate_segment(i);
		if (is_valid != 0) {
			std::cout << "oops!" << std::endl;
		}
	}

	int index = start2;
	while (index >= 0) {
		seg1.index.insert(index);
		int next = p[index].next;
		int prev = p[index].prev;

		p[index].next = prev;
		p[index].prev = next;

		index = next;
	}
	p[start1].prev = start2;
	p[start2].next = start1;
	seg1.start = seg2.end;

	{
		int is_valid = validate_segment(i);
		if (is_valid != 0) {
			std::cout << "oops!" << std::endl;
		}
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament2D::sow_ee(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int end1 = seg1.end;
	int end2 = seg2.end;

	int index = end2;
	while (index >= 0) {
		seg1.index.insert(index);
		int next = p[index].next;
		int prev = p[index].prev;

		p[index].next = prev;
		p[index].prev = next;

		index = prev;
	}
	p[end1].next = end2;
	p[end2].prev = end1;
	seg1.end = seg2.start;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament2D::sow_es(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int start2 = seg2.start;
	int end1 = seg1.end;

	int index = start2;
	while (index >= 0) {
		seg1.index.insert(index);
		index = p[index].next;
	}
	p[end1].next = start2;
	p[start2].prev = end1;
	seg1.end = seg2.end;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament2D::sow_se(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int start1 = seg1.start;
	int end2 = seg2.end;

	int index = end2;
	while (index >= 0) {
		seg1.index.insert(index);
		index = p[index].prev;
	}
	seg1.start = seg2.start;
	p[start1].prev = end2;
	p[end2].next = start1;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament3D::sow_ss(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];
	int size1 = seg1.index.size();
	int size2 = seg1.index.size();

	int start1 = seg1.start;
	int start2 = seg2.start;

	{
		int is_valid = validate_segment(i);
		if (is_valid != 0) {
			std::cout << "oops!" << std::endl;
		}
	}

	int index = start2;
	while (index >= 0) {
		seg1.index.insert(index);
		int next = p[index].next;
		int prev = p[index].prev;

		p[index].next = prev;
		p[index].prev = next;

		index = next;
	}
	p[start1].prev = start2;
	p[start2].next = start1;
	seg1.start = seg2.end;

	{
		int is_valid = validate_segment(i);
		if (is_valid != 0) {
			std::cout << "oops!" << std::endl;
		}
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament3D::sow_ee(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int end1 = seg1.end;
	int end2 = seg2.end;

	int index = end2;
	while (index >= 0) {
		seg1.index.insert(index);
		int next = p[index].next;
		int prev = p[index].prev;

		p[index].next = prev;
		p[index].prev = next;

		index = prev;
	}
	p[end1].next = end2;
	p[end2].prev = end1;
	seg1.end = seg2.start;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament3D::sow_es(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int start2 = seg2.start;
	int end1 = seg1.end;

	int index = start2;
	while (index >= 0) {
		seg1.index.insert(index);
		index = p[index].next;
	}
	p[end1].next = start2;
	p[start2].prev = end1;
	seg1.end = seg2.end;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament3D::sow_se(int i, int j) {
	auto& seg1 = segments[i];
	auto& seg2 = segments[j];

	int start1 = seg1.start;
	int end2 = seg2.end;

	int index = end2;
	while (index >= 0) {
		seg1.index.insert(index);
		index = p[index].prev;
	}
	seg1.start = seg2.start;
	p[start1].prev = end2;
	p[end2].next = start1;

	int is_valid = validate_segment(i);
	if (is_valid != 0) {
		std::cout << "oops!" << std::endl;
	}

	seg2.index.clear();
	seg2.start = -1; seg2.end = -1;
	segments[j] = segments[Nseg - 1]; Nseg--;
}

void Filament2D::fuse_segments() {
	bool sown = false; int pass = 0;
	while (!sown) {
		sown = true;
		std::cout << "passes: " << pass << std::endl; pass++;
		for (int i = 0; i < Nseg - 1; i++) {
			if (segments[i].index.size() > 0) {
				auto& seg1 = segments[i];
				for (int j = i + 1; j < Nseg; j++) {
					if ((segments[j].index.size() > 0) && (i != j)) {
						auto& seg2 = segments[j];

						int start1 = seg1.start;
						int start2 = seg2.start;
						int end1 = seg1.end;
						int end2 = seg2.end;

						int d2ss = dist2(start1, start2);
						int d2se = dist2(start1, end2);
						int d2es = dist2(end1, start2);
						int d2ee = dist2(end1, end2);

						//std::cout << i << " (" << seg1.index.size() << ") " << j << " (" << seg2.index.size() << ") " << Nseg << " " << d2ss << " " << d2se << " " << d2es << " " << d2ee << std::endl;

						if (d2ss <= threshold) { /* segments touch start to start, but the smaller is flipped - flip it back, attach! */
							if (seg1.index.size() >= seg2.index.size()) sow_ss(i, j);
							else sow_ss(j, i);
							sown = false;
							break;
						}
						else if (d2ee <= threshold) { /* segments touch end to end, but the smaller is flipped - flip it back, attach! */
							if (seg1.index.size() >= seg2.index.size()) sow_ee(i, j);
							else sow_ee(j, i);
							sown = false;
							break;
						}
						else if (d2se <= threshold) {
							if (seg1.index.size() >= seg2.index.size()) sow_se(i, j);
							else sow_es(j, i);
							sown = false;
							break;
						}
						else if (d2es <= threshold) {
							if (seg1.index.size() >= seg2.index.size()) sow_es(i, j);
							else sow_se(j, i);
							sown = false;
							break;
						}
					}
				}
			}
		}
	}
	start = start;
	end = end;
}

void Filament3D::fuse_segments() {
	bool sown = false; int pass = 0;
	while (threshold < 10) {
		sown = false;
		while (!sown) {
			sown = true;
			//std::cout << "passes: " << pass << std::endl; pass++;
			for (int i = 0; i < Nseg - 1; i++) {
				if (segments[i].index.size() > 0) {
					auto& seg1 = segments[i];
					for (int j = i + 1; j < Nseg; j++) {
						if ((segments[j].index.size() > 0) && (i != j)) {
							auto& seg2 = segments[j];

							int start1 = seg1.start;
							int start2 = seg2.start;
							int end1 = seg1.end;
							int end2 = seg2.end;

							int d2ss = dist2(start1, start2);
							int d2se = dist2(start1, end2);
							int d2es = dist2(end1, start2);
							int d2ee = dist2(end1, end2);

							//std::cout << i << " (" << seg1.index.size() << ") " << j << " (" << seg2.index.size() << ") " << Nseg << " " << d2ss << " " << d2se << " " << d2es << " " << d2ee << std::endl;

							if (d2ss <= threshold) { /* segments touch start to start, but the smaller is flipped - flip it back, attach! */
								if (seg1.index.size() >= seg2.index.size()) sow_ss(i, j);
								else sow_ss(j, i);
								sown = false;
								break;
							}
							else if (d2ee <= threshold) { /* segments touch end to end, but the smaller is flipped - flip it back, attach! */
								if (seg1.index.size() >= seg2.index.size()) sow_ee(i, j);
								else sow_ee(j, i);
								sown = false;
								break;
							}
							else if (d2se <= threshold) {
								if (seg1.index.size() >= seg2.index.size()) sow_se(i, j);
								else sow_es(j, i);
								sown = false;
								break;
							}
							else if (d2es <= threshold) {
								if (seg1.index.size() >= seg2.index.size()) sow_es(i, j);
								else sow_se(j, i);
								sown = false;
								break;
							}
						}
					}
				}
			}
		}
		threshold++;
	}
	start = start;
	end = end;
}

void Filament2D::connect() {
	flip();
	//find_connections();
	//continuous_segments();
	find_segments();
	fuse_segments();
	std::cout << Nseg << std::endl;
	for (int i = 0; i < Nseg; i++) {
		std::cout << segments[i].index.size() << " ";
	}
	std::cout << std::endl;
}

void Filament3D::connect() {
	flip();
	//find_connections();
	//continuous_segments();
	find_segments();
	fuse_segments();
	//std::cout << Nseg << std::endl;
	for (int i = 0; i < Nseg; i++) {
		//std::cout << segments[i].index.size() << " ";
	}
	//std::cout << std::endl;
}

std::vector<std::vector<Point3D>> load_filaments(int period, std::string filename) {
	auto filaments = load3D(period, filename);
	std::vector<std::vector<Point3D>> res;
	for (int i = 0; i < period; i++) {
		for (int j = 0; j < period; j++) {
			for (int k = 0; k < period; k++) {
				int fil_index = filaments[i][j][k];
				if (fil_index > 0) {
					while (res.size() < fil_index) res.push_back({});
					Point3D p; p.x = i; p.y = j; p.z = k;
					res[fil_index - 1].push_back(p);
				}
			}
			delete filaments[i][j];
		}
		delete filaments[i];
	}
	delete[] filaments;
	return res;
}