#include "fitting.h"

vec3D glob_dir = vec3D(0.70465703391667, 0.6097087120454031, 0.3629238914809779);

std::random_device rd;
std::mt19937 gen(rd());

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_rms = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	double res = 0;
	for (auto& point : data) {
		double x = point.first;
		double y = point.second;
		double y_theory = model(x, param);
		res += (y - y_theory) * (y - y_theory);
	}
	res = std::sqrt(res / data.size());
	return res;
	};

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_rms_perc = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	return 100 * dist_rms(data, model, param);
	};

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_abs = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	double res = 0;
	for (auto& point : data) {
		double x = point.first;
		double y = point.second;
		double y_theory = model(x, param);
		res += std::abs(y - y_theory);
	}
	res /= data.size();
	return res;
	};

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	double res = 0;
	for (auto& point : data) {
		double x = point.first;
		double y = point.second;
		double y_theory = model(x, param);
		res = std::max(std::abs(y - y_theory), res);
	}
	return res;
	};

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max_rel = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	double res = 0;
	for (auto& point : data) {
		double x = point.first;
		double y = point.second;
		double y_theory = model(x, param);
		res = std::max(std::abs(y - y_theory) / y, res);
	}
	return res;
	};

std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)> dist_max_rel_perc = [](const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, const std::vector<double>& param) {
	return 100 * dist_max_rel(data, model, param);
	};

double line_dist2(const vec3D& P1, const vec3D& P2, const Point3D& p) {
	double denom = dist2(P1, P2); // find distance from p to the line defined by P1->P2, P1=/=P2
	double t = -(P1 - p) * (P2 - P1) / denom; // find the value of parameter along the line which corresponds to the closest point
	if (t < 0) return dist2(p, P1); // p hangs out past point P1, so P1 is closest
	else if (t > 1) return dist2(p, P2); // p hangs out past point P2, so P2 is closest
	else { // p is closer to some point between P1 and P2, we find it
		double x = P1.x + (P2.x - P1.x) * t; // perpendicular projection from p onto the line P1->P2, this is the closest point
		double y = P1.y + (P2.y - P1.y) * t;
		double z = P1.z + (P2.z - P1.z) * t;
		return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z);
	}
}

double line_dist2_unconstrained(const vec3D& P1, const vec3D& P2, const Point3D& p) {
	double denom = dist2(P1, P2); // find distance from p to the line defined by P1->P2, P1=/=P2
	double t = -(P1 - p) * (P2 - P1) / denom; // find the value of parameter along the line which corresponds to the closest point
	double x = P1.x + (P2.x - P1.x) * t; // perpendicular projection from p onto the line P1->P2, this is the closest point
	double y = P1.y + (P2.y - P1.y) * t;
	double z = P1.z + (P2.z - P1.z) * t;
	return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z);
}

double multilinear3D_avg(const Filament3D& fil, int segment_index, std::vector<vec3D>& points) {
	// points.size() >= 2 !!!
	double dist = 0;
	int length = fil.segments[segment_index].index.size();
	for (auto& point_index : fil.segments[segment_index].index) {
		const Point3D& p = fil.p[point_index];
		double mindist2 = 1e10;
		for (int lineseg = 0; lineseg < points.size() - 1; lineseg++) {
			const vec3D& P1 = points[lineseg];
			const vec3D& P2 = points[lineseg + 1];
			double dist2 = line_dist2(P1, P2, p);
			if (dist2 < mindist2) mindist2 = dist2;
		}
		dist += std::sqrt(mindist2);
	}
	return dist / length;
}

double multilinear3D_max(const Filament3D& fil, int segment_index, std::vector<vec3D>& points) {
	// points.size() >= 2 !!!
	double dist = 0;
	for (auto& point_index : fil.segments[segment_index].index) {
		const Point3D& p = fil.p[point_index];
		double mindist2 = 1e10;
		for (int lineseg = 0; lineseg < points.size() - 1; lineseg++) {
			const vec3D& P1 = points[lineseg];
			const vec3D& P2 = points[lineseg + 1];
			double dist2 = line_dist2(P1, P2, p);
			if (dist2 < mindist2) mindist2 = dist2;
		}
		dist = std::max(dist, std::sqrt(mindist2));
	}
	return dist;
}

double multilinear3D_max(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2) {
	double dist = 0;
	int index = start_index;
	while (true) {
		const Point3D& p = fil.p[index];
		double dist2 = line_dist2(P1, P2, p);
		dist = std::max(dist, std::sqrt(dist2));
		if (index == end_index) break;
		index = p.next;
	}
	return dist;
}

double multilinear3D_max_unconstrained(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2) {
	double dist = 0;
	int index = start_index;
	while (true) {
		const Point3D& p = fil.p[index];
		double dist2 = line_dist2_unconstrained(P1, P2, p); // only the line given by P1, P2 matters
		dist = std::max(dist, std::sqrt(dist2));
		if (index == end_index) break;
		index = p.next;
	}
	return dist;
}

double multilinear3D_rms_unconstrained(const Filament3D& fil, int start_index, int end_index, const vec3D& P1, const vec3D& P2) {
	double dist2 = 0;
	int index = start_index;
	int length = 0;
	while (true) {
		const Point3D& p = fil.p[index];
		dist2 += line_dist2_unconstrained(P1, P2, p);
		length++;
		if (index == end_index) break;
		index = p.next;
	}
	return std::sqrt(dist2 / length);
}

std::map<double, int> arclen(const Filament3D& fil, int segment_index) {
	std::map<double, int> foo, res;
	int index = fil.segments[segment_index].start;
	foo.insert({ (double)0, index });
	double total_dist = 0;
	while (index >= 0) {
		auto p_this = fil.p[index];
		int index_next = p_this.next;
		if (index_next >= 0) {
			Point3D p_next = fil.p[index_next];
			int d2 = ::dist2(p_this, p_next);
			total_dist += std::sqrt((double)d2);
			index = index_next;
			foo.insert({ total_dist, index });
		}
		else break;
	}
	for (auto& f : foo) {
		res.insert({ f.first / total_dist, f.second });
	}
	return res;
}

void initial_multi(const Filament3D& fil, int segment_index, std::map<double, int> arc_length, std::vector<vec3D>& points, int n_seg) {
	points.clear();
	vec3D start = fil.p[fil.segments[segment_index].start];
	vec3D end = fil.p[fil.segments[segment_index].end];
	for (int i = 0; i < n_seg + 1; i++) {
		if (i == 0) {
			points.push_back(start);
		}
		else if (i == n_seg) {
			points.push_back(end);
		}
		else {
			double len_along = (double)i / ((double)n_seg);
			auto low = arc_length.lower_bound(len_along);
			vec3D p = fil.p[low->second];
			points.push_back(p);
		}
	}
}

void shuffle_point(vec3D& p, const double& sigma) {
	std::normal_distribution<double> dx(p.x, sigma);
	std::normal_distribution<double> dy(p.y, sigma);
	std::normal_distribution<double> dz(p.z, sigma);
	p.x = dx(gen);
	p.y = dy(gen);
	p.z = dz(gen);
}

bool box_condition(const vec3D& p, const double& box_size) {
	return (p.x >= 0) && (p.y >= 0) && (p.z >= 0) && (p.x <= box_size) && (p.y <= box_size) && (p.z <= box_size);
}

void shuffle_point(vec3D& p, const double& sigma, const double& box_size) {
	std::normal_distribution<double> dx(p.x, sigma);
	std::normal_distribution<double> dy(p.y, sigma);
	std::normal_distribution<double> dz(p.z, sigma);
	/*p.x = -1;
	p.y = -1;
	p.z = -1;
	while (!box_condition(p, box_size)) {
		p.x = dx(gen);
		p.y = dy(gen);
		p.z = dz(gen);
	}*/
	p.x = dx(gen);
	p.y = dy(gen);
	p.z = dz(gen);
}

double iterate_distance(const Filament3D& fil, int segment_index, std::vector<vec3D>& points, int n_seg) {
	double dist = multilinear3D_avg(fil, segment_index, points);
	double mindist = dist;
	double sigma = (double)fil.period / 10;
	int repeated = 0;
	while (sigma > 1e-10) {
		auto points_new = points;
		for (int i = 1; i < points_new.size() - 1; i++) {
			shuffle_point(points_new[i], sigma, fil.period - 1);
		}
		dist = multilinear3D_avg(fil, segment_index, points_new);
		if (dist < mindist) {
			mindist = dist;
			points = points_new;
			repeated = 0;
			sigma *= 10;
			if (sigma > fil.period / 2) sigma = fil.period / 2;
		}
		else repeated++;
		if (repeated > 100) {
			sigma *= 0.99;
		}
	}
	return mindist;
}

void save_points(std::vector<vec3D>& points, std::string filename) {
	std::ofstream out; out.open(filename, std::ios_base::out);
	for (int i = 0; i < points.size(); i++) {
		out << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
	}
	out.close();
}

void find_fit_global(const Filament3D& fil, int segment_index, std::string filename) {
	auto len = arclen(fil, segment_index);
	std::vector<vec3D> points, points_old;
	int n_seg = 2;
	while (true) {
		initial_multi(fil, segment_index, len, points, n_seg);

		save_points(points, filename + "_" + std::to_string(segment_index) + ".txt");

		double dist_init = multilinear3D_avg(fil, segment_index, points);
		double dist_min = iterate_distance(fil, segment_index, points, n_seg);

		save_points(points, filename + "_" + std::to_string(segment_index) + ".txt");

		std::cout << n_seg << " " << dist_init << " " << dist_min << " " << n_seg * dist_min << std::endl;
		if ((dist_init - dist_min) / dist_init < 1e-3) {
			points = points_old;
			break;
		}
		else points_old = points;
		n_seg++;
	}
	std::ofstream out; out.open(filename + "_" + std::to_string(segment_index) + ".txt", std::ios_base::out);
	for (int i = 0; i < points.size(); i++) {
		out << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
	}
	out.close();
}

void extend_line(vec3D& P1, vec3D& P2, const Point3D& p1, const Point3D& p2) {
	double denom = dist2(P1, P2);
	double t1 = -(P1 - p1) * (P2 - P1) / denom;
	double t2 = -(P1 - p2) * (P2 - P1) / denom;
	double x1 = P1.x + (P2.x - P1.x) * t1;
	double y1 = P1.y + (P2.y - P1.y) * t1;
	double z1 = P1.z + (P2.z - P1.z) * t1;
	double x2 = P1.x + (P2.x - P1.x) * t2;
	double y2 = P1.y + (P2.y - P1.y) * t2;
	double z2 = P1.z + (P2.z - P1.z) * t2;
	P1.x = x1; P1.y = y1; P1.z = z1;
	P2.x = x2; P2.y = y2; P2.z = z2;
}

vec3D find_avg(const Filament3D& fil, int start_index, int end_index) {
	vec3D Pmid;
	int index = start_index; int count = 0;
	while (true) {
		Pmid = Pmid + fil.p[index];
		count++;
		if (index == end_index) break;
		else index = fil.p[index].next;
	}
	Pmid = Pmid / count;
	return Pmid;
}

double iterate_distance_rms(const Filament3D& fil, int start_index, int end_index, vec3D& P1, vec3D& P2) {
	vec3D Pmid = find_avg(fil, start_index, end_index);
	P2 = fil.p[end_index];
	double dist = multilinear3D_rms_unconstrained(fil, start_index, end_index, Pmid, P2);
	double mindist = dist;
	double sigma = (double)fil.period / 10;
	int repeated = 0;
	while ((sigma > 1e-10) && (mindist > 1e-10)) {
		auto P_new = P2;
		shuffle_point(P_new, sigma, fil.period - 1);
		dist = multilinear3D_rms_unconstrained(fil, start_index, end_index, Pmid, P_new);
		if (dist < mindist) {
			mindist = dist;
			P2 = P_new;
			repeated = 0;
			sigma *= 2;
			if (sigma > fil.period / 2) sigma = fil.period / 2;
		}
		else repeated++;
		if (repeated > 100) {
			sigma *= 0.99;
		}
	}
	P1 = Pmid;
	extend_line(P1, P2, fil.p[start_index], fil.p[end_index]);
	return mindist;
}

double iterate_distance_max(const Filament3D& fil, int start_index, int end_index, vec3D& P1, vec3D& P2) {
	vec3D Pmid = find_avg(fil, start_index, end_index);
	P1 = Pmid;
	P2 = fil.p[end_index];
	double dist = multilinear3D_max_unconstrained(fil, start_index, end_index, P1, P2);
	double mindist = dist;
	double sigma = (double)fil.period / 10;
	int repeated = 0;
	while ((sigma > 1e-10) && (mindist > 1e-10)) {
		auto P1_new = P1;
		auto P2_new = P2;
		shuffle_point(P1_new, sigma, fil.period - 1);
		shuffle_point(P2_new, sigma, fil.period - 1);
		dist = multilinear3D_max_unconstrained(fil, start_index, end_index, P1_new, P2_new);
		if (dist < mindist) {
			mindist = dist;
			P1 = P1_new;
			P2 = P2_new;
			repeated = 0;
			sigma *= 2;
			if (sigma > fil.period / 2) sigma = fil.period / 2;
		}
		else repeated++;
		if (repeated > 100) {
			sigma *= 0.99;
		}
	}
	extend_line(P1, P2, fil.p[start_index], fil.p[end_index]);
	return mindist;
}

void find_fit_progressive(const Filament3D& fil, int segment_index, std::ofstream& out, double line_rms_threshold, double line_length_threshold) {
	//double dist_threshold = 0.9; // for maximum metric
	double dist_threshold = 0.65; // for rms metric
	double len_threshold = line_length_threshold * fil.period;
	std::vector<std::pair<vec3D, vec3D>> lines;
	if (fil.segments[segment_index].index.size() > 2) {
		int start = fil.segments[segment_index].start;
		int end = fil.p[start].next;
		vec3D P1 = fil.p[start];
		vec3D P2 = fil.p[end];
		int len = 2;
		while (true) {
			double dist = iterate_distance_rms(fil, start, end, P1, P2);
			if (fil.p[end].next < 0) {
				//if (len >= len_threshold) lines.push_back({ P1, P2 });
				if (::dist2(P1, P2) >= len_threshold * len_threshold) lines.push_back({ P1, P2 });
				break;
			}
			else if (dist > line_rms_threshold) {
				end = fil.p[end].prev; len--; // go back one point
				dist = iterate_distance_rms(fil, start, end, P1, P2);

				//if (len >= len_threshold) lines.push_back({ P1, P2 });
				if (::dist2(P1, P2) >= len_threshold * len_threshold) lines.push_back({ P1, P2 });
				//std::cout << lines.size() << " " << start << " " << end << " " << len << std::endl;
				start = fil.p[end].next;
				end = fil.p[start].next;
				len = 2;
			}
			else {
				end = fil.p[end].next;
				P2 = fil.p[end]; len++;
			}
		}
		if (lines.size() > 0) {
			for (auto& line : lines) {
				out << line.first.x << " " << line.first.y << " " << line.first.z << std::endl;
				out << line.second.x << " " << line.second.y << " " << line.second.z << std::endl;
			}
		}
	}
	/*else {
		int start = fil.segments[segment_index].start;
		int end = fil.segments[segment_index].end;
		out << fil.p[start].x << " " << fil.p[start].y << " " << fil.p[start].z << std::endl;
		out << fil.p[end].x << " " << fil.p[end].y << " " << fil.p[end].z << std::endl;
	}*/
}

std::vector<std::pair<vec3D, vec3D>> find_fit_progressive(const Filament3D& fil, int segment_index, double line_rms_threshold, double line_length_threshold) {
	//double dist_threshold = 0.9; // for maximum metric
	double dist_threshold = 0.65; // for rms metric
	double len_threshold = line_length_threshold * fil.period;
	std::vector<std::pair<vec3D, vec3D>> lines;
	if (fil.segments[segment_index].index.size() > 2) {
		int start = fil.segments[segment_index].start;
		int end = fil.p[start].next;
		vec3D P1 = fil.p[start];
		vec3D P2 = fil.p[end];
		int len = 2;
		while (true) {
			double dist = iterate_distance_rms(fil, start, end, P1, P2);
			if (fil.p[end].next < 0) {
				//if (len >= len_threshold) lines.push_back({ P1, P2 });
				if (::dist2(P1, P2) >= len_threshold * len_threshold) lines.push_back({ P1, P2 });
				break;
			}
			else if (dist > line_rms_threshold) {
				end = fil.p[end].prev; len--; // go back one point
				dist = iterate_distance_rms(fil, start, end, P1, P2);

				//if (len >= len_threshold) lines.push_back({ P1, P2 });
				if (::dist2(P1, P2) >= len_threshold * len_threshold) lines.push_back({ P1, P2 });
				//std::cout << lines.size() << " " << start << " " << end << " " << len << std::endl;
				start = fil.p[end].next;
				end = fil.p[start].next;
				len = 2;
			}
			else {
				end = fil.p[end].next;
				P2 = fil.p[end]; len++;
			}
		}
	}
	return lines;
}

void sanitize_point(vec3D& p, int Nx, int Ny, int Nz) {
	if (p.x < 0) p.x += Nx; if (p.x >= Nx) p.x -= Nx;
	if (p.y < 0) p.y += Ny; if (p.y >= Ny) p.y -= Ny;
	if (p.z < 0) p.z += Nz; if (p.z >= Nz) p.z -= Nz;
}

void sanitize_point(vec3D& p, int period) {
	sanitize_point(p, period, period, period);
}

int get_nseg(const vec3D& P1, const vec3D& P2, double dL) {
	return std::ceil(std::sqrt(::dist2(P1, P2)) / dL);
}

std::vector<double> spine(const vec3D& P1, const vec3D& P2, const std::vector<std::vector<std::vector<double>>>& field, int period, double dL) {
	int N = get_nseg(P1, P2, dL);
	std::vector<double> res;
	for (int i = 0; i <= N; i++) {
		double u = (double)i / ((double)N);
		vec3D P = P1 * (1 - u) + P2 * u;
		sanitize_point(P, period);
		res.push_back(get_value(field, P));
	}
	return res;
}

std::vector<double> spine(const vec3D& P1, const vec3D& P2, double*** field, int Nx, int Ny, int Nz, double dL) {
	int N = get_nseg(P1, P2, dL);
	std::vector<double> res;
	for (int i = 0; i <= N; i++) {
		double u = (double)i / ((double)N);
		vec3D P = P1 * (1 - u) + P2 * u;
		sanitize_point(P, Nx, Ny, Nz);
		res.push_back(get_value(field, Nx, Ny, Nz, P));
	}
	return res;
}

std::vector<double> spine(const vec3D& P1, const vec3D& P2, double* field, int Nx, int Ny, int Nz, double dL) {
	int N = get_nseg(P1, P2, dL);
	std::vector<double> res;
	for (int i = 0; i <= N; i++) {
		double u = (double)i / ((double)N);
		vec3D P = P1 * (1 - u) + P2 * u;
		sanitize_point(P, Nx, Ny, Nz);
		res.push_back(get_value(field, Nx, Ny, Nz, P));
	}
	return res;
}

std::vector<double> spine(const vec3D& P1, const vec3D& P2, double*** field, int N, double dL) {
	return spine(P1, P2, field, N, N, N, dL);
}

std::vector<double> spine(const vec3D& P1, const vec3D& P2, double* field, int N, double dL) {
	return spine(P1, P2, field, N, N, N, dL);
}

vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, double r, double phi) {
	vec3D P = P1 * (1 - u) + P2 * u;
	vec3D n = (P2 - P1) / std::sqrt(::dist2(P1, P2));
	vec3D perp1 = cross(n, glob_dir); normalize(perp1); // one direction perpendicular to n
	vec3D perp2 = cross(perp1, n); normalize(perp2); // another direction perpendicular to both n and the first one
	vec3D res = P + r * (perp1 * std::cos(phi) + perp2 * std::sin(phi));
	return res;
}

vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, double r, double phi, int period) {
	vec3D P = P1 * (1 - u) + P2 * u;
	vec3D n = (P2 - P1) / std::sqrt(::dist2(P1, P2));
	vec3D perp1 = cross(n, glob_dir); normalize(perp1); // one direction perpendicular to n
	vec3D perp2 = cross(perp1, n); normalize(perp2); // another direction perpendicular to both n and the first one
	vec3D res = P + r * (perp1 * std::cos(phi) + perp2 * std::sin(phi));
	sanitize_point(res, period);
	return res;
}

vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u) {
	return P1 * (1 - u) + P2 * u;
}

vec3D cylinder_coord(const vec3D& P1, const vec3D& P2, double u, int period) {
	vec3D res = P1 * (1 - u) + P2 * u; sanitize_point(res, period);
	return res;
}

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double normalizing_factor, double*** field, int period) {
	double val = 0;
	for (int k = 0; k < Nphi; k++) {
		double phi = 2 * M_PI * k / ((double)Nphi);
		vec3D P = cylinder_coord(P1, P2, u, r, phi, period);
		val += get_value(field, period, P) / normalizing_factor;
	}
	return val;
}

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double normalizing_factor, double* field, int period) {
	double val = 0;
	for (int k = 0; k < Nphi; k++) {
		double phi = 2 * M_PI * k / ((double)Nphi);
		vec3D P = cylinder_coord(P1, P2, u, r, phi, period);
		val += get_value(field, period, P) / normalizing_factor;
	}
	return val;
}

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double*** field, int period) {
	return angular_average(P1, P2, u, r, Nphi, 1, field, period);
}

double angular_average(const vec3D& P1, const vec3D& P2, double u, double r, int Nphi, double* field, int period) {
	return angular_average(P1, P2, u, r, Nphi, 1, field, period);
}

std::vector<std::vector<double>> cylinder_accumulant(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period) {
	if (dL <= 0) dL = std::sqrt(::dist2(P1, P2)) / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (Rmax <= 0) Rmax = -Rmax * std::sqrt(::dist2(P1, P2));
	if (dr <= 0) dr = Rmax / Nu;
	int Nr = std::ceil(Rmax / dr);
	int Nphi = 10;
	std::vector<std::vector<double>> res;
	for (int i = 0; i <= Nu; i++) {
		double u = (double)i / ((double)Nu);
		std::vector<double> rvalues;
		vec3D P0 = cylinder_coord(P1, P2, u, period);
		double central_value = get_value(field, period, P0);
		rvalues.push_back(central_value);
		double value = 0;
		for (int j = 1; j <= Nr; j++) {
			double r = Rmax * j / ((double)Nr);
			double val = 0;
			for (int k = 0; k < Nphi; k++) {
				double phi = 2 * M_PI * k / ((double)Nphi);
				vec3D P = cylinder_coord(P1, P2, u, r, phi, period);
				val += get_value(field, period, P);
			}
			val /= Nphi;
			value += val;
			rvalues.push_back(value);
		}
		res.push_back(rvalues);
	}
	return res;
}

std::vector<std::vector<double>> cylinder_accumulant_normalized(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period) {
	if (dL <= 0) dL = std::sqrt(::dist2(P1, P2)) / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (Rmax <= 0) Rmax = -Rmax * std::sqrt(::dist2(P1, P2));
	if (dr <= 0) dr = Rmax / Nu;
	int Nr = std::ceil(Rmax / dr);
	int Nphi = 10;
	std::vector<std::vector<double>> res;
	for (int i = 0; i <= Nu; i++) {
		double u = (double)i / ((double)Nu);
		std::vector<double> rvalues;
		vec3D P0 = cylinder_coord(P1, P2, u, period);
		double central_value = get_value(field, period, P0);
		rvalues.push_back(1);
		double value = 0;
		for (int j = 1; j <= Nr; j++) {
			double r = Rmax * j / ((double)Nr);
			double val = angular_average(P1, P2, u, r, Nphi, field, period) / (central_value * Nphi);
			value += val;
			rvalues.push_back(value);
		}
		res.push_back(rvalues);
	}
	return res;
}

std::vector<std::pair<double, double>> cylinder_accumulant_avg(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period) {
	if (dL <= 0) dL = std::sqrt(::dist2(P1, P2)) / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (Rmax <= 0) Rmax = -Rmax * std::sqrt(::dist2(P1, P2));
	if (dr <= 0) dr = Rmax / Nu;
	int Nr = std::ceil(Rmax / dr);
	int Nphi = 10;
	std::vector<std::pair<double, double>> res;

	res.push_back({ 0,0 });
	double val, value;
	value = 0;
	for (int j = 1; j <= Nr; j++) {
		double r = Rmax * j / ((double)Nr);
		val = 0;
		for (int i = 0; i <= Nu; i++) {
			double u = (double)i / ((double)Nu);
			val += angular_average(P1, P2, u, r, Nphi, field, period);
		}
		val *= 2 * M_PI * r / ((Nu + 1) * Nphi);
		value += val;
		res.push_back({ r,value });
	}
	return res;
}

std::vector<std::pair<double, double>> cylinder_accumulant_avg_norm(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period) {
	if (dL <= 0) dL = std::sqrt(::dist2(P1, P2)) / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (Rmax <= 0) Rmax = -Rmax * std::sqrt(::dist2(P1, P2));
	if (dr <= 0) dr = Rmax / Nu;
	int Nr = std::ceil(Rmax / dr);
	int Nphi = 10;
	std::vector<std::pair<double, double>> res;

	std::vector<double> central_values;
	for (int i = 0; i <= Nu; i++) {
		double u = (double)i / ((double)Nu);
		vec3D P = cylinder_coord(P1, P2, u, period);
		central_values.push_back(get_value(field, period, P));
	}

	res.push_back({ 0,0 });
	double val, value;
	value = 0;
	for (int j = 1; j <= Nr; j++) {
		double r = Rmax * j / ((double)Nr);
		double jacobian = r;
		val = 0;
		for (int i = 0; i <= Nu; i++) {
			double u = (double)i / ((double)Nu);
			val += angular_average(P1, P2, u, r, Nphi, central_values[i], field, period);
		}
		val *= jacobian / ((Nu + 1) * Nphi); // the proper integration including jacobian
		value += val;
		res.push_back({ r, value / jacobian }); // to normalize we divide back by the jacobian
	}
	return res;
}

std::vector<std::pair<double, double>> cylinder_radial_avg_norm(const vec3D& P1, const vec3D& P2, double dL, double dr, double Rmax, double*** field, int period) {
	if (dL <= 0) dL = std::sqrt(::dist2(P1, P2)) / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (Rmax <= 0) Rmax = -Rmax * std::sqrt(::dist2(P1, P2));
	if (dr <= 0) dr = Rmax / Nu;
	int Nr = std::ceil(Rmax / dr);
	int Nphi = 10;
	std::vector<std::pair<double, double>> res;

	auto central_values = spine(P1, P2, field, period, dL);

	double value, val;
	res.push_back({ 0,1 });
	double jacobian = 2 * M_PI;
	for (int j = 1; j <= Nr; j++) {
		double r = Rmax * j / ((double)Nr);
		val = 0;
		for (int i = 0; i <= Nu; i++) {
			double u = (double)i / ((double)Nu);
			val += angular_average(P1, P2, u, r, Nphi, central_values[i], field, period);
		}
		val /= (Nu + 1) * Nphi;
		res.push_back({ r, val });
	}
	return res;
}

std::vector<std::pair<double, double>> cylinder_radial_avg_norm_constrained(const vec3D& P1, const vec3D& P2, double dL, double dr, double threshold, double*** field, int period) {
	double L = std::sqrt(::dist2(P1, P2));
	if (dL <= 0) dL = L / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (dr <= 0) dr = - dr * dL;
	int Nphi = 10;
	std::vector<std::pair<double, double>> res;

	auto central_values = spine(P1, P2, field, period, dL);

	double value, val;
	val = 1;
	res.push_back({ 0,val });
	double jacobian = 2 * M_PI;
	int j = 1;
	double val_min = 1;
	while ((val >= threshold) && (val <= 1.001) && (val <= val_min)) {
		if (val_min > val) val_min = val;
		double r = j * dr;
		if (r > L) {
			return res;
		}
		val = 0;
		for (int i = 0; i <= Nu; i++) {
			double u = (double)i / ((double)Nu);
			val += angular_average(P1, P2, u, r, Nphi, central_values[i], field, period);
		}
		val /= (Nu + 1) * Nphi;
		res.push_back({ r, val });
		j++;
	}
	return res;
}

std::vector<std::pair<double, double>> cylinder_radial_avg_norm_constrained(const vec3D& P1, const vec3D& P2, double dL, double dr, double threshold, double* field, int period) {
	double L = std::sqrt(::dist2(P1, P2));
	if (dL <= 0) dL = L / 10;
	int Nu = get_nseg(P1, P2, dL);
	if (dr <= 0) dr = - dr * dL;
	int Nphi = 10;
	std::vector<std::pair<double, double>> res;

	auto central_values = spine(P1, P2, field, period, dL);

	double value, val;
	val = 1;
	res.push_back({ 0,val });
	double jacobian = 2 * M_PI;
	int j = 1;
	double val_min = 1;
	while ((val >= threshold) && (val <= 1.001) && (val <= val_min)) {
		if (val_min > val) val_min = val;
		double r = j * dr;
		if (r > L) {
			return res;
		}
		val = 0;
		for (int i = 0; i <= Nu; i++) {
			double u = (double)i / ((double)Nu);
			val += angular_average(P1, P2, u, r, Nphi, central_values[i], field, period);
		}
		val /= (Nu + 1) * Nphi;
		res.push_back({ r, val });
		j++;
	}
	return res;
}

void generate_params(const std::vector<double>& param, std::vector<double>& new_param, const std::vector<double>& sigma, const std::vector<std::pair<double, double>>& constraint) {
	int Npar = param.size();
	for (int i = 0; i < Npar; i++) {
		double parmin = constraint[i].first;
		double parmax = constraint[i].second;
		std::normal_distribution<double> par_dist(param[i], sigma[i]);
		new_param[i] = parmin - 1;
		while ((new_param[i] <= parmin) || (new_param[i] >= parmax)) {
			new_param[i] = par_dist(gen);
		}
	}
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<double>& starting_values, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun) {
	int Npar = param.size();
	for (int i = 0; i < Npar; i++) {
		param[i] = starting_values[i];
	}
	double dist = dist_fun(data, model, param);
	double new_dist = 1;
	auto new_param = param;
	int repeated = 0;
	double max_sigma = 100;
	while (max_sigma > 1e-6) {
		generate_params(param, new_param, sigma, constraint);
		new_dist = dist_fun(data, model, new_param);
		if (new_dist < dist) {
			param = new_param;
			dist = new_dist;
			repeated = 0;
			max_sigma = 100;
			for (int i = 0; i < Npar; i++) {
				sigma[i] = param[i];
				max_sigma = std::min(sigma[i], max_sigma);
			}
		}
		else {
			repeated++;
		}
		if (repeated > 100) {
			max_sigma = 100;
			for (int i = 0; i < Npar; i++) {
				sigma[i] *= 0.99;
				max_sigma = std::min(sigma[i], max_sigma);
			}
		}
	}
	return dist;
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun) {
	auto starting_values = param;
	return fit(data, model, param, starting_values, sigma, constraint, dist_fun);
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun) {
	auto sigma = param;
	for (auto& s : sigma) {
		s = std::max(1.0, s);
	}
	return fit(data, model, param, sigma, constraint, dist_fun);
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun) {
	std::vector<std::pair<double, double>> constraint;
	for (int i = 0; i < param.size(); i++) {
		constraint.push_back({ -1e10,1e10 });
	}
	return fit(data, model, param, constraint, dist_fun);
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param) {
	auto dist_fun = dist_rms;
	return fit(data, model, param, dist_fun);
}

double fit(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint) {
	auto dist_fun = dist_rms;
	return fit(data, model, param, constraint, dist_fun);
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<double>& starting_values, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N) {
	int Npar = data.size();
	if (Npar <= 3) {
		for (int i = 0; i < param.size(); i++) param[i] = constraint[i].first - 1;
		return -1;
	}
	else {
		N = 3;
		double dist = 1;
		while (N <= Npar) {
			std::vector<std::pair<double, double>> data_chopped; data_chopped.resize(N);
			for (int i = 0; i < N; i++) {
				data_chopped[i] = data[i];
			}
			dist = fit(data_chopped, model, param, starting_values, sigma, constraint, dist_fun);
			if (dist > threshold) {
				N--; data_chopped.resize(N);
				for (int i = 0; i < N; i++) {
					data_chopped[i] = data[i];
				}
				dist = fit(data_chopped, model, param, starting_values, sigma, constraint, dist_fun);
				return dist;
			}
			else N++;
		}
		return dist;
	}
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, std::vector<double> sigma, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N) {
	auto starting_values = param;
	return fit_progressive(data, model, param, starting_values, sigma, constraint, dist_fun, threshold, N);
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N) {
	auto sigma = param;
	for (auto& s : sigma) {
		s = std::max(1.0, s);
	}
	return fit_progressive(data, model, param, sigma, constraint, dist_fun, threshold, N);
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::function<double(const std::vector<std::pair<double, double>>&, const std::function<double(double, std::vector<double>)>&, const std::vector<double>&)>& dist_fun, double threshold, int& N) {
	std::vector<std::pair<double, double>> constraint;
	for (int i = 0; i < param.size(); i++) {
		constraint.push_back({ -1e10,1e10 });
	}
	return fit_progressive(data, model, param, constraint, dist_fun, threshold, N);
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, double threshold, int& N) {
	auto dist_fun = dist_rms;
	return fit_progressive(data, model, param, dist_fun, threshold, N);
}

double fit_progressive(const std::vector<std::pair<double, double>>& data, const std::function<double(double, std::vector<double>)>& model, std::vector<double>& param, const std::vector<std::pair<double, double>>& constraint, double threshold, int& N) {
	auto dist_fun = dist_rms;
	return fit_progressive(data, model, param, constraint, threshold, N);
}

double mean_field(const vec3D& P1, const vec3D& P2, double*** field, double period, double dL) {
	double res = 0;
	int Nu = get_nseg(P1, P2, dL);
	for (int i = 0; i <= Nu; i++) {
		double u = i * dL / Nu;
		vec3D P = cylinder_coord(P1, P2, u, period);
		res += get_value(field, period, P);
	}
	res /= Nu + 1;
	return res;
}

double mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double*** field, double period, double dL, double dr, double Nphi) {
	double res = 0;
	int Nu = get_nseg(P1, P2, dL); dL = 1.0 / ((double)Nu);
	if (dr < 0) dr = -dr * dL;
	double dphi = 2 * M_PI / ((double)Nphi);
	double norm = 0; double dA; vec3D P;
	for (int i = 0; i < Nu; i++) { // this is the integral over a cylinder with weight given by the profile function
		double u = (i + 0.5) * dL; // midpoint for u variable
		int j = 1;
		while (true) { // radial loop starts at the first circle
			double r = (j + 0.5) * dr;
			for (int k = 0; k < Nphi; k++) {
				double phi = (k + 0.5) * dphi;
				dA = r * dL * dr * dphi;
				P = cylinder_coord(P1, P2, u, r, phi, period);
				res += get_value(field, period, P) * dA;
				norm += profile(r, param) * dA;
			}
		}
		dA = M_PI * dr * dr;
		P = cylinder_coord(P1, P2, u, period); // the central line has radius 0 and thus no angle needed
		res += get_value(field, period, P) * dA; norm += profile(0, param) * dA;
	}
	res /= norm;
	return res;
}

vec3D mean_field(const vec3D& P1, const vec3D& P2, double*** Fx, double*** Fy, double*** Fz, double period, double dL) {
	vec3D res = vec3D(0, 0, 0);
	int Nu = get_nseg(P1, P2, dL);
	for (int i = 0; i <= Nu; i++) {
		double u = i * dL / Nu;
		vec3D P = cylinder_coord(P1, P2, u, period);
		res += get_value(Fx, Fy, Fz, period, P);
	}
	res /= Nu + 1;
	return res;
}

vec3D mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double*** Fx, double*** Fy, double*** Fz, double period, double dL, double dr, double Nphi) {
	vec3D res = vec3D(0, 0, 0);
	int Nu = get_nseg(P1, P2, dL); dL = 1.0 / ((double)Nu);
	if (dr < 0) dr = -dr * dL;
	double dphi = 2 * M_PI / ((double)Nphi);
	double norm = 0; double dV, weight; vec3D P;
	for (int i = 0; i < Nu; i++) { // this is the integral over a cylinder with weight given by the profile function
		double u = (i + 0.5) * dL; // midpoint for u variable
		int j = 1;
		while (true) { // radial loop starts at the first circle
			double r = (j + 0.5) * dr;
			weight = profile(r, param);
			for (int k = 0; k < Nphi; k++) {
				double phi = (k + 0.5) * dphi;
				dV = r * dL * dr * dphi;
				P = cylinder_coord(P1, P2, u, r, phi, period);
				res += get_value(Fx, Fy, Fz, period, P) * weight * dV;
				norm += weight * dV;
			}
			if (weight < 1e-3) break;
			j++;
		}
		weight = profile(0, param);
		dV = M_PI * dr * dr * dL; // the area is that of a circle of radius dr
		P = cylinder_coord(P1, P2, u, period); // central line has radius 0 and thus no angle needed
		res += get_value(Fx, Fy, Fz, period, P) * weight * dV;
		norm += weight * dV;
	}
	res /= norm;
	return res;
}

double mean_field(const vec3D& P1, const vec3D& P2, double* field, double period, double dL) {
	double res = 0;
	int Nu = get_nseg(P1, P2, dL);
	for (int i = 0; i <= Nu; i++) {
		double u = i * dL / Nu;
		vec3D P = cylinder_coord(P1, P2, u, period);
		res += get_value(field, period, P);
	}
	res /= Nu + 1;
	return res;
}

double mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double* field, double period, double dL, double dr, double Nphi) {
	double res = 0;
	int Nu = get_nseg(P1, P2, dL); dL = 1.0 / ((double)Nu);
	if (dr < 0) dr = -dr * dL;
	double dphi = 2 * M_PI / ((double)Nphi);
	double norm = 0; double dA; vec3D P;
	for (int i = 0; i < Nu; i++) { // this is the integral over a cylinder with weight given by the profile function
		double u = (i + 0.5) * dL; // midpoint for u variable
		int j = 1;
		while (true) { // radial loop starts at the first circle
			double r = (j + 0.5) * dr;
			for (int k = 0; k < Nphi; k++) {
				double phi = (k + 0.5) * dphi;
				dA = r * dL * dr * dphi;
				P = cylinder_coord(P1, P2, u, r, phi, period);
				res += get_value(field, period, P) * dA;
				norm += profile(r, param) * dA;
			}
		}
		dA = M_PI * dr * dr;
		P = cylinder_coord(P1, P2, u, period); // the central line has radius 0 and thus no angle needed
		res += get_value(field, period, P) * dA; norm += profile(0, param) * dA;
	}
	res /= norm;
	return res;
}

vec3D mean_field(const vec3D& P1, const vec3D& P2, double* Fx, double* Fy, double* Fz, double period, double dL) {
	vec3D res = vec3D(0, 0, 0);
	int Nu = get_nseg(P1, P2, dL);
	for (int i = 0; i <= Nu; i++) {
		double u = i * dL / Nu;
		vec3D P = cylinder_coord(P1, P2, u, period);
		res += get_value(Fx, Fy, Fz, period, P);
	}
	res /= Nu + 1;
	return res;
}

vec3D mean_field(const vec3D& P1, const vec3D& P2, const std::function<double(double, std::vector<double>)>& profile, std::vector<double>& param, double* Fx, double* Fy, double* Fz, double period, double dL, double dr, double Nphi) {
	vec3D res = vec3D(0, 0, 0);
	int Nu = get_nseg(P1, P2, dL); dL = 1.0 / ((double)Nu);
	if (dr < 0) dr = -dr * dL;
	double dphi = 2 * M_PI / ((double)Nphi);
	double norm = 0; double dV, weight; vec3D P;
	for (int i = 0; i < Nu; i++) { // this is the integral over a cylinder with weight given by the profile function
		double u = (i + 0.5) * dL; // midpoint for u variable
		int j = 1;
		while (true) { // radial loop starts at the first circle
			double r = (j + 0.5) * dr;
			weight = profile(r, param);
			for (int k = 0; k < Nphi; k++) {
				double phi = (k + 0.5) * dphi;
				dV = r * dL * dr * dphi;
				P = cylinder_coord(P1, P2, u, r, phi, period);
				res += get_value(Fx, Fy, Fz, period, P) * weight * dV;
				norm += weight * dV;
			}
			if (weight < 1e-3) break;
			j++;
		}
		weight = profile(0, param);
		dV = M_PI * dr * dr * dL; // the area is that of a circle of radius dr
		P = cylinder_coord(P1, P2, u, period); // central line has radius 0 and thus no angle needed
		res += get_value(Fx, Fy, Fz, period, P) * weight * dV;
		norm += weight * dV;
	}
	res /= norm;
	return res;
}