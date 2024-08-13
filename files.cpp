#include "files.h"

void load1D(double* field, int count, std::string filename) {
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < count; i++) {
		file_in >> field[i];
	}
	file_in.close();
}

void load1D(std::vector<double>& field, std::string filename) {
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < field.size(); i++) {
		file_in >> field[i];
	}
	file_in.close();
}

void load2D(double** field, int countX, int countY, std::string filename) {
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		for (int j = 0; j < countY; j++) {
			file_in >> field[i][j];
		}
	}
	file_in.close();
}

void load2D(double** field, int count, std::string filename) {
	load2D(field, count, count, filename);
}

void load3D(double*** field, int countX, int countY, int countZ, std::string filename) {
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		for (int j = 0; j < countY; j++) {
			for (int k = 0; k < countZ; k++) {
				file_in >> field[i][j][k];
			}
		}
	}
	file_in.close();
}

void load3D(double*** field, int countX, int countY, int countZ, std::string filename, double (*transform)(double)) {
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		for (int j = 0; j < countY; j++) {
			for (int k = 0; k < countZ; k++) {
				file_in >> field[i][j][k];
				field[i][j][k] = transform(field[i][j][k]);
			}
		}
	}
	file_in.close();
}

void load3D(double*** field, int count, std::string filename) {
	load3D(field, count, count, count, filename);
}

void load3D(double*** field, int count, std::string filename, double (*transform)(double)) {
	load3D(field, count, count, count, filename, transform);
}

double* load1D(int count, std::string filename) {
	double* field = new double[count];
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < count; i++) {
		file_in >> field[i];
	}
	file_in.close();
	return field;
}

double** load2D(int countX, int countY, std::string filename) {
	double** field = new double* [countX];
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		field[i] = new double[countY];
		for (int j = 0; j < countY; j++) {
			file_in >> field[i][j];
		}
	}
	file_in.close();
	return field;
}

double** load2D(int count, std::string filename) {
	return load2D(count, count, filename);
}

double*** load3D(int countX, int countY, int countZ, std::string filename, double (*transform)(double)) {
	double*** field = new double** [countX];
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		field[i] = new double* [countY];
		for (int j = 0; j < countY; j++) {
			field[i][j] = new double[countZ];
			for (int k = 0; k < countZ; k++) {
				file_in >> field[i][j][k];
				field[i][j][k] = transform(field[i][j][k]);
			}
		}
	}
	file_in.close();
	return field;
}

double*** load3D(int countX, int countY, int countZ, std::string filename) {
	double*** field = new double** [countX];
	std::ifstream file_in; file_in.open(filename);
	for (int i = 0; i < countX; i++) {
		field[i] = new double* [countY];
		for (int j = 0; j < countY; j++) {
			field[i][j] = new double[countZ];
			for (int k = 0; k < countZ; k++) {
				file_in >> field[i][j][k];
				field[i][j][k] = field[i][j][k];
			}
		}
	}
	file_in.close();
	return field;
}

double*** load3D(int count, std::string filename, double (*transform)(double)) {
	return load3D(count, count, count, filename, transform);
}

double*** load3D(int count, std::string filename) {
	return load3D(count, count, count, filename);
}

void load_lines(std::vector<std::pair<vec3D, vec3D>>& lines, std::string filename) {
	std::ifstream file_in; file_in.open(filename);
	lines.clear();
	double x, y, z;
	while (file_in >> x >> y >> z) {
		vec3D first(x, y, z);
		file_in >> x >> y >> z;
		vec3D second(x, y, z);
		lines.push_back({ first, second });
	}
	file_in.close();
}

std::vector<std::pair<vec3D, vec3D>> load_lines(std::string filename) {
	std::vector<std::pair<vec3D, vec3D>> res;
	load_lines(res, filename);
	return res;
}

void save_spine(const std::vector<double>& spine, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	double dL = 1.0f / ((double)spine.size() - 1);
	for (int i = 0; i < spine.size(); i++) {
		file_out << i * dL << " " << spine[i] << std::endl;
	}
	file_out.close();
}

void save_num(double x, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x;
	file_out.close();
}

void save_num(double x, double y, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x << " " << y;
	file_out.close();
}

void save_num(double x, double y, double z, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x << " " << y << " " << z;
	file_out.close();
}

void save_num(double x1, double x2, double x3, double x4, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x1 << " " << x2 << " " << x3 << " " << x4;
	file_out.close();
}

void save_num(double x1, double x2, double x3, double x4, double x5, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5;
	file_out.close();
}

void save_num(double x1, double x2, double x3, double x4, double x5, double x6, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	file_out << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6;
	file_out.close();
}

void save_num(const std::vector<double>& x, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	for (auto& x : x) file_out << x << " ";
	file_out.close();
}

void save_num(const std::vector<double>& x, std::ofstream& file_out) {
	for (auto& x : x) file_out << x << " "; file_out << std::endl;
}

void save1D(const std::vector<double>& field, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	for (int i = 0; i < field.size(); i++) {
		file_out << field[i] << std::endl;
	}
	file_out.close();
}

void save2D(const std::vector<std::vector<double>>& field, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	int Nx = field.size(); int Ny = field[0].size(); // assume square
	if (Ny <= 1e6) {
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				file_out << field[i][j] << " ";
			}
			file_out << std::endl;
		}
	}
	else {
		file_out << Nx << " " << Ny << std::endl;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				file_out << field[i][j] << std::endl;
			}
		}
	}
	file_out.close();
}

void save2D(const std::vector<std::pair<double, double>>& field, std::string filename) {
	std::ofstream file_out; file_out.open(filename);
	int N = field.size();
	for (int i = 0; i < N; i++) {
		file_out << field[i].first << " " << field[i].second << std::endl;
	}
	file_out.close();
}

int getSize(std::string filename, std::string groupname) {
	std::string fieldname = "Density";
	try {
		H5File file(filename.c_str(), H5F_ACC_RDONLY);
		//std::cout << "foo1" << std::endl;
		Group group = file.openGroup(groupname.c_str());
		//std::cout << "foo2" << std::endl;
		DataSet dataset = group.openDataSet(fieldname);
		//std::cout << "foo3" << std::endl;
		DataSpace dataspace = dataset.getSpace();
		//std::cout << "foo4" << std::endl;
		int rank = dataspace.getSimpleExtentNdims();
		//std::cout << "foo5 " << rank << std::endl;
		hsize_t dims[rank];
		//std::cout << "foo6" << std::endl;
		int ndims = dataspace.getSimpleExtentDims(dims, NULL);
		//std::cout << "foo7 " << ndims << " " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
		//hsize_t offset[rank];
		//for (auto & o: offset) o = 0;
		//std::cout << "foo8 " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
		//dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
		//std::cout << "foo9" << std::endl;
		//DataSpace memspace(rank, dims);
		int res_size = 1;
		for (int i = 0; i < rank; i++) {
			res_size *= dims[i];
		}
		return res_size;
	}
	
	// catch failure caused by the H5File operations
	catch( FileIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSet operations
	catch( DataSetIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ) {
		error.printErrorStack();
		return -1;
	}
	return 0;
}

int read_from_hdf5(std::string filename, std::string groupname, std::string fieldname, long double res[], int size) {
	if ((fieldname == "density") || (fieldname == "rho")) fieldname = "Density";
	else if (fieldname == "vx") fieldname = "x-velocity";
	else if (fieldname == "vy") fieldname = "y-velocity";
	else if (fieldname == "vz") fieldname = "z-velocity";
	//std::cout << "poop " << filename << std::endl;
	//std::cout << "poop " << groupname << std::endl;
	//std::cout << "poop " << fieldname << std::endl;
	try {
		if ((fieldname == "absv") || (fieldname == "velocity-magnitude") || (fieldname == "v")) {
			long double vx[size];
			long double vy[size];
			long double vz[size];
			read_from_hdf5(filename, groupname, "x-velocity", vx, size);
			read_from_hdf5(filename, groupname, "y-velocity", vy, size);
			read_from_hdf5(filename, groupname, "z-velocity", vz, size);
			for (int i = 0; i < size; i++) res[i] = std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
		}
		else {
			H5File file(filename.c_str(), H5F_ACC_RDONLY);
			//std::cout << "foo1" << std::endl;
			Group group = file.openGroup(groupname.c_str());
			//std::cout << "foo2" << std::endl;
			DataSet dataset = group.openDataSet(fieldname);
			//std::cout << "foo3" << std::endl;
			DataSpace dataspace = dataset.getSpace();
			//std::cout << "foo4" << std::endl;
			int rank = dataspace.getSimpleExtentNdims();
			//std::cout << "foo5 " << rank << std::endl;
			hsize_t dims[rank];
			//std::cout << "foo6" << std::endl;
			int ndims = dataspace.getSimpleExtentDims(dims, NULL);
			//std::cout << "foo7 " << ndims << " " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
			hsize_t offset[rank];
			for (auto & o: offset) o = 0;
			//std::cout << "foo8 " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
			dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
			//std::cout << "foo9" << std::endl;
			DataSpace memspace(rank, dims);
			dataset.read(res, PredType::NATIVE_LDOUBLE, memspace, dataspace);
		}
	}
	
	// catch failure caused by the H5File operations
	catch( FileIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSet operations
	catch( DataSetIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ) {
		error.printErrorStack();
		return -1;
	}
	return 0;
}
