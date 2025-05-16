#include "hist.h"

std::string pad(int number, int places, char fill) {
    std::ostringstream oss;
    oss << std::setw(places) << std::setfill(fill) << number;
    return oss.str();
}

std::string pad(int number) {
    return pad(number, 4, '0');
}

bool is_int(const std::string& s) {
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

bool is_float(const std::string& string) {
    std::string::const_iterator it = string.begin();
    bool decimalPoint = false;
    int minSize = 0;
    if (string.size() > 0 && (string[0] == '-' || string[0] == '+')) {
        it++;
        minSize++;
    }
    while(it != string.end()) {
        if(*it == '.') {
            if(!decimalPoint) decimalPoint = true;
            else break;
        }
            else if(!std::isdigit(*it) && ((*it!='f') || it + 1 != string.end() || !decimalPoint)){
            break;
        }
        ++it;
    }
    return string.size()>minSize && it == string.end();
}

bool is_path_valid(const std::string& path) {
    fs::path p(path);
    return fs::exists(p);
}

std::vector<int> get_frames(const std::string& path) {
    std::regex pattern(R"(^DD(\d{4})$)");

    std::vector<int> frames;

    for (const auto& entry : fs::directory_iterator(path)) {
        if (entry.is_directory()) {  // Ensure it's a directory
            std::string dirname = entry.path().filename().string();
            std::smatch match;
            if (std::regex_match(dirname, match, pattern)) {
                int number = std::stoi(match[1].str());  // Extract and convert xxxx to int
                frames.push_back(number);
            }
        }
    }

    std::sort(frames.begin(), frames.end());

    return frames;
}

bool validate_input(int argc, char *argv[], std::string& sim_path, std::string& path_out, std::vector<int>& frames) {
    int frame_start = -1;
    int frame_end = -1;
    frames.resize(0);
    if (argc < 3) {
        std::cout << "Give me a path to filaments and output path!" << std::endl;
        return false;
    }
    else {
        sim_path = argv[1];
        path_out = argv[2];
        while (!sim_path.empty() && sim_path.back() == '/') {
            sim_path.pop_back();
        }
        sim_path.push_back('/');
        while (!path_out.empty() && path_out.back() == '/') {
            path_out.pop_back();
        }
        path_out.push_back('/');
        if (!is_path_valid(sim_path)) {
            std::cout << sim_path << " is not a valid path, aborting..." << std::endl;
            return false;
        }
        fs::create_directory(path_out);
        if (argc == 4) {
            if (is_int(argv[3])) {
                frame_start = std::stoi(argv[3]);
                frame_end = frame_start;
            }
            else {
                std::cout << "Invalid parameter " << argv[3] << ", aborting..." << std::endl;
                return false;
            }
        }
        if (argc == 5) {
            if (is_int(argv[3]) && is_int(argv[4])) {
                frame_start = std::stoi(argv[3]);
                frame_end = std::stoi(argv[4]);
            }
            else {
                std::cout << "Invalid parameters " << argv[3] << " " << argv[4] << ", aborting..." << std::endl;
                return false;
            }
        }
        if (frame_start > frame_end) std::swap(frame_start, frame_end);
        std::vector<int> available_frames = get_frames(sim_path);
        for (auto & f : available_frames) {
            if ((f >= frame_start) && ((f <= frame_end) || (frame_end < 0))) {
                frames.push_back(f);
            }
        }
        if (frames.empty()) {
            std::cout << "No available frames found in range [" << frame_start << ", " << frame_end << "], aborting..." << std::endl;
            return false;
        }
        else {
            for (auto & f: frames) {
                fs::create_directory(path_out + "DD" + pad(f, 4, '0'));
            }
        }
        return true;
    }
}

std::pair<double, double> find_bounds(const std::vector<double>& data) {
    double xmin = 1e10;
    double xmax = -1e10;
    for (const double& d: data) {
        if (xmin > d) xmin = d;
        if (xmax < d) xmax = d;
    }
    return std::make_pair(xmin, xmax);
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbins, std::pair<double, double> bounds) {
    std::pair<std::vector<double>, std::vector<double>> res;
    if (nbins > 1) {
        int n = data.size();
        double npt = (double)n / nbins;
        int remainder = n - static_cast<int>(npt) * nbins;
        std::multiset<double> sorted;
        double xmin = bounds.first;
        double xmax = bounds.second;

        // Sort input data within bounds
        for (const double& d : data) {
            if ((d >= xmin) && (d <= xmax)) sorted.insert(d);
        }

        std::vector<double> bin_dividers, weights;
        bin_dividers.push_back(xmin);
        int count = 0;
        int bin_index = 0;
        std::multiset<double>::iterator it;

        for (it = sorted.begin(); it != sorted.end(); ++it) {
            if (count < npt) {
                count++;
            } else {
                double val_old = *std::prev(it);
                double val = *it;
                double new_bin = 0.5 * (val_old + val);

                if (new_bin > bin_dividers.back()) {
                    bin_dividers.push_back(new_bin);
                    weights.push_back(count);
                    bin_index++;
                    count = 1;
                }
            }
        }
        if (!bin_dividers.empty() && xmax > bin_dividers.back()) {
            bin_dividers.push_back(xmax);
            weights.push_back(count + remainder);
        } else if (!weights.empty()) {
            weights.back() += remainder;
        }

        // Normalize weights to sum to 1
        for (auto& w : weights) w /= n;

        res.first = bin_dividers;
        res.second = weights;
        return res;
    } else {
        std::cout << "Number of bins must be > 1!" << std::endl;
        return res;
    }
}


std::vector<std::vector<double>> make_joint_histogram(
    const std::vector<double>& data_x,
    const std::vector<double>& data_y,
    const std::vector<double>& x_bins,
    const std::vector<double>& y_bins
) {
    int nbins_x = x_bins.size() - 1;
    int nbins_y = y_bins.size() - 1;

    // Initialize 2D histogram with zeros
    std::vector<std::vector<double>> hist(nbins_x, std::vector<double>(nbins_y, 0.0));

    int total_count = 0;

    // Loop over data points and distribute into bins
    for (size_t i = 0; i < data_x.size(); i++) {
        double x = data_x[i];
        double y = data_y[i];

        // Find the bin index for x
        int bin_x = -1;
        for (int j = 0; j < nbins_x; j++) {
            if (x >= x_bins[j] && x < x_bins[j + 1]) {
                bin_x = j;
                break;
            }
        }

        // Find the bin index for y
        int bin_y = -1;
        for (int j = 0; j < nbins_y; j++) {
            if (y >= y_bins[j] && y < y_bins[j + 1]) {
                bin_y = j;
                break;
            }
        }

        // Skip out-of-bounds points
        if (bin_x >= 0 && bin_y >= 0) {
            hist[bin_x][bin_y] += 1.0;
            total_count++;
        }
    }

    // Normalize so that sum of all bins = 1
    if (total_count > 0) {
        for (int i = 0; i < nbins_x; i++) {
            for (int j = 0; j < nbins_y; j++) {
                hist[i][j] /= total_count;
            }
        }
    }

    return hist;
}

std::vector<std::vector<std::vector<double>>> make_joint_histogram(
    const std::vector<double>& data_x,
    const std::vector<double>& data_y,
    const std::vector<double>& data_z,
    const std::vector<double>& x_bins,
    const std::vector<double>& y_bins,
    const std::vector<double>& z_bins
) {
    int nbins_x = x_bins.size() - 1;
    int nbins_y = y_bins.size() - 1;
    int nbins_z = z_bins.size() - 1;

    // Initialize 2D histogram with zeros
    std::vector<std::vector<std::vector<double>>> hist(
        nbins_x, std::vector<std::vector<double>>(nbins_y, std::vector<double>(nbins_z, 0.0))
    );

    int total_count = 0;

    for (size_t i = 0; i < data_x.size(); i++) {
        double x = data_x[i];
        double y = data_y[i];
        double z = data_z[i];

        int bin_x = -1;
        for (int j = 0; j < nbins_x; j++) {
            if (x >= x_bins[j] && x < x_bins[j + 1]) {
                bin_x = j;
                break;
            }
        }

        int bin_y = -1;
        for (int j = 0; j < nbins_y; j++) {
            if (y >= y_bins[j] && y < y_bins[j + 1]) {
                bin_y = j;
                break;
            }
        }

        int bin_z = -1;
        for (int j = 0; j < nbins_z; j++) {
            if (z >= z_bins[j] && z < z_bins[j + 1]) {
                bin_z = j;
                break;
            }
        }

        if (bin_x >= 0 && bin_y >= 0 && bin_z >= 0) {
            hist[bin_x][bin_y][bin_z] += 1;
            total_count += 1;
        }
    }

    if (total_count > 0) {
        for (int i = 0; i < nbins_x; i++) {
            for (int j = 0; j < nbins_y; j++) {
                for (int k = 0; k < nbins_z; k++) {
                    hist[i][j][k] /= total_count;
                }
            }
        }
    }

    return hist;
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, int nbins) {
    std::pair<double, double> bounds = find_bounds(data);
    return make_histogram(data, nbins, bounds);
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins, std::pair<double, double> bounds) {
    std::pair<std::vector<double>, std::vector<double>> res;
    if (nbins > 1) {
        int n = data.size();
        double total_weight = 0;
        for (const double& d : data) {
            total_weight+= weight_fun(d);
        }
        //std::cout << total_weight << std::endl;
        double bin_weight = total_weight / nbins;
        std::multiset<double> sorted;
        double xmin = bounds.first;
        double xmax = bounds.second;
        for (const double& d: data) {
            if ((d >= xmin) && (d <= xmax)) sorted.insert(d);
        }
        std::vector<double> bin_dividers, weights;
        bin_dividers.push_back(xmin);
        double accu_weight = 0; int bin_index = 0;
        std::multiset<double>::iterator it;
        for (it = sorted.begin(); it != sorted.end(); ++it) {
            if (accu_weight < bin_weight) {
                accu_weight += weight_fun(*it);
            }
            else {
                double val_old = *std::prev(it);
                double val = *it;
                double new_bin = 0.5 * (val_old + val);

                // Ensure bin width is nonzero before adding
                if (new_bin > bin_dividers.back()) {
                    bin_dividers.push_back(new_bin);
                    weights.push_back(accu_weight);
                    bin_index++;
                    accu_weight = weight_fun(val);
                }
            }
        }
        if (accu_weight < 0.5 * bin_weight) {
            bin_dividers.back() = xmax; weights.back() += accu_weight;
        }
        else {
            bin_dividers.push_back(xmax); weights.push_back(accu_weight);
        }
        for (auto& w: weights) w /= total_weight;
        res.first = bin_dividers; res.second = weights;
        return res;
    }
    else {
        std::cout << "Number of bins must be > 1!" << std::endl;
        return res;
    }
}

std::vector<std::vector<double>> make_joint_histogram(
    const std::vector<double>& data_x,
    const std::vector<double>& data_y,
    const std::vector<double>& x_bins,
    const std::vector<double>& y_bins,
    std::function<double(double, double)> weight_fun
) {
    int nbins_x = x_bins.size() - 1;
    int nbins_y = y_bins.size() - 1;
    
    // Initialize 2D histogram with zeros
    std::vector<std::vector<double>> hist(nbins_x, std::vector<double>(nbins_y, 0.0));

    double total_weight = 0.0;

    // Loop over data points and distribute into bins
    for (size_t i = 0; i < data_x.size(); i++) {
        double x = data_x[i];
        double y = data_y[i];
        double w = weight_fun(x, y);
        
        // Find the bin index for x
        int bin_x = -1;
        for (int j = 0; j < nbins_x; j++) {
            if (x >= x_bins[j] && x < x_bins[j + 1]) {
                bin_x = j;
                break;
            }
        }

        // Find the bin index for y
        int bin_y = -1;
        for (int j = 0; j < nbins_y; j++) {
            if (y >= y_bins[j] && y < y_bins[j + 1]) {
                bin_y = j;
                break;
            }
        }

        // Skip out-of-bounds points
        if (bin_x >= 0 && bin_y >= 0) {
            hist[bin_x][bin_y] += w;
            total_weight += w;
        }
    }

    // Normalize so that sum of all bins = 1
    if (total_weight > 0) {
        for (int i = 0; i < nbins_x; i++) {
            for (int j = 0; j < nbins_y; j++) {
                hist[i][j] /= total_weight;
            }
        }
    }

    return hist;
}

std::vector<std::vector<std::vector<double>>> make_joint_histogram(
    const std::vector<double>& data_x,
    const std::vector<double>& data_y,
    const std::vector<double>& data_z,
    const std::vector<double>& x_bins,
    const std::vector<double>& y_bins,
    const std::vector<double>& z_bins,
    std::function<double(double, double, double)> weight_fun
) {
    int nbins_x = x_bins.size() - 1;
    int nbins_y = y_bins.size() - 1;
    int nbins_z = z_bins.size() - 1;
    
    std::vector<std::vector<std::vector<double>>> hist(
        nbins_x, std::vector<std::vector<double>>(nbins_y, std::vector<double>(nbins_z, 0.0))
    );

    double total_weight = 0.0;

    for (size_t i = 0; i < data_x.size(); i++) {
        double x = data_x[i];
        double y = data_y[i];
        double z = data_z[i];
        double w = weight_fun(x, y, z);
        
        int bin_x = -1;
        for (int j = 0; j < nbins_x; j++) {
            if (x >= x_bins[j] && x < x_bins[j + 1]) {
                bin_x = j;
                break;
            }
        }

        int bin_y = -1;
        for (int j = 0; j < nbins_y; j++) {
            if (y >= y_bins[j] && y < y_bins[j + 1]) {
                bin_y = j;
                break;
            }
        }

        int bin_z = -1;
        for (int j = 0; j < nbins_z; j++) {
            if (z >= z_bins[j] && z < z_bins[j + 1]) {
                bin_z = j;
                break;
            }
        }

        if (bin_x >= 0 && bin_y >= 0 && bin_z >= 0) {
            hist[bin_x][bin_y][bin_z] += w;
            total_weight += w;
        }
    }

    // Normalize so that sum of all bins = 1
    if (total_weight > 0) {
        for (int i = 0; i < nbins_x; i++) {
            for (int j = 0; j < nbins_y; j++) {
                for (int k = 0; k < nbins_z; k++) {
                    hist[i][j][k] /= total_weight;
                }
            }
        }
    }

    return hist;
}

std::pair<std::vector<double>, std::vector<double>> make_histogram(const std::vector<double>& data, std::function<double(double)> weight_fun, int nbins) {
    std::pair<double, double> bounds = find_bounds(data);
    return make_histogram(data, weight_fun, nbins, bounds);
}

void make_histograms(const std::string& path, const std::vector<int>& frames, const std::string& path_out) {

}