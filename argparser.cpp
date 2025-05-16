#include "argparser.h"

std::map<std::string, std::string> default_pars = {{"-samplingDistance", "0.25"},
{"-radialSamplingDistance", "-0.25"},
{"-linearRMSThreshold", "0.65"},
{"-minLineLength", "0.05"},
{"-radialFraction", "0.75"},
{"-radialFitThreshold", "0.1"},
{"-resolution", "128"},
{"-saveFailed", "0"}};

std::vector<std::string> parse(int argc, char *argv[], std::vector<std::string> options) {
    std::vector<std::string> values; values.resize(options.size());
    for (int j = 0; j < options.size(); j++) {
        std::string option = options[j];
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == option) {
                if (i + 1 >= argc) {
                    std::cout << "Option " << option << " not followed by a valid argument!" << std::endl;
                }
                values[j] = std::string(argv[i + 1]);
                break;
            }
        }
    }
    return values;
}

void sanitize_filename(std::string full_name, std::string& main, std::string& postfix, std::string& extension) {
    std::string filename = full_name.substr(full_name.find_last_of('/') + 1);
    auto p = filename.find_last_of(".");
    extension = filename.substr(p + 1);
    std::string file_without_extension = filename.substr(0, p);
    p = file_without_extension.find_first_of("_");
    main = file_without_extension.substr(0, p);
    postfix = file_without_extension.substr(p + 1);
}

void handle_pars(int argc, char *argv[], std::string mainName, std::string& outDir, std::string& outName, double& dL, double& dR, double& line_rms_threshold, double& line_length_threshold, double& radial_frac, double& radial_fit_threshold, int& resolution, bool& saveFailed) {
    std::vector<std::string> parsed = parse(argc, argv,
    {"-outName", "-outDir", "-samplingDistance", "-radialSamplingDistance", "-linearRMSThreshold",
    "-minLineLength", "-radialFraction", "-radialFitThreshold", "-resolution", "-saveFailed"});
    outName = parsed[0];
    outDir = parsed[1];
    if (outName == "") outName = mainName;
    if ((outDir != "") && (outDir.back() != '/')) outDir = outDir + "/";
    if (parsed[2] == "") dL = std::stod(default_pars["-samplingDistance"]);
    else dL = std::stod(parsed[2]);
    if (parsed[3] == "") dR = std::stod(default_pars["-radialSamplingDistance"]);
    else dR = std::stod(parsed[3]);
    if (parsed[4] == "") line_rms_threshold = std::stod(default_pars["-linearRMSThreshold"]);
    else line_rms_threshold = std::stod(parsed[4]);
    if (parsed[5] == "") line_length_threshold = std::stod(default_pars["-minLineLength"]);
    else line_length_threshold = std::stod(parsed[5]);
    if (parsed[6] == "") radial_frac = std::stod(default_pars["-radialFraction"]);
    else radial_frac = std::stod(parsed[6]);
    if (parsed[7] == "") radial_fit_threshold = std::stod(default_pars["-radialFitThreshold"]);
    else radial_fit_threshold = std::stod(parsed[7]);
    if (parsed[8] == "") resolution = std::stoi(default_pars["-resolution"]);
    else resolution = std::stoi(parsed[8]);
    if (parsed[9] == "") saveFailed = (default_pars["-saveFailed"] != "0");
    else {
        if ((parsed[9] == "true") || (parsed[9] == "True") || (parsed[9] == "TRUE") || (parsed[9] == "yes") || (parsed[9] == "Yes") || (parsed[9] == "YES") || (parsed[9] == "1")) saveFailed = true;
        else saveFailed = false;
    }
}

//char* path = std::getenv("PATH");
//file_in = std::ifstream("$PATH/default_pars.txt");