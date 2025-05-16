#include "options.h"

void show_options() {
    std::cout << std::endl;
    std::cout << "Available arguments:" << std::endl;

    std::cout << "-inFilaments:\t\t path to the file containing filaments (mandatory)" << std::endl;

    std::cout << "-inDensity:\t\t path to the file containing density" << std::endl;

    std::cout << "-magnetized:\t\t if set true, magnetic fields required to run. If set false, magnetic field optional." << std::endl;

    std::cout << "-inHx:\t\t\t path to the file containing x-component of the magnetic field" << std::endl;

    std::cout << "-inHy:\t\t\t path to the file containing y-component of the magnetic field" << std::endl;

    std::cout << "-inHz:\t\t\t path to the file containing z-component of the magnetic field" << std::endl;

    std::cout << "-fieldsDir:\t\t path to the folder containing simulation fields (density, magnetic field)." << std::endl;
    std::cout << "\t\t\t -filename inferred from the naming scheme of the filament file." << std::endl;
    std::cout << "\t\t\t -do not combine with -inDensity, inHx, inHy, inHz" << std::endl;

    std::cout << "-outName:\t\t base output name (otherwise mimics the filament file naming scheme)." << std::endl;

    std::cout << "-outDir:\t\t output directory (default: same directory as this application)." << std::endl;

    std::cout << "-samplingDistance:\t sampling distance along the filament (default: quarter of a cell)" << std::endl;
    
    std::cout << "-radialSamplingDistance: radial sampling distance perpendicular to the filament, taken as a fraction" << std::endl;
    std::cout << "\t\t\t of -samplingDistance if negative. (default: -0.25 => 1/4th of the sampling Distance)" << std::endl;

    std::cout << "-linearRMSThreshold:\t how well we want to fit a filament with a line - smaller is stricter (default: 0.65)" << std::endl;

    std::cout << "-minLineLength:\t\t minimum acceptable linear filament length as a fraction of the box size (default: 0.05)" << std::endl;

    std::cout << "-radialFraction:\t fraction of the central value of the filament density in the radial direction" << std::endl;
    std::cout << "\t\t\t at which the fitter stops looking (default: 0.75)" << std::endl;

    std::cout << "-radialFitThreshold:\t radial profile RMS fit metric in %. Can be small (default: 0.1)" << std::endl;
    std::cout << std::endl;
}

bool check_pars(int argc, char *argv[], std::string& inFilaments, std::string& inDensity, std::string& inHx, std::string& inHy, std::string& inHz, std::string& mainName, bool& magnetized, std::string& message) {
    message = "";
    magnetized = true;
    if (argc == 1) {
        show_options();
        return false;
    }
    else {
        std::set<std::string> pars;
        std::string postfix, extension;
        for (int i = 1; i < argc; i++) { // first we load all arguments into a set for easy checking.
            pars.insert(argv[i]);
        }
        if (pars.find("-inFilaments") == pars.end()) { // if inFilaments option doesn't exist, we abort - no filaments could be loaded!
            message = "No filament file provided! Use option -inFilaments to provide the filament file! Aborting...";
            return false;
        }
        else { // even if filaments are specified, the file might no exist...
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-inFilaments") { // the next argv must be a valid file. Otherwise abort.
                    if ((i == argc - 1) || !std::filesystem::is_regular_file(argv[i + 1])) {
                        message = "-inFilaments must be followed by a valid file location! Aborting...";
                        return false;
                    }
                    else {
                        inFilaments = argv[i + 1];
                        // we break the full inFilaments path+filename into the main name (like cubeXXXX), the postfix (filaments_cut_1_...) and extension (.fits, .h5, ...)
                        sanitize_filename(inFilaments, mainName, postfix, extension);
                    }
                }
                
            }
        } // if we got past this point, we have the filaments file location saved in variable inFilaments
        if (pars.find("-magnetized") != pars.end()) { // by default the simulation is magnetized
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-magnetized") {
                    std::string val = argv[i + 1];
                    if ((val == "true") || (val == "True") || (val == "TRUE") || (val == "yes") || (val == "Yes") || (val == "YES") || val == "1") magnetized = true;
                    else if ((val == "false") || (val == "False") || (val == "FALSE") || (val == "no") || (val == "No") || (val == "NO") || val == "0") magnetized = false;
                    else {
                        message = "Unrecognized option \"" + val + "\" following -magnetized. Use keywords true or false to specify presence of magnetic field.";
                        return false;
                    }
                }
            }
        }
        if (pars.find("-fieldsDir") != pars.end()) { // if fieldsDir is specified...
            if ((pars.find("-inDensity") != pars.end()) || (pars.find("-inHx") != pars.end()) || (pars.find("-inHy") != pars.end()) || (pars.find("-inHz") != pars.end())) {
                message = "Directory for fields specified via -fieldsDir, options -inDensity, -inHx, -inHy, -inHz ignored.";
                // these options will be deduced automatically, so are ignored. This is just a warning
            }
            std::string fieldsDir;
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-fieldsDir") { // the next argv must be a valid path. Otherwise abort.
                    if ((i == argc - 1) || !std::filesystem::exists(argv[i + 1])) {
                        message = "-fieldsDir must be followed by a valid path! Aborting...";
                        return false;
                    }
                    else {
                        fieldsDir = argv[i + 1];
                        if (fieldsDir.back() != '/') fieldsDir = fieldsDir + "/"; // we found the correct fields directory, we need to check it contains all the necessary fields.
                    }
                }
            }
            // mainName from the sanitized filename is crucial here.
            inDensity = fieldsDir + mainName + "_density.fits"; // Collins naming scheme.
            inHx = fieldsDir + mainName + "_magnetic_field_x.fits"; // Collins naming scheme.
            inHy = fieldsDir + mainName + "_magnetic_field_y.fits"; // Collins naming scheme.
            inHz = fieldsDir + mainName + "_magnetic_field_z.fits"; // Collins naming scheme.
            if (!std::filesystem::exists(inDensity)) {
                message = "File " + mainName + "_density.fits does not exist in\n" + fieldsDir + ", use option -inDensity to specify file location! Aborting...";
                return false;
            }
            if (!std::filesystem::exists(inHx) && magnetized) {
                magnetized = false;
                message = "File " + mainName + "_magnetic_field_x.fits does not exist in\n" + fieldsDir + ", use option -inHx to specify file location! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
            if (!std::filesystem::exists(inHy) && magnetized) {
                message = "File " + mainName + "_magnetic_field_y.fits does not exist in\n" + fieldsDir + ", use option -inHy to specify file location! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
            if (!std::filesystem::exists(inHz) && magnetized) {
                message = "File " + mainName + "_magnetic_field_z.fits does not exist in\n" + fieldsDir + ", use option -inHz to specify file location! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
        }
        else { // if -fieldsDir isn't specified...
            if (pars.find("-inDensity") == pars.end()) { // first we ensure all fields have been provided
                message = "No -inDensity or -fieldsDir specified! Aborting...";
                return false;
            }
            if (pars.find("-inHx") == pars.end() && magnetized) {
                message = "No -inHx or -fieldsDir specified! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
            if (pars.find("-inHy") == pars.end() && magnetized) {
                message = "No -inHy or -fieldsDir specified! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
            if (pars.find("-inHz") == pars.end() && magnetized) {
                std::cout << "No -inHz or -fieldsDir specified! Alternatively, set -magnetized to false to skip magnetic fields. Aborting...";
                return false;
            }
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-inDensity") { // the next argv must be a valid file.
                    if ((i == argc - 1) || !std::filesystem::exists(argv[i + 1])) {
                        message = "-inDensity must be followed by a valid file! Aborting...";
                        return false;
                    }
                    else inDensity = argv[i + 1];
                    break;
                }
            }
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-inHx") {
                    if (((i == argc - 1) || !std::filesystem::exists(argv[i + 1])) && magnetized) {
                        message = "-inHx must be followed by a valid file! Aborting...";
                        return false;
                    }
                    else inHx = argv[i + 1];
                    break;
                }
            }
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-inHy") {
                    if (((i == argc - 1) || !std::filesystem::exists(argv[i + 1])) && magnetized) {
                        message = "-inHy must be followed by a valid file! Aborting...";
                        return false;
                    }
                    else inHy = argv[i + 1];
                    break;
                }
            }
            for (int i = 1; i < argc; i++) {
                if (std::string(argv[i]) == "-inHz") {
                    if (((i == argc - 1) || !std::filesystem::exists(argv[i + 1])) && magnetized) {
                        message = "-inHz must be followed by a valid file! Aborting...";
                        return false;
                    }
                    else inHz = argv[i + 1];
                    break;
                }
            }
        }
        return true;
    }
}