#include "cli/CLIController.hpp"
#include <iostream>
#include <string>
#include <vector>

void print_usage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [options]\n\n"
              << "Options:\n"
              << "  --help, -h          Display this help message and exit\n"
              << "  --demo [1|2|3|all]  Run reproducible benchmark demonstration without interactive prompts:\n"
              << "                        1: Methane at 300 K, 5.0 MPa (Single-root gas / supercritical)\n"
              << "                        2: CO2 at 300 K, 7.0 MPa (Dense compressed fluid)\n"
              << "                        3: Propane at 320 K, 1.0 MPa (Three physical roots & phase preference)\n"
              << "                        all: Run all three demonstrations in sequence (default if omitted)\n\n"
              << "Interactive Mode:\n"
              << "  Run without arguments to launch the full interactive menu.\n";
}

int main(int argc, char* argv[]) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }

    thermo::cli::CLIController controller(std::cin, std::cout);

    if (args.empty()) {
        return controller.run();
    }

    if (args[0] == "--help" || args[0] == "-h") {
        print_usage(argv[0]);
        return 0;
    }

    if (args[0] == "--demo") {
        if (args.size() == 1 || args[1] == "all") {
            controller.run_demonstration(0);
            return 0;
        }
        if (args[1] == "1") {
            controller.run_demonstration(1);
            return 0;
        }
        if (args[1] == "2") {
            controller.run_demonstration(2);
            return 0;
        }
        if (args[1] == "3") {
            controller.run_demonstration(3);
            return 0;
        }
        std::cerr << "Error: Unknown demonstration target '" << args[1] << "'. Use 1, 2, 3, or all.\n";
        print_usage(argv[0]);
        return 1;
    }

    std::cerr << "Error: Unrecognized command-line argument '" << args[0] << "'.\n";
    print_usage(argv[0]);
    return 1;
}
