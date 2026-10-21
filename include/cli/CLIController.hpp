#ifndef THERMO_CLI_CONTROLLER_HPP
#define THERMO_CLI_CONTROLLER_HPP

#include "thermo/ComponentDatabase.hpp"
#include "thermo/PengRobinsonEOS.hpp"
#include "thermo/StateProperties.hpp"
#include <iostream>
#include <string>
#include <optional>
#include <vector>

namespace thermo::cli {

/**
 * @brief Interactive command-line interface controller for the Peng-Robinson EOS software.
 * Decoupled from std::cin/std::cout to allow deterministic automated testing via streams.
 */
class CLIController {
public:
    /**
     * @brief Constructs a CLI controller with provided input and output streams.
     * @param in Input stream (defaults to std::cin).
     * @param out Output stream (defaults to std::cout).
     */
    explicit CLIController(std::istream& in = std::cin, std::ostream& out = std::cout);

    /**
     * @brief Runs the interactive menu loop until Exit (5) is chosen or EOF is reached.
     * @return Exit code (0 for normal exit).
     */
    int run();

    /**
     * @brief Executes a specific demonstration benchmark directly without interactive prompts.
     * @param case_num 1 (Methane 300K 5MPa), 2 (CO2 300K 7MPa), 3 (Propane 320K 1MPa), or 0 for all.
     */
    void run_demonstration(int case_num);

    // Menu handlers (public for direct unit testing):
    void handle_pvt();
    void handle_thermo();
    void handle_fugacity();
    void handle_all_properties();

    // Display & Formatting utilities:
    void print_header(const std::string& title) const;
    void print_component_summary(const PureComponent& comp) const;
    void print_state_conditions(double T, double P, const PureComponent& comp) const;
    void print_eos_parameters(const EOSParameters& params, double A, double B, double discr) const;
    void print_root_pvt(size_t index, const PVTRoot& root, double B) const;
    void print_root_thermo(size_t index, const PVTRoot& root, const PureComponent& comp, double T) const;
    void print_root_fugacity(size_t index, const PVTRoot& root) const;
    void print_root_all(size_t index, const PVTRoot& root, const PureComponent& comp, double T, double B) const;
    void print_equilibrium_summary(const PVTResult& pvt) const;

    // Interactive Input Prompters:
    std::optional<PureComponent> prompt_component();
    std::optional<double> prompt_temperature();
    std::optional<double> prompt_pressure();

private:
    std::istream& in_;
    std::ostream& out_;
    ComponentDatabase db_;

    bool read_line(std::string& line);
    std::optional<int> prompt_int(const std::string& prompt_text,
                                  int min_val, int max_val,
                                  std::optional<int> default_val = std::nullopt);
    std::optional<double> prompt_positive_double(const std::string& prompt_text,
                                                 std::optional<double> default_val = std::nullopt);
};

} // namespace thermo::cli

#endif // THERMO_CLI_CONTROLLER_HPP
