# CHL2002 Thermodynamics Project: Implementation Roadmap & Architecture

**Document Version:** 1.0  
**Phase:** Engineering Architecture & Milestone Planning  
**Target Completion:** Prior to November 10, 2026 Submission  

---

## 1. Modular C++ Architecture & Directory Layout

The application follows a modular, layered C++17 design adhering to the Single Responsibility Principle:

```text
Thermo Software/
├── CMakeLists.txt              # Top-level reproducible build configuration
├── README.md                   # User guide, build instructions, and verification summary
├── .gitignore                  # Git ignore rules for build artifacts and IDE cache
├── docs/
│   ├── PROJECT_SPECIFICATION.md # Complete equations, thermodynamic foundations, references
│   ├── IMPLEMENTATION_ROADMAP.md# Phased milestone schedule, dependencies, test plan
│   └── REPORT_DRAFT.md         # 4-page technical report draft mapped to assignment criteria
├── include/
│   └── thermo/
│       ├── Constants.hpp       # Universal constants (R, standard conditions, precision epsilons)
│       ├── Component.hpp       # PureComponent data model and Cp(T) integration methods
│       ├── ComponentDatabase.hpp# Curated catalog (8 gases) and custom component registry
│       ├── CubicSolver.hpp     # Analytical cubic equation solver (Cardano + Trigonometric)
│       ├── PengRobinsonEOS.hpp # PR parameters (a, b, da/dT), PVT calculations, departure functions
│       ├── StateProperties.hpp # Comprehensive thermodynamic state container (Z, v, H, S, phi, f)
│       └── CommandLineInterface.hpp # Interactive menu, input validation, formatted reporting
├── src/
│   ├── Component.cpp
│   ├── ComponentDatabase.cpp
│   ├── CubicSolver.cpp
│   ├── PengRobinsonEOS.cpp
│   ├── CommandLineInterface.cpp
│   └── main.cpp                # Application entry point
├── tests/
│   ├── CMakeLists.txt          # Test build target configuration
│   ├── test_cubic_solver.cpp   # Cubic polynomial solver unit tests
│   ├── test_pvt_eos.cpp        # PVT consistency, inverse PVT, ideal-gas limits
│   ├── test_thermo_derivatives.cpp # Finite-difference validation of da/dT
│   ├── test_departure_properties.cpp# Residual enthalpy, entropy, and fugacity benchmarks
│   └── test_literature_benchmarks.cpp# End-to-end validation against NIST and textbook cases
└── data/                       # Optional structured component definition files if required
```

---

## 2. Module Responsibilities & C++ Interface Contracts

```mermaid
graph TD
    Constants["Constants.hpp<br/>(Universal Constants & Reference States)"] --> Component["Component.hpp<br/>(Critical Properties & Ideal Cp)"]
    Component --> Database["ComponentDatabase.hpp<br/>(Curated Catalog & Custom Components)"]
    Constants --> CubicSolver["CubicSolver.hpp<br/>(Cardano / Trigonometric Root Engine)"]
    Component --> PR["PengRobinsonEOS.hpp<br/>(PVT, da/dT, Residual Functions)"]
    CubicSolver --> PR
    PR --> State["StateProperties.hpp<br/>(Full Thermodynamic State Record)"]
    Database --> CLI["CommandLineInterface.hpp<br/>(User Menus & Input Validation)"]
    PR --> CLI
    CLI --> Main["main.cpp<br/>(Application Entrypoint)"]
```

### 2.1 Interface Declarations

#### A. `Component.hpp`
```cpp
namespace thermo {

struct HeatCapacityPolynomial {
    double c0{0.0}, c1{0.0}, c2{0.0}, c3{0.0}, c4{0.0}; // J/(mol*K)
    double T_min{200.0}, T_max{1500.0};                 // K

    [[nodiscard]] double evaluate_cp(double T) const noexcept;
    [[nodiscard]] double integrate_cp(double T_ref, double T) const noexcept;
    [[nodiscard]] double integrate_cp_over_T(double T_ref, double T) const noexcept;
};

struct PureComponent {
    std::string name;
    std::string formula;
    double molar_mass;     // kg/mol
    double Tc;             // K
    double Pc;             // Pa
    double omega;          // dimensionless
    HeatCapacityPolynomial cp_model;
    std::string citation;  // Data provenance
};

} // namespace thermo
```

#### B. `CubicSolver.hpp`
```cpp
namespace thermo {

enum class RootType {
    SINGLE_PHYSICAL,
    VAPOR_LIKE,
    LIQUID_LIKE,
    UNSTABLE_INTERMEDIATE,
    NON_PHYSICAL
};

struct CubicRoot {
    double value{0.0};
    RootType type{RootType::NON_PHYSICAL};
    bool is_physical{false};
};

struct CubicSolution {
    int num_real_roots{0};
    std::vector<CubicRoot> roots;
    double discriminant{0.0};
};

class CubicSolver {
public:
    // Solves Z^3 + c2*Z^2 + c1*Z + c0 = 0
    static CubicSolution solve(double c2, double c1, double c0, double B_bound);
};

} // namespace thermo
```

#### C. `PengRobinsonEOS.hpp`
```cpp
namespace thermo {

struct ThermodynamicState {
    double T{0.0};                   // K
    double P{0.0};                   // Pa
    double Z{1.0};                   // dimensionless
    double molar_volume{0.0};        // m^3/mol
    double H_res{0.0};               // J/mol
    double S_res{0.0};               // J/(mol*K)
    double ln_phi{0.0};              // dimensionless
    double phi{1.0};                 // dimensionless
    double fugacity{0.0};            // Pa
    double H_total{0.0};             // J/mol (relative to T_ref, P_ref)
    double S_total{0.0};             // J/(mol*K) (relative to T_ref, P_ref)
    RootType phase_identification;
    bool is_stable_phase{true};
};

class PengRobinsonEOS {
public:
    explicit PengRobinsonEOS(PureComponent comp);

    [[nodiscard]] double calculate_P(double T, double v) const;
    [[nodiscard]] std::vector<ThermodynamicState> calculate_state_TP(double T, double P) const;
    [[nodiscard]] double a_parameter(double T) const;
    [[nodiscard]] double da_dT(double T) const;
    [[nodiscard]] double b_parameter() const;

private:
    PureComponent component_;
    double a_c_{0.0};
    double b_{0.0};
    double kappa_{0.0};
};

} // namespace thermo
```

---

## 3. Phased Implementation Roadmap

```mermaid
gantt
    title CHL2002 Development Roadmap
    dateFormat  YYYY-MM-DD
    section Specification & Foundations
    Session 1 : Spec, Equations & Architecture :done, 2026-10-20, 1d
    section Core Infrastructure
    Session 2 : Components, Database & Cubic Solver :2026-10-21, 2d
    section Thermodynamic Engine
    Session 3 : PR EOS Forward/Inverse PVT :2026-10-23, 2d
    Session 4 : Departure Enthalpy, Entropy & Fugacity :2026-10-25, 2d
    section Testing & Benchmarks
    Session 5 : Automated Test Suite & NIST Validation :2026-10-27, 2d
    section UI & Reporting
    Session 6 : Interactive CLI & Error Handling :2026-10-29, 2d
    Session 7 : 4-Page Technical Report & Viva Prep :2026-10-31, 3d
```

### Phase Breakdown & Dependencies

| Phase / Session | Objective | Key Deliverables | Pre-requisites | Measurable Completion Criteria |
| :--- | :--- | :--- | :--- | :--- |
| **Session 1** (Current) | Specification, Foundations & Git Repo | `PROJECT_SPECIFICATION.md`, `IMPLEMENTATION_ROADMAP.md`, `.gitignore`, git init | Course Prompt | Approved engineering specification; zero unverified equations. |
| **Session 2** | Core Data Types & Cubic Solver | `Constants.hpp`, `Component.hpp`, `ComponentDatabase.cpp`, `CubicSolver.cpp` | Session 1 | Solves known cubic polynomials with residuals $< 10^{-12}$; 8 gases loaded. |
| **Session 3** | PR EOS Engine & PVT Calculations | `PengRobinsonEOS.cpp` (PVT forward & inverse), CMake infrastructure | Session 2 | $P(T, v(T, P)) == P$ within machine precision; root identification verified. |
| **Session 4** | Departure Functions & Fugacity | $H^{\text{res}}, S^{\text{res}}, \phi, f$, small-$B$ regularization, ideal $C_p$ integrals | Session 3 | $\lim_{P\to 0}\phi = 1.0, H^{\text{res}} \to 0$; $da/dT$ finite-difference error $< 10^{-7}$. |
| **Session 5** | Test Suite & Literature Validation | Complete test suite in `tests/`, benchmark validation tables | Session 4 | 100% test pass rate across unit tests and NIST benchmark cases. |
| **Session 6** | CLI & Usability Polish | Interactive command-line interface with formatted outputs, input guardrails | Session 5 | Robust user workflow, handles bad input without crashing, viva ready. |
| **Session 7** | Final Technical Report & Archive | 4-page LaTeX/Markdown report, final GitHub release | Session 6 | Technical report fully aligned with actual program outputs. |

---

## 4. Test & Validation Plan

Every calculation will be validated using four independent layers of verification:

### 4.1 Layer 1: Mathematical Limits & Asymptotic Checks
1. **Ideal Gas Limit ($P \to 0$):**
   - At $T = 300\text{ K}$, evaluate at $P = 10\text{ Pa}, 1\text{ Pa}, 10^{-2}\text{ Pa}$:
     - $Z \to 1.00000000$ (residual $< 10^{-7}$)
     - $\phi \to 1.00000000$ (residual $< 10^{-7}$)
     - $H^{\text{res}} \to 0.00\text{ J/mol}$
     - $S^{\text{res}} \to 0.00\text{ J/(mol*K)}$
2. **Cubic Root Inversion:**
   - Every evaluated root $Z$ must satisfy $|Z^3 + c_2 Z^2 + c_1 Z + c_0| < 10^{-12}$.
3. **EOS Reciprocity Check:**
   - Compute $v = Z R T / P$. Evaluate $P_{\text{calc}} = P(T, v)$.
   - Acceptance: $|P_{\text{calc}} - P| / P < 10^{-11}$.

### 4.2 Layer 2: Analytical Derivative Cross-Validation
The analytical derivative $\frac{da}{dT}$ is verified against central finite differences:
$$\left( \frac{da}{dT} \right)_{\text{num}} = \frac{a(T + \delta) - a(T - \delta)}{2\delta}, \quad \delta = 10^{-5} T$$
- Acceptance criteria: Relative discrepancy $\left| \frac{(da/dT)_{\text{ana}} - (da/dT)_{\text{num}}}{(da/dT)_{\text{ana}}} \right| < 10^{-8}$.

### 4.3 Layer 3: Traceable Literature Benchmarks

| Case ID | Substance | State ($T, P$) | Expected Property | Published / Reference Value | Source / Benchmark Reference |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **BM-01** | **Methane** | $300\text{ K}, 5.0\text{ MPa}$ | $Z$<br>$v$<br>$\phi$<br>$H^{\text{res}}$ | $0.9329$<br>$4.654 \times 10^{-4}\text{ m}^3/\text{mol}$<br>$0.9338$<br>$-756.2\text{ J/mol}$ | Smith, Van Ness, Abbott (App. C) / PR-1976 baseline |
| **BM-02** | **Carbon Dioxide** | $300\text{ K}, 4.0\text{ MPa}$ | Phase: Vapor<br>$Z$<br>$\phi$ | Vapor<br>$0.6974$<br>$0.7670$ | Poling et al. Table 3-3 / NIST Chemistry WebBook |
| **BM-03** | **Carbon Dioxide** | $300\text{ K}, 7.0\text{ MPa}$ | Phase: 3 Roots<br>$Z_{\text{liq}}, Z_{\text{vap}}$<br>Stable Phase | $Z_{\text{liq}} \approx 0.088$<br>$Z_{\text{vap}} \approx 0.384$<br>Liquid ($\phi_{\text{liq}} < \phi_{\text{vap}}$) | Subcritical two-phase region near saturation |
| **BM-04** | **Propane** | $320\text{ K}, 1.0\text{ MPa}$ | $Z$<br>$\phi$<br>$H^{\text{res}}$ | $0.8524$<br>$0.8576$<br>$-1685.4\text{ J/mol}$ | Perry's Chemical Engineers' Handbook (9th Ed.) |

### 4.4 Layer 4: Edge Cases & Numerical Robustness
- **Supercritical Region:** $T \gg T_c, P \gg P_c$ (single supercritical root, no unphysical branches).
- **Near-Critical Point:** $T \to T_c, P \to P_c$ (triple root coalescence, stable evaluation without discriminant sign oscillation).
- **Extreme Low-Density:** $P \to 0$ without numeric overflow in $1/B$.
- **Input Guardrails:** Prevent non-positive temperatures ($T \le 0$), pressures ($P \le 0$), or negative acentric factors ($\omega < -1$).

---

## 5. Report Structure & Page Budget Plan (Max 4 Pages)

The technical report will be structured to fit strictly within the 4-page limit:

1. **Page 1: Introduction, Objectives, and Thermodynamic Foundations**
   - Scope of project and problem statement.
   - Governing Peng–Robinson EOS formulations and parameter derivations ($a_c, b, \alpha(T), \kappa, da/dT$).
   - Dimensionless polynomial formulation in $Z$.
2. **Page 2: Residual Properties, Pure Fugacity, and Reference States**
   - Exact derivation and analytical equations for $H^{\text{res}}$, $S^{\text{res}}$, and $\ln \phi$.
   - Ideal gas heat capacity integration and reference state definition ($T_{\text{ref}} = 298.15\text{ K}, P_{\text{ref}} = 1.0\text{ bar}$).
   - Low-pressure asymptotic regularization.
3. **Page 3: Software Architecture, Numerical Algorithm, and Stability Criterion**
   - C++ modular design and execution workflow.
   - Depressed cubic Cardano and Vieta trigonometric root solver.
   - Physical root filtering ($Z > B$) and equal-fugacity phase stability criterion.
4. **Page 4: Verification Results, Discussion, Limitations, and Conclusions**
   - Tabulated verification against NIST / textbook reference cases.
   - Derivative finite-difference validation and reciprocity checks.
   - Model limitations (pure components only, polar fluids, supercritical liquid density accuracy).
   - Traceable bibliography and references.
