# CHL2002 Thermodynamics Project — Viva Voce Preparation Guide

**Document Version:** 1.0  
**Author:** Harshith N  
**Course:** CHL2002 Thermodynamics  
**Project:** Pure-Component Peng–Robinson (PR-1976) Equation of State Software  
**Repository:** `https://github.com/harsh-bitverse/thermo-eos-solver.git`  
**Target Audience:** Student preparation for oral examination / defense with course professor or external evaluators.

---

## Table of Contents

1. [Core Thermodynamics Questions and Answers](#1-core-thermodynamics-questions-and-answers)
   - [Q1. Why do real gases deviate from ideal behavior, and how does PR-1976 capture that?](#q1-why-do-real-gases-deviate-from-ideal-behavior-and-how-does-pr-1976-capture-that)
   - [Q2. What do parameters $a(T)$ and $b$ represent physically? Why is $b$ constant while $a$ depends on $T$?](#q2-what-do-parameters-at-and-b-represent-physically-why-is-b-constant-while-a-depends-on-t)
   - [Q3. What are the physical meanings of $T_c$, $P_c$, and acentric factor $\omega$?](#q3-what-are-the-physical-meanings-of-t_c-p_c-and-acentric-factor-omega)
   - [Q4. Why does the Peng–Robinson EOS yield a cubic polynomial in $Z$?](#q4-why-does-the-peng-robinson-eos-yield-a-cubic-polynomial-in-z)
   - [Q5. What are residual properties, and why are departure functions useful?](#q5-what-are-residual-properties-and-why-are-departure-functions-useful)
   - [Q6. What is fugacity, and why is fugacity equality the equilibrium condition instead of pressure alone?](#q6-what-is-fugacity-and-why-is-fugacity-equality-the-equilibrium-condition-instead-of-pressure-alone)
   - [Q7. Why are there three roots under subcritical conditions, and why is the middle root unphysical?](#q7-why-are-there-three-roots-under-subcritical-conditions-and-why-is-the-middle-root-unphysical)
   - [Q8. How is phase preference determined, and why is the phase with lower fugacity coefficient preferred?](#q8-how-is-phase-preference-determined-and-why-is-the-phase-with-lower-fugacity-coefficient-preferred)
   - [Q9. How does phase preference differ from saturation pressure calculation?](#q9-how-does-phase-preference-differ-from-saturation-pressure-calculation)
   - [Q10. Why cannot total $H$ and $S$ be computed for steam relative to 298.15 K using NIST Shomate?](#q10-why-cannot-total-h-and-s-be-computed-for-steam-relative-to-29815-k-using-nist-shomate)
2. [Numerical and Implementation Questions](#2-numerical-and-implementation-questions)
   - [Q11. How does the cubic solver work? Why Cardano + Vieta? What determines the number of roots?](#q11-how-does-the-cubic-solver-work-why-cardano--vieta-what-determines-the-number-of-roots)
   - [Q12. How are unphysical roots ($Z \le B$) filtered?](#q12-how-are-unphysical-roots-z-le-b-filtered)
   - [Q13. How is mechanical stability evaluated in code?](#q13-how-is-mechanical-stability-evaluated-in-code)
   - [Q14. What is the role of analytical $da/dT$, and what would happen if evaluated numerically?](#q14-what-is-the-role-of-analytical-dadt-and-what-would-happen-if-evaluated-numerically)
   - [Q15. How does the low-pressure Taylor series expansion work, and why is it necessary?](#q15-how-does-the-low-pressure-taylor-series-expansion-work-and-why-is-it-necessary)
   - [Q16. What do the automated tests prove, and what do they NOT prove?](#q16-what-do-the-automated-tests-prove-and-what-do-they-not-prove)
3. [Code Walkthrough & Equation-to-Code Mapping](#3-code-walkthrough--equation-to-code-mapping)
   - [Thermodynamic Equations to Code Architecture Mapping Table](#thermodynamic-equations-to-code-architecture-mapping-table)
   - [Step-by-Step Execution Trace of a Calculation](#step-by-step-execution-trace-of-a-calculation)

---

## 1. Core Thermodynamics Questions and Answers

### Q1. Why do real gases deviate from ideal behavior, and how does PR-1976 capture that?

**Answer:**
The ideal gas law ($P v = R T$) relies on two fundamental assumptions:
1. Gas molecules are point masses that occupy zero physical volume ($v_{\text{mol}} = 0$).
2. There are zero intermolecular attractive or repulsive forces between molecules ($U_{\text{inter}} = 0$).

At low pressure and high temperature, molecules are widely separated, rendering these assumptions valid. However, at elevated pressures or reduced temperatures:
- **Repulsive effects:** Molecules occupy non-zero volume (finite electron cloud exclusion). The volume accessible to molecular motion is reduced from the total molar volume $v$ to the free volume $(v - b)$, where $b$ is the repulsive co-volume.
- **Attractive effects:** Attractive forces (London dispersion, dipole interactions) pull molecules inward, reducing the momentum transfer (pressure) imparted to the container walls.

The Peng–Robinson (1976) cubic equation of state incorporates both effects:
$$P = \underbrace{\frac{R T}{v - b}}_{\text{Repulsive pressure}} - \underbrace{\frac{a(T)}{v^2 + 2bv - b^2}}_{\text{Attractive pressure}}$$

Peng and Robinson selected the quadratic denominator $(v^2 + 2bv - b^2)$ specifically to improve upon the Soave–Redlich–Kwong (SRK) EOS by yielding a universal critical compressibility factor $Z_c \approx 0.3074$ (closer to experimental values for non-polar hydrocarbons $\approx 0.27\text{--}0.29$ than SRK's $0.333$) and significantly superior liquid molar volume predictions.

---

### Q2. What do parameters $a(T)$ and $b$ represent physically? Why is $b$ constant while $a$ depends on $T$?

**Answer:**
- **Co-volume parameter $b$:** Represents the hard-core repulsive volume physically occupied by one mole of molecules (the excluded volume, roughly 4 times the actual van der Waals volume of the molecules). Because atomic and electronic radii are largely insensitive to thermal kinetic energies across normal engineering temperature ranges, $b$ is treated as a temperature-independent constant for a given fluid:
  $$b = \Omega_b \frac{R T_c}{P_c}, \quad \text{where } \Omega_b \approx 0.0777960739$$
- **Attractive parameter $a(T)$:** Represents the magnitude of intermolecular cohesive attraction. At elevated temperatures, molecular thermal velocity increases significantly. Faster collisions and increased rotational/vibrational agitation disrupt directional attractive alignments (London dispersion and dipole interactions), diminishing their net pressure-reduction effect. Consequently, the effective attractive parameter decreases monotonically as temperature rises:
  $$a(T) = a_c \cdot \alpha(T_r, \omega)$$
  $$a_c = \Omega_a \frac{R^2 T_c^2}{P_c}, \quad \Omega_a \approx 0.4572355289$$
  $$\alpha(T_r, \omega) = \left[1 + \kappa \left(1 - \sqrt{T_r}\right)\right]^2$$

---

### Q3. What are the physical meanings of $T_c$, $P_c$, and acentric factor $\omega$?

**Answer:**
- **Critical Temperature ($T_c$):** The maximum temperature at which a distinct liquid phase can exist. Above $T_c$, thermal kinetic energy exceeds all intermolecular cohesion; gas cannot be liquefied regardless of applied pressure.
- **Critical Pressure ($P_c$):** The vapor pressure at the critical temperature. At $(T_c, P_c)$, the properties of coexisting liquid and vapor phases (density, refractive index, enthalpy) become completely identical, and the meniscus disappears.
- **Acentric Factor ($\omega$):** Introduced by Kenneth Pitzer (1955) to quantify the non-sphericity (acentricity) of a molecule's force field. It is defined from the reduced vapor pressure at a reduced temperature of $T_r = 0.7$:
  $$\omega = -\log_{10}\left(\frac{P^{\text{sat}}(T_r = 0.7)}{P_c}\right) - 1.0$$
  - For simple, spherical monoatomic gases (Argon, Krypton, Xenon), $P^{\text{sat}}/P_c \approx 0.1$ at $T_r = 0.7$, so $\omega \approx 0$.
  - For elongated hydrocarbons, branched isomers, or polar fluids, intermolecular forces are anisotropic, increasing molecular cohesion and lowering the vapor pressure at $T_r = 0.7$, which produces $\omega > 0$ (e.g., Propane $\omega = 0.152$, Water $\omega = 0.344$, n-Heptane $\omega = 0.349$).

---

### Q4. Why does the Peng–Robinson EOS yield a cubic polynomial in $Z$?

**Answer:**
Starting with the pressure equation:
$$P = \frac{R T}{v - b} - \frac{a}{v^2 + 2bv - b^2}$$
Clearing denominators by multiplying both sides by $(v - b)(v^2 + 2bv - b^2)$:
$$P (v - b)(v^2 + 2bv - b^2) = R T (v^2 + 2bv - b^2) - a (v - b)$$
Expanding this expression produces a 3rd-degree (cubic) polynomial in molar volume $v$:
$$P v^3 + (b P - R T) v^2 + (a - 3 b^2 P - 2 b R T) v + (b^3 P + b^2 R T - a b) = 0$$

Substituting $v = \frac{Z R T}{P}$ and non-dimensionalizing using $A = \frac{a P}{(R T)^2}$ and $B = \frac{b P}{R T}$ converts the equation to the standard dimensionless cubic polynomial in compressibility factor $Z$:
$$Z^3 - (1 - B) Z^2 + (A - 3B^2 - 2B) Z - (A B - B^2 - B^3) = 0$$

**Thermodynamic necessity of degree 3:** A cubic polynomial is the lowest-order algebraic equation capable of capturing an inflection point with horizontal tangents ($(\partial P/\partial v)_T = 0$ and $(\partial^2 P/\partial v^2)_T = 0$ at the critical point) and supporting two stable phases (liquid and vapor) connected by an unstable transition loop.

---

### Q5. What are residual properties, and why are departure functions useful?

**Answer:**
- **Definition:** A residual property $M^{\text{res}}(T, P)$ is defined as the difference between the actual molar thermodynamic property of the real fluid and that of an ideal gas at the *exact same temperature and pressure*:
  $$M^{\text{res}}(T, P) \equiv M(T, P) - M^{\text{ig}}(T, P)$$
- **Why departure functions are useful:**
  Thermodynamic properties cannot be measured directly across all $(T, P)$ from first principles. However:
  1. **Ideal-gas properties ($M^{\text{ig}}$):** Depend purely on intramolecular vibrations, rotations, and ideal translations. For an ideal gas, enthalpy $H^{\text{ig}}$ and heat capacity $C_p^\circ$ depend strictly on temperature, and entropy has a simple pressure dependence $-R \ln(P / P_{\text{ref}})$. These are readily obtained from spectroscopic data and empirical heat-capacity correlations (such as NIST Shomate polynomials).
  2. **Equation of State:** The EOS provides the residual departure $M^{\text{res}}(T, P)$ by integrating volumetric departures from infinite volume ($P \to 0$, where real fluids behave ideally) to the real density:
     $$H^{\text{res}} = R T (Z - 1) + \frac{T \frac{da}{dT} - a}{2\sqrt{2}b} \ln\left(\frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B}\right)$$
  Combining the two layers allows evaluation of total state properties at any pressure and temperature:
  $$M(T, P) = M^{\text{ig}}(T, P) + M^{\text{res}}(T, P)$$

---

### Q6. What is fugacity, and why is fugacity equality the equilibrium condition instead of pressure alone?

**Answer:**
- **Chemical Potential:** The rigorous criterion for phase equilibrium between two phases at the same temperature and pressure is the equality of chemical potentials:
  $$\mu_{\text{liq}}(T, P) = \mu_{\text{vap}}(T, P)$$
- **Fugacity definition:** Chemical potential approaches $-\infty$ as pressure approaches zero ($\lim_{P \to 0} \mu = -\infty$), making it inconvenient for numerical computation. Gilbert N. Lewis introduced **fugacity** $f$ as an "effective thermodynamic pressure" (units of Pa) defined by:
  $$d\mu = R T \, d\ln f \quad (\text{at constant } T)$$
  with the ideal-gas boundary condition:
  $$\lim_{P \to 0} \frac{f}{P} = 1$$
- **Fugacity coefficient ($\phi$):** The dimensionless ratio $\phi = f / P$.
- **Why pressure equality alone is insufficient:**
  Pressure equality ($P_{\text{liq}} = P_{\text{vap}}$) guarantees **mechanical equilibrium** (the interface does not accelerate). Temperature equality ($T_{\text{liq}} = T_{\text{vap}}$) guarantees **thermal equilibrium** (no net heat transfer).
  However, **phase equilibrium** also requires **diffusive/chemical equilibrium** (no net transfer of mass across the interface). At the same $T$ and $P$:
  $$\mu_{\text{liq}} - \mu_{\text{vap}} = R T \ln\left(\frac{f_{\text{liq}}}{f_{\text{vap}}}\right) = R T \ln\left(\frac{\phi_{\text{liq}}}{\phi_{\text{vap}}}\right) = 0 \iff f_{\text{liq}} = f_{\text{vap}} \iff \phi_{\text{liq}} = \phi_{\text{vap}}$$
  If $\phi_{\text{liq}} \neq \phi_{\text{vap}}$, molecules will spontaneously evaporate or condense until the chemical potentials equalize.

---

### Q7. Why are there three roots under subcritical conditions, and why is the middle root unphysical?

**Answer:**
- **Origin of three roots:** For temperatures below the critical temperature ($T < T_c$) and within the subcritical spinodal envelope, the cubic polynomial exhibits three real roots ($\Delta < 0$):
  - Root 1 ($Z_{\text{liq}}$): Smallest compressibility factor $\to$ high density $\to$ liquid-like state.
  - Root 2 ($Z_{\text{mid}}$): Intermediate compressibility factor $\to$ unphysical spinodal state.
  - Root 3 ($Z_{\text{vap}}$): Largest compressibility factor $\to$ low density $\to$ vapor-like state.
- **Why the intermediate root is unphysical (Mechanical Instability):**
  The intrinsic thermodynamic criterion for mechanical stability is that the fluid must resist compression: increasing pressure must reduce volume, meaning the isothermal compressibility must be positive:
  $$\kappa_T = -\frac{1}{v} \left(\frac{\partial v}{\partial P}\right)_T > 0 \iff \left(\frac{\partial P}{\partial v}\right)_T < 0$$
  At the intermediate root, the slope of the van der Waals isotherm loop is positive:
  $$\left(\frac{\partial P}{\partial v}\right)_T > 0$$
  If a substance had $(\partial P/\partial v)_T > 0$, an infinitesimal compressive fluctuation would reduce the opposing internal pressure, causing runaway collapse into liquid; conversely, an infinitesimal expansion would increase internal pressure, causing runaway expansion into vapor (spinodal decomposition). The intermediate root is therefore physically unobservable in static equilibrium.

---

### Q8. How is phase preference determined, and why is the phase with lower fugacity coefficient preferred?

**Answer:**
At fixed temperature and pressure, the second law of thermodynamics dictates that a closed system reaches equilibrium at the global minimum of the Gibbs free energy $G$.

For any phase root $k$:
$$G_k(T, P) = G^{\text{ig}}(T, P) + G^{\text{res}}_k(T, P) = G^{\text{ig}}(T, P) + R T \ln \phi_k$$

Comparing the liquid and vapor roots at the identical $(T, P)$:
$$\Delta G^{\text{res}} = G_{\text{liq}} - G_{\text{vap}} = R T \ln\left(\frac{\phi_{\text{liq}}}{\phi_{\text{vap}}}\right)$$

- **Case 1: $\phi_{\text{vap}} < \phi_{\text{liq}}$ ($\Delta G^{\text{res}} > 0$):**
  The vapor phase has lower Gibbs free energy ($G_{\text{vap}} < G_{\text{liq}}$) and lower chemical potential. The vapor phase is thermodynamically preferred (superheated vapor).
- **Case 2: $\phi_{\text{liq}} < \phi_{\text{vap}}$ ($\Delta G^{\text{res}} < 0$):**
  The liquid phase has lower Gibbs free energy ($G_{\text{liq}} < G_{\text{vap}}$). The liquid phase is thermodynamically preferred (subcooled liquid).
- **Case 3: $\phi_{\text{liq}} = \phi_{\text{vap}}$ within tolerance ($\Delta G^{\text{res}} = 0$):**
  Chemical potentials are equal. Both phases coexist in thermodynamic equilibrium at the saturation pressure $P^{\text{sat}}(T)$.

---

### Q9. How does phase preference differ from saturation pressure calculation?

**Answer:**
- **Phase preference (State Evaluation):** Evaluates the stability of existing cubic roots at a *single specified state point* $(T, P)$. The software solves the cubic at that given $(T, P)$, extracts the physical roots, and compares $\phi_{\text{liq}}$ and $\phi_{\text{vap}}$ to report whether the state is subcooled liquid, superheated vapor, or saturated.
- **Saturation pressure calculation (Iterative Boundary Search):** Requires finding the unknown pressure $P^{\text{sat}}$ at temperature $T$ that satisfies the non-linear equation:
  $$f_{\text{liq}}(T, P^{\text{sat}}) - f_{\text{vap}}(T, P^{\text{sat}}) = 0$$
  This requires an iterative numerical method (such as Newton–Raphson or successive substitution) iterating across multiple pressure evaluations. The current project evaluates states directly; when provided with the literature saturation pressure (e.g. BM-06: Propane at 320 K, 1.6036 MPa), it validates that $|\phi_L - \phi_V|/\phi_V \le 10^{-5}$.

---

### Q10. Why cannot total $H$ and $S$ be computed for steam relative to 298.15 K using NIST Shomate?

**Answer:**
- The NIST Shomate heat-capacity polynomial for gas-phase Water Vapor is parameterized and experimentally validated only over the temperature range:
  $$T \in [500.0, 1700.0]\text{ K}$$
- The project's universal reference state is defined at $T_{\text{ref}} = 298.15\text{ K}$ and $P_{\text{ref}} = 100\,000\text{ Pa}$ (1.0 bar).
- At 298.15 K and 1 bar, water is physically a liquid. Evaluating an ideal gas heat capacity polynomial below its valid domain ($298.15\text{ K} < 500.0\text{ K}$) would constitute unphysical extrapolation of high-order polynomial coefficients ($t = T/1000$).
- Rather than printing an unvalidated or misleading number, the C++ code enforces domain safety: `integrate_cp` validates temperature boundaries and throws `std::out_of_range`. The CLI catches this and clearly reports `H: N/A` and `S: N/A` with an explanatory diagnostic message, while all EOS-level properties ($Z, v, (\partial P/\partial v)_T, H^{\text{res}}, S^{\text{res}}, \phi, f$) remain exact.

---

## 2. Numerical and Implementation Questions

### Q11. How does the cubic solver work? Why Cardano + Vieta? What determines the number of roots?

**Answer:**
Given the cubic equation $Z^3 + a_2 Z^2 + a_1 Z + a_0 = 0$:
1. **Tschirnhaus Transformation:** Eliminates the quadratic term using $Z = y - \frac{a_2}{3}$, converting the cubic into depressed form:
   $$y^3 + p y + q = 0$$
   where:
   $$p = a_1 - \frac{a_2^2}{3}, \quad q = a_0 - \frac{a_1 a_2}{3} + \frac{2 a_2^3}{27}$$
2. **Depressed Discriminant:**
   $$\Delta = \left(\frac{q}{2}\right)^2 + \left(\frac{p}{3}\right)^3$$
3. **Cardano's Method ($\Delta > 0$):**
   One real root and two complex conjugate roots. The single real root is:
   $$y_1 = \sqrt[3]{-\frac{q}{2} + \sqrt{\Delta}} + \sqrt[3]{-\frac{q}{2} - \sqrt{\Delta}}$$
   Solved using `std::cbrt` to avoid negative root domain issues.
4. **Vieta's Trigonometric Method ($\Delta \le 0$):**
   When $\Delta \le 0$, three real roots exist (*casus irreducibilis*). Cardano's formula would require taking cube roots of complex numbers, causing rounding errors and floating-point ambiguity. Vieta's trigonometric substitution ($y = 2 \sqrt{-p/3} \cos\theta$) uses the identity $4\cos^3\theta - 3\cos\theta = \cos(3\theta)$:
   $$\cos(3\theta) = -\frac{q / 2}{\sqrt{(-p/3)^3}}$$
   $$\theta = \frac{1}{3} \arccos\left(\text{clamp}\left(-\frac{q / 2}{\sqrt{(-p/3)^3}}, -1, 1\right)\right)$$
   The three real roots are:
   $$y_k = 2 \sqrt{-\frac{p}{3}} \cos\left(\theta + \frac{2\pi k}{3}\right), \quad k \in \{0, 1, 2\}$$
5. **Back-shift & Sort:** $Z_k = y_k - \frac{a_2}{3}$, followed by ascending sorting ($Z_0 \le Z_1 \le Z_2$).

---

### Q12. How are unphysical roots ($Z \le B$) filtered?

**Answer:**
- In the Peng–Robinson EOS, $B = \frac{b P}{R T}$ and molar volume is $v = \frac{Z R T}{P}$.
- The condition $Z > B$ is mathematically identical to $v > b$.
- Because $b$ represents the impenetrable co-volume of the molecules, any root with $v \le b$ corresponds to negative free volume $(v - b \le 0)$, violating the physical impossibility of interpenetrating matter.
- Furthermore, the term $\ln(Z - B)$ in residual entropy and fugacity:
  $$\ln\phi = (Z - 1) - \ln(Z - B) - \frac{A}{2\sqrt{2}B} \ln\left(\frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B}\right)$$
  is undefined ($\ln(\le 0)$) for $Z \le B$.
- In `PengRobinsonEOS.cpp`, roots returned by `CubicSolver::solve` are strictly filtered:
  ```cpp
  if (root_val > B + 1e-12) {
      // Physically admissible root retained
  }
  ```

---

### Q13. How is mechanical stability evaluated in code?

**Answer:**
In `PengRobinsonEOS::calculate_dP_dv(double T, double v)`:
$$\left(\frac{\partial P}{\partial v}\right)_T = -\frac{R T}{(v - b)^2} + \frac{2 a(T) (v + b)}{(v^2 + 2bv - b^2)^2}$$
- If $\left(\frac{\partial P}{\partial v}\right)_T < 0$: Classified as `MechanicalStability::STABLE`. Compressing the fluid increases pressure, resisting displacement.
- If $\left(\frac{\partial P}{\partial v}\right)_T \ge 0$: Classified as `MechanicalStability::UNSTABLE`. The root lies on the unphysical spinodal branch where compression reduces pressure.

---

### Q14. What is the role of analytical $da/dT$, and what would happen if evaluated numerically?

**Answer:**
The temperature derivative of the attractive parameter appears directly in residual enthalpy and entropy:
$$H^{\text{res}} = R T (Z - 1) + \frac{T \frac{da}{dT} - a(T)}{2\sqrt{2}b} L$$
$$S^{\text{res}} = R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} L$$

In PR-1976, $a(T) = a_c \left[1 + \kappa(1 - \sqrt{T/T_c})\right]^2$. Its exact analytical derivative is:
$$\frac{da}{dT} = -\frac{a_c \, \kappa \sqrt{\alpha(T)}}{\sqrt{T \cdot T_c}}$$

- **Why analytical is superior:** It executes in ~5 CPU cycles and has zero truncation error.
- **Failure of numerical differentiation:** A finite-difference approximation $\frac{a(T+h) - a(T-h)}{2h}$ requires tuning step size $h$. If $h$ is too large, $O(h^2)$ truncation error corrupts $H^{\text{res}}$ and $S^{\text{res}}$. If $h$ is too small, subtractive floating-point cancellation destroys precision. Near critical temperatures ($T \approx T_c$), numerical derivatives can introduce significant relative errors.

---

### Q15. How does the low-pressure Taylor series expansion work, and why is it necessary?

**Answer:**
The PR residual departure functions contain the quotient:
$$I_L = \frac{L}{2\sqrt{2}B} = \frac{1}{2\sqrt{2}B} \ln\left(\frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B}\right)$$

As pressure approaches zero ($P \to 0$):
- $B = \frac{b P}{R T} \to 0$
- $Z \to 1$
- The argument of the natural logarithm approaches $1$, causing the numerator $L \to 0$ and the denominator $2\sqrt{2}B \to 0$.
- Direct evaluation yields $\frac{0}{0}$ indeterminate form, causing severe subtractive cancellation and catastrophic precision loss.

To resolve this, we define:
$$u \equiv \frac{2\sqrt{2}B}{Z + (1-\sqrt{2})B}$$
Then the logarithmic argument is exactly $1 + u$, and:
$$I_L = \frac{1}{Z + (1-\sqrt{2})B} \cdot \frac{\ln(1 + u)}{u}$$
For $|u| < 10^{-5}$, we evaluate $\frac{\ln(1+u)}{u}$ via its 4-term Taylor expansion:
$$\frac{\ln(1+u)}{u} = 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4} + \mathcal{O}(u^4)$$
This guarantees smooth, machine-precision continuity down to $P = 10^{-5}\text{ Pa}$, eliminating any floating-point singularity.

---

### Q16. What do the automated tests prove, and what do they NOT prove?

**Answer:**
- **What the tests prove:**
  1. **Algorithmic Correctness:** The C++ code implements the PR-1976 equations, Shomate integrations, and cubic solver with complete mathematical fidelity.
  2. **Internal Thermodynamic Consistency:** Differential Maxwell relations $\left(\frac{\partial \ln\phi}{\partial P}\right)_T = \frac{Z-1}{P}$ and Gibbs–Helmholtz departures $H^{\text{res}} = -R T^2 \left(\frac{\partial \ln\phi}{\partial T}\right)_P$ hold within relative tolerances $\le 10^{-6}$.
  3. **Reference State Invariance:** At $(T_{\text{ref}}, P_{\text{ref}})$, total enthalpy and entropy evaluate to $0.000\text{ J/mol}$ within machine precision.
  4. **Domain Robustness:** The program never crashes or returns silent NaNs when fed non-physical inputs, negative temperatures, or low pressures.
- **What the tests DO NOT prove:**
  1. **Real-world accuracy for all fluids:** The tests do not prove that PR-1976 predicts real fluid behavior with zero error. PR is a two-parameter cubic model; liquid volumes can deviate by $5\text{--}15\%$ from experimental values, and polar/hydrogen-bonding fluids (like liquid water) require advanced association models (such as CPA or SAFT).
  2. **Multi-component mixture VLE:** The tests only validate pure-component calculations (mixture mixing rules $a_{\text{mix}} = \sum \sum x_i x_j a_{ij}$ were excluded by design from CHL2002 scope).

---

## 3. Code Walkthrough & Equation-to-Code Mapping

### Thermodynamic Equations to Code Architecture Mapping Table

| Thermodynamic Quantity / Equation | Mathematical Expression | C++ Source File & Method |
| :--- | :--- | :--- |
| **Universal Gas Constant** | $R = 8.31446261815324\text{ J/(mol}\cdot\text{K)}$ | [`Constants.hpp:16`](file:///c:/Workspace/Projects/Thermo%20Software/include/thermo/Constants.hpp#L16) (`constants::R`) |
| **PR Parameters ($a_c, b, \kappa$)** | $a_c = \Omega_a \frac{R^2 T_c^2}{P_c}, b = \Omega_b \frac{R T_c}{P_c}$ | [`PengRobinsonEOS.cpp:32-38`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L32-L38) (`PengRobinsonEOS::PengRobinsonEOS`) |
| **Reduced Temperature & Alpha** | $T_r = T/T_c, \alpha = [1 + \kappa(1 - \sqrt{T_r})]^2$ | [`PengRobinsonEOS.cpp:58-75`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L58-L75) (`Tr()`, `alpha()`, `a()`) |
| **Analytical Temperature Derivative** | $\frac{da}{dT} = -\frac{a_c \kappa \sqrt{\alpha}}{\sqrt{T \cdot T_c}}$ | [`PengRobinsonEOS.cpp:92-99`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L92-L99) (`da_dT()`) |
| **Inverse Pressure** | $P = \frac{RT}{v-b} - \frac{a(T)}{v^2 + 2bv - b^2}$ | [`PengRobinsonEOS.cpp:115-123`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L115-L123) (`calculate_pressure()`) |
| **Mechanical Stability Derivative** | $\left(\frac{\partial P}{\partial v}\right)_T = -\frac{RT}{(v-b)^2} + \frac{a(2v+2b)}{(v^2+2bv-b^2)^2}$ | [`PengRobinsonEOS.cpp:125-139`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L125-L139) (`calculate_dP_dv()`) |
| **Cubic Equation in $Z$** | $Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0$ | [`PengRobinsonEOS.cpp:177-187`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L177-L187) (`calculate_pvt()`) |
| **Cubic Polynomial Solver** | Depressed cubic reduction, Cardano & Vieta | [`CubicSolver.cpp:18-85`](file:///c:/Workspace/Projects/Thermo%20Software/src/CubicSolver.cpp#L18-L85) (`CubicSolver::solve()`) |
| **Physical Root Filtering** | $Z > B + 10^{-12} \iff v > b$ | [`PengRobinsonEOS.cpp:196-198`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L196-L198) (`calculate_pvt()`) |
| **Taylor Expansion of Departure Factor** | $u = \frac{2\sqrt{2}B}{Z+(1-\sqrt{2})B}, \frac{\ln(1+u)}{u} \approx 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4}$ | [`PengRobinsonEOS.cpp:306-324`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L306-L324) (`calculate_I_L()`) |
| **Residual Enthalpy** | $H^{\text{res}} = RT(Z-1) + \frac{T \frac{da}{dT} - a}{2\sqrt{2}b} L$ | [`PengRobinsonEOS.cpp:333-353`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L333-L353) (`calculate_H_res()`) |
| **Residual Entropy** | $S^{\text{res}} = R\ln(Z-B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} L$ | [`PengRobinsonEOS.cpp:360-377`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L360-L377) (`calculate_S_res()`) |
| **Fugacity Coefficient & Fugacity** | $\ln\phi = (Z-1) - \ln(Z-B) - \frac{A}{2\sqrt{2}B} L, f = \phi P$ | [`PengRobinsonEOS.cpp:385-424`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L385-L424) (`calculate_ln_phi()`, `calculate_fugacity()`) |
| **Phase Preference & $\Delta G^{\text{res}}$** | $\Delta G^{\text{res}} = RT\ln(\phi_{\text{liq}} / \phi_{\text{vap}})$ | [`PengRobinsonEOS.cpp:468-522`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L468-L522) (`evaluate_phase_preference()`) |
| **Shomate Ideal-Gas Integrals** | $\Delta H^{\text{ig}} = \int C_p^\circ dT, \Delta S^{\text{ig}} = \int \frac{C_p^\circ}{T} dT$ | [`Component.cpp:78-124`](file:///c:/Workspace/Projects/Thermo%20Software/src/Component.cpp#L78-L124) (`integrate_cp()`, `integrate_cp_over_T()`) |
| **Reference-State Total $H$ and $S$** | $H = \Delta H^{\text{ig}} + H^{\text{res}} - H^{\text{res}}_{\text{ref}}$ | [`PengRobinsonEOS.cpp:428-464`](file:///c:/Workspace/Projects/Thermo%20Software/src/PengRobinsonEOS.cpp#L428-L464) (`calculate_H_total()`, `calculate_S_total()`) |
| **Interactive CLI Controller** | Menu loop, input prompts, formatted output | [`CLIController.cpp:218-315`](file:///c:/Workspace/Projects/Thermo%20Software/src/cli/CLIController.cpp#L218-L315) (`handle_pvt_calculation()`) |

---

### Step-by-Step Execution Trace of a Calculation

When a user runs the CLI binary `./thermo_cli` and enters a calculation:

```
[1] User selects Option 1: "PVT Calculation (All Properties)"
    --> CLIController::run() dispatches to CLIController::handle_pvt_calculation()
```

```
[2] Fluid Selection & Validation:
    --> User selects fluid (e.g., Propane, Index 3).
    --> CLIController::prompt_component() queries ComponentDatabase::get_by_index(2).
    --> Validates presence of critical properties (Tc = 369.83 K, Pc = 4.248e6 Pa, omega = 0.152).
```

```
[3] Temperature & Pressure Prompting:
    --> CLIController::prompt_temperature() prompts for T (e.g., 320.0 K).
    --> CLIController::prompt_pressure() prompts for P (e.g., 1000000.0 Pa).
    --> Domain check: strictly requires T > 0, P > 0, finite IEEE 754 doubles.
```

```
[4] Model Initialization:
    --> Constructs thermo::PengRobinsonEOS eos(component).
    --> Computes a_c = 0.4572355289 * R^2 * Tc^2 / Pc, b = 0.0777960739 * R * Tc / Pc.
    --> Computes kappa = 0.37464 + 1.54226*omega - 0.26992*omega^2.
    --> Calls init_reference_departures(): solves EOS at T_ref = 298.15 K, P_ref = 100,000 Pa
        and caches (H_res_ref, S_res_ref).
```

```
[5] Forward PVT & Cubic Solving:
    --> Calls eos.calculate_pvt(320.0, 1.0e6).
    --> Computes reduced parameters A = a(T)*P / (RT)^2 and B = b*P / (RT).
    --> Calls CubicSolver::solve(1.0, -(1 - B), A - 3B^2 - 2B, -(AB - B^2 - B^3)).
    --> Solver reduces equation to depressed cubic, computes discriminant Delta.
        For Propane at 320 K, 1 MPa: Delta < 0 --> invokes Vieta trigonometric solver.
        Yields 3 real roots: Z_0 = 0.036453, Z_1 = 0.091826, Z_2 = 0.850555.
    --> Filters unphysical roots (Z > B + 1e-12). All 3 roots satisfy Z > B.
```

```
[6] Root Characterization & Thermodynamic Departure Evaluation:
    --> For each physical root:
        * Computes molar volume: v = Z * R * T / P.
        * Evaluates isothermal derivative (dP/dv)_T via calculate_dP_dv(T, v).
          - Root 0: (dP/dv)_T < 0 --> STABLE (Liquid-like)
          - Root 1: (dP/dv)_T > 0 --> UNSTABLE (Intermediate spinodal branch)
          - Root 2: (dP/dv)_T < 0 --> STABLE (Vapor-like)
        * Evaluates departure functions:
          - calculate_I_L(Z, B) computes log factor with small-B Taylor safety.
          - calculate_H_res(T, P, Z) computes residual enthalpy.
          - calculate_S_res(T, P, Z) computes residual entropy.
          - calculate_ln_phi(T, P, Z) computes fugacity coefficient.
          - calculate_fugacity() computes f = phi * P.
          - Integrates Shomate polynomials for total H and S relative to (298.15 K, 1 bar).
```

```
[7] Phase Equilibrium & Preference Resolution:
    --> eos.evaluate_phase_preference() compares liquid (Root 0) and vapor (Root 2) roots.
    --> Computes Delta G_res = R * T * ln(phi_liq / phi_vap) = +947.42 J/mol.
    --> Because Delta G_res > 0 (phi_vap = 0.868 < phi_liq = 1.239), the vapor phase has
        lower chemical potential and is thermodynamically preferred.
    --> Flags preference as VAPOR_PREFERRED.
```

```
[8] Display Results:
    --> CLIController::display_pvt_results() formats output into a clear diagnostic table.
    --> Displays Z, molar volume, mechanical stability, fugacity, and total H/S.
    --> Displays thermodynamic phase preference and explanation.
```
