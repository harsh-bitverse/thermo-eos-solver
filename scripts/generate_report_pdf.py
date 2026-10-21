import os
import sys
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether, HRFlowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfgen import canvas

class NumberedCanvas(canvas.Canvas):
    """
    Two-pass canvas to dynamically compute and render exact total page numbers,
    running headers, and running footers.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_decorations(num_pages)
            super().showPage()
        super().save()

    def draw_page_decorations(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#4A5568"))

        # Running header (pages 2+)
        if self._pageNumber > 1:
            self.drawString(36, 762, "CHL2002 Thermodynamics — Technical Report | Peng–Robinson EOS Solver")
            self.setStrokeColor(colors.HexColor("#CBD5E0"))
            self.setLineWidth(0.5)
            self.line(36, 756, 576, 756)

        # Running footer (all pages)
        self.setStrokeColor(colors.HexColor("#CBD5E0"))
        self.setLineWidth(0.5)
        self.line(36, 42, 576, 42)

        footer_text = f"Page {self._pageNumber} of {page_count}"
        self.drawRightString(576, 30, footer_text)
        self.drawString(36, 30, "Confidential — CHL2002 Course Submission — November 2026")
        self.restoreState()


def build_pdf(filename="docs/FINAL_PROJECT_REPORT.pdf"):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    doc = SimpleDocTemplate(
        filename,
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=40,
        bottomMargin=46
    )

    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=15,
        leading=18,
        textColor=colors.HexColor("#1A365D"),
        spaceAfter=4,
        alignment=1 # Center
    )

    subtitle_style = ParagraphStyle(
        'DocSubTitle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=9.5,
        leading=12,
        textColor=colors.HexColor("#2B6CB0"),
        spaceAfter=6,
        alignment=1
    )

    meta_style = ParagraphStyle(
        'DocMeta',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.5,
        leading=10,
        textColor=colors.HexColor("#4A5568"),
        alignment=1
    )

    h1_style = ParagraphStyle(
        'SectionH1',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10.5,
        leading=13,
        textColor=colors.HexColor("#1A365D"),
        spaceBefore=6,
        spaceAfter=3,
        keepWithNext=True
    )

    h2_style = ParagraphStyle(
        'SectionH2',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8.5,
        leading=11,
        textColor=colors.HexColor("#2B6CB0"),
        spaceBefore=4,
        spaceAfter=2,
        keepWithNext=True
    )

    body_style = ParagraphStyle(
        'BodyDark',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.8,
        leading=10.2,
        textColor=colors.HexColor("#2D3748"),
        spaceAfter=3,
        alignment=4 # Justified
    )

    bullet_style = ParagraphStyle(
        'BulletText',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.5,
        leading=9.8,
        textColor=colors.HexColor("#2D3748"),
        leftIndent=10,
        spaceAfter=2
    )

    eq_style = ParagraphStyle(
        'EquationText',
        parent=styles['Normal'],
        fontName='Helvetica-Oblique',
        fontSize=7.5,
        leading=9.8,
        textColor=colors.HexColor("#1A202C"),
        leftIndent=15,
        spaceBefore=1,
        spaceAfter=2
    )

    table_cell_style = ParagraphStyle(
        'TableCell',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=6.8,
        leading=8.6,
        textColor=colors.HexColor("#1A202C")
    )

    table_header_style = ParagraphStyle(
        'TableHeader',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=7.0,
        leading=9.0,
        textColor=colors.HexColor("#FFFFFF"),
        alignment=1
    )

    ref_style = ParagraphStyle(
        'RefText',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=6.8,
        leading=8.8,
        textColor=colors.HexColor("#2D3748"),
        leftIndent=12,
        firstLineIndent=-12,
        spaceAfter=2
    )

    story = []

    # =========================================================================
    # PAGE 1: TITLE, OBJECTIVES, AND MATHEMATICAL FOUNDATIONS
    # =========================================================================
    story.append(Paragraph("CHL2002 THERMODYNAMICS COURSE PROJECT", subtitle_style))
    story.append(Paragraph("Prediction of Real-Gas PVT Behaviour, Thermodynamic Departure Functions, and Pure-Component Fugacity Using the Peng–Robinson Equation of State", title_style))
    story.append(Paragraph("<b>Author:</b> Harsh Agarwal &nbsp;|&nbsp; <b>Language:</b> C++17 &nbsp;|&nbsp; <b>Build:</b> CMake 3.20+ (MinGW GCC 13.2) &nbsp;|&nbsp; <b>Date:</b> November 2026", meta_style))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#1A365D"), spaceBefore=4, spaceAfter=5))

    story.append(Paragraph("1. Introduction, Objectives, and Physical Scope", h1_style))
    story.append(Paragraph(
        "Accurate prediction of phase equilibria, volumetric behaviour, and thermodynamic properties of real fluids is a cornerstone "
        "of chemical process design. Ideal-gas models fail drastically at elevated pressures and near the vapour–liquid coexistence envelope. "
        "The objective of this project is to develop a scientifically rigorous, fully audited C++17 software application implementing the "
        "canonical Peng–Robinson (PR-1976) cubic equation of state for pure fluids. The software is engineered to predict: "
        "(1) Forward and inverse pressure–volume–temperature (PVT) behaviour; "
        "(2) Analytical compressibility factor roots with complete mechanical stability classification; "
        "(3) Residual enthalpy (<i>H</i><sup>res</sup>) and residual entropy (<i>S</i><sup>res</sup>); "
        "(4) Pure-component fugacity coefficients (ln &phi;, &phi;) and fugacities (<i>f</i>); "
        "(5) Reference-state total enthalpy and entropy referenced to the ideal gas at 298.15 K and 1.0 bar; and "
        "(6) Thermodynamic phase preference and saturation coexistence via rigorous equality of chemical potentials.",
        body_style
    ))
    story.append(Paragraph(
        "<b>Physical Scope:</b> The software strictly treats pure chemical components. Curated datasets are provided for eight industrial "
        "fluids: Methane (CH<sub>4</sub>), Ethane (C<sub>2</sub>H<sub>6</sub>), Propane (C<sub>3</sub>H<sub>8</sub>), n-Butane (C<sub>4</sub>H<sub>10</sub>), "
        "Nitrogen (N<sub>2</sub>), Oxygen (O<sub>2</sub>), Carbon Dioxide (CO<sub>2</sub>), and Water Vapor (H<sub>2</sub>O). All calculations "
        "are executed internally in strict SI units (Pa, K, m<sup>3</sup>/mol, J/mol, J/(mol·K)).",
        body_style
    ))

    story.append(Paragraph("2. Governing Mathematical Model and Thermodynamic Equations", h1_style))
    story.append(Paragraph("2.1 The Peng–Robinson Cubic Equation of State (PR-1976)", h2_style))
    story.append(Paragraph(
        "The canonical pressure-explicit formulation proposed by Peng and Robinson (1976) expresses pressure <i>P</i> as the sum of a "
        "repulsive hard-sphere term and an attractive intermolecular dispersion term:",
        body_style
    ))
    story.append(Paragraph("<i>P</i>(<i>T</i>, <i>v</i>) = [<i>R</i><i>T</i> / (<i>v</i> &minus; <i>b</i>)] &minus; [<i>a</i>(<i>T</i>) / (<i>v</i><sup>2</sup> + 2<i>b</i><i>v</i> &minus; <i>b</i><sup>2</sup>)]", eq_style))
    story.append(Paragraph(
        "where <i>R</i> = 8.314462618 J/(mol·K) is the universal gas constant (CODATA 2018), <i>v</i> is molar volume [m<sup>3</sup>/mol], "
        "<i>b</i> is the co-volume parameter representing finite molecular volume, and <i>a</i>(<i>T</i>) is the temperature-dependent attractive parameter. "
        "Enforcing the critical point inflection conditions (&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i><sub><i>c</i></sub></sub> = 0 and "
        "(&part;<sup>2</sup><i>P</i>/&part;<i>v</i><sup>2</sup>)<sub><i>T</i><sub><i>c</i></sub></sub> = 0 yields the universal critical parameters:",
        body_style
    ))
    story.append(Paragraph(
        "<i>a</i><sub><i>c</i></sub> = &Omega;<sub><i>a</i></sub> (<i>R</i><sup>2</sup> <i>T</i><sub><i>c</i></sub><sup>2</sup> / <i>P</i><sub><i>c</i></sub>), &nbsp;&nbsp; "
        "<i>b</i> = &Omega;<sub><i>b</i></sub> (<i>R</i> <i>T</i><sub><i>c</i></sub> / <i>P</i><sub><i>c</i></sub>), &nbsp;&nbsp; "
        "where &Omega;<sub><i>a</i></sub> &approx; 0.4572355289, &nbsp;&Omega;<sub><i>b</i></sub> &approx; 0.0777960739",
        eq_style
    ))
    story.append(Paragraph("2.2 Temperature Dependence and Analytical Derivative", h2_style))
    story.append(Paragraph(
        "The attractive parameter accounts for temperature effects through the dimensionless scaling function &alpha;(<i>T</i><sub><i>r</i></sub>, &omega;):",
        body_style
    ))
    story.append(Paragraph(
        "<i>a</i>(<i>T</i>) = <i>a</i><sub><i>c</i></sub> &alpha;(<i>T</i>), &nbsp;&nbsp; "
        "&alpha;(<i>T</i>) = [ 1 + &kappa; (1 &minus; (<i>T</i> / <i>T</i><sub><i>c</i></sub>)<sup>1/2</sup>) ]<sup>2</sup>, &nbsp;&nbsp; "
        "&kappa; = 0.37464 + 1.54226 &omega; &minus; 0.26992 &omega;<sup>2</sup>",
        eq_style
    ))
    story.append(Paragraph(
        "where <i>T</i><sub><i>r</i></sub> = <i>T</i> / <i>T</i><sub><i>c</i></sub> is the reduced temperature and &omega; is the Pitzer acentric factor. "
        "The exact analytical temperature derivative (<i>da</i>/<i>dT</i>), vital for departure functions, is given without empirical approximation by:",
        body_style
    ))
    story.append(Paragraph("<i>da</i>/<i>dT</i> = &minus; [ <i>a</i><sub><i>c</i></sub> &kappa; &alpha;<sup>1/2</sup>(<i>T</i>) ] / [ (<i>T</i> <i>T</i><sub><i>c</i></sub>)<sup>1/2</sup> ] &nbsp;&nbsp; [Pa·m<sup>6</sup>/(mol<sup>2</sup>·K)]", eq_style))

    story.append(Paragraph("2.3 Dimensionless Compressibility Polynomial", h2_style))
    story.append(Paragraph(
        "Defining the compressibility factor <i>Z</i> &equiv; <i>P</i><i>v</i> / (<i>R</i><i>T</i>) and the dimensionless parameters "
        "<i>A</i> &equiv; <i>a</i><i>P</i> / (<i>R</i><i>T</i>)<sup>2</sup> and <i>B</i> &equiv; <i>b</i><i>P</i> / (<i>R</i><i>T</i>), "
        "the Peng–Robinson EOS transforms into a monic cubic polynomial:",
        body_style
    ))
    story.append(Paragraph("<i>Z</i><sup>3</sup> &minus; (1 &minus; <i>B</i>) <i>Z</i><sup>2</sup> + (<i>A</i> &minus; 3<i>B</i><sup>2</sup> &minus; 2<i>B</i>) <i>Z</i> &minus; (<i>A</i><i>B</i> &minus; <i>B</i><sup>2</sup> &minus; <i>B</i><sup>3</sup>) = 0", eq_style))
    story.append(Paragraph(
        "At physical states, only roots satisfying <i>Z</i> &gt; <i>B</i> correspond to positive molar volumes exceeding the co-volume (<i>v</i> &gt; <i>b</i>). "
        "Roots violating this boundary are strictly discarded as mathematical artefacts.",
        body_style
    ))

    story.append(PageBreak())

    # =========================================================================
    # PAGE 2: DEPARTURE FUNCTIONS, FUGACITY, AND REFERENCE STATES
    # =========================================================================
    story.append(Paragraph("2.4 Residual Thermodynamic Departure Functions", h2_style))
    story.append(Paragraph(
        "Residual departure properties describe the difference between the actual thermodynamic state of a fluid and an ideal gas at the same temperature and pressure: "
        "<i>M</i><sup>res</sup>(<i>T</i>, <i>P</i>) &equiv; <i>M</i>(<i>T</i>, <i>P</i>) &minus; <i>M</i><sup>ig</sup>(<i>T</i>, <i>P</i>). "
        "Using the fundamental Maxwell relations and integrating along an isothermal path from infinite volume to molar volume <i>v</i>:",
        body_style
    ))
    story.append(Paragraph(
        "<i>H</i><sup>res</sup>(<i>T</i>, <i>P</i>, <i>Z</i>) = <i>R</i><i>T</i> (<i>Z</i> &minus; 1) + "
        "[ (<i>T</i> <i>da</i>/<i>dT</i> &minus; <i>a</i>) / (2&radic;2 <i>b</i>) ] &middot; ln [ (<i>Z</i> + (1 + &radic;2)<i>B</i>) / (<i>Z</i> + (1 &minus; &radic;2)<i>B</i>) ] &nbsp;&nbsp; [J/mol]",
        eq_style
    ))
    story.append(Paragraph(
        "<i>S</i><sup>res</sup>(<i>T</i>, <i>P</i>, <i>Z</i>) = <i>R</i> ln(<i>Z</i> &minus; <i>B</i>) + "
        "[ (<i>da</i>/<i>dT</i>) / (2&radic;2 <i>b</i>) ] &middot; ln [ (<i>Z</i> + (1 + &radic;2)<i>B</i>) / (<i>Z</i> + (1 &minus; &radic;2)<i>B</i>) ] &nbsp;&nbsp; [J/(mol·K)]",
        eq_style
    ))
    story.append(Paragraph(
        "These closed-form equations (Peng & Robinson 1976, Eqs. 17–18) satisfy the rigorous thermodynamic identity "
        "<i>G</i><sup>res</sup> = <i>H</i><sup>res</sup> &minus; <i>T</i> <i>S</i><sup>res</sup> = <i>R</i><i>T</i> ln &phi; identically.",
        body_style
    ))

    story.append(Paragraph("2.5 Pure-Component Fugacity Coefficient and Fugacity", h2_style))
    story.append(Paragraph(
        "The fugacity coefficient &phi; &equiv; <i>f</i> / <i>P</i> represents the dimensionless activity of the fluid relative to an ideal gas. "
        "From the fundamental thermodynamic relation ln &phi; = &int;<sub>0</sub><sup><i>P</i></sup> [(<i>Z</i> &minus; 1) / <i>P</i>] <i>dP</i>, the exact PR expression is:",
        body_style
    ))
    story.append(Paragraph(
        "ln &phi; = (<i>Z</i> &minus; 1) &minus; ln(<i>Z</i> &minus; <i>B</i>) &minus; "
        "[ <i>A</i> / (2&radic;2 <i>B</i>) ] &middot; ln [ (<i>Z</i> + (1 + &radic;2)<i>B</i>) / (<i>Z</i> + (1 &minus; &radic;2)<i>B</i>) ], &nbsp;&nbsp; "
        "<i>f</i> = &phi; &middot; <i>P</i> &nbsp;&nbsp; [Pa]",
        eq_style
    ))
    story.append(Paragraph(
        "In the ideal gas limit (<i>P</i> &rarr; 0), <i>Z</i> &rarr; 1, <i>B</i> &rarr; 0, ln &phi; &rarr; 0, and &phi; &rarr; 1.0 identically.",
        body_style
    ))

    story.append(Paragraph("2.6 NIST Shomate Heat Capacity and Reference-State Total H and S", h2_style))
    story.append(Paragraph(
        "To calculate absolute state enthalpy <i>H</i>(<i>T</i>, <i>P</i>) and entropy <i>S</i>(<i>T</i>, <i>P</i>), the program integrates "
        "the authentic 5-parameter NIST Shomate isobaric heat capacity correlation for the ideal gas state (<i>t</i> = <i>T</i> / 1000):",
        body_style
    ))
    story.append(Paragraph("<i>C</i><sub><i>p</i></sub><sup>&deg;</sup>(<i>t</i>) = <i>A</i> + <i>B</i><i>t</i> + <i>C</i><i>t</i><sup>2</sup> + <i>D</i><i>t</i><sup>3</sup> + <i>E</i> / <i>t</i><sup>2</sup> &nbsp;&nbsp; [J/(mol·K)]", eq_style))
    story.append(Paragraph(
        "Exact analytical integration yields the temperature departures &Delta;<i>H</i><sup>ig</sup> and &Delta;<i>S</i><sup>ig</sup>:",
        body_style
    ))
    story.append(Paragraph(
        "&Delta;<i>H</i><sup>ig</sup>(<i>T</i><sub>ref</sub> &rarr; <i>T</i>) = 1000 [ <i>A</i><i>t</i> + <i>B</i><i>t</i><sup>2</sup>/2 + <i>C</i><i>t</i><sup>3</sup>/3 + <i>D</i><i>t</i><sup>4</sup>/4 &minus; <i>E</i>/<i>t</i> ]<sub><i>t</i><sub>ref</sub></sub><sup><i>t</i></sup> &nbsp;&nbsp; [J/mol]",
        eq_style
    ))
    story.append(Paragraph(
        "&Delta;<i>S</i><sup>ig</sup>(<i>T</i><sub>ref</sub> &rarr; <i>T</i>) = [ <i>A</i> ln <i>t</i> + <i>B</i><i>t</i> + <i>C</i><i>t</i><sup>2</sup>/2 + <i>D</i><i>t</i><sup>3</sup>/3 &minus; <i>E</i> / (2<i>t</i><sup>2</sup>) ]<sub><i>t</i><sub>ref</sub></sub><sup><i>t</i></sup> &nbsp;&nbsp; [J/(mol·K)]",
        eq_style
    ))
    story.append(Paragraph(
        "<b>Reference State Convention:</b> Total enthalpy and entropy are defined relative to the ideal gas reference state at "
        "<i>T</i><sub>ref</sub> = 298.15 K (25.0 &deg;C) and <i>P</i><sub>ref</sub> = 100,000 Pa (1.0 bar), where "
        "<i>H</i><sub>total</sub>(<i>T</i><sub>ref</sub>, <i>P</i><sub>ref</sub>) &equiv; 0.0 J/mol and "
        "<i>S</i><sub>total</sub>(<i>T</i><sub>ref</sub>, <i>P</i><sub>ref</sub>) &equiv; 0.0 J/(mol·K). The total properties at arbitrary (<i>T</i>, <i>P</i>, <i>Z</i>) are:",
        body_style
    ))
    story.append(Paragraph(
        "<i>H</i>(<i>T</i>, <i>P</i>, <i>Z</i>) = &Delta;<i>H</i><sup>ig</sup>(<i>T</i><sub>ref</sub> &rarr; <i>T</i>) + "
        "<i>H</i><sup>res</sup>(<i>T</i>, <i>P</i>, <i>Z</i>) &minus; <i>H</i><sup>res</sup>(<i>T</i><sub>ref</sub>, <i>P</i><sub>ref</sub>, <i>Z</i><sub>ref</sub>)",
        eq_style
    ))
    story.append(Paragraph(
        "<i>S</i>(<i>T</i>, <i>P</i>, <i>Z</i>) = &Delta;<i>S</i><sup>ig</sup>(<i>T</i><sub>ref</sub> &rarr; <i>T</i>) &minus; <i>R</i> ln(<i>P</i> / <i>P</i><sub>ref</sub>) + "
        "<i>S</i><sup>res</sup>(<i>T</i>, <i>P</i>, <i>Z</i>) &minus; <i>S</i><sup>res</sup>(<i>T</i><sub>ref</sub>, <i>P</i><sub>ref</sub>, <i>Z</i><sub>ref</sub>)",
        eq_style
    ))

    story.append(Paragraph("2.7 Low-Pressure Taylor Series Regularization", h2_style))
    story.append(Paragraph(
        "At near-vacuum pressures (<i>P</i> &lt; 10 Pa), <i>B</i> &rarr; 0, causing the quotient ln[(<i>Z</i> + (1+&radic;2)<i>B</i>)/(<i>Z</i> + (1&minus;&radic;2)<i>B</i>)] / (2&radic;2 <i>B</i>) "
        "to approach 0 / 0. To prevent catastrophic floating-point cancellation, the software switches to a 4-term analytical Taylor series when <i>B</i>/<i>Z</i> &lt; 10<sup>&minus;5</sup>:",
        body_style
    ))
    story.append(Paragraph("<i>I</i><sub><i>L</i></sub> = (1 / <i>Z</i>) &middot; [ 1 &minus; (<i>B</i>/<i>Z</i>) + (5/3)(<i>B</i>/<i>Z</i>)<sup>2</sup> &minus; 3(<i>B</i>/<i>Z</i>)<sup>3</sup> + (33/5)(<i>B</i>/<i>Z</i>)<sup>4</sup> ]", eq_style))
    story.append(Paragraph("This maintains double-precision continuity to 10<sup>&minus;16</sup> down to <i>P</i> = 10<sup>&minus;5</sup> Pa.", body_style))

    story.append(PageBreak())

    # =========================================================================
    # PAGE 3: NUMERICAL METHODS, STABILITY, AND ARCHITECTURE
    # =========================================================================
    story.append(Paragraph("2.8 Analytical Cubic Root Solver and Stability Classification", h2_style))
    story.append(Paragraph(
        "The cubic equation in <i>Z</i> is solved analytically using Cardano's formulation transformed into depressed form "
        "<i>y</i><sup>3</sup> + <i>p</i><i>y</i> + <i>q</i> = 0 via the Tschirnhaus substitution <i>Z</i> = <i>y</i> &minus; <i>c</i><sub>2</sub>/3. "
        "The discriminant &Delta; = <i>q</i><sup>2</sup>/4 + <i>p</i><sup>3</sup>/27 governs the root structure: "
        "(1) If &Delta; &gt; 0, one real root exists (Cardano's real radical formula); "
        "(2) If &Delta; &lt; 0, three distinct real roots exist (Vieta's trigonometric formulation); and "
        "(3) If &Delta; = 0, multiple real roots coalesce.",
        body_style
    ))
    story.append(Paragraph(
        "<b>Mechanical Stability Criterion:</b> For every physical root (<i>Z</i> &gt; <i>B</i>), the isothermal pressure–volume derivative "
        "(&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub> is evaluated analytically: "
        "(&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub> = &minus;<i>R</i><i>T</i>/(<i>v</i> &minus; <i>b</i>)<sup>2</sup> + <i>a</i>(<i>T</i>)(2<i>v</i> + 2<i>b</i>)/(<i>v</i><sup>2</sup> + 2<i>b</i><i>v</i> &minus; <i>b</i><sup>2</sup>)<sup>2</sup>. "
        "States with (&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub> &lt; 0 are <b>STABLE</b> (compression increases pressure). "
        "States with (&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub> &ge; 0 are <b>UNSTABLE</b> (spinodal states; physically unobservable).",
        body_style
    ))
    story.append(Paragraph(
        "<b>Phase Preference Criterion:</b> In the three-root subcritical region, the smallest root (<i>Z</i><sub>min</sub>) is liquid-like, "
        "the largest root (<i>Z</i><sub>max</sub>) is vapour-like, and the intermediate root is mechanically unstable. "
        "Phase stability is decided by the minimum Gibbs free energy, which corresponds to the lower fugacity coefficient &phi;: "
        "&Delta;<i>G</i><sup>res</sup> = <i>G</i><sub>liq</sub><sup>res</sup> &minus; <i>G</i><sub>vap</sub><sup>res</sup> = <i>R</i><i>T</i> ln(&phi;<sub>liq</sub> / &phi;<sub>vap</sub>). "
        "If &phi;<sub>vap</sub> &lt; &phi;<sub>liq</sub> (&Delta;<i>G</i><sup>res</sup> &gt; 0), the vapour phase is preferred. "
        "If &phi;<sub>liq</sub> &lt; &phi;<sub>vap</sub> (&Delta;<i>G</i><sup>res</sup> &lt; 0), the liquid phase is preferred. "
        "If &phi;<sub>liq</sub> &approx; &phi;<sub>vap</sub> within 10<sup>&minus;5</sup> relative tolerance, saturated two-phase equilibrium exists.",
        body_style
    ))

    story.append(Paragraph("3. Software Architecture and Implementation", h1_style))
    story.append(Paragraph(
        "The software is designed as a modular, header/source C++17 library (<code>thermo_core</code>) decoupled from user-interface logic. "
        "It relies exclusively on the C++ standard library (<code>&lt;cmath&gt;</code>, <code>&lt;vector&gt;</code>, <code>&lt;string&gt;</code>, <code>&lt;optional&gt;</code>, <code>&lt;iostream&gt;</code>) with zero external dependencies.",
        body_style
    ))

    # Architecture Table
    arch_data = [
        [Paragraph("Module / Component", table_header_style), Paragraph("Source Files", table_header_style), Paragraph("Architectural Responsibility & Function", table_header_style)],
        [Paragraph("<b>Constants</b>", table_cell_style), Paragraph("include/thermo/Constants.hpp", table_cell_style), Paragraph("CODATA universal gas constant <i>R</i>, standard reference state constants, and analytical PR-1976 parameters.", table_cell_style)],
        [Paragraph("<b>Component & Database</b>", table_cell_style), Paragraph("src/Component.cpp<br/>src/ComponentDatabase.cpp", table_cell_style), Paragraph("Validates physical parameters; stores curated database for 8 industrial gases; evaluates analytical Shomate integrals.", table_cell_style)],
        [Paragraph("<b>Cubic Solver</b>", table_cell_style), Paragraph("src/CubicSolver.cpp", table_cell_style), Paragraph("Analytical Cardano/Vieta polynomial solver; filters duplicate roots; guarantees residuals &le; 10<sup>&minus;12</sup>.", table_cell_style)],
        [Paragraph("<b>PR EOS Engine</b>", table_cell_style), Paragraph("src/PengRobinsonEOS.cpp", table_cell_style), Paragraph("Forward/inverse PVT, (&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub>, departure functions <i>H</i><sup>res</sup>, <i>S</i><sup>res</sup>, fugacity &phi;, and phase preference.", table_cell_style)],
        [Paragraph("<b>Interactive CLI</b>", table_cell_style), Paragraph("src/cli/CLIController.cpp<br/>src/cli_main.cpp", table_cell_style), Paragraph("Stream-decoupled interactive menu (`thermo_cli`), unit conversions (K/&deg;C, Pa/bar/MPa), benchmark demo runner.", table_cell_style)],
        [Paragraph("<b>Automated Test Suite</b>", table_cell_style), Paragraph("tests/test_*.cpp (5 targets)", table_cell_style), Paragraph("373 automated assertions covering cubic roots, Shomate quadrature, PVT reciprocity, Maxwell derivatives, and CLI.", table_cell_style)]
    ]
    arch_table = Table(arch_data, colWidths=[105, 135, 300])
    arch_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#1A365D")),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ('TOPPADDING', (0, 0), (-1, -1), 3),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E0")),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor("#F7FAFC"), colors.HexColor("#FFFFFF")])
    ]))
    story.append(arch_table)
    story.append(Spacer(1, 4))

    story.append(Paragraph("3.1 Build System and Execution Instructions", h2_style))
    story.append(Paragraph(
        "The project is built using CMake 3.20+ with strict compiler warnings (`-Wall -Wextra -Wpedantic`):",
        body_style
    ))
    story.append(Paragraph("<code>cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build</code>", eq_style))
    story.append(Paragraph("To run the full automated CTest suite (all 5 targets, 373 assertions):", body_style))
    story.append(Paragraph("<code>ctest --test-dir build --output-on-failure</code>", eq_style))
    story.append(Paragraph("To run the interactive CLI or automated demonstration benchmarks:", body_style))
    story.append(Paragraph("<code>./build/thermo_cli &nbsp;&nbsp;# Interactive mode<br/>./build/thermo_cli --demo all &nbsp;&nbsp;# Automated demonstration benchmarks</code>", eq_style))

    story.append(PageBreak())

    # =========================================================================
    # PAGE 4: RESULTS, BENCHMARKS, LIMITATIONS, AND REFERENCES
    # =========================================================================
    story.append(Paragraph("4. Verification Results and Audited Numerical Benchmarks", h1_style))
    story.append(Paragraph(
        "All thermodynamic and volumetric calculations were independently verified against literature benchmarks and internal mathematical consistency checks. "
        "The table below reports values confirmed by direct execution of `thermo_cli`:",
        body_style
    ))

    # Benchmark Table
    bench_data = [
        [Paragraph("Case / State", table_header_style), Paragraph("Fluid & Phase", table_header_style), Paragraph("<i>Z</i> (PR-1976)", table_header_style), Paragraph("<i>v</i> [m<sup>3</sup>/mol]", table_header_style), Paragraph("&phi; (Fugacity Coeff)", table_header_style), Paragraph("<i>H</i><sup>res</sup> [J/mol]", table_header_style), Paragraph("<i>S</i><sup>res</sup> [J/(mol·K)]", table_header_style), Paragraph("Status / Preference", table_header_style)],
        [Paragraph("<b>Demo 1</b><br/>300 K, 5.0 MPa", table_cell_style), Paragraph("Methane<br/>(Gas)", table_cell_style), Paragraph("0.901844", table_cell_style), Paragraph("4.4990&times;10<sup>&minus;4</sup>", table_cell_style), Paragraph("0.901385<br/>(<i>f</i>=4.507 MPa)", table_cell_style), Paragraph("&minus;902.34", table_cell_style), Paragraph("&minus;2.1446", table_cell_style), Paragraph("Single Phase<br/>(Stable)", table_cell_style)],
        [Paragraph("<b>Demo 2</b><br/>300 K, 7.0 MPa", table_cell_style), Paragraph("CO<sub>2</sub><br/>(Dense Fluid)", table_cell_style), Paragraph("0.198057", table_cell_style), Paragraph("7.0575&times;10<sup>&minus;5</sup>", table_cell_style), Paragraph("0.638285<br/>(<i>f</i>=4.468 MPa)", table_cell_style), Paragraph("&minus;9374.37", table_cell_style), Paragraph("&minus;27.5150", table_cell_style), Paragraph("Single Phase<br/>(Stable Liquid-like)", table_cell_style)],
        [Paragraph("<b>Demo 3 (Root 1)</b><br/>320 K, 1.0 MPa", table_cell_style), Paragraph("Propane<br/>(Liquid-like)", table_cell_style), Paragraph("0.036453", table_cell_style), Paragraph("9.6988&times;10<sup>&minus;5</sup>", table_cell_style), Paragraph("1.239280<br/>(<i>f</i>=1.239 MPa)", table_cell_style), Paragraph("&minus;14,854.34", table_cell_style), Paragraph("&minus;48.2035", table_cell_style), Paragraph("Liquid Branch<br/>(Metastable)", table_cell_style)],
        [Paragraph("<b>Demo 3 (Root 2)</b><br/>320 K, 1.0 MPa", table_cell_style), Paragraph("Propane<br/>(Intermediate)", table_cell_style), Paragraph("0.091826", table_cell_style), Paragraph("2.4432&times;10<sup>&minus;4</sup>", table_cell_style), Paragraph("1.385615<br/>(<i>f</i>=1.386 MPa)", table_cell_style), Paragraph("&minus;8,210.60", table_cell_style), Paragraph("&minus;28.3698", table_cell_style), Paragraph("<b>UNSTABLE</b><br/>(&part;<i>P</i>/&part;<i>v</i> &gt; 0, Spinodal)", table_cell_style)],
        [Paragraph("<b>Demo 3 (Root 3)</b><br/>320 K, 1.0 MPa", table_cell_style), Paragraph("Propane<br/>(Vapour-like)", table_cell_style), Paragraph("0.850555", table_cell_style), Paragraph("2.2630&times;10<sup>&minus;3</sup>", table_cell_style), Paragraph("0.868006<br/>(<i>f</i>=0.868 MPa)", table_cell_style), Paragraph("&minus;1,131.04", table_cell_style), Paragraph("&minus;2.3575", table_cell_style), Paragraph("<b>VAPOUR PREFERRED</b><br/>(&Delta;<i>G</i><sup>res</sup> = +947.42 J/mol)", table_cell_style)]
    ]
    bench_table = Table(bench_data, colWidths=[65, 65, 55, 65, 75, 65, 65, 85])
    bench_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#1A365D")),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
        ('TOPPADDING', (0, 0), (-1, -1), 2),
        ('LEFTPADDING', (0, 0), (-1, -1), 2),
        ('RIGHTPADDING', (0, 0), (-1, -1), 2),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E0")),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor("#F7FAFC"), colors.HexColor("#FFFFFF")])
    ]))
    story.append(bench_table)
    story.append(Spacer(1, 4))

    story.append(Paragraph("4.1 In-Depth Analysis of Multi-Root Propane Behaviour (Demo 3)", h2_style))
    story.append(Paragraph(
        "At 320 K and 1.0 MPa, Propane lies in the subcritical two-phase region (&Delta; = &minus;1.083&times;10<sup>&minus;5</sup> &lt; 0), generating three real roots: "
        "Root 1 (liquid-like, <i>Z</i> = 0.036453), Root 2 (intermediate, <i>Z</i> = 0.091826), and Root 3 (vapour-like, <i>Z</i> = 0.850555). "
        "The intermediate root exhibits (&part;<i>P</i>/&part;<i>v</i>)<sub><i>T</i></sub> = +1.8825&times;10<sup>10</sup> Pa/(m<sup>3</sup>/mol) &gt; 0, "
        "confirming it is mechanically unstable and located inside the spinodal zone. For the two mechanically stable roots, "
        "&phi;<sub>vap</sub> = 0.868006 &lt; &phi;<sub>liq</sub> = 1.239280, producing &Delta;<i>G</i><sup>res</sup> = <i>R</i><i>T</i> ln(&phi;<sub>liq</sub>/&phi;<sub>vap</sub>) = +947.42 J/mol &gt; 0. "
        "The vapour phase possesses strictly lower chemical potential (&mu;<sub>vap</sub> &lt; &mu;<sub>liq</sub>) and is selected as the stable equilibrium state.",
        body_style
    ))

    story.append(Paragraph("4.2 Automated Testing and Independent Differential Maxwell Cross-Checks", h2_style))
    story.append(Paragraph(
        "The repository contains an automated test suite of <b>373 assertions across 5 CTest executables (100% pass rate)</b>. "
        "Circularity in validation was strictly avoided by implementing true independent differential checks: "
        "(1) <b>Gibbs–Helmholtz Cross-Check:</b> <i>H</i><sup>res</sup> was independently verified against numerical central differences "
        "&minus;<i>R</i><i>T</i><sup>2</sup> (&part; ln &phi;/&part;<i>T</i>)<sub><i>P</i></sub> with &delta;<i>T</i> = 0.1 mK, matching to a relative error of 3.8&times;10<sup>&minus;10</sup>; "
        "(2) <b>Pressure Fugacity Differential:</b> (&part; ln &phi;/&part;<i>P</i>)<sub><i>T</i></sub> was matched to (<i>Z</i> &minus; 1)/<i>P</i> to within 1.4&times;10<sup>&minus;9</sup>; and "
        "(3) <b>Numerical Quadrature:</b> Analytical Shomate integrals for &Delta;<i>H</i><sup>ig</sup> and &Delta;<i>S</i><sup>ig</sup> were cross-checked "
        "against an independent 1000-subinterval composite Simpson's 1/3 quadrature to within 10<sup>&minus;6</sup>.",
        body_style
    ))

    story.append(Paragraph("5. Assumptions, Limitations, and Conclusion", h1_style))
    story.append(Paragraph(
        "<b>Model Scope and Limitations:</b> "
        "(1) <i>Pure Fluids Only:</i> Multicomponent mixtures and mixing rules are out of scope. "
        "(2) <i>Steam Reference State Limitation:</i> NIST Shomate gas-phase parameters for Water Vapor are valid only in [500.0, 1700.0] K. "
        "Because <i>T</i><sub>ref</sub> = 298.15 K lies outside this interval, total <i>H</i> and <i>S</i> are reported as N/A with clear diagnostics; "
        "PVT, fugacity, and residual departure functions remain fully operational. "
        "(3) <i>No Unphysical Extrapolation:</i> Calculating outside documented Shomate intervals throws exceptions rather than generating erroneous values. "
        "(4) <i>Cubic Liquid Volume Accuracy:</i> Near the critical point, two-parameter cubic EOSs typically underestimate liquid densities by 5–10% without volume translation (Peneloux).",
        body_style
    ))
    story.append(Paragraph(
        "<b>Conclusion:</b> The software provides a mathematically rigorous, verified, and robust engineering tool for pure-fluid thermodynamics. "
        "It successfully models complex non-ideal fluid behaviour, multi-root phase preference, and departure functions while maintaining high numerical precision.",
        body_style
    ))

    story.append(Paragraph("References", h2_style))
    story.append(Paragraph("[1] Peng, D.-Y., & Robinson, D. B. (1976). A New Two-Constant Equation of State. <i>Ind. Eng. Chem. Fundam.</i>, 15(1), 59–64.", ref_style))
    story.append(Paragraph("[2] Linstrom, P. J., & Mallard, W. G., Eds. (2023). <i>NIST Chemistry WebBook, NIST Standard Reference Database Number 69</i>, National Institute of Standards and Technology, Gaithersburg MD. https://webbook.nist.gov", ref_style))
    story.append(Paragraph("[3] Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001). <i>The Properties of Gases and Liquids</i>, 5th Ed., McGraw-Hill, New York.", ref_style))

    doc.build(story, canvasmaker=NumberedCanvas)
    print(f"[SUCCESS] Generated PDF report: {filename}")

if __name__ == "__main__":
    out_pdf = "docs/FINAL_PROJECT_REPORT.pdf"
    if len(sys.argv) > 1:
        out_pdf = sys.argv[1]
    build_pdf(out_pdf)
