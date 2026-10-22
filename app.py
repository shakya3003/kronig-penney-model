"""
app.py: Interactive Python Interface for the Kronig-Penney Model

A Streamlit web dashboard for investigating:
- 1D periodic potential V(x) across N wells
- Wavefunction psi(x) and probability density |psi(x)|^2
- Kronig-Penney dispersion relation f(E) and band structures
- Numerical differential equation solution (Finite Difference)

Run locally with:
    streamlit run app.py
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

import kp_engine as kp

# Page configuration and custom CSS styling
st.set_page_config(
    page_title="Kronig-Penney Model Explorer",
    page_icon="KP",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #4B5563;
        margin-bottom: 1.5rem;
    }
    .status-allowed {
        background-color: #DEF7EC;
        color: #03543F;
        padding: 0.75rem 1rem;
        border-radius: 8px;
        font-weight: 600;
        border-left: 5px solid #31C48D;
        margin-bottom: 1rem;
    }
    .status-forbidden {
        background-color: #FDE8E8;
        color: #9B1C1C;
        padding: 0.75rem 1rem;
        border-radius: 8px;
        font-weight: 600;
        border-left: 5px solid #F05252;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)


# Sidebar: user-adjustable model parameters (Question 2 inputs)
st.sidebar.title("Model Parameters")
st.sidebar.markdown("Configure the 1D crystalline superlattice:")

# Quick presets for common physical scenarios
preset = st.sidebar.selectbox(
    "Choose a Physical Preset:",
    [
        "Custom Parameters",
        "Benchmark Semiconductor (Allowed)",
        "Forbidden Gap State",
        "Nearly Free Electron (Weak Barriers)",
        "Tight-Binding (Deep Wells)",
    ]
)

if preset == "Benchmark Semiconductor (Allowed)":
    def_N, def_E, def_a, def_b, def_V0 = 4, 1.20, 0.50, 0.10, 10.00
elif preset == "Forbidden Gap State":
    def_N, def_E, def_a, def_b, def_V0 = 4, 3.20, 0.50, 0.10, 10.00
elif preset == "Nearly Free Electron (Weak Barriers)":
    def_N, def_E, def_a, def_b, def_V0 = 5, 2.50, 0.80, 0.05, 3.00
elif preset == "Tight-Binding (Deep Wells)":
    def_N, def_E, def_a, def_b, def_V0 = 4, 1.50, 0.40, 0.30, 20.00
else:
    def_N, def_E, def_a, def_b, def_V0 = 4, 1.20, 0.50, 0.10, 10.00

# Slider controls for all five physical parameters
N = st.sidebar.slider(
    "Number of Potential Wells (N)",
    min_value=1,
    max_value=12,
    value=def_N,
    step=1,
    help="Number of repeating atomic wells in the finite superlattice."
)

V0 = st.sidebar.slider(
    "Barrier Height V0 (eV)",
    min_value=1.0,
    max_value=30.0,
    value=float(def_V0),
    step=0.5,
    help="Potential energy height of the repulsive barriers."
)

E = st.sidebar.slider(
    "Electron Energy E (eV)",
    min_value=0.1,
    max_value=float(max(25.0, V0 * 1.5)),
    value=float(def_E),
    step=0.05,
    help="Total energy of the incident electron."
)

col_a, col_b = st.sidebar.columns(2)
with col_a:
    a = st.slider(
        "Well Width a (nm)",
        min_value=0.10,
        max_value=2.00,
        value=float(def_a),
        step=0.05,
        help="Spatial width of region where V = 0."
    )
with col_b:
    b = st.slider(
        "Barrier Width b (nm)",
        min_value=0.02,
        max_value=1.00,
        value=float(def_b),
        step=0.02,
        help="Thickness of barrier where V = V0."
    )

d = a + b
st.sidebar.markdown(f"**Lattice Constant (d = a + b):** `{d:.3f} nm`")
st.sidebar.markdown(f"**Quantum Constant $\\frac{{\\hbar^2}}{{2m_e}}$:** `0.0762 eV nm^2`")


# Run the dispersion and wavefunction calculations using the engine
disp_result = kp.compute_dispersion(E, a, b, V0)
f_E = disp_result["f_E"]
is_allowed = disp_result["is_allowed"]
k_val = disp_result["k"]
kappa_val = disp_result["kappa"]

# Compute continuous wavefunction and probability density via TMM
x_arr, psi_arr, prob_arr, V_arr = kp.compute_wavefunction_tmm(N, E, a, b, V0)

# Compute allowed bands for context display
bands, gaps = kp.find_allowed_bands(a, b, V0, E_max=max(20.0, V0 * 1.3))


# Page header and state classification banner
st.markdown('<div class="main-header">Kronig-Penney Model & Wavefunction Explorer</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Interactive investigation of electron wavefunctions and energy band structures in a 1D periodic crystalline lattice.</div>', unsafe_allow_html=True)

# Display whether the current energy falls in an allowed band or forbidden gap
if is_allowed:
    st.markdown(
        f"""
        <div class="status-allowed">
            <b>ALLOWED ELECTRONIC STATE</b> (Inside an Energy Band)<br>
            Dispersion Function: <code>f(E) = {f_E:.4f}</code> (within [-1, 1]) &nbsp;|&nbsp;
            Bloch Wavevector: <code>k = {k_val:.4f} nm^-1</code> &nbsp;|&nbsp;
            Wavelength: <code>lambda = {2 * np.pi / k_val:.3f} nm</code>
        </div>
        """,
        unsafe_allow_html=True
    )
else:
    st.markdown(
        f"""
        <div class="status-forbidden">
            <b>FORBIDDEN BAND GAP</b> (Bragg Reflection / Attenuation)<br>
            Dispersion Function: <code>f(E) = {f_E:.4f}</code> (outside [-1, 1]) &nbsp;|&nbsp;
            Decay Factor: <code>kappa = {kappa_val:.4f} nm^-1</code> &nbsp;|&nbsp;
            Penetration Depth: <code>1/kappa = {1 / kappa_val:.3f} nm</code>
        </div>
        """,
        unsafe_allow_html=True
    )

# Summary metric cards at the top of the page
col1, col2, col3, col4 = st.columns(4)
col1.metric("Selected Energy E", f"{E:.2f} eV")
col2.metric("Dispersion f(E)", f"{f_E:.3f}")
col3.metric("State Classification", "Allowed" if is_allowed else "Forbidden")
col4.metric("Lattice Period d", f"{d:.2f} nm")


# Tabbed layout for organized investigation
tab1, tab2, tab3, tab4 = st.tabs([
    "Wavefunction & Probability Density",
    "Dispersion Relation f(E) & Bands",
    "Differential Equation Solver (TISE)",
    "Theory, Physical Meaning & References"
])


# Tab 1: Potential landscape, wavefunction psi(x), and probability density (Questions 4 & 5)
with tab1:
    st.subheader("Spatial Profiles across the Periodic Lattice")
    st.markdown("Visualizing the periodic potential, the continuous wavefunction, and the probability density.")

    fig1, (ax_v, ax_psi, ax_prob) = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)

    # Plot 1: Periodic potential V(x) with the electron energy level overlaid
    ax_v.plot(x_arr, V_arr, "k-", lw=2, label="Periodic Potential $V(x)$")
    ax_v.axhline(E, color="crimson", linestyle="--", lw=1.8, label=f"Electron Energy $E = {E:.2f}$ eV")
    ax_v.fill_between(x_arr, 0, V_arr, color="lightgray", alpha=0.5, label="Barrier Regions ($V_0$)")
    ax_v.set_ylabel("Energy (eV)", fontsize=11)
    ax_v.set_title(f"1D Periodic Superlattice ($N = {N}$ Wells, $a = {a}$ nm, $b = {b}$ nm, $V_0 = {V0}$ eV)", fontsize=12, fontweight="bold")
    ax_v.legend(loc="upper right", framealpha=0.9)
    ax_v.grid(True, linestyle=":", alpha=0.6)

    # Plot 2: Wavefunction psi(x) computed via Transfer Matrix Method
    ax_psi.plot(x_arr, psi_arr, color="royalblue", lw=2, label=r"Wavefunction $\psi(x)$ (Transfer Matrix Method)")
    ax_psi.axhline(0, color="gray", linestyle=":", alpha=0.7)
    ax_psi.set_ylabel(r"$\psi(x)$ Amplitude", fontsize=11)
    ax_psi.set_title(r"Continuous Wavefunction $\psi(x)$ across Wells & Barriers", fontsize=12, fontweight="bold")
    ax_psi.legend(loc="upper right", framealpha=0.9)
    ax_psi.grid(True, linestyle=":", alpha=0.6)

    # Plot 3: Probability density |psi(x)|^2, normalized to integrate to 1
    ax_prob.plot(x_arr, prob_arr, color="darkmagenta", lw=2, label=r"Probability Density $|\psi(x)|^2$")
    ax_prob.fill_between(x_arr, prob_arr, color="magenta", alpha=0.25)
    ax_prob.set_xlabel("Spatial Position $x$ (nm)", fontsize=11)
    ax_prob.set_ylabel(r"$|\psi(x)|^2$ (nm$^{-1}$)", fontsize=11)
    ax_prob.set_title(r"Probability Density Distribution $|\psi(x)|^2$ (Normalized to 1.0)", fontsize=12, fontweight="bold")
    ax_prob.legend(loc="upper right", framealpha=0.9)
    ax_prob.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    st.pyplot(fig1)

    # Physical explanation of the wavefunction behavior (Question 5 requirement)
    st.info(r"""
    **Physical Significance of the Plots:**
    - **Inside Potential Wells ($V = 0$):** The wavefunction $\psi(x)$ oscillates sinusoidally with spatial frequency $\alpha = \sqrt{2mE}/\hbar$. Peaks in $|\psi(x)|^2$ represent regions where the electron has the highest probability of being located.
    - **Inside Barriers ($V = V_0$):** When $E < V_0$, the electron classically lacks the energy to cross. However, quantum mechanics permits **barrier penetration (tunneling)**; the wavefunction decays exponentially with decay parameter $\beta = \sqrt{2m(V_0 - E)}/\hbar$.
    - **Allowed vs. Forbidden Behavior:** In an allowed state, the wave amplitude remains balanced across the whole periodic crystal. In a forbidden band gap, destructive interference (Bragg reflection) attenuates the wave across successive wells.
    """)


# Tab 2: Dispersion relation f(E) vs E and band structure table (Question 3)
with tab2:
    st.subheader("Kronig-Penney Dispersion Relation & Band Spectrum")

    col_disp1, col_disp2 = st.columns([2, 1])

    with col_disp1:
        # Compute f(E) over a continuous energy range for the dispersion plot
        E_max_plot = max(20.0, V0 * 1.4)
        energies_plot, f_plot, allowed_plot, k_plot = kp.compute_dispersion_curve(a, b, V0, E_min=0.01, E_max=E_max_plot, num_points=1200)

        fig2, ax_disp = plt.subplots(figsize=(8, 5))
        ax_disp.plot(energies_plot, f_plot, "b-", lw=1.8, label=r"$f(E) = \cos(kd)$")
        ax_disp.axhline(1.0, color="red", linestyle="--", lw=1.2, label=r"$f(E) = +1$ Boundary")
        ax_disp.axhline(-1.0, color="red", linestyle="--", lw=1.2, label=r"$f(E) = -1$ Boundary")
        ax_disp.axhspan(-1.0, 1.0, color="lightgreen", alpha=0.25, label=r"Allowed Bands ($|f(E)| \leq 1$)")

        # Mark the currently selected energy on the dispersion curve
        ax_disp.plot(E, f_E, "ro", markersize=8, label=f"Current E = {E:.2f} eV")

        ax_disp.set_ylim(-3.5, 3.5)
        ax_disp.set_xlabel("Electron Energy $E$ (eV)", fontsize=11)
        ax_disp.set_ylabel(r"$f(E)$", fontsize=11)
        ax_disp.set_title("Kronig-Penney Dispersion Curve $f(E)$ vs. $E$", fontsize=12, fontweight="bold")
        ax_disp.legend(loc="upper right", framealpha=0.9, fontsize=9)
        ax_disp.grid(True, linestyle=":", alpha=0.6)

        st.pyplot(fig2)

    with col_disp2:
        # Table of all detected allowed energy bands
        st.markdown("#### Allowed Energy Bands")
        if bands:
            band_data = []
            for b_info in bands:
                band_data.append({
                    "Band #": f"Band {b_info['band']}",
                    "E_min (eV)": f"{b_info['E_min']:.3f}",
                    "E_max (eV)": f"{b_info['E_max']:.3f}",
                    "Width (eV)": f"{b_info['width']:.3f}"
                })
            st.dataframe(band_data, use_container_width=True)
        else:
            st.write("No bands found in the given range.")

        # Table of all detected forbidden band gaps
        st.markdown("#### Forbidden Band Gaps")
        if gaps:
            gap_data = []
            for g_info in gaps:
                gap_data.append({
                    "Gap #": f"Gap {g_info['gap']}",
                    "Interval (eV)": f"{g_info['E_min']:.3f} -> {g_info['E_max']:.3f}",
                    "Gap Width (eV)": f"{g_info['width']:.3f}"
                })
            st.dataframe(gap_data, use_container_width=True)


# Tab 3: Numerical TISE solver via finite-difference method (Question 6)
with tab3:
    st.subheader("Direct Numerical Solution of the Schrodinger Differential Equation")
    st.markdown("""
    Here, the 1D Time-Independent Schrodinger Differential Equation:
    $$-\\frac{\\hbar^2}{2m} \\frac{d^2\\psi(x)}{dx^2} + V(x)\\psi(x) = E\\psi(x)$$
    is solved **directly on a discrete spatial grid** using the Finite-Difference Method (tridiagonal Hamiltonian matrix diagonalization).
    """)

    num_evals = min(N * 3, 15)
    evals_num, evecs_num, x_grid_num, V_grid_num = kp.solve_tise_numerical(N, a, b, V0, num_states=num_evals)

    fig3, (ax_spec, ax_modes) = plt.subplots(1, 2, figsize=(12, 5))

    # Left subplot: discrete eigenvalues overlaid on the potential and analytical bands
    ax_spec.plot(x_grid_num, V_grid_num, "k-", lw=1.5, alpha=0.5, label="Potential $V(x)$")
    for idx, e_val in enumerate(evals_num):
        ax_spec.axhline(e_val, color="royalblue", linestyle="-", lw=1.2, alpha=0.8)
    ax_spec.axhline(evals_num[0], color="royalblue", linestyle="-", lw=1.2, label=f"Numerical Eigenvalues ({num_evals} states)")

    # Shade allowed band regions from the analytical solution for comparison
    for b_info in bands:
        ax_spec.axhspan(b_info["E_min"], b_info["E_max"], color="lightgreen", alpha=0.3)

    ax_spec.set_ylabel("Energy (eV)")
    ax_spec.set_xlabel("Position $x$ (nm)")
    ax_spec.set_title("Discrete Eigenvalues in Potential Wells", fontsize=11, fontweight="bold")
    ax_spec.legend(loc="upper right")
    ax_spec.grid(True, linestyle=":", alpha=0.5)

    # Right subplot: plot the first 3 eigenstates
    colors = ["crimson", "darkblue", "darkgreen"]
    for i in range(min(3, len(evals_num))):
        ax_modes.plot(x_grid_num, evecs_num[:, i], color=colors[i], lw=1.8, label=f"State {i+1} ($E = {evals_num[i]:.2f}$ eV)")
    ax_modes.axhline(0, color="gray", linestyle=":", alpha=0.7)
    ax_modes.set_xlabel("Position $x$ (nm)")
    ax_modes.set_ylabel(r"$\psi_n(x)$ Amplitude")
    ax_modes.set_title("Lowest Numerical Eigen-wavefunctions", fontsize=11, fontweight="bold")
    ax_modes.legend(loc="upper right")
    ax_modes.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    st.pyplot(fig3)

    # Explanation of the physical significance of the eigenvalue spectrum
    st.success(f"""
    **Key Physical Insight:**
    For a finite chain of **$N = {N}$ wells**, each atomic energy level splits into **exactly $N$ discrete eigen-energies**.
    Notice how these $N$ states cluster exclusively inside the green shaded **Allowed Bands** calculated from the analytical Kronig-Penney model!
    In the limit as $N \\to \\infty$ (infinite bulk crystal), these discrete levels become infinitely dense, forming the continuous energy bands of crystalline solids.
    """)

    st.markdown("---")
    st.markdown("### Question 6: Parameter Study - Varying Well Width $a$")
    st.markdown("""
    Explore how the energy spectrum and band gaps evolve as the potential well width $a$ is varied while keeping barrier width $b$ and height $V_0$ fixed.
    """)

    if st.checkbox("Show Parameter Study: Energy Levels vs. Well Width $a$", value=True):
        st.image("figures/figure4_parameter_sweep_a.png", caption="Numerical TISE Solution: Energy Eigenvalues and Bandwidths vs. Well Width a", use_container_width=True)
        st.info(r"""
        **Physical Explanation of the Behavior (Question 6):**
        1. **Quantum Confinement Effect ($E \propto 1/a^2$):** As the well width $a$ increases, the electron has more room to spread out. Its kinetic energy drops according to $\sim 1/a^2$, shifting all energy levels downward.
        2. **Proliferation of Bound States:** Wider wells accommodate more discrete bound states below the barrier height $V_0 = 10\text{ eV}$.
        3. **Bandwidth Narrowing (Tight-Binding Limit):** At small $a$, states sit near the barrier top with large tunneling overlap ($\beta$ small), creating broad bands. At large $a$, states sink deep into the wells, barrier penetration weakens, and bands compress into razor-thin levels resembling isolated atoms.
        """)


# Tab 4: Theoretical background and academic references (Questions 1 & 7)
with tab4:
    st.subheader("Theoretical Context, Physical Interpretation & References")

    st.markdown("""
    ### 1. The Kronig-Penney Model & Its Significance (Question 1)
    - **Origin:** Formulated in 1931 by Ralph de Laer Kronig and William George Penney to investigate how an electron behaves in an idealized one-dimensional periodic crystal lattice.
    - **Physical Significance:**
      1. **Energy Band Formation:** Proves mathematically that a periodic potential splits the continuous energy spectrum into alternating **allowed bands** and **forbidden band gaps**.
      2. **Classification of Solids:** Provides the quantum basis for explaining why some crystalline solids are **metals/conductors** (partially filled bands or overlapping bands), **semiconductors** (narrow band gap $\\sim 1$ eV), or **insulators** (wide band gap $> 3$ eV).
      3. **Bragg Reflection:** Demonstrates that electron waves undergo constructive back-reflection at Brillouin zone boundaries ($k = \\pm \\pi/d$), forming standing waves that open energy gaps.
      4. **Effective Mass ($m^*$):** The curvature of the $E(k)$ dispersion relation determines electron dynamics inside the lattice: $m^* = \\hbar^2 \\left( \\frac{d^2E}{dk^2} \\right)^{-1}$.

    ### 2. Analytical Formulation & Bloch's Theorem
    The 1D Time-Independent Schrodinger Equation is:
    $$-\\frac{\\hbar^2}{2m} \\frac{d^2\\psi(x)}{dx^2} + V(x)\\psi(x) = E\\psi(x)$$

    By **Bloch's Theorem**, the stationary states in a periodic lattice satisfy:
    $$\\psi(x + d) = e^{ikd} \\psi(x)$$

    Matching boundary conditions of continuity for $\\psi(x)$ and $\\frac{d\\psi}{dx}$ across well ($a$) and barrier ($b$) interfaces yields:
    $$f(E) = \\cos(\\alpha a)\\cosh(\\beta b) + \\left(\\frac{\\beta^2 - \\alpha^2}{2\\alpha\\beta}\\right)\\sin(\\alpha a)\\sinh(\\beta b) = \\cos(kd)$$
    where $\\alpha = \\sqrt{2mE}/\\hbar$ and $\\beta = \\sqrt{2m(V_0 - E)}/\\hbar$.

    Since $-1 \\le \\cos(kd) \\le 1$:
    - **Allowed Bands:** $|f(E)| \\le 1 \\implies$ Real crystal momentum $k$.
    - **Forbidden Gaps:** $|f(E)| > 1 \\implies$ Evanescent wave, exponential decay.

    ---
    ### Academic References (Question 7)
    1. **Kronig, R. de L., & Penney, W. G.** (1931). *A Quantum Mechanics of Electrons in Crystal Lattices*. Proceedings of the Royal Society of London. Series A, 130(814), 499-513.
    2. **Kittel, C.** (2004). *Introduction to Solid State Physics* (8th ed.). John Wiley & Sons. (Chapters 7 & 8: Energy Bands and Kronig-Penney Model).
    3. **Ashcroft, N. W., & Mermin, N. D.** (1976). *Solid State Physics*. Saunders College Publishing. (Chapters 8 & 9: Electron Levels in a Periodic Potential).
    4. **Griffiths, D. J.** (2018). *Introduction to Quantum Mechanics* (3rd ed.). Cambridge University Press. (Section on Periodic Potentials & Band Theory).
    """)
