"""
numerical_analysis.py: Numerical Differential Equation Solver & Parameter Sweep (Question 6)
===========================================================================================
This script fulfills Question 6 of the assignment:
- Numerically solves the 1D Schrödinger Differential Equation:
    -[hbar^2 / (2m)] * d^2 psi / dx^2 + V(x)*psi(x) = E * psi(x)
- Sweeps the well width 'a' from 0.20 nm to 1.20 nm with fixed b = 0.10 nm, V0 = 10.0 eV.
- Compares the numerical eigenvalues against analytical Kronig-Penney bands.
- Investigates the quantum confinement effect (E ~ 1/a^2) and band narrowing.
- Generates publication figure: 'figure4_parameter_sweep_a.png'.

Run with:
    python3 numerical_analysis.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive publication backend
import matplotlib.pyplot as plt

import kp_engine as kp


def perform_parameter_sweep(b=0.10, V0=10.0, num_a_points=25):
    """
    Sweeps well width 'a' and computes both numerical eigenvalues and analytical band edges.
    """
    a_values = np.linspace(0.20, 1.20, num_a_points)
    num_eigenstates = 8

    # Store numerical eigenvalues for N=4 wells
    numerical_energies = np.zeros((len(a_values), num_eigenstates))

    # Store analytical band edges for bands 1, 2, 3
    band1_min = np.zeros(len(a_values))
    band1_max = np.zeros(len(a_values))
    band2_min = np.zeros(len(a_values))
    band2_max = np.zeros(len(a_values))
    band3_min = np.zeros(len(a_values))
    band3_max = np.zeros(len(a_values))

    print("Running numerical differential equation solver across parameter sweep...")
    for idx, a_val in enumerate(a_values):
        # 1. Numerical Differential Equation solution (N=4 wells)
        evals, _, _, _ = kp.solve_tise_numerical(
            N=4,
            a=a_val,
            b=b,
            V0=V0,
            num_states=num_eigenstates,
            points_per_cell=100
        )
        numerical_energies[idx, :len(evals)] = evals

        # 2. Analytical Band Edges
        bands, _ = kp.find_allowed_bands(a=a_val, b=b, V0=V0, E_max=22.0)
        if len(bands) >= 1:
            band1_min[idx] = bands[0]["E_min"]
            band1_max[idx] = bands[0]["E_max"]
        if len(bands) >= 2:
            band2_min[idx] = bands[1]["E_min"]
            band2_max[idx] = bands[1]["E_max"]
        if len(bands) >= 3:
            band3_min[idx] = bands[2]["E_min"]
            band3_max[idx] = bands[2]["E_max"]

    return a_values, numerical_energies, (band1_min, band1_max, band2_min, band2_max, band3_min, band3_max)


def generate_parameter_sweep_figure(a_values, numerical_energies, band_edges, b=0.10, V0=10.0, filename="figure4_parameter_sweep_a.png"):
    """
    Generates a comprehensive 3-panel figure visualizing the E vs. a behavior.
    """
    b1_min, b1_max, b2_min, b2_max, b3_min, b3_max = band_edges

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0])

    # -----------------------------------------------------------------------
    # PANEL 1: Energy Eigenvalues & Bands vs. Well Width 'a'
    # -----------------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, :])

    # Shaded analytical bands
    ax1.fill_between(a_values, b1_min, b1_max, color="lightgreen", alpha=0.35, label="Analytical Band 1 (KP)")
    ax1.fill_between(a_values, b2_min, b2_max, color="lightblue", alpha=0.35, label="Analytical Band 2 (KP)")
    ax1.fill_between(a_values, b3_min, b3_max, color="navajowhite", alpha=0.35, label="Analytical Band 3 (KP)")

    # Barrier height threshold V0
    ax1.axhline(V0, color="crimson", linestyle="--", lw=1.5, label=r"Barrier Height $V_0 = 10.0$ eV")

    # Plot numerical eigenvalues (discrete dots / curves)
    for state_i in range(numerical_energies.shape[1]):
        label_text = "Numerical TISE Eigenvalues ($N=4$)" if state_i == 0 else ""
        ax1.plot(a_values, numerical_energies[:, state_i], "o-", color="darkblue", markersize=3.5, lw=1.2, alpha=0.8, label=label_text)

    # Quantum confinement 1/a^2 trend guide
    a_fit = np.linspace(0.20, 1.20, 100)
    fit_curve = numerical_energies[0, 0] * (a_values[0] / a_fit)**1.8
    ax1.plot(a_fit, fit_curve, "r:", lw=1.5, label=r"Quantum Confinement Scaling Trend ($\sim 1/a^2$)")

    ax1.set_ylim(0, 16)
    ax1.set_xlim(0.18, 1.22)
    ax1.set_xlabel("Potential Well Width $a$ (nm)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Energy (eV)", fontsize=11, fontweight="bold")
    ax1.set_title(r"Numerical Solution of TISE: Energy States vs. Well Width $a$ (with $b = 0.1$ nm, $V_0 = 10$ eV)", fontsize=13, fontweight="bold")
    ax1.legend(loc="upper right", framealpha=0.9, fontsize=9.5)
    ax1.grid(True, linestyle=":", alpha=0.6)

    # -----------------------------------------------------------------------
    # PANEL 2: Bandwidth Evolution (Band 1 & Band 2)
    # -----------------------------------------------------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    w1 = b1_max - b1_min
    w2 = b2_max - b2_min

    ax2.plot(a_values, w1, "g-o", lw=2, markersize=4, label="Band 1 Width")
    ax2.plot(a_values, w2, "b-s", lw=2, markersize=4, label="Band 2 Width")
    ax2.set_xlabel("Well Width $a$ (nm)", fontsize=11)
    ax2.set_ylabel(r"Band Width $\Delta E$ (eV)", fontsize=11)
    ax2.set_title("Bandwidth Narrowing (Transition to Tight-Binding)", fontsize=11, fontweight="bold")
    ax2.legend(loc="upper right", framealpha=0.9)
    ax2.grid(True, linestyle=":", alpha=0.6)

    # -----------------------------------------------------------------------
    # PANEL 3: Wavefunction Comparison (Narrow well vs. Wide well)
    # -----------------------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 1])

    # Narrow well: a = 0.25 nm
    _, psi_narrow, _, V_narrow = kp.compute_wavefunction_tmm(N=3, E=3.50, a=0.25, b=0.10, V0=10.0)
    x_narrow = np.linspace(0, 0.10 + 3 * (0.25 + 0.10), len(psi_narrow))

    # Wide well: a = 0.90 nm
    _, psi_wide, _, V_wide = kp.compute_wavefunction_tmm(N=3, E=0.70, a=0.90, b=0.10, V0=10.0)
    x_wide = np.linspace(0, 0.10 + 3 * (0.90 + 0.10), len(psi_wide))

    ax3.plot(x_narrow, psi_narrow + 2.0, color="crimson", lw=1.8, label=r"Narrow Well ($a = 0.25$ nm): Strong Tunneling Overlap")
    ax3.plot(x_wide, psi_wide - 1.0, color="navy", lw=1.8, label=r"Wide Well ($a = 0.90$ nm): Deep Atomic Localization")
    ax3.axhline(2.0, color="gray", linestyle=":", alpha=0.5)
    ax3.axhline(-1.0, color="gray", linestyle=":", alpha=0.5)

    ax3.set_xlabel("Spatial Position $x$ (nm)", fontsize=11)
    ax3.set_ylabel(r"Wavefunction Amplitude $\psi(x)$ (Offset)", fontsize=11)
    ax3.set_title(r"Wavefunction Spatial Profile: Narrow vs. Wide Wells", fontsize=11, fontweight="bold")
    ax3.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax3.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(filename, dpi=250)
    plt.close()
    print(f"[OK] Saved: {filename}")


def main():
    print("=" * 75)
    print("PHASE 5: NUMERICAL DIFFERENTIAL EQUATION SOLVER & PARAMETER SWEEP")
    print("=" * 75)

    b = 0.10   # fixed barrier width (nm)
    V0 = 10.0  # fixed barrier height (eV)

    a_values, numerical_energies, band_edges = perform_parameter_sweep(b=b, V0=V0, num_a_points=21)

    print("\nSample Output: Lowest 4 Eigenvalues (eV) across Selected Well Widths:")
    print(f"{'a (nm)':<10} | {'State 1 (eV)':<14} | {'State 2 (eV)':<14} | {'State 3 (eV)':<14} | {'State 4 (eV)':<14}")
    print("-" * 72)
    sample_indices = [0, 4, 8, 12, 16, 20]
    for idx in sample_indices:
        a_val = a_values[idx]
        evals = numerical_energies[idx, :4]
        print(f"{a_val:<10.2f} | {evals[0]:<14.4f} | {evals[1]:<14.4f} | {evals[2]:<14.4f} | {evals[3]:<14.4f}")

    print("-" * 75)
    print("Generating comprehensive publication figure 'figure4_parameter_sweep_a.png'...")
    generate_parameter_sweep_figure(a_values, numerical_energies, band_edges, b=b, V0=V0)

    print("=" * 75)
    print("Phase 5 calculations completed successfully!")


if __name__ == "__main__":
    main()
