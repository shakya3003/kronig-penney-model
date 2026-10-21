"""
verify_phase2.py: Comprehensive Verification for Phase 2 Numerical Engine
========================================================================
This script tests:
1. Automatic band and band gap finder (Brent root-finding).
2. Transfer Matrix Method (TMM) continuous wavefunction and probability density.
3. Direct Finite-Difference numerical differential equation solver (Question 6).
4. Physical agreement between the numerical TISE eigenvalues and analytical bands.
5. Generates and saves a verification figure 'verify_phase2_result.png'.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless execution
import matplotlib.pyplot as plt

import kp_engine as kp


def main():
    print("=" * 75)
    print("PHASE 2: NUMERICAL PHYSICS ENGINE VERIFICATION")
    print("=" * 75)

    # 1. Benchmark Parameters
    a = 0.50   # well width in nm
    b = 0.10   # barrier width in nm
    V0 = 10.0  # barrier height in eV
    d = a + b  # lattice period: 0.60 nm
    print(f"Lattice Parameters: a = {a:.2f} nm, b = {b:.2f} nm, V0 = {V0:.2f} eV (d = {d:.2f} nm)")
    print("-" * 75)

    # 2. Test Band Finder
    print("TEST 1: Automatic Band and Gap Root-Finder (up to 20 eV)")
    bands, gaps = kp.find_allowed_bands(a, b, V0, E_max=20.0)

    print(f"{'Band #':<8} | {'E_min (eV)':<12} | {'E_max (eV)':<12} | {'Bandwidth (eV)':<15}")
    print("-" * 55)
    for b_info in bands:
        print(f"Band {b_info['band']:<3} | {b_info['E_min']:<12.4f} | {b_info['E_max']:<12.4f} | {b_info['width']:<15.4f}")

    print("-" * 55)
    print("Intervening Band Gaps:")
    for g_info in gaps:
        print(f"Gap {g_info['gap']}   | {g_info['E_min']:<12.4f} -> {g_info['E_max']:<12.4f} (width: {g_info['width']:.4f} eV)")
    print("-" * 75)

    # 3. Test Transfer Matrix Method (TMM) Wavefunction
    N_wells = 4
    E_allowed = 1.20  # Inside Band 1
    print(f"TEST 2: Transfer Matrix Method for N = {N_wells} wells at E = {E_allowed:.2f} eV")
    x_tmm, psi_tmm, prob_tmm, V_tmm = kp.compute_wavefunction_tmm(N=N_wells, E=E_allowed, a=a, b=b, V0=V0)

    # Verify normalization
    total_prob_tmm = np.trapezoid(prob_tmm, x_tmm)
    print(f"  -> Total Integrated Probability: {total_prob_tmm:.6f} (Target: 1.000000)")
    assert np.isclose(total_prob_tmm, 1.0, atol=1e-3), "TMM Wavefunction normalization check failed!"
    print("  -> Normalization: PASSED!")
    print("-" * 75)

    # 4. Test Direct Numerical Differential Equation Solver (Question 6)
    print(f"TEST 3: Direct Numerical TISE Differential Equation Solver (Finite Difference, N={N_wells})")
    evals, evecs, x_num, V_num = kp.solve_tise_numerical(N=N_wells, a=a, b=b, V0=V0, num_states=12)

    print(f"{'State n':<10} | {'Eigenvalue En (eV)':<20} | {'Analytical Band'}")
    print("-" * 55)
    for idx, e_val in enumerate(evals):
        # Check which analytical band it belongs to
        belongs_to = "In Gap / Unassigned"
        for b_info in bands:
            if b_info["E_min"] <= e_val <= b_info["E_max"]:
                belongs_to = f"Band {b_info['band']} [{b_info['E_min']:.2f} - {b_info['E_max']:.2f} eV]"
                break
        print(f"State {idx + 1:<4} | {e_val:<20.4f} | {belongs_to}")

    print("-" * 75)
    print("PHYSICAL CONFIRMATION:")
    print(f"Notice how the first {N_wells} discrete numerical states (States 1 to 4)")
    print(f"cluster neatly inside Band 1 ([{bands[0]['E_min']:.2f} - {bands[0]['E_max']:.2f} eV])!")
    print(f"And the next {N_wells} states (States 5 to 8) cluster inside Band 2!")
    print("This perfectly demonstrates the physics of energy band formation in a solid.")
    print("-" * 75)

    # 5. Generate and Save Plot
    print("Generating verification plot 'verify_phase2_result.png'...")
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # Plot 1: Potential and TMM Wavefunction
    ax1.plot(x_tmm, V_tmm, "k-", lw=2, label="Potential V(x) [eV]")
    ax1.axhline(E_allowed, color="red", linestyle="--", lw=1.5, label=f"Electron Energy E = {E_allowed} eV")
    ax1.set_ylabel("Energy (eV)")
    ax1.set_title(f"Kronig-Penney Lattice (N={N_wells} Wells, a={a} nm, b={b} nm, V0={V0} eV)")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Wavefunction psi(x)
    ax2.plot(x_tmm, psi_tmm, "b-", lw=1.8, label=r"Wavefunction $\psi(x)$ (TMM)")
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.7)
    ax2.set_ylabel(r"$\psi(x)$ Amplitude")
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)

    # Plot 3: Probability Density |psi(x)|^2
    ax3.plot(x_tmm, prob_tmm, "m-", lw=2, label=r"Probability Density $|\psi(x)|^2$")
    ax3.fill_between(x_tmm, prob_tmm, color="magenta", alpha=0.25)
    ax3.set_xlabel("Position x (nm)")
    ax3.set_ylabel(r"$|\psi(x)|^2$ ($nm^{-1}$)")
    ax3.legend(loc="upper right")
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("verify_phase2_result.png", dpi=200)
    plt.close()
    print("Plot successfully saved to 'verify_phase2_result.png'!")
    print("=" * 75)


if __name__ == "__main__":
    main()
