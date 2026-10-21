"""
plot_generator.py: Standalone Publication-Quality Plot Generator

Generates high-resolution figures for the assignment deliverables:
1. figure1_allowed_state.png  : Wavefunction and probability density for an allowed state.
2. figure2_forbidden_gap.png  : Spatial decay of wavefunction in a forbidden band gap.
3. figure3_ek_dispersion.png  : E(k) band structure in the first Brillouin zone.

Run with:
    python3 plot_generator.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for saving figures without a display
import matplotlib.pyplot as plt

import kp_engine as kp


def plot_allowed_state(N=5, E=1.20, a=0.50, b=0.10, V0=10.0, filename="figures/figure1_allowed_state.png"):
    """Generates a 3-panel figure showing an allowed propagating Bloch state."""
    d = a + b
    res = kp.compute_dispersion(E, a, b, V0)
    x, psi, prob, V = kp.compute_wavefunction_tmm(N, E, a, b, V0)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    # Panel 1: Periodic potential landscape with electron energy level
    ax1.plot(x, V, "k-", lw=2, label=r"Periodic Potential $V(x)$")
    ax1.axhline(E, color="crimson", linestyle="--", lw=1.8, label=f"Electron Energy $E = {E:.2f}$ eV")
    ax1.fill_between(x, 0, V, color="lightgray", alpha=0.5, label=f"Barriers ($V_0 = {V0:.1f}$ eV)")
    ax1.set_ylabel("Energy (eV)", fontsize=11)
    ax1.set_title(f"ALLOWED STATE: $E = {E:.2f}$ eV, $f(E) = {res['f_E']:.3f}$, $k = {res['k']:.3f}$ nm$^{{-1}}$\n"
                  f"($N = {N}$ Wells, $a = {a}$ nm, $b = {b}$ nm, $V_0 = {V0}$ eV)", fontsize=12, fontweight="bold")
    ax1.legend(loc="upper right", framealpha=0.9)
    ax1.grid(True, linestyle=":", alpha=0.6)

    # Panel 2: Wavefunction psi(x) computed via TMM
    ax2.plot(x, psi, color="#1f77b4", lw=2, label=r"Wavefunction $\psi(x)$ (TMM)")
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.7)
    ax2.set_ylabel(r"$\psi(x)$ Amplitude", fontsize=11)
    ax2.set_title(r"Continuous Wavefunction $\psi(x)$ (Oscillatory in Wells, Tunneling in Barriers)", fontsize=11, fontweight="bold")
    ax2.legend(loc="upper right", framealpha=0.9)
    ax2.grid(True, linestyle=":", alpha=0.6)

    # Panel 3: Probability density |psi(x)|^2
    ax3.plot(x, prob, color="#9467bd", lw=2, label=r"Probability Density $|\psi(x)|^2$")
    ax3.fill_between(x, prob, color="#9467bd", alpha=0.3)
    ax3.set_xlabel("Spatial Position $x$ (nm)", fontsize=11)
    ax3.set_ylabel(r"$|\psi(x)|^2$ (nm$^{-1}$)", fontsize=11)
    ax3.set_title(r"Probability Distribution (Balanced Extended Bloch Wave across All Wells)", fontsize=11, fontweight="bold")
    ax3.legend(loc="upper right", framealpha=0.9)
    ax3.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(filename, dpi=250)
    plt.close()
    print(f"Saved: {filename}")


def plot_forbidden_state(N=5, E=3.20, a=0.50, b=0.10, V0=10.0, filename="figures/figure2_forbidden_gap.png"):
    """Generates a 3-panel figure showing wave attenuation inside a forbidden band gap."""
    d = a + b
    res = kp.compute_dispersion(E, a, b, V0)
    x, psi, prob, V = kp.compute_wavefunction_tmm(N, E, a, b, V0)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    # Panel 1: Potential landscape with energy in a forbidden gap
    ax1.plot(x, V, "k-", lw=2, label=r"Periodic Potential $V(x)$")
    ax1.axhline(E, color="crimson", linestyle="--", lw=1.8, label=f"Electron Energy $E = {E:.2f}$ eV")
    ax1.fill_between(x, 0, V, color="lightgray", alpha=0.5, label=f"Barriers ($V_0 = {V0:.1f}$ eV)")
    ax1.set_ylabel("Energy (eV)", fontsize=11)
    ax1.set_title(f"FORBIDDEN BAND GAP: $E = {E:.2f}$ eV, $f(E) = {res['f_E']:.3f}$, Decay $\\kappa = {res['kappa']:.3f}$ nm$^{{-1}}$\n"
                  f"($N = {N}$ Wells, $a = {a}$ nm, $b = {b}$ nm, $V_0 = {V0}$ eV)", fontsize=12, fontweight="bold", color="darkred")
    ax1.legend(loc="upper right", framealpha=0.9)
    ax1.grid(True, linestyle=":", alpha=0.6)

    # Panel 2: Attenuated wavefunction showing evanescent decay
    ax2.plot(x, psi, color="#d62728", lw=2, label=r"Wavefunction $\psi(x)$ (Evanescent Attenuation)")
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.7)
    ax2.set_ylabel(r"$\psi(x)$ Amplitude", fontsize=11)
    ax2.set_title(r"Wavefunction Attenuation (Evanescent Decay via Bragg Reflection)", fontsize=11, fontweight="bold")
    ax2.legend(loc="upper right", framealpha=0.9)
    ax2.grid(True, linestyle=":", alpha=0.6)

    # Panel 3: Probability density showing rapid exponential decay
    ax3.plot(x, prob, color="#e377c2", lw=2, label=r"Probability Density $|\psi(x)|^2$")
    ax3.fill_between(x, prob, color="#e377c2", alpha=0.3)
    ax3.set_xlabel("Spatial Position $x$ (nm)", fontsize=11)
    ax3.set_ylabel(r"$|\psi(x)|^2$ (nm$^{-1}$)", fontsize=11)
    ax3.set_title(r"Probability Distribution (Rapid Exponential Decay through the Crystal)", fontsize=11, fontweight="bold")
    ax3.legend(loc="upper right", framealpha=0.9)
    ax3.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(filename, dpi=250)
    plt.close()
    print(f"Saved: {filename}")


def plot_ek_band_structure(a=0.50, b=0.10, V0=10.0, E_max=20.0, filename="figures/figure3_ek_dispersion.png"):
    """Plots E vs k dispersion in the Reduced Brillouin Zone Scheme."""
    d = a + b
    bands, gaps = kp.find_allowed_bands(a, b, V0, E_max=E_max)

    fig, ax = plt.subplots(figsize=(7, 7))

    k_bz = np.pi / d  # Brillouin zone boundary wavevector

    # Compute E(k) curve within each allowed band
    energies = np.linspace(0.01, E_max, 2500)
    for b_info in bands:
        e_mask = (energies >= b_info["E_min"]) & (energies <= b_info["E_max"])
        e_subset = energies[e_mask]
        k_subset = []
        for e_val in e_subset:
            res = kp.compute_dispersion(e_val, a, b, V0)
            k_subset.append(res["k"])
        k_subset = np.array(k_subset)

        # Plot both positive and negative k branches (reduced zone scheme)
        ax.plot(k_subset, e_subset, color="darkblue", lw=2.2)
        ax.plot(-k_subset, e_subset, color="darkblue", lw=2.2)

        # Shade allowed band region and label it
        ax.axhspan(b_info["E_min"], b_info["E_max"], color="lightgreen", alpha=0.2)
        ax.text(k_bz * 0.95, (b_info["E_min"] + b_info["E_max"]) / 2, f"Band {b_info['band']}",
                va="center", ha="right", fontsize=10, fontweight="bold", color="darkgreen")

    # Shade forbidden gap regions
    for g_info in gaps:
        ax.axhspan(g_info["E_min"], g_info["E_max"], color="lightcoral", alpha=0.15)
        ax.text(0.0, (g_info["E_min"] + g_info["E_max"]) / 2, f"Gap {g_info['gap']} ({g_info['width']:.2f} eV)",
                va="center", ha="center", fontsize=9, color="darkred", fontstyle="italic")

    # Brillouin zone boundary markers
    ax.axvline(k_bz, color="gray", linestyle="--", lw=1.2, label=r"Brillouin Zone Edge ($k = +\pi/d$)")
    ax.axvline(-k_bz, color="gray", linestyle="--", lw=1.2, label=r"Brillouin Zone Edge ($k = -\pi/d$)")
    ax.axvline(0, color="lightgray", linestyle=":", lw=1)

    ax.set_xlim(-k_bz * 1.15, k_bz * 1.15)
    ax.set_ylim(0, E_max)
    ax.set_xlabel(r"Bloch Wavevector $k$ (nm$^{-1}$)", fontsize=11)
    ax.set_ylabel(r"Electron Energy $E$ (eV)", fontsize=11)
    ax.set_title(f"Electronic Band Structure $E(k)$ (Reduced Zone Scheme)\n"
                 f"($a = {a}$ nm, $b = {b}$ nm, $V_0 = {V0}$ eV, $d = {d}$ nm)", fontsize=12, fontweight="bold")
    ax.set_xticks([-k_bz, 0, k_bz])
    ax.set_xticklabels([r"$-\pi/d$", r"$0$", r"$+\pi/d$"])
    ax.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(filename, dpi=250)
    plt.close()
    print(f"Saved: {filename}")


def main():
    print("Generating publication-quality figures...")
    print("-" * 50)

    os.makedirs("figures", exist_ok=True)
    a, b, V0 = 0.50, 0.10, 10.00
    plot_allowed_state(N=5, E=1.20, a=a, b=b, V0=V0, filename="figures/figure1_allowed_state.png")
    plot_forbidden_state(N=5, E=3.20, a=a, b=b, V0=V0, filename="figures/figure2_forbidden_gap.png")
    plot_ek_band_structure(a=a, b=b, V0=V0, E_max=20.0, filename="figures/figure3_ek_dispersion.png")

    print("-" * 50)
    print("All 3 publication figures generated successfully.")


if __name__ == "__main__":
    main()
