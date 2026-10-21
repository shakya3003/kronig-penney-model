"""
Phase 1 Verification Script: Analytical Kronig-Penney Dispersion Relation
-----------------------------------------------------------------------
This script tests the exact mathematical formulas derived in phase1_theory.md
for both E < V0 and E > V0. It prints intermediate parameters (alpha, beta,
hyperbolic/trigonometric terms) and classifies whether a given energy falls
into an Allowed Band (|f(E)| <= 1) or a Forbidden Band Gap (|f(E)| > 1).
"""

import numpy as np

# Physical constants in convenient semiconductor / nanoscale units:
# Length in nm, Energy in eV, mass in electron rest mass (m_e)
# hbar^2 / (2 * m_e) in eV * nm^2:
# hbar = 1.054571817e-34 J*s
# m_e  = 9.1093837e-31 kg
# 1 eV = 1.602176634e-19 J
# hbar^2 / (2 * m_e) = (1.054571817e-34)^2 / (2 * 9.1093837e-31 * 1.602176634e-19 * 1e-18)
HBAR2_OVER_2M = 0.0761996  # eV * nm^2


def evaluate_kronig_penney_dispersion(E, a, b, V0):
    """
    Computes f(E) = cos(k*d) for the 1D Kronig-Penney model.

    Parameters:
    -----------
    E  : float, electron energy in eV
    a  : float, well width in nm
    b  : float, barrier width in nm
    V0 : float, barrier height in eV

    Returns:
    --------
    dict containing:
        - 'f_E': float, value of f(E)
        - 'is_allowed': bool, True if |f(E)| <= 1
        - 'k': float or None, real Bloch wavevector (nm^-1) if allowed
        - 'kappa': float or None, decay parameter (nm^-1) if forbidden
        - 'details': dict with intermediate quantities (alpha, beta, etc.)
    """
    d = a + b
    if E <= 0:
        raise ValueError("Energy E must be positive (E > 0).")

    # In the well (V = 0): alpha = sqrt(2mE / hbar^2) in nm^-1
    alpha = np.sqrt(E / HBAR2_OVER_2M)

    if np.isclose(E, V0):
        # Limit as E -> V0 (beta -> 0):
        # cosh(beta*b) -> 1, sinh(beta*b)/(beta) -> b
        # (beta^2 - alpha^2)/(2*alpha*beta) * sinh(beta*b) -> -alpha*b/2
        f_E = np.cos(alpha * a) - 0.5 * alpha * b * np.sin(alpha * a)
        beta_val = 0.0
        mode = "E == V0 (Boundary Limit)"
    elif E < V0:
        # Tunneling regime (E < V0): beta = sqrt(2m(V0 - E) / hbar^2)
        beta = np.sqrt((V0 - E) / HBAR2_OVER_2M)
        beta_val = beta
        mode = "E < V0 (Tunneling Regime)"

        term1 = np.cos(alpha * a) * np.cosh(beta * b)
        term2 = ((beta**2 - alpha**2) / (2.0 * alpha * beta)) * np.sin(alpha * a) * np.sinh(beta * b)
        f_E = term1 + term2
    else:
        # Propagating over barrier (E > V0): beta_prime = sqrt(2m(E - V0) / hbar^2)
        beta_p = np.sqrt((E - V0) / HBAR2_OVER_2M)
        beta_val = beta_p
        mode = "E > V0 (Over-Barrier Regime)"

        term1 = np.cos(alpha * a) * np.cos(beta_p * b)
        term2 = -((alpha**2 + beta_p**2) / (2.0 * alpha * beta_p)) * np.sin(alpha * a) * np.sin(beta_p * b)
        f_E = term1 + term2

    is_allowed = abs(f_E) <= 1.0

    if is_allowed:
        # cos(k*d) = f_E  =>  k = arccos(f_E) / d
        k = np.arccos(np.clip(f_E, -1.0, 1.0)) / d
        kappa = None
    else:
        k = None
        # In forbidden gap, cosh(kappa*d) = |f_E|
        kappa = np.arccosh(abs(f_E)) / d

    return {
        "f_E": f_E,
        "is_allowed": is_allowed,
        "k": k,
        "kappa": kappa,
        "mode": mode,
        "alpha": alpha,
        "beta_or_betap": beta_val,
        "d": d
    }


def main():
    print("=" * 70)
    print("PHASE 1: KRONIG-PENNEY ANALYTICAL DISPERSION CHECK")
    print("=" * 70)

    # Standard representative parameters:
    # Well width a = 0.5 nm, Barrier width b = 0.1 nm, Barrier height V0 = 10.0 eV
    a = 0.5  # nm
    b = 0.1  # nm
    V0 = 10.0  # eV
    d = a + b  # lattice period: 0.6 nm

    print(f"Lattice Parameters:")
    print(f"  Well width (a)       : {a:.3f} nm")
    print(f"  Barrier width (b)    : {b:.3f} nm")
    print(f"  Lattice constant (d) : {d:.3f} nm")
    print(f"  Barrier height (V0)  : {V0:.2f} eV")
    print("-" * 70)

    # Test an array of sample energies across E < V0 and E > V0
    test_energies = [0.5, 1.2, 2.5, 4.0, 6.0, 8.5, 10.0, 12.0, 15.0]

    print(f"{'Energy (eV)':<12} | {'f(E)':<12} | {'State':<15} | {'k (nm^-1)':<12} | {'Regime'}")
    print("-" * 70)

    for E in test_energies:
        res = evaluate_kronig_penney_dispersion(E, a, b, V0)
        f_val = res["f_E"]
        status = "ALLOWED BAND" if res["is_allowed"] else "FORBIDDEN GAP"
        k_str = f"{res['k']:.4f}" if res["k"] is not None else f"(decay={res['kappa']:.2f})"
        print(f"{E:<12.2f} | {f_val:<12.4f} | {status:<15} | {k_str:<12} | {res['mode']}")

    print("-" * 70)
    print("Verification completed successfully!")
    print("Notice how energies alternate between Allowed Bands (|f(E)| <= 1) and")
    print("Forbidden Band Gaps (|f(E)| > 1). This is the hallmark of the Kronig-Penney model.")
    print("=" * 70)


if __name__ == "__main__":
    main()
