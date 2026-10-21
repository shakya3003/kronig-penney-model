"""
kp_engine.py: Core Numerical Physics Engine for the Kronig-Penney Model

This module implements both analytical and numerical solutions for the 1D
Kronig-Penney model and the 1D Time-Independent Schrodinger Equation (TISE).

Key Capabilities:
1. Physical units normalized for nanostructures (eV, nm, m_e).
2. Analytical dispersion relation f(E) and exact band edge root-finding.
3. Transfer Matrix Method (TMM) to construct continuous wavefunctions psi(x).
4. Direct Finite-Difference differential equation solver for the 1D TISE.
5. Probability density |psi(x)|^2 computation and normalization.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal

# Physical constant: hbar^2 / (2 * m_e) in eV * nm^2
# Using nanoscale units (eV, nm, m_e) to prevent numerical overflow/underflow.
HBAR2_OVER_2M = 0.0761996  # eV * nm^2


# Dispersion relation: computes f(E) = cos(kd) for a given electron energy
def compute_dispersion(E, a, b, V0):
    """
    Computes the Kronig-Penney dispersion function f(E) = cos(k*d).

    Parameters:
        E  : float, electron energy in eV (must be > 0)
        a  : float, well width in nm
        b  : float, barrier width in nm
        V0 : float, barrier potential height in eV

    Returns:
        dict containing:
            - 'f_E': float, value of f(E)
            - 'is_allowed': bool, True if |f(E)| <= 1
            - 'k': float or None, real Bloch wavevector in nm^-1 (first Brillouin zone)
            - 'kappa': float or None, decay parameter in nm^-1 if in a band gap
            - 'alpha': float, wavevector inside the well (nm^-1)
            - 'beta': float, decay/wavevector inside the barrier (nm^-1)
            - 'd': float, lattice constant a + b (nm)
    """
    d = a + b
    if E <= 0:
        raise ValueError("Electron energy E must be strictly positive (E > 0).")

    # Wavevector inside the well (V = 0): alpha = sqrt(E / gamma)
    alpha = np.sqrt(E / HBAR2_OVER_2M)

    if np.isclose(E, V0, atol=1e-7):
        # Limiting case as E approaches V0 (beta approaches 0)
        f_E = np.cos(alpha * a) - 0.5 * alpha * b * np.sin(alpha * a)
        beta_val = 0.0
    elif E < V0:
        # Tunneling regime (E < V0): wavefunction decays exponentially in barriers
        beta = np.sqrt((V0 - E) / HBAR2_OVER_2M)
        beta_val = beta
        term1 = np.cos(alpha * a) * np.cosh(beta * b)
        term2 = ((beta**2 - alpha**2) / (2.0 * alpha * beta)) * np.sin(alpha * a) * np.sinh(beta * b)
        f_E = term1 + term2
    else:
        # Over-barrier regime (E > V0): wavefunction oscillates in both regions
        beta_p = np.sqrt((E - V0) / HBAR2_OVER_2M)
        beta_val = beta_p
        term1 = np.cos(alpha * a) * np.cos(beta_p * b)
        term2 = -((alpha**2 + beta_p**2) / (2.0 * alpha * beta_p)) * np.sin(alpha * a) * np.sin(beta_p * b)
        f_E = term1 + term2

    is_allowed = abs(f_E) <= 1.0

    if is_allowed:
        # Allowed state: extract real Bloch wavevector from cos(k*d) = f_E
        k = np.arccos(np.clip(f_E, -1.0, 1.0)) / d
        kappa = None
    else:
        # Forbidden gap: extract imaginary wavevector (decay rate) from cosh(kappa*d) = |f_E|
        k = None
        kappa = np.arccosh(abs(f_E)) / d

    return {
        "f_E": float(f_E),
        "is_allowed": is_allowed,
        "k": float(k) if k is not None else None,
        "kappa": float(kappa) if kappa is not None else None,
        "alpha": float(alpha),
        "beta": float(beta_val),
        "d": float(d),
    }


# Evaluates f(E) over a continuous energy range for dispersion curve plotting
def compute_dispersion_curve(a, b, V0, E_min=0.01, E_max=20.0, num_points=1000):
    """
    Evaluates f(E) over a continuous range of energies for plotting.

    Returns:
        energies : np.ndarray, array of energy values (eV)
        f_vals   : np.ndarray, array of f(E) values
        allowed  : np.ndarray (bool), True where |f(E)| <= 1
        k_vals   : np.ndarray, Bloch wavevector k (nm^-1) where allowed, NaN elsewhere
    """
    energies = np.linspace(E_min, E_max, num_points)
    f_vals = np.empty_like(energies)
    allowed = np.empty(len(energies), dtype=bool)
    k_vals = np.full(len(energies), np.nan)

    d = a + b
    for i, E in enumerate(energies):
        res = compute_dispersion(E, a, b, V0)
        f_vals[i] = res["f_E"]
        allowed[i] = res["is_allowed"]
        if res["is_allowed"]:
            k_vals[i] = res["k"]

    return energies, f_vals, allowed, k_vals


# Band edge root-finder: locates all allowed [E_min, E_max] intervals using Brent's method
def find_allowed_bands(a, b, V0, E_max=20.0, num_samples=1500):
    """
    Locates all allowed energy bands [E_min, E_max] up to E_max using
    bracketed Brent root-finding on f(E) = +1 and f(E) = -1.

    Returns:
        bands : list of dicts, each with keys:
                'band': int (1-based index)
                'E_min': float (eV)
                'E_max': float (eV)
                'width': float (eV)
        gaps  : list of dicts, each with keys:
                'gap': int
                'E_min': float
                'E_max': float
                'width': float
    """
    E_grid = np.linspace(0.01, E_max, num_samples)
    f_grid = np.array([compute_dispersion(E, a, b, V0)["f_E"] for E in E_grid])

    # Helper functions for root-finding at f(E) = +1 and f(E) = -1
    def f_minus_1(E):
        return compute_dispersion(E, a, b, V0)["f_E"] - 1.0

    def f_plus_1(E):
        return compute_dispersion(E, a, b, V0)["f_E"] + 1.0

    # Identify continuous segments where |f(E)| <= 1 (allowed regions)
    is_allowed = np.abs(f_grid) <= 1.0
    transitions = np.diff(is_allowed.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]

    # Handle cases where the grid starts or ends inside an allowed region
    if is_allowed[0]:
        starts = np.insert(starts, 0, 0)
    if is_allowed[-1]:
        ends = np.append(ends, len(E_grid) - 2)

    bands = []
    band_idx = 1
    for s_idx, e_idx in zip(starts, ends):
        # Refine the lower band edge using Brent's method
        E_start_bracket = (E_grid[s_idx], E_grid[min(s_idx + 1, len(E_grid) - 1)])
        f_s0 = f_grid[s_idx]
        f_s1 = f_grid[min(s_idx + 1, len(E_grid) - 1)]

        try:
            if (f_s0 - 1.0) * (f_s1 - 1.0) <= 0:
                E_min = brentq(f_minus_1, E_start_bracket[0], E_start_bracket[1])
            elif (f_s0 + 1.0) * (f_s1 + 1.0) <= 0:
                E_min = brentq(f_plus_1, E_start_bracket[0], E_start_bracket[1])
            else:
                E_min = E_grid[s_idx]
        except ValueError:
            E_min = E_grid[s_idx]

        # Refine the upper band edge using Brent's method
        E_end_bracket = (E_grid[e_idx], E_grid[min(e_idx + 1, len(E_grid) - 1)])
        f_e0 = f_grid[e_idx]
        f_e1 = f_grid[min(e_idx + 1, len(E_grid) - 1)]

        try:
            if (f_e0 - 1.0) * (f_e1 - 1.0) <= 0:
                E_max_band = brentq(f_minus_1, E_end_bracket[0], E_end_bracket[1])
            elif (f_e0 + 1.0) * (f_e1 + 1.0) <= 0:
                E_max_band = brentq(f_plus_1, E_end_bracket[0], E_end_bracket[1])
            else:
                E_max_band = E_grid[e_idx]
        except ValueError:
            E_max_band = E_grid[e_idx]

        if E_max_band > E_min:
            bands.append({
                "band": band_idx,
                "E_min": float(E_min),
                "E_max": float(E_max_band),
                "width": float(E_max_band - E_min),
            })
            band_idx += 1

    # Derive band gaps as the intervals between successive allowed bands
    gaps = []
    for i in range(len(bands) - 1):
        gap_min = bands[i]["E_max"]
        gap_max = bands[i + 1]["E_min"]
        if gap_max > gap_min:
            gaps.append({
                "gap": i + 1,
                "E_min": float(gap_min),
                "E_max": float(gap_max),
                "width": float(gap_max - gap_min),
            })

    return bands, gaps


# Transfer matrix for propagation through a well region (V = 0)
def _well_matrix(dx, alpha):
    """Transfer matrix across distance dx inside a well (V = 0)."""
    cos_val = np.cos(alpha * dx)
    sin_val = np.sin(alpha * dx)
    return np.array([
        [cos_val, sin_val / alpha],
        [-alpha * sin_val, cos_val]
    ], dtype=float)


# Transfer matrix for propagation through a barrier region (V = V0)
def _barrier_matrix(dx, E, V0):
    """Transfer matrix across distance dx inside a barrier (V = V0)."""
    if E < V0:
        beta = np.sqrt((V0 - E) / HBAR2_OVER_2M)
        cosh_val = np.cosh(beta * dx)
        sinh_val = np.sinh(beta * dx)
        return np.array([
            [cosh_val, sinh_val / beta],
            [beta * sinh_val, cosh_val]
        ], dtype=float)
    elif np.isclose(E, V0):
        return np.array([
            [1.0, dx],
            [0.0, 1.0]
        ], dtype=float)
    else:
        beta_p = np.sqrt((E - V0) / HBAR2_OVER_2M)
        cos_val = np.cos(beta_p * dx)
        sin_val = np.sin(beta_p * dx)
        return np.array([
            [cos_val, sin_val / beta_p],
            [-beta_p * sin_val, cos_val]
        ], dtype=float)


# TMM wavefunction builder: constructs psi(x) across N wells + barriers
def compute_wavefunction_tmm(N, E, a, b, V0, points_per_segment=80):
    """
    Computes the continuous, smooth wavefunction psi(x) across a finite
    chain of N potential wells separated by barriers using the Transfer Matrix Method.

    Lattice layout:
        [Barrier b] + N * ([Well a] + [Barrier b])
    Total segments = 2*N + 1.

    Parameters:
        N  : int, number of potential wells (e.g. 1 to 20)
        E  : float, electron energy in eV
        a  : float, well width in nm
        b  : float, barrier width in nm
        V0 : float, barrier height in eV
        points_per_segment : int, spatial resolution points per well/barrier segment

    Returns:
        x_total     : np.ndarray, spatial coordinates (nm)
        psi_total   : np.ndarray, real-valued normalized wavefunction psi(x)
        prob_total  : np.ndarray, probability density |psi(x)|^2 (nm^-1)
        V_total     : np.ndarray, potential profile V(x) in eV
    """
    alpha = np.sqrt(E / HBAR2_OVER_2M)

    # Initial state vector at x = 0: [psi(0), psi'(0)]
    state = np.array([1.0, 0.0], dtype=float)

    # Build the segment list: outer barrier, then alternating well-barrier pairs
    segments = [("barrier", b)]
    for _ in range(N):
        segments.append(("well", a))
        segments.append(("barrier", b))

    x_list = []
    psi_list = []
    V_list = []

    current_x = 0.0

    for seg_type, seg_width in segments:
        dx_grid = np.linspace(0.0, seg_width, points_per_segment, endpoint=False)
        V_seg = V0 if seg_type == "barrier" else 0.0

        # Compute psi at each grid point within this segment
        for dx in dx_grid:
            if seg_type == "well":
                M_dx = _well_matrix(dx, alpha)
            else:
                M_dx = _barrier_matrix(dx, E, V0)

            local_state = M_dx @ state
            x_list.append(current_x + dx)
            psi_list.append(local_state[0])
            V_list.append(V_seg)

        # Propagate the state vector to the end of this segment
        if seg_type == "well":
            M_end = _well_matrix(seg_width, alpha)
        else:
            M_end = _barrier_matrix(seg_width, E, V0)

        state = M_end @ state
        current_x += seg_width

    # Append the final endpoint
    x_list.append(current_x)
    psi_list.append(state[0])
    V_list.append(segments[-1][0] == "barrier" and V0 or 0.0)

    x_arr = np.array(x_list)
    psi_arr = np.array(psi_list)
    V_arr = np.array(V_list)

    # Normalize the wavefunction so that integral(|psi|^2 dx) = 1
    dx = np.diff(x_arr)
    prob_unnorm = psi_arr**2
    norm_factor = np.sqrt(np.trapezoid(prob_unnorm, x_arr))
    if norm_factor > 1e-12:
        psi_norm = psi_arr / norm_factor
    else:
        psi_norm = psi_arr

    prob_density = psi_norm**2

    return x_arr, psi_norm, prob_density, V_arr


# Finite-difference TISE solver: discretizes the Schrodinger equation on a spatial grid
def solve_tise_numerical(N, a, b, V0, num_states=10, points_per_cell=100):
    """
    Directly solves the 1D Time-Independent Schrodinger Differential Equation
    using the Finite-Difference Method and tridiagonal matrix diagonalization:

        -[hbar^2 / (2m)] * d^2 psi / dx^2 + V(x)*psi(x) = E * psi(x)

    Parameters:
        N : int, number of potential wells
        a : float, well width in nm
        b : float, barrier width in nm
        V0: float, barrier height in eV
        num_states: int, number of lowest energy eigenstates to return
        points_per_cell: int, spatial resolution grid points per period (d = a + b)

    Returns:
        eigenvalues  : np.ndarray, lowest eigen-energies E_n in eV
        eigenvectors : np.ndarray, normalized wavefunctions psi_n(x) [shape: (M, num_states)]
        x_grid       : np.ndarray, spatial coordinate grid in nm
        V_grid       : np.ndarray, potential profile V(x) in eV
    """
    d = a + b
    # Total lattice length: outer barrier (b) + N * (well a + barrier b)
    total_length = b + N * (a + b)
    dx = d / points_per_cell
    M = int(np.round(total_length / dx))
    x_grid = np.linspace(0.0, total_length, M)

    # Construct the periodic potential V(x) on the spatial grid
    V_grid = np.zeros(M, dtype=float)
    for i, x in enumerate(x_grid):
        if x < b:
            V_grid[i] = V0  # initial outer barrier
        else:
            rel_x = (x - b) % d
            if rel_x < a:
                V_grid[i] = 0.0  # inside well
            else:
                V_grid[i] = V0   # inside barrier

    # Build finite-difference Hamiltonian matrix diagonals:
    # Main diagonal: 2*gamma/dx^2 + V(x_i)
    # Off-diagonal:  -gamma/dx^2
    diag_kinetic = 2.0 * HBAR2_OVER_2M / (dx**2)
    offdiag_kinetic = -HBAR2_OVER_2M / (dx**2)

    diag = diag_kinetic + V_grid
    offdiag = np.full(M - 1, offdiag_kinetic)

    # Solve the tridiagonal symmetric eigenvalue problem: H * psi = E * psi
    num_evals = min(num_states, M - 2)
    evals, evecs = eigh_tridiagonal(diag, offdiag, select="i", select_range=(0, num_evals - 1))

    # Normalize each eigenvector so that integral(|psi(x)|^2 dx) = 1
    for n in range(num_evals):
        norm = np.sqrt(np.trapezoid(evecs[:, n]**2, x_grid))
        if norm > 1e-12:
            evecs[:, n] /= norm

    return evals, evecs, x_grid, V_grid


# Parameter sweep engine: varies well width 'a' and tracks eigenvalue evolution
def sweep_well_width_a(a_values, b, V0, num_states=5):
    """
    Sweeps the well width 'a' to observe how energy eigenvalues and band
    broadening evolve with quantum confinement.

    Parameters:
        a_values   : np.ndarray or list of float, well widths to test (in nm)
        b          : float, fixed barrier width (nm)
        V0         : float, fixed barrier height (eV)
        num_states : int, number of lowest energy levels to track

    Returns:
        a_arr      : np.ndarray, array of tested well widths
        energies   : np.ndarray, shape (len(a_values), num_states)
    """
    a_arr = np.array(a_values)
    energies = np.zeros((len(a_arr), num_states))

    for idx, a_val in enumerate(a_arr):
        evals, _, _, _ = solve_tise_numerical(
            N=4,
            a=a_val,
            b=b,
            V0=V0,
            num_states=num_states,
            points_per_cell=80
        )
        energies[idx, :len(evals)] = evals

    return a_arr, energies
