# Phase 2: Numerical Engine & Computational Physics (Clean & Intuitive Guide)

Welcome to **Phase 2**! In Phase 1, we established the analytical pen-and-paper formula for the Kronig–Penney model.

Now, in Phase 2, we answer: **How do we make a computer solve this equation, compute the wavefunctions across multiple wells, and solve the differential equation numerically?**

This document explains the 3 core computational algorithms we will implement in `kp_engine.py`:
1. **The Choice of Physical Units** (Why computers hate SI units and love nanoscale units).
2. **The Band Edge Root-Finder** (How we automatically find all allowed energy bands $[E_{\min}, E_{\max}]$).
3. **The Transfer Matrix Method (TMM)** (How we compute the exact continuous wavefunction $\psi(x)$ across $N$ wells).
4. **The Finite-Difference Differential Equation Solver** (How we solve the 1D Schrödinger equation directly as a matrix eigenvalue problem, directly fulfilling Question 6).

---

## 1. The Choice of Units: Eliminating Underflow & Overflow

In standard SI units:
- Planck's constant: $\hbar \approx 1.055 \times 10^{-34}\text{ J}\cdot\text{s}$
- Electron mass: $m_e \approx 9.109 \times 10^{-31}\text{ kg}$
- Electron charge: $e \approx 1.602 \times 10^{-19}\text{ C}$

Notice the exponents: $10^{-34}$ squared is $10^{-68}$! When a 64-bit floating-point number tries to divide $10^{-68}$ by $10^{-31}$, floating-point roundoff errors and underflows ruin numerical precision.

### The Solution: Nanoscale Solid-State Units
Instead of meters and Joules, we choose units natural to solid-state physics:
- **Length**: nanometers ($\text{nm}$), where $1\text{ nm} = 10^{-9}\text{ m}$.
- **Energy**: electronvolts ($\text{eV}$), where $1\text{ eV} = 1.6021766 \times 10^{-19}\text{ J}$.
- **Mass**: electron rest mass ($m_e = 1.0$).

In these friendly units, the quantum kinetic constant becomes a clean, well-behaved number:
$$\gamma \equiv \frac{\hbar^2}{2m_e} \approx 0.0761996 \text{ eV}\cdot\text{nm}^2$$

Now:
- In the well: $\alpha = \sqrt{\frac{E}{\gamma}} = \sqrt{\frac{E}{0.0762}} \quad (\text{units: }\text{nm}^{-1})$
- In the barrier: $\beta = \sqrt{\frac{V_0 - E}{\gamma}} = \sqrt{\frac{V_0 - E}{0.0762}} \quad (\text{units: }\text{nm}^{-1})$

All numbers are now between $0.1$ and $100$—perfect for floating-point calculations!

---

## 2. Algorithm 1: Finding Allowed Energy Bands Automatically

In Phase 1, we learned that an energy $E$ is **allowed** if and only if:
$$-1 \le f(E) \le 1$$

And the band edges (the exact transition energies between an allowed band and a forbidden gap) occur where:
$$f(E) = +1 \quad \text{or} \quad f(E) = -1$$

```text
       f(E) ^
            |      Band 1 Gap         Band 2 Gap
         +1 +----+ - - - - - - +----+ - - - - - - -
            |   / \           /      \
            |  /   \         /        \
          0 +-+-----+-------+----------+----------> Energy E
            |/       \     /            \
         -1 + - - - - +---+ - - - - - - -+--------
            |
            | [Band 1]     [Band 2]      [Band 3]
            | E_min E_max  E_min  E_max
```

### How the Engine Finds Bands:
1. We sample $E$ over a fine energy grid from $E = 0.01\text{ eV}$ to $E_{\max}$.
2. For each point, we evaluate $f(E)$ and check whether $|f(E)| \le 1$.
3. Whenever $f(E)$ crosses $+1$ or $-1$, a band edge has been crossed!
4. We use **Brent's method (bisection root-finding)** on $g(E) = f(E) \pm 1 = 0$ to find the exact boundary energy to 8 decimal places.
5. The result is a list of clean band intervals:
   $$\text{Band 1: } [E_{\min, 1}, E_{\max, 1}], \quad \text{Band 2: } [E_{\min, 2}, E_{\max, 2}], \quad \dots$$

---

## 3. Algorithm 2: The Transfer Matrix Method (TMM) for Wavefunctions

How do we calculate the actual continuous wavefunction $\psi(x)$ across a finite structure of $N$ potential wells?

### The State Vector
At any point $x$, the quantum state of the electron is completely determined by two numbers:
$$\mathbf{w}(x) = \begin{pmatrix} \psi(x) \\ \psi'(x) \end{pmatrix}$$
where $\psi'(x) = \frac{d\psi}{dx}$ is the spatial derivative.

### The Propagation Matrix
Inside any region of constant potential $V$:
- **Inside a Well ($V=0$, wavenumber $\alpha$):**
  If we start at $x_0$ with state $\mathbf{w}(x_0)$, after traveling a distance $\Delta x = x - x_0$, the state at $x$ is:
  $$\begin{pmatrix} \psi(x) \\ \psi'(x) \end{pmatrix} = \begin{pmatrix} \cos(\alpha \Delta x) & \frac{1}{\alpha}\sin(\alpha \Delta x) \\ -\alpha \sin(\alpha \Delta x) & \cos(\alpha \Delta x) \end{pmatrix} \begin{pmatrix} \psi(x_0) \\ \psi'(x_0) \end{pmatrix}$$

- **Inside a Barrier ($V=V_0$, decay parameter $\beta$, for $E < V_0$):**
  $$\begin{pmatrix} \psi(x) \\ \psi'(x) \end{pmatrix} = \begin{pmatrix} \cosh(\beta \Delta x) & \frac{1}{\beta}\sinh(\beta \Delta x) \\ \beta \sinh(\beta \Delta x) & \cosh(\beta \Delta x) \end{pmatrix} \begin{pmatrix} \psi(x_0) \\ \psi'(x_0) \end{pmatrix}$$

- **Inside a Barrier ($E > V_0$, with $\beta' = \sqrt{(E - V_0)/\gamma}$):**
  $$\begin{pmatrix} \psi(x) \\ \psi'(x) \end{pmatrix} = \begin{pmatrix} \cos(\beta' \Delta x) & \frac{1}{\beta'}\sin(\beta' \Delta x) \\ -\beta' \sin(\beta' \Delta x) & \cos(\beta' \Delta x) \end{pmatrix} \begin{pmatrix} \psi(x_0) \\ \psi'(x_0) \end{pmatrix}$$

### Why This Method Is Brilliant:
1. **Guaranteed Smoothness**: Because matrix multiplication inherently preserves the continuity of both $\psi(x)$ and $\psi'(x)$, the wavefunction is mathematically guaranteed to have no tears and no kinks across all interfaces!
2. **Any Number of Wells $N$**: We simply start at the left boundary $x = 0$ with an initial state $\mathbf{w}(0)$, and multiply through each successive well and barrier across the entire multi-well chain.
3. **Normalization**: Once the full profile $\psi(x)$ is computed across the structure of total length $L$, we normalize it so the total probability equals 1:
   $$\int_0^L |\psi(x)|^2 dx = 1 \implies \psi_{\text{norm}}(x) = \frac{\psi(x)}{\sqrt{\sum_i |\psi(x_i)|^2 \Delta x}}$$

---

## 4. Algorithm 3: Numerical Differential Equation Solver (Question 6)

Question 6 specifically asks:
> *"Solve the relevant differential equation numerically and plot the required function with respect to the specified variable a, explaining the obtained behaviour."*

Instead of assuming analytical solutions, we solve the **1D Time-Independent Schrödinger Differential Equation** directly on a discrete spatial grid:

$$-\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$

### Step 4A: Spatial Discretization
We divide the entire length $L$ into $M$ discrete grid points:
$$x_0, x_1, x_2, \dots, x_{M-1}$$
with spacing $\Delta x = x_{i+1} - x_i$.

At each grid point, the potential is known: $V_i = V(x_i)$.

### Step 4B: Finite Difference Approximation
Using Taylor expansion, the second derivative of $\psi(x)$ at point $i$ is:
$$\frac{d^2\psi}{dx^2}\Big|_{x_i} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{(\Delta x)^2}$$

Substitute this into the Schrödinger equation:
$$-\gamma \left( \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{(\Delta x)^2} \right) + V_i \psi_i = E \psi_i$$

Rearranging term-by-term:
$$\left(-\frac{\gamma}{(\Delta x)^2}\right)\psi_{i-1} + \left( \frac{2\gamma}{(\Delta x)^2} + V_i \right)\psi_i + \left(-\frac{\gamma}{(\Delta x)^2}\right)\psi_{i+1} = E \psi_i$$

### Step 4C: The Hamiltonian Matrix Equation
This set of equations is simply a **matrix eigenvalue problem**:
$$\mathbf{H} \vec{\psi} = E \vec{\psi}$$

where $\mathbf{H}$ is a symmetric, tridiagonal matrix:
$$\mathbf{H} = \begin{pmatrix}
\frac{2\gamma}{\Delta x^2} + V_0 & -\frac{\gamma}{\Delta x^2} & 0 & \dots & 0 \\
-\frac{\gamma}{\Delta x^2} & \frac{2\gamma}{\Delta x^2} + V_1 & -\frac{\gamma}{\Delta x^2} & \dots & 0 \\
0 & -\frac{\gamma}{\Delta x^2} & \frac{2\gamma}{\Delta x^2} + V_2 & \dots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & -\frac{\gamma}{\Delta x^2} & \frac{2\gamma}{\Delta x^2} + V_{M-1}
\end{pmatrix}$$

### Step 4D: Solving with `scipy.linalg.eigh_tridiagonal`
Because $\mathbf{H}$ is tridiagonal, standard numerical linear algebra algorithms can find its eigenvalues (energies $E_n$) and eigenvectors (wavefunctions $\vec{\psi}_n$) in **milliseconds** with $O(M)$ time complexity!

This provides an independent, purely numerical solution of the differential equation, giving us a perfect direct comparison against the analytical Kronig–Penney model.

---

## 5. Summary of What `kp_engine.py` Will Provide

| Function Name | Inputs | Outputs | Description |
| :--- | :--- | :--- | :--- |
| `compute_dispersion(E, a, b, V0)` | $E, a, b, V_0$ | $f(E), \text{is\_allowed}, k, \kappa$ | Evaluates master equation at energy $E$. |
| `find_allowed_bands(a, b, V0, E_max)` | $a, b, V_0, E_{\max}$ | List of `(E_min, E_max, band_index)` | Automatically locates all allowed band intervals. |
| `generate_potential_profile(N, a, b, V0)` | $N, a, b, V_0$ | $x, V(x)$ | Creates the spatial grid and potential landscape. |
| `compute_wavefunction_tmm(N, E, a, b, V0)` | $N, E, a, b, V_0$ | $x, \psi(x), |\psi(x)|^2$ | Transfer Matrix Method for continuous $\psi(x)$. |
| `solve_tise_numerical(N, a, b, V0, num_states)` | $N, a, b, V_0$ | $E_n, \psi_n(x)$ | Solves differential equation via tridiagonal finite-difference. |

With these four pillars explained, we will now code `kp_engine.py` and test it with a dedicated verification script `verify_phase2.py`.
