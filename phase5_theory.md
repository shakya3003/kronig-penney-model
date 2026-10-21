# Phase 5: Numerical Solution of the Differential Equation & Parameter Study vs. Well Width $a$

Welcome to **Phase 5**! This document directly addresses **Question 6** of the assignment:
> *"Solve the relevant differential equation numerically and plot the required function with respect to the specified variable a, explaining the obtained behaviour."*

---

## 1. The Relevant Differential Equation

The equation governing the electron is the 1D **Time-Independent Schrödinger Equation (TISE)**:

$$-\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$

where:
- $\hbar = \frac{h}{2\pi}$ is the reduced Planck constant,
- $m = m_e$ is the electron mass ($\frac{\hbar^2}{2m_e} \approx 0.0762\text{ eV}\cdot\text{nm}^2$),
- $V(x)$ is the periodic potential landscape defined by well width $a$, barrier width $b$, and barrier height $V_0$:
  $$V(x) = \begin{cases} 0, & \text{inside wells (width } a \text{)} \\ V_0, & \text{inside barriers (width } b \text{)} \end{cases}$$
- $E$ is the energy eigenvalue, and $\psi(x)$ is the eigen-wavefunction.

---

## 2. Numerical Solution via Matrix Diagonalization (Finite Difference)

Rather than assuming any analytical ansatz (like $\cos$ or $\sin$), we solve the differential equation **purely numerically**:

1. **Discretize Space**: We define a discrete spatial grid $x_i = i \cdot \Delta x$ with fine step size $\Delta x \approx 0.005\text{ nm}$.
2. **Second Derivative Discretization**:
   $$\frac{d^2\psi}{dx^2}\Big|_{x_i} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\Delta x^2}$$
3. **The Matrix Equation**:
   Substituting into the TISE yields:
   $$\left(-\frac{\hbar^2}{2m \Delta x^2}\right)\psi_{i-1} + \left( \frac{\hbar^2}{m \Delta x^2} + V(x_i) \right)\psi_i + \left(-\frac{\hbar^2}{2m \Delta x^2}\right)\psi_{i+1} = E \psi_i$$
   This forms a tridiagonal eigenvalue equation:
   $$\mathbf{H} \vec{\psi} = E \vec{\psi}$$
4. **Diagonalization**: We compute the eigenvalues $E_n$ and eigen-vectors $\vec{\psi}_n$ using the high-performance symmetric tridiagonal solver `scipy.linalg.eigh_tridiagonal`.

---

## 3. Parameter Study: What Happens as Well Width $a$ Varies?

We hold the barrier parameters fixed ($b = 0.10\text{ nm}, V_0 = 10.0\text{ eV}$) and systematically vary the well width $a$ from narrow quantum wells ($a = 0.15\text{ nm}$) to wide wells ($a = 1.50\text{ nm}$).

### Key Physical Behaviors Observed:

#### 1. Quantum Confinement Effect ($E \propto 1/a^2$)
- In an infinite potential well, energy levels scale strictly as:
  $$E_n = \frac{n^2 \pi^2 \hbar^2}{2m a^2}$$
- In our periodic Kronig–Penney lattice, as $a$ increases:
  - The electron has more spatial room to spread out.
  - Its kinetic energy $\langle T \rangle \sim \frac{\hbar^2}{2m} \left(\frac{\pi}{a}\right)^2$ decreases rapidly.
  - Consequently, **all energy eigenvalues shift downward toward lower energies as $1/a^2$**.

#### 2. Proliferation of Allowed Bound States
- When the well is narrow ($a = 0.15\text{ nm}$), only 1 allowed band fits below the barrier height $V_0 = 10\text{ eV}$.
- As $a$ widens to $1.0\text{ nm}$ and beyond, the well can accommodate 4, 5, or more allowed bands below $V_0$.
- The energy spacing between adjacent bands ($\Delta E = E_{n+1} - E_n$) shrinks substantially.

#### 3. Band Broadening vs. Band Narrowing
- **At small $a$**:
  - The ground state energy $E_1$ sits high up, close to the barrier top $V_0$.
  - The decay parameter $\beta = \sqrt{2m(V_0 - E)}/\hbar$ is small, meaning the tunneling penetration depth $\delta = 1/\beta$ is large.
  - Neighboring wells strongly overlap (strong inter-well coupling).
  - **Result**: The energy bands broaden into wide intervals.
- **At large $a$**:
  - The lowest energy states drop deep to the bottom of the well ($E \ll V_0$).
  - The barrier decay parameter $\beta$ becomes very large, and the penetration depth shrinks to a tiny fraction of a nanometer.
  - The electron wavefunctions are tightly trapped inside each individual well and cannot easily tunnel to neighbors.
  - **Result**: The bands become extremely narrow, approaching the discrete, razor-sharp delta-like levels of isolated atoms (the **Tight-Binding Limit**).

#### 4. Behavior of the Forbidden Band Gaps
- At large $a$, because the bands sink and compress, the band gaps also become narrower in absolute energy spacing.
- In the limit $a \to \infty$, the discrete bands merge into a continuous semi-classical spectrum (Bohr's correspondence principle).

---

## 4. Summary Table of the $a$-Dependence

| Well Width $a$ | Quantum Confinement | Level Spacing | Inter-well Tunneling | Band Width | Physical Regime |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Small ($< 0.3\text{ nm}$)** | Extremely strong | Very large | High (large overlap) | Broad bands | Nearly Free Electron / Quantum Dot |
| **Medium ($0.4 - 0.7\text{ nm}$)** | Moderate | Medium | Balanced | Moderate bands & gaps | Typical Semiconductor Superlattice |
| **Large ($> 1.0\text{ nm}$)** | Weak | Small | Very low (isolated wells) | Razor-thin bands | Tight-Binding / Isolated Atoms |

---

We will now execute `numerical_analysis.py` to calculate these eigenvalues numerically, plot $E$ vs. $a$, and generate the visual proof.
