# Phase 4: Wavefunction Physics, Probability Density & Visualization Theory

Welcome to **Phase 4**! In this phase, we dive deep into the physical interpretation of the **wavefunction $\psi(x)$** and the **probability density $|\psi(x)|^2$**, addressing **Question 4** and **Question 5** of the assignment.

---

## 1. Physical Anatomy of the Wavefunction $\psi(x)$

The 1D Time-Independent Schrödinger Equation can be rewritten as:
$$\frac{d^2\psi(x)}{dx^2} = -\frac{2m}{\hbar^2}\Big[ E - V(x) \Big] \psi(x)$$

The second derivative $\frac{d^2\psi}{dx^2}$ represents the **spatial curvature** of the wavefunction. This yields two completely different physical regimes:

```text
                        Curvature of psi(x)
       =========================================================
       Regime              E vs V(x)       Curvature Behavior
       ---------------------------------------------------------
       Inside Wells        E > V(x) = 0    d²psi/dx² has OPPOSITE sign to psi
                                           --> Curves toward axis (Oscillatory)

       Inside Barriers     E < V(x) = V0   d²psi/dx² has SAME sign to psi
                                           --> Curves away from axis (Exponential/Tunneling)
       =========================================================
```

### 1.1 Inside the Wells ($V = 0$): Kinetic Energy & Wavelength
- In the well, the electron has positive kinetic energy $T = E - 0 = E$.
- The wavevector is $\alpha = \sqrt{\frac{2mE}{\hbar^2}}$.
- The de Broglie wavelength is:
  $$\lambda = \frac{2\pi}{\alpha} = \frac{h}{\sqrt{2mE}}$$
- **Visual feature**: The wave behaves like a sinusoidal ripple. Higher energy $E$ means shorter wavelength $\lambda$ and more rapid oscillations (more nodes).

### 1.2 Inside the Barriers ($V = V_0$): Quantum Tunneling & Penetration Depth
- In classical mechanics, an electron with $E < V_0$ hitting a barrier of height $V_0$ bounces off with 100% certainty; it can **never** enter the barrier.
- In quantum mechanics, the wavefunction does **not** drop to zero instantly at the barrier edge. Instead, it enters the barrier as an **evanescent wave**:
  $$\psi(x) \sim e^{-\beta x}$$
  where $\beta = \sqrt{\frac{2m(V_0 - E)}{\hbar^2}}$.
- **Penetration Depth ($\delta$)**:
  The characteristic distance over which the wave amplitude drops by a factor of $1/e \approx 37\%$ is:
  $$\delta = \frac{1}{\beta} = \frac{\hbar}{\sqrt{2m(V_0 - E)}}$$
- If the barrier thickness $b$ is comparable to or smaller than $\delta$ ($b \lesssim \delta$), the decaying tail reaches the next well with a non-zero amplitude! This is the physical mechanism of **quantum tunneling** between adjacent atomic sites.

---

## 2. Physical Significance of Probability Density $|\psi(x)|^2$

According to Max Born's statistical interpretation of quantum mechanics:
$$P(x_1 \le x \le x_2) = \int_{x_1}^{x_2} |\psi(x)|^2 dx$$

$|\psi(x)|^2$ represents the **probability per unit length** of finding the electron at position $x$.

### 2.1 Normalization Condition
Since the electron must exist somewhere in the total physical system of length $L$:
$$\int_0^L |\psi(x)|^2 dx = 1.0 \quad (100\%)$$

### 2.2 What $|\psi(x)|^2$ Tells Us in a Crystal:
1. **Charge Distribution**: In a real solid, multiplying $|\psi(x)|^2$ by the elementary charge $(-e)$ gives the electronic charge density $\rho(x) = -e |\psi(x)|^2$.
2. **Localization vs. Delocalization**:
   - In tight-binding states (deep wells / thick barriers), $|\psi(x)|^2$ is strongly localized inside the wells with almost zero probability inside the barriers.
   - In nearly free states (shallow wells or energies near $V_0$), $|\psi(x)|^2$ spreads smoothly across the whole lattice.
3. **Tunneling Leakage**: The non-zero shaded area of $|\psi(x)|^2$ inside the barrier regions directly proves that electrons spend a measurable fraction of their time *inside* the classically forbidden barriers!

---

## 3. Allowed States vs. Forbidden Gap States: Visual Comparison

| Visual Characteristic | Allowed Band State ($|f(E)| \le 1$) | Forbidden Gap State ($|f(E)| > 1$) |
| :--- | :--- | :--- |
| **Propagation** | Extended Bloch wave; travels across all $N$ wells. | Evanescent wave; strongly reflected. |
| **Amplitude Envelope** | Uniform / balanced across the entire multi-well chain. | Decays exponentially across successive wells ($e^{-\kappa x}$). |
| **Phase Shift per Cell** | Pure phase shift: $\psi(x + d) = e^{ikd}\psi(x)$ ($k \in \mathbb{R}$). | Growth/decay factor: $\psi(x + d) = e^{-\kappa d}\psi(x)$. |
| **Transmission** | Electron freely transmits through the crystal wire. | Electron is reflected back (Bragg backscattering). |

---

## 4. Origin of the Band Gap at Brillouin Zone Boundaries ($k = \pm \pi/d$)

Why does a band gap open right at $k = \pi/d$?
- When $k = \pi/d$, the electron de Broglie wavelength satisfies $\lambda = 2d$, exactly matching the Bragg diffraction condition!
- Forward-traveling waves $e^{ikx}$ and backward-traveling waves $e^{-ikx}$ interfere to form two distinct standing waves:
  $$\psi_+(x) \sim \cos\left(\frac{\pi x}{d}\right) \quad \text{(Symmetric state)}$$
  $$\psi_-(x) \sim \sin\left(\frac{\pi x}{d}\right) \quad \text{(Antisymmetric state)}$$
- **Energy Difference ($\Delta E_g$)**:
  - In $\psi_+(x)$, the electron probability peaks **inside the wells** (closer to the positive ions $\implies$ lower potential energy).
  - In $\psi_-(x)$, the probability peaks **inside the barriers** (farther from the wells $\implies$ higher potential energy).
- This energy difference $\Delta E_g = \langle V \rangle_- - \langle V \rangle_+$ creates the **forbidden band gap** where no propagating states can exist!

---

## 5. Summary of Phase 4 Goals
To provide a complete visual demonstration for your project, we will now create `plot_generator.py` to generate publication-quality figures comparing:
1. An **Allowed State** ($E = 1.20\text{ eV}$) showing extended Bloch transmission across all $N$ wells.
2. A **Forbidden State** ($E = 3.20\text{ eV}$) showing exponential decay and quantum attenuation.
3. The **$E$-$k$ Band Diagram** within the first Brillouin zone.
