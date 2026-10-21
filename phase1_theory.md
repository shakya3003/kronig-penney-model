# Phase 1: Physics Foundations & Mathematical Theory (Clean & Intuitive Guide)

Welcome to **Phase 1**! The goal of this document is to explain the **Kronig–Penney Model** in plain English, step-by-step, with crystal-clear math and physical intuition.

---

## 1. The Big Picture: What Are We Trying to Solve?

Imagine an electron traveling through a crystal wire (like copper or silicon):
- The wire is not empty space. It is made of atomic nuclei arranged in a neat, repeating row with spacing $d$.
- Near each nucleus, the electron feels an attractive pull (a **potential well**, where potential energy $V = 0$).
- Between the nuclei, there is a repulsive barrier (a **potential barrier**, where potential energy $V = V_0$).
- As the electron flies down this 1D "egg carton" of repeating wells and barriers, how does it behave?

**The Big Discovery:**
Unlike a free electron (which can have *any* energy) or an isolated single atom (which only has *sharp, discrete* energy levels), an electron in a periodic crystal can only travel in continuous ranges of energies called **Allowed Energy Bands**, separated by forbidden zones called **Band Gaps** where electrons cannot exist.

---

## 2. Table of Symbols (Your Cheat Sheet)

Before looking at any equation, here is what every letter represents:

| Symbol | Name | What it Physically Means | Units |
| :--- | :--- | :--- | :--- |
| $x$ | Position | Coordinate along the 1D crystal wire | $\text{nm}$ (nanometers) |
| $a$ | Well Width | Width of the region near the atom where $V = 0$ | $\text{nm}$ |
| $b$ | Barrier Width | Thickness of the barrier separating atoms | $\text{nm}$ |
| $d = a + b$ | Lattice Constant | Repeat distance (period) of the crystal unit cell | $\text{nm}$ |
| $V_0$ | Barrier Height | Potential energy of the barrier | $\text{eV}$ (electronvolts) |
| $E$ | Electron Energy | Total kinetic + potential energy of the electron | $\text{eV}$ |
| $m$ | Electron Mass | Rest mass of the electron ($m_e$) | $\text{kg}$ |
| $\hbar$ | Reduced Planck Constant | Fundamental quantum constant ($h / 2\pi$) | $\text{J}\cdot\text{s}$ or $\text{eV}\cdot\text{s}$ |
| $\alpha$ | Well Wavevector | Spatial frequency of wave oscillation inside the well | $\text{nm}^{-1}$ |
| $\beta$ | Barrier Decay Rate | How fast the wave decays/tunnels inside the barrier | $\text{nm}^{-1}$ |
| $k$ | Bloch Wavevector | The crystal momentum ($p = \hbar k$) of the traveling electron | $\text{nm}^{-1}$ |
| $\psi(x)$ | Wavefunction | Quantum wave amplitude; its square $|\psi(x)|^2$ is probability | dimensionless |

---

## 3. The Potential Landscape $V(x)$

We divide space into repeating cells of length $d = a + b$:

```text
       Potential Energy V(x)
          ^
          |      Barrier Region               Barrier Region
          |      (V = V0, width b)            (V = V0, width b)
      V0  +---------+                 +---------+                 +---------+
          |         |                 |         |                 |         |
          |         |    Well Region  |         |    Well Region  |         |
          |         |  (V = 0, width a)         |  (V = 0, width a)         |
       0  +---------+--------+--------+---------+--------+--------+---------+--> Position x
          |<---b--->|<-------a------->|
          |<---------- d ------------>|
```

In words:
- **Region I (The Well)**: From $x = 0$ to $x = a$, potential energy is zero: $V(x) = 0$.
- **Region II (The Barrier)**: From $x = -b$ to $x = 0$, potential energy is high: $V(x) = V_0$.
- **Periodicity**: Every $d$ nanometers, the pattern repeats: $V(x + d) = V(x)$.

---

## 4. Step-by-Step Schrödinger Equation

The fundamental rule of quantum mechanics is the **Time-Independent Schrödinger Equation (TISE)**:

$$\text{(Kinetic Energy)} + \text{(Potential Energy)} = \text{Total Energy}$$

Mathematically:
$$-\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$

Let's solve this in the two regions separately:

### Step 4A: Inside the Well ($0 < x < a$, where $V = 0$)
Since $V(x) = 0$:
$$-\frac{\hbar^2}{2m} \frac{d^2\psi_I}{dx^2} = E\psi_I$$

Rearranging:
$$\frac{d^2\psi_I}{dx^2} + \left( \frac{2mE}{\hbar^2} \right) \psi_I = 0$$

To make this clean, we define:
$$\alpha \equiv \sqrt{\frac{2mE}{\hbar^2}}$$

Now the equation looks just like a simple vibrating guitar string:
$$\frac{d^2\psi_I}{dx^2} + \alpha^2 \psi_I = 0$$

**Solution in the well:**
$$\psi_I(x) = A \cos(\alpha x) + B \sin(\alpha x) \quad \text{or} \quad A e^{i\alpha x} + B e^{-i\alpha x}$$
> **Physical meaning**: Inside the well, the electron is an oscillating ripple with wavelength $\lambda = \frac{2\pi}{\alpha}$.

---

### Step 4B: Inside the Barrier ($-b < x < 0$, where $V = V_0$)
First consider when the electron's energy is below the barrier ($E < V_0$):
$$-\frac{\hbar^2}{2m} \frac{d^2\psi_{II}}{dx^2} + V_0 \psi_{II} = E\psi_{II}$$

Rearranging:
$$\frac{d^2\psi_{II}}{dx^2} - \left( \frac{2m(V_0 - E)}{\hbar^2} \right) \psi_{II} = 0$$

To make this clean, we define:
$$\beta \equiv \sqrt{\frac{2m(V_0 - E)}{\hbar^2}}$$

Now the equation becomes:
$$\frac{d^2\psi_{II}}{dx^2} - \beta^2 \psi_{II} = 0$$

**Solution in the barrier ($E < V_0$):**
$$\psi_{II}(x) = C e^{\beta x} + D e^{-\beta x} \quad \text{or} \quad C \cosh(\beta x) + D \sinh(\beta x)$$
> **Physical meaning**: Inside the barrier, classical physics says the electron is forbidden. But in quantum mechanics, the electron can **tunnel** into the barrier; its wavefunction decays exponentially at rate $\beta$.

*(Note: If $E > V_0$, the electron has enough energy to fly over the barrier! In that case, we define $\beta' = \sqrt{\frac{2m(E - V_0)}{\hbar^2}}$, and the wave in the barrier oscillates as $\cos(\beta' x)$ and $\sin(\beta' x)$).*

---

## 5. What is Bloch's Theorem? (Plain English Intuition)

Because the crystal repeats identically every distance $d$, the probability of finding an electron must also repeat identically:
$$|\psi(x + d)|^2 = |\psi(x)|^2$$

Felix Bloch showed that this means the wavefunction at $(x + d)$ is simply the wavefunction at $x$ multiplied by a phase angle factor $e^{ikd}$:
$$\psi(x + d) = e^{ikd} \psi(x)$$

Here, $k$ is called the **Bloch wavevector**. It represents the crystal momentum of the electron traveling down the lattice.

---

## 6. The 4 Boundary Conditions (The "Handshakes")

A real electron wave cannot have sudden tears or sharp kinks:
1. **No tears (Smoothness)**: The value of $\psi(x)$ must match at interfaces.
2. **No sharp kinks (Continuous Slope)**: The derivative $\frac{d\psi}{dx}$ must match at interfaces.

Let's write out the 4 matching points between the well and barrier:

1. **At boundary $x = 0$ (Wave value matches)**:
   $$\psi_I(0) = \psi_{II}(0) \implies A + B = C + D$$

2. **At boundary $x = 0$ (Slope matches)**:
   $$\psi'_I(0) = \psi'_{II}(0) \implies i\alpha (A - B) = \beta (C - D)$$

3. **At boundary $x = a$ connecting to next cell at $x = -b$ (Wave matches with Bloch phase)**:
   $$\psi_I(a) = e^{ikd} \psi_{II}(-b)$$

4. **At boundary $x = a$ connecting to next cell at $x = -b$ (Slope matches with Bloch phase)**:
   $$\psi'_I(a) = e^{ikd} \psi'_{II}(-b)$$

These 4 equations form a system of 4 linear equations for the 4 unknown coefficients $(A, B, C, D)$.

For a solution to exist (other than the trivial $\psi = 0$), the determinant of the coefficient matrix must equal zero:
$$\det \mathbf{M} = 0$$

---

## 7. The Master Kronig–Penney Equation

When you evaluate that determinant, all the algebra condenses into **one single formula**:

### When $E < V_0$ (Tunneling Regime):
$$f(E) = \cos(\alpha a) \cosh(\beta b) + \left( \frac{\beta^2 - \alpha^2}{2\alpha\beta} \right) \sin(\alpha a) \sinh(\beta b) = \cos(kd)$$

### When $E > V_0$ (Over-Barrier Regime):
$$f(E) = \cos(\alpha a) \cos(\beta' b) - \left( \frac{\alpha^2 + \beta'^2}{2\alpha\beta'} \right) \sin(\alpha a) \sin(\beta' b) = \cos(kd)$$

---

## 8. Why Do Energy Bands and Band Gaps Appear?

Look at the right side of the master equation:
$$f(E) = \cos(kd)$$

The cosine of any real angle can **only** produce numbers between $-1$ and $+1$:
$$-1 \le \cos(kd) \le 1$$

Therefore:

### 🟢 Case 1: Allowed Band (When $|f(E)| \le 1$)
- Because $-1 \le f(E) \le 1$, we can find a **real** wavevector:
  $$k = \frac{1}{d} \arccos(f(E))$$
- The electron travels freely across the entire crystal as a continuous Bloch wave.
- **Physical result**: The electron is allowed to have this energy!

### 🔴 Case 2: Forbidden Band Gap (When $|f(E)| > 1$)
- If $f(E) > +1$ or $f(E) < -1$, **no real angle $kd$ can have that cosine!**
- Mathematically, $k$ becomes a complex number: $k = \text{real} + i\kappa$.
- The wave factor becomes $e^{ikx} = e^{i k_{\text{real}} x} \cdot e^{-\kappa x}$, which **dies out exponentially**.
- **Physical result**: The electron wave reflects back destructively (Bragg reflection). It cannot travel through the crystal at this energy!

```text
       f(E)
        ^
        |              FORBIDDEN REGION (f(E) > 1)
     +1 + - - - - - - - - - - - - - - - - - - - - - - - - - - -
        |        /\                  /\                  /\
        |       /  \                /  \                /  \
      0 +------+----+--------------+----+--------------+----+-------> Energy E
        |     /      \            /      \            /      \
        |    /        \          /        \          /        \
     -1 + - - - - - - - - - - - - - - - - - - - - - - - - - - -
        |              FORBIDDEN REGION (f(E) < -1)
        |
        |  [ALLOWED]   [BAND GAP]   [ALLOWED]   [BAND GAP]
        |    BAND 1                    BAND 2
```

---

## 9. Concrete Numerical Walkthrough

Let's plug in real numbers so you see exactly how the formula works.

Take:
- Well width: $a = 0.50\text{ nm}$
- Barrier width: $b = 0.10\text{ nm}$
- Lattice constant: $d = a + b = 0.60\text{ nm}$
- Barrier height: $V_0 = 10.0\text{ eV}$
- Fundamental constant: $\frac{\hbar^2}{2m_e} \approx 0.0762\text{ eV}\cdot\text{nm}^2$

### Test 1: Let's test $E = 0.50\text{ eV}$
1. Calculate well frequency $\alpha$:
   $$\alpha = \sqrt{\frac{0.50}{0.0762}} = 2.56\text{ nm}^{-1}$$
2. Calculate barrier decay $\beta$:
   $$\beta = \sqrt{\frac{10.0 - 0.50}{0.0762}} = 11.16\text{ nm}^{-1}$$
3. Compute the terms:
   - $\cos(\alpha a) = \cos(2.56 \times 0.50) = \cos(1.28) \approx 0.287$
   - $\cosh(\beta b) = \cosh(11.16 \times 0.10) = \cosh(1.116) \approx 1.688$
   - Impedance factor $\frac{\beta^2 - \alpha^2}{2\alpha\beta} = \frac{124.5 - 6.55}{2 \times 2.56 \times 11.16} \approx 2.06$
   - $\sin(\alpha a) \sinh(\beta b) = \sin(1.28) \times \sinh(1.116) \approx 0.958 \times 1.362 \approx 1.305$
4. Plug into $f(E)$:
   $$f(E) = (0.287 \times 1.688) + (2.06 \times 1.305) = 0.484 + 2.69 = \mathbf{3.18}$$
5. **Verdict**: $3.18 > 1$. Therefore, **$E = 0.50\text{ eV}$ is FORBIDDEN (in a band gap)!**

---

### Test 2: Let's test $E = 1.20\text{ eV}$
1. Calculate $\alpha = \sqrt{1.20 / 0.0762} = 3.97\text{ nm}^{-1}$
2. Calculate $\beta = \sqrt{(10.0 - 1.20) / 0.0762} = 10.75\text{ nm}^{-1}$
3. Evaluating $f(E)$ gives:
   $$f(E) = \mathbf{0.729}$$
4. **Verdict**: Since $|0.729| \le 1$, this energy is **ALLOWED**!
5. The electron travels with crystal wavevector:
   $$k = \frac{1}{0.60} \arccos(0.729) = \mathbf{1.26\text{ nm}^{-1}}$$

---

## 10. Why Does This Matter in the Real World?

| Phenomenon | How the Kronig–Penney Model Explains It |
| :--- | :--- |
| **Metals vs. Insulators** | If the highest electrons fill an allowed band partially, they can accelerate easily $\implies$ **Metal**. If a band is totally full and separated from the next empty band by a wide forbidden gap, electrons cannot jump $\implies$ **Insulator**. |
| **Semiconductors (Silicon)** | The gap between the valence band and conduction band is small ($\sim 1.1\text{ eV}$ in Silicon), so thermal energy or photons can bump electrons across $\implies$ **Transistors, Solar Cells, LEDs**. |
| **Effective Mass ($m^*$)** | The curve of $E(k)$ tells us how heavy an electron feels inside the crystal: $m^* = \hbar^2 / (d^2E/dk^2)$. In some materials like GaAs, electrons move faster than in vacuum because $m^* \approx 0.067 m_0$! |

---

## Summary of Phase 1
You now know:
1. What the potential landscape $V(x)$ looks like.
2. Why the wave oscillates in the well ($\alpha$) and tunnels in the barrier ($\beta$).
3. How Bloch's theorem provides the periodic boundary phase $e^{ikd}$.
4. How the master equation $f(E) = \cos(kd)$ creates allowed bands and forbidden gaps.

With the theory crystal clear, we are ready to move to **Phase 2** (building the Python physics engine `kp_engine.py`)!
