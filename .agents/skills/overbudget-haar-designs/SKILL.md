---
name: overbudget-haar-designs
description: >-
  Canonical mathematical definitions, theorems, Delsarte LP certificates, and
  zero-dependency exact-arithmetic verification workflow for Overbudget
  Haar-Purity Subsystem Designs Beyond the Quadratic Threshold.
---

# Overbudget Haar-Purity Subsystem Designs Skill

Author: Thijs Missiaen (ORCID: `0009-0007-3678-0810`)

This skill equips an AI agent with the exact mathematical setup, theorem structure, and standalone verification workflow for complex-Haar subsystem-purity variance minimization on $\overline{J}(2k,k)$ beyond the Hadamard budget ($m > n$).

## 1. Canonical Objects and Normalizations

Let $n = 2k \ge 6$ be even, $d \ge 2$ be an integer local dimension, and $D = d^n$. Balanced bipartitions $\lbrace A, A^c \rbrace$ of $\lbrack n \rbrack$ ($\lvert A \rvert = k$) form the folded Johnson scheme $\overline{J}(2k, k)$ of size $M_n = \frac{1}{2}\binom{2k}{k}$, canonically indexed by $k$-subsets containing $0$. For canonical representatives $A_a, A_b$, define $q _ {ab} := \lvert A_a \cap A_b \rvert$ and overlap discrepancy $g _ {ab} := 4 q _ {ab} - n$.

The Haar purity kernel and objective over a family $\mathcal{F} \subseteq \overline{J}(2k, k)$ of size $m$ are given by:

$$
h_{n,d}(g) = d^{(n-g)/2} + d^{(n+g)/2}, \qquad S_{n,d}(\mathcal{F}) = \sum_{1 \le a < b \le m} h_{n,d}(g_{ab}), \qquad Q(\mathcal{F}) = \sum_{1 \le a < b \le m} g_{ab}^2.
$$

Lemma 1 establishes the exact linear relation between complex-Haar empirical variance and $S _ {n,d}(\mathcal{F})$:

$$
\mathrm{Var}_{\psi \sim \mathrm{Haar}}\bigl(\overline{P}_{\mathcal{F}}(\psi)\bigr) = \frac{2}{m^2 (D+1)(D+2)(D+3)} \Bigl( m(D+1) + 2 S_{n,d}(\mathcal{F}) - \frac{4 m^2 D}{D+1} \Bigr).
$$

## 2. Core Theorems and Structural Dichotomy

**1. Alphabet Dichotomy and Separation.** For $n \in \lbrace 6, 8, 10 \rbrace$, the off-diagonal alphabet $\lbrace \lvert g _ {ab} \rvert : a \ne b \rbrace$ has size at most $2$, so $h _ {n,d}(g) = \alpha _ {n,d} + \beta _ {n,d} g^2$ holds exactly ($\beta _ {n,d} > 0$) and $E(s^2)$ agrees identically with Haar optimization. For every even $n \ge 12$, there are $\lfloor k/2 \rfloor \ge 3$ off-diagonal levels and $x \mapsto h _ {n,d}(\sqrt{x})$ is strictly convex, breaking quadratic minorization. Bitmask pairs $(\mathcal{F} _ Q, \mathcal{F} _ H)$ at $n=12$ ($m \in \lbrace 20, 40, 50, 60, 140, 220 \rbrace{}$) and $n=14$ ($m \in \lbrace 40, 80, 140 \rbrace{}$) achieve $Q(\mathcal{F} _ Q) \le Q(\mathcal{F} _ H)$ alongside $S _ {n,d}(\mathcal{F} _ Q) > S _ {n,d}(\mathcal{F} _ H)$ for every integer $d \ge 2$ (peaking at ratio $8.4373$ for $n=12, m=60, d=10$).

**2. Folded Johnson Eigenvalue Formula.** For all even $n = 2k \ge 4$, real $d > 1$, and even $0 \le i \le k$, the eigenvalue $\Lambda_i(n,d)$ of $H _ {\mathrm{fold}}$ on eigenspace $V_i$ satisfies:

$$
\Lambda_i(n,d) = \sum_{r=i}^{k} (d^2 - 1)^r \binom{k-i}{k-r} \binom{2k-r-i}{k-r},
$$

giving strict monotonicity $\Lambda_0 > \Lambda_2 > \dots > \Lambda _ {2\lfloor k/2 \rfloor} > 0$ and recovering the $E(s^2)$ bound $Q(\mathcal{F}) \ge \frac{n^2}{2(n-1)} m (m - n + 1)$ from $Q_2(\delta) = \frac{(n-1)g _ {\delta}^2 - n^2}{n^2(n-2)}$.

**3. Witt Phase Transition at n = 12.** In $\overline{J}(12,6)$ ($M _ {12} = 462$), any family $\mathcal{F}$ of size $m$ obeys $c_8(\mathcal{F}) \ge \max\bigl(0, \frac{m(m-66)}{22}\bigr)$, saturated at $m=66$ uniquely by $W_1 = S(5,6,12)/\mathbb{Z}_2$ derived from the extended ternary Golay code $G _ {12}$. Explicit polynomial dual certificates in $t = d^2 - 1 > 0$ certify dimension-uniform global optimality for all integers $d \ge 2$ at $m \in \lbrace 11, 22, 55, 66, 132, 330, 396, 462 \rbrace$.

**4. Overbudget Parseval Sandwich.** Defining the row sum $R _ {n,d} = \sum _ {q=0}^k \binom{k}{q}^2 d^{2q} = \Lambda_0(n,d)$, exact bounds bracket $B _ {n,d}(m) = \min _ {\lvert\mathcal{F}\rvert = m} S _ {n,d}(\mathcal{F})$ for every $n = 2k \ge 4$, $d \ge 2$, and $1 \le m \le M_n$:

$$
\max\Bigl(0, \frac{m}{2}\Bigl(\frac{m}{M_n} R_{n,d} - (d^n + 1)\Bigr)\Bigr) \le B_{n,d}(m) \le \frac{m(m-1)}{2(M_n - 1)} (R_{n,d} - d^n - 1).
$$

Whenever $m / M_n \gg \sqrt{n} (d/(d+1))^n$, the relative sandwich gap vanishes as $o(1)$.

## 3. Standalone Verification and Build Commands

Run the standard-library exact-arithmetic verifier from the repository root:

```powershell
python scripts/verify_overbudget_designs.py
```

Regenerate exact figures (`fig1_es2_vs_haar_n12.{pdf,png}` and `fig2_asymptotic_collapse.{pdf,png}`):

```powershell
python scripts/generate_figures.py
```

Compile the AMS-LaTeX manuscript cleanly via Tectonic/pdfTeX checks and check academic rigor:

```powershell
python C:\Users\Thijs\.gemini\config\skills\scientific-publisher\scripts\build_academic_pdf.py paper\overbudget_haar_designs.tex
python C:\Users\Thijs\.gemini\config\skills\scientific-publisher\scripts\anti_crank_linter.py paper\overbudget_haar_designs.tex
```
