# Overbudget Haar-Purity Subsystem Designs Beyond the Quadratic Threshold

**Author:** Thijs Missiaen (Independent Researcher, Bruges, Belgium)  
**ORCID:** [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)  
**DOI:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24244470.svg)](https://doi.org/10.5281/zenodo.24244470)  
**Manuscript:** [paper/overbudget_haar_designs.pdf](paper/overbudget_haar_designs.pdf) (AMS-LaTeX source: [paper/overbudget_haar_designs.tex](paper/overbudget_haar_designs.tex), 11 pages)

## Overview

For an even number $n = 2k$ of $d$-dimensional quantum subsystems, minimizing the complex-Haar variance of the empirical cut-averaged subsystem purity $\overline{P} _ {\mathcal{F}}(\psi) = \frac{1}{m}\sum _ {a=1}^m \mathrm{Tr}(\rho _ {A_a}^2)$ over a family $\mathcal{F}$ of $m$ distinct balanced bipartitions reduces identically to minimizing the exponential overlap sum $S _ {n,d}(\mathcal{F}) = \sum _ {a < b} \bigl(d^{(n-g _ {ab})/2} + d^{(n+g _ {ab})/2}\bigr)$ on the folded Johnson association scheme $\overline{J}(2k, k)$ of size $M_n = \frac{1}{2}\binom{2k}{k}$.

Previous work by Shiroshita (`10.5281/zenodo.22773593`, 2026) resolved the square Hadamard budget $m = n$ for $n \in \lbrace 6, 8, 10, 12, 14, 16 \rbrace$ using a $\lbrace 0, 4 \rbrace$ quadratic minorant. This repository resolves the **overbudget regime** ($m > n$ across arbitrary cut counts $1 \le m \le M_n$ and all integer local dimensions $d \ge 2$), where orthogonal cuts are exhausted and the square-budget minorant fails.

## Main Results

**1. Alphabet Dichotomy and E(s^2)-versus-Haar Separation (Theorem 1).** For $n \in \lbrace 6, 8, 10 \rbrace$, $\overline{J}(2k, k)$ has at most $2$ off-diagonal overlap levels $\lvert g \rvert$, so minimizing $S _ {n,d}(\mathcal{F})$ is identically equivalent to minimizing the classical quadratic objective $Q(\mathcal{F}) = \sum _ {a<b} g _ {ab}^2$ for all $1 \le m \le M_n$ and all $d \ge 2$. For every even $n \ge 12$, there are $\lfloor k/2 \rfloor \ge 3$ off-diagonal levels and strict convexity of $x \mapsto h _ {n,d}(\sqrt{x})$ breaks quadratic minorization. Explicit bitmask certificates at $n=12$ ($m \in \lbrace 20, 40, 50, 60, 140, 220 \rbrace{}$) and $n=14$ ($m \in \lbrace 40, 80, 140 \rbrace{}$) satisfy $Q(\mathcal{F} _ Q) \le Q(\mathcal{F} _ H)$ alongside $S _ {n,d}(\mathcal{F} _ Q) > S _ {n,d}(\mathcal{F} _ H)$ for every integer $d \ge 2$, reaching ratio $S _ {12,10}(\mathcal{F} _ Q) / S _ {12,10}(\mathcal{F} _ H) = 8.4373$ at $m=60$.

**2. Folded Johnson Spectral Formula (Theorem 2).** For every even $n = 2k \ge 4$, real $d > 1$, and even $0 \le i \le k$, the eigenvalue of the folded Haar kernel $H _ {\mathrm{fold}}$ on eigenspace $V_i$ is:

$$
\Lambda_i(n,d) = \sum_{r=i}^{k} (d^2 - 1)^r \binom{k-i}{k-r} \binom{2k-r-i}{k-r}.
$$

This closed form establishes strict spectral monotonicity $\Lambda_0 > \Lambda_2 > \dots > \Lambda _ {2\lfloor k/2 \rfloor} > 0$ and derives the folded $E(s^2)$ floor $Q(\mathcal{F}) \ge \frac{n^2}{2(n-1)} m (m - n + 1)$ from $Q_2(\delta) = \frac{(n-1)g _ {\delta}^2 - n^2}{n^2(n-2)}$.

**3. The Witt S(5,6,12)/Z_2 Phase Transition at n = 12 (Theorem 3).** Delsarte linear programming yields $c_8(\mathcal{F}) \ge \max\bigl(0, \frac{m(m-66)}{22}\bigr)$, uniquely saturated at $m=66$ by the self-complementary Witt block system $W_1 = S(5,6,12)/\mathbb{Z}_2$. Explicit $2$-class and $3$-class dual polynomial certificates in $t = d^2 - 1 > 0$ establish dimension-uniform global optimality across all integers $d \ge 2$ at $m \in \lbrace 11, 22, 55, 66, 132, 330, 396, 462 \rbrace$.

**4. Overbudget Parseval Sandwich and Asymptotic Collapse (Theorem 4).** For all $n = 2k \ge 4$, integers $d \ge 2$, and $1 \le m \le M_n$, the exact Parseval row sum $R _ {n,d} = \Lambda_0(n,d)$ brackets $B _ {n,d}(m) = \min _ {\lvert\mathcal{F}\rvert = m} S _ {n,d}(\mathcal{F})$ inside:

$$
\max\Bigl(0, \frac{m}{2}\Bigl(\frac{m}{M_n} R_{n,d} - (d^n + 1)\Bigr)\Bigr) \le B_{n,d}(m) \le \frac{m(m-1)}{2(M_n - 1)} (R_{n,d} - d^n - 1).
$$

Whenever $m / M_n \gg \sqrt{n} (d/(d+1))^n$, $B _ {n,d}(m) = \frac{d+1}{2\sqrt{d}} m^2\bigl(\frac{d+1}{2}\bigr)^n (1 + o(1))$ with exponentially decaying relative gap at any fixed cut density $\rho \in (0, 1]$.

## Repository Structure

```text
HaarPurityDesigns/
├── README.md                                  # Manuscript overview and reproduction instructions
├── CITATION.cff                               # Citation metadata (ORCID: 0009-0007-3678-0810)
├── .zenodo.json                               # Zenodo archive metadata
├── AGENT_SKILL_MANIFEST.json                  # Machine-readable theorem and artifact manifest
├── LICENSE                                    # MIT License
├── .agents/skills/overbudget-haar-designs/
│   └── SKILL.md                               # AI-agent skill for independent verification
├── paper/
│   ├── overbudget_haar_designs.tex            # AMS-LaTeX manuscript source
│   ├── overbudget_haar_designs.pdf            # Compiled 11-page PDF manuscript
│   ├── references.bib                         # Verified BibTeX references
│   └── figures/
│       ├── fig1_es2_vs_haar_n12.pdf           # Separation ratios and Witt phase transition (vector PDF)
│       ├── fig1_es2_vs_haar_n12.png           # Separation ratios and Witt phase transition (300 DPI PNG)
│       ├── fig2_asymptotic_collapse.pdf       # Folded Johnson spectrum and O(1/n) Parseval collapse (PDF)
│       └── fig2_asymptotic_collapse.png       # Folded Johnson spectrum and O(1/n) Parseval collapse (PNG)
└── scripts/
    ├── verify_overbudget_designs.py           # Zero-dependency exact-arithmetic verifier (< 0.9s)
    ├── build_certificates.py                  # Deterministic builder for canonical bitmask certificates
    ├── generate_figures.py                    # Figure generator reading solely from verified_results.json
    ├── certificates_data.json                 # Canonical bitmask witnesses for Tables 1 and 2
    └── verified_results.json                  # Exact integer/rational output consumed by tables and figures
```

## Independent Verification

All finite exact certificates, $S_4$ moment counts, folded Johnson spectra ($6 \le n \le 20$), Delsarte polynomial identities ($n \in \lbrace 12, 14 \rbrace{}$), Parseval convergence values ($n \le 320$), and negative mutation checks run in under one second using Python's standard library alone (`fractions.Fraction` and exact integer arithmetic, zero floating-point relaxations):

```powershell
python scripts/verify_overbudget_designs.py
```

Measured execution: `ALL_CHECKS_PASSED in 0.870s (peak mem 0.18 MB)`.

To regenerate Figures 1 and 2 directly from `scripts/verified_results.json` (requires `matplotlib`):

```powershell
python scripts/generate_figures.py
```
