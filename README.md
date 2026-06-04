# 6N Prime Triplet: a Three-Member Two-Centre CRT Pattern (Part XV)

The prime triplet (6N-1, 6N+1, 6N+5) has a **non-monotone** omega-dependence —
a peak at omega=3, then a fall — reproduced by the three-member extension of the
Paper X two-centre CRT mechanism, no fitted parameters.

**The triplet straddles two centres.** The only admissible triplet patterns within
a span of six are (0,2,6) and its mirror (0,4,6); (0,2,4) is inadmissible (occupies
all residues mod 3). On the 6N skeleton, (p,p+2,p+6) with p=6N-1 is

```
    (6N-1, 6N+1, 6N+5) = (6N-1, 6N+1, 6(N+1)-1),
```

so its three members sit at (centre-offset, wing) = (0,-1), (0,+1), (1,-1): two
wings on centre N, one on the neighbour N+1. Like cousins and sexy pairs (Parts
IX-X), and unlike the single-centre twin, the triplet straddles two consecutive
centres — so its omega-dependence is a two-centre distortion, not the monotone
single-centre enrichment of the twin.

**A non-monotone distortion.** Binned by omega_{>3}(N), normalised to omega=1, the
triplet rate rises to a peak near omega=3 then falls — the sexy-B signature of two
competing effects: the twin component (6N±1) on the conditioned centre N enriches
with omega, while the third member on N+1 (like the cousin's partner) suppresses
the rate at high omega. The balance gives the peak.

**Three-member two-centre CRT model.** For a prime q>3, with dead(q,s) = -s·6⁻¹
mod q, the triplet CRT factor over the three members {(0,-1),(0,+1),(1,-1)} is

```
    q | N : deterministic — each member q-safe unless its offset hits dead(q,s)
    q ∤ N : averaged over admissible nonzero residues
```

The rate model is prod_{q in POOL} f_q, averaged within each omega stratum. No
fitted parameters; only the wing geometry enters.

**Result (S9 and S10).** Error ≤ 2.4% for omega ≤ 6 on S10 (~2.3×10⁶ triplets).

| omega | triplets (S10) | measured | model |   | shell | omega=6 residual |
|------:|---------------:|---------:|------:|---|------:|-----------------:|
| 1 | 356,239 | 1.000 | 1.000 |   | S9  | 6.5% (56 triplets) |
| 2 | 875,673 | 1.025 | 1.021 |   | S10 | 0.5% (1832 triplets) |
| 3 | 768,157 | 1.048 | 1.036 |   |     |                  |
| 4 | 288,800 | 1.039 | 1.020 |   |     |                  |
| 5 |  43,127 | 0.950 | 0.927 |   |     |                  |
| 6 |   1,832 | 0.691 | 0.687 |   |     |                  |
| 7 |      11 | 0.327 | 0.254 |   |     |                  |

> **Small-sample check.** The omega=6 residual falls from 6.5% (56 triplets, S9)
> to 0.5% (1832, S10) — the S9 excess was sampling noise. The omega=7 stratum has
> only 11 triplets on S10; the 22% residual there is not interpretable. The result
> is stated for omega ≤ 6.
>
> **Attribution.** The Hardy-Littlewood triplet singular series S(0,2,6)=2.858 is
> classical (our sieve reproduces it); this paper proposes no new constant. The
> contribution is identifying the triplet's omega-dependence as a non-monotone
> two-centre distortion and reproducing it with the parameter-free CRT model.
>
> **Scope.** No statement about the infinitude of prime triplets or any k-tuple
> conjecture. This is NOT a "divide by S and collapse to a constant" result (that
> is Goldbach/Polignac): the triplet straddles two centres and is a distortion,
> like cousins/sexy.

Part IX: doi:10.5281/zenodo.20520492 · X: doi:10.5281/zenodo.20525642 ·
XIV: doi:10.5281/zenodo.20533747

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── triplet_data.csv    shell, omega, triplets, measured, model, err_pct (S9+S10)
├── code/
│   ├── triplet.py          counts triplets (6N-1,6N+1,6N+5), evaluates the
│   │                       three-member two-centre CRT model vs measured by omega;
│   │                       emits triplet_S{K}.csv. Streaming, low memory. Default S10.
│   └── make_triplet_fig.py builds the 2-panel figure from ../data
├── figures/                fig_paper15_triplet.{pdf,png}
└── paper/                  Chen_6N_Paper15.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Triplet model vs measured. Default S10 (~28 min). Memory-light streaming.
#    Emits triplet_S{K}.csv.
python code/triplet.py            # S10
MAXK=9 python code/triplet.py     # S9 (faster; shows omega=6 small-sample)

# 2. Figure (reads ../data/triplet_data.csv).
cd code && python make_triplet_fig.py
```

### Conventions (same as Parts I-XIV)

- Triplet (6N-1, 6N+1, 6N+5); members at (offset, wing) = (0,-1),(0,+1),(1,-1).
- omega_{>3}(N) = number of distinct prime factors >3 of the left centre N.
- dead(q,s) = -s·6⁻¹ mod q : residue of N making wing 6(N+off)+s divisible by q.
- Rates normalised to omega=1 for shape comparison.
- Engine: complete segmented-sieve factorisation + interval-sieve primality;
  S10 twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
