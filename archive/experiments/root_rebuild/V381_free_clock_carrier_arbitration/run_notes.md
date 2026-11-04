# V381 - the free-clock carrier arbitration: run notes

Verdict: **CARRIER-ARBITRATION-BREAKS** (tier 2, O8/O9; the
booked fallback branch FIRED - the carrier does not dissolve the
clock arbitration). Wall 3.8 s. Suite 294.

## The per-beta table (the whole result)

| beta | closure (cap 24) | gapM | R signature | 0-point mod 2pi |
|---|---|---|---|---|
| pi/4 | 2 | 0.0 | [32, 48] | 0 |
| pi/3 | 6 | 1.8e-02 | [16, 48] | 0 |
| pi/5 | 10 | 3.1e-01 | [16, 48] | 0 |
| 2pi/5 | 10 | 7.5e-01 | [16, 48] | 0 |
| pi/6 | 6 | 5.2e-01 | [16, 48] | 0 |
| 1/5 rad | none | 6.1e-01 | [16, 48] | 0 |
| g_obs | none | 6.1e-01 | [16, 48] | 0 |
| 0.3 rad | none | 1.3e-01 | [16, 48] | 0 |

Unitarity <= 2e-15 everywhere; internal mirror exact.

## What it means

- The CLOSURE family is genuinely freed: every rational
  multiple of pi closes (orders 2, 6, 10, 10, 6), far beyond the
  booked forcing set {pi/2, pi/4, 3pi/4}.
- The PHYSICS rides only at pi/4 - the published clock: every
  other closed beta splits the R 32-cluster exactly dyadically
  ([16, 48] - the same split the V379 near-miss map found at the
  2I points) and gaps the M node.
- The rational-RADIANS pins (1/5, g_obs, 0.3) never close to cap
  24 - the V378 transport-incommensurability lesson verbatim.
- The 0-point survives on the carrier at every beta (the pi/4
  sum wraps to exactly 2 pi - zero mod 2 pi).

So V334's rigidity re-asserts itself where it matters: closure
freedom without spectrum contact is arithmetic only. F2's clock
arbitration stands; its question transfers to whether the
carrier LADDER re-reads the coefficient - V382.

## Declared corrections (same day)

Two arithmetic self-test bugs in my check lines (mirror line
computed |ph| + ph; zero-point summed raw phases instead of mod
2 pi) failed the first run; both fixed and declared in the
pre-registration file before re-issue. No tolerance on any
physics gate was touched; the T2 contact gate fails on real
data and the verdict branch is the booked one.
