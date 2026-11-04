# V380 - the cubic-axis clock-free core: run notes

Verdict: **DILEMMA-BROKEN** (tier 2, O8/O9; the booked solution
branch FIRED). Wall 4.5 s (host numpy). Suite 293.

## What was found

The one carrier class every tribunal missed: finite-order twists
of the SPATIAL-plane Spin(3) (Sx, Sy, Sz of the V337-certified
triple, S_a^2 = I/4 exactly, anticommutators exactly 0) about
CUBIC axes, on the published crystallographic transport
(beta = G10, the order-8 root clock). Of the 10 booked points,
3 carry the full three-leg conjunction: u_2i(pi/2, e_x),
u_2i(pi/2, e_y), u_2i(pi/2, e_z) -

| gate | reading |
|---|---|
| internal closure | order 4 (clock-free, independent of g) |
| FULL-step closure at GAMMA | order 2 |
| gapM over all 8 M-points | 0.0 exactly |
| R signature | [32, 48] - the published multiplets exactly |
| unitarity | 2.2e-16 |
| internal mirror | 0.0 |

The order-2 points (pi/2 about coordinate axes) are the carriers;
the pi-order points gap M (gapM 5.24e-01) and split R to
[64, 96]; the body-diagonal points shatter R ([24..40 x plus
96]) - the axis choice is everything.

## The mechanism (measured, then used in V381)

The staggered transport INVERTS the spatial-plane twist about a
coordinate axis: (P R)(P R) = +-I at beta = pi/4 (step^2 defect
1.33e-15) but NOT at beta = 1/5 rad (defect 0.66) - the
involution is a rational-pi phenomenon. The closure axiom no
longer forces {pi/2, pi/4, 3pi/4} on this carrier: V381 measures
the freed family.

## Standing

The crystalline dilemma (V377/V378: the three certified legs
never co-occur) is BROKEN as an existence statement - they do
co-occur, on the published engine, at the published clock, in
three points related by cubic symmetry. What it does NOT break
is the clock arbitration: V381 (booked follow-up obligation)
re-runs F2 on the carrier.
