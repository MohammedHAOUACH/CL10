# V374 - full lattice-sum threshold scan: run notes

Verdict: **NO-THRESHOLD-LATTICE-SUM** (tier 1, residual). Wall
15.2 s, CPU (GPU-free). Suite 287.

## Machine
The V371 contact-cut window machine lifted to full lattice sums:
composite Q = beta_A + beta_B + V251 omega rows (full sum over the
N^3 momentum grid) and continuum P = beta_A + beta_B + iota + iota
(the V368 machinery sampled 24-of-216 grid points per K at seed
0x51; here every grid point enters - the booked "full lattice-sum
version"). Dressed by the MEASURED spectra only: composite shifts
{-96,-64,0} (V373 dk sets, ell=1/2/3 value sets agree), continuum
2-quanta shifts = endpoint sums of the V371 J.L fiber band
[-64, 8]: {-128,-56,16} (declared endpoint bookkeeping). N in
{4, 6}; K = the V373 gap_order top-12 interpreted on each grid
(declared); eta in [0,1] step 0.02. Fibers: comp 27648 / cont
14155776 per K at N=6 - inside the booked 2e7 gate; ell=3 value
sets coincide with ell<=2, long-string sector declared (contingency
never silently taken).

## What was measured
- g(eta=0) max = 0.0312 (N=6), -0.0405 (N=4): the composite ladder
  sits at the continuum edge, no gap - consistent with the free
  V368 picture.
- The dressing DOES open a nominal phase window: max g = 0.1930
  rad/step at (eta 0.48, K = (0,0,0)) on the resolved N=6 sums;
  the N=4 coarse-lattice maximum 0.8874 at eta 0.32 collapses
  under refinement (0.19) - coarse-continuum under-resolution.
- The booked finite-size discriminator is decisive and fails
  cleanly: the set-model composite states are plane-wave pairs,
  <d_geo(r)^2> = 4.500 (N=4) -> 9.500 (N=6), ratio 2.111 = the
  exact uniform-pair value (closed form gated to 0.0e+00) -
  continuum-like spreading, no bound-state localization.
- Booked conjunction: (g >= 0.05) AND (size-stable) is FALSE on the
  second clause at every scan point: NO-THRESHOLD-LATTICE-SUM.
  Door n-2 stays closed at full lattice-sum scope.

## Anchors
V368 free min gap 2.017083e-05 to 1e-9 (verbatim seed 0x51);
V371 contact-cut window scan reproduced point-by-point (best
2.0171e-05 at eta 0.00, artifact match <= 1e-9); V373 dk sets from
artifact (2-link exact, 3-link value set); V371 band/smear 72;
V251 omega_rows pair +-8.418126; edge function (min of mod-2pi
over set x shifts) self-tested against brute force, dev 0.0.

## Machinery lessons
1. min mod-2pi over (set + shift): the discontinuity of the saw
   is at t = mod2pi(-shift) and the minimum sits at the first set
   point >= t, NOT the last one <= t - the first implementation
   got the side wrong and the self-test caught it with a 2 pi
   deviation. Brute-force self-test on random sets, always.
2. A "below the edge" claim on a phase circle needs the
   localization discriminator; the set model trivially puts plane-
   wave ladders "below" dense-continuum edges. The booking's
   two-clause criterion is what makes the number meaningful -
   report both clauses, never the g-crossing alone.
3. Coarse-lattice (N=4) continuum sets under-resolve the edge by
   ~1 rad; always report the resolved-N max separately.

## What this licenses
V375 cross-checks the dressed eigenphase differences (this scan's
shifted spectra + V373 kernel spectra) against the five empirical
ladder contacts; the V371 "full lattice-sum version stays booked
unrun" clause is retired.
