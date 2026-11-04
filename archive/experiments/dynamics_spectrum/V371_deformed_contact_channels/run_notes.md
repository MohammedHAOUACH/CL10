# V371 - deformed contact channels: run notes

Verdict: **NO-THRESHOLD-CONTACT-CUT** (tier 2, residual). Wall 421 s, GPU 1.

## What was measured
- ell=1 (E4 fiber 46080): kernel 2 taken from the V369 SM-400 cache
  (V370's own source; fold order guarantees it). J.L = sum_a (T1+T2).L
  acts on it as **-64 I_2** exactly (block diag [-64, -64], leak
  2.1e-15 residual / ~9e-14 invariance leak reproduces V370's 8.8e-14).
- ell=2 (fiber 32x32x45x45 = 2073600): kernel of the diagonal Casimir
  C = sum_a (J_a^tot)^2 = **6** (V372's count re-measured), Ritz-polished
  residuals ~5e-15, Gram defect ~5e-14, JL-invariance leak (rel) 3.0e-14.
- Deformation spectra D = W^dag (J.L) W:
  ell=1: -64, -64 (pure phase - V370).
  ell=2: **d_k in {-96, -96, -64, -64, 0, 0}** - three two-fold bands.
  The deformation is NOT a phase: two channel bands sit at -64 (same as
  ell=1) and two at -96; two channels are untouched (0).
- Continuum deformation band on E4: J.L spectrum [-64, 8], smear unit 72.
- Window scan at V368 kinematics (verbatim, seed 0x51): free-step min gap
  2.0171e-05 reproduces V368 to 1e-9. For every eta > 0 with
  eta * 72 >= 2 * free_min the channel sets (27648 points per K) sit in a
  smeared continuum of 1572864 points per K: no window. GATE (max_eta
  window >= 0.05) fails at every grid eta: best window 0.0000 at eta 0.00.
  Hard cutoff: eta >= 2*pi/72 = 0.087 zeroes every window by smear alone.

## Verdict meaning
The published bond-current deformation splits the six ell=2 singlet
channels into three bands (-96/-64/0), i.e. it genuinely acts on the
pair's algebraic channel structure (it is not a pure phase), but at the
contact cut (pair at contact bond, one flux line of length <= 2) it
NEVER re-isolates a singlet channel against the two-quanta continuum:
the deformation displaces the continuum by up to 72 rad/step, vastly
more than the 2e-5 rad/step free isolation. Door n°2 does not open at
this cut. The full lattice-sum version of V371 stays booked unrun.

## Machinery corrections (all booked, none changed a gate threshold)
1. The kernel operator is the diagonal **Casimir sum_a (J_a^tot)^2**
   (V366/V370/V372 convention). A mid-run "linear collapse" to
   (sum_a J_a)^2 is a DIFFERENT operator (dense continuum at zero -
   plain-Lanczos probe bg_39 showed Ritz values falling like j^-2);
   it failed the booked ker-count gate and was reverted.
2. Fast evaluation via precomputed leg-pair couplings (K_M, K_F,
   C_MM, C_FF, C_MF with pair factor 2 for AB+BA). One inverse
   permutation typo for the (1,2) matter-flux pair silently changed
   the operator - caught by the new IN-SCRIPT operator self-test gate
   (coupling form == 45-loop anchored form, 5.7e-16) that now guards
   this class of bug permanently.
3. ell=1 kernel span: a fresh SA eigsh gave residuals 3e-15 but a
   JL-invariance leak 3.2e-3 - NOT the published anchor span where
   J.L acts as -64 I (leak 8.8e-14). The run anchors on the V369
   cache span (literal V370 reproduction) and a JL-leak gate
   (ell=2: rel 3.0e-14 <= 1e-3) now certifies fiber-span invariance.
4. First full-fiber eigsh attempt found 4/6 kernel vectors; after the
   operator fix the SA-endpoint Lanczos + Rayleigh-Ritz polish
   (tol=1e-9 above the roundoff floor, ncv=30, partial-Ritz fallback)
   lands all 6.
5. eta scan early-exit (running-min window is exactly 0 once hit);
   a-fortiori zero only for eta*smear >= 2*pi (free-step gap is NOT a
   bound for shifted channel sets).
