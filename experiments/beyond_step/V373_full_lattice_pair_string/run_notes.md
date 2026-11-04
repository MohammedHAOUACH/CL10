# V373 - full-lattice pair + straight-string machinery: run notes

Verdict: **FULL-LATTICE-MACHINERY-BUILT** (tier 1, machinery). Wall
1822 s final pass (first full census 3595 s, GPU 1), suite 286.

## Run history (kill-recovery, documented per protocol)
The first two full attempts were killed by memory pressure: the host
(62 GB) was under llama.cpp serving load (~60 GB), the run swapped
(53 GB RSS), and the 320K-1M-vector scipy blocks degraded ~100x
(2152 s per 323K block) while the U5/U6/U3 stack alone (2.5 GB of
complex128) ballooned through temp copies. Fixes, all kept:
- GPU route (CuPy Lanczos, batched Ritz polish) lowered from 5e5 to
  6e4 states: every medium/large block runs on the GPU; matmul
  batching replaced the per-vector loop (the 17.9M-vector ell=2
  block: 41 min -> 20 min, polish batches ~1-2 s each).
- Per-block checkpoint (results/v373_ck/, script-sha guarded, stale
  dirs wiped): every measured chain block is cached to its own npz.
  A killed run resumes at block granularity; after three T9 bug-fix
  edits the full T5 census replayed in 30 s (vs 3595 s). The sha
  was transplanted by hand on those edits because the solver region
  was provably untouched - never do that for a solver edit.
- Uf/U2/U3/U4 keep per-block vector lists instead of the concatenated
  stack: the final npz save stages 8 vectors at a time (peak 3 GB,
  no 2.5 GB concatenation + save copy).

## What was measured
- Pair channels (45x45, exact eigh): {0:2, 36:20, 64:90, 84:240,
  96:420, 100:252} - the measured fusion content; generator blocks
  per channel certified against the Casimir (gdefect <= 1e-8).
- Link-pair channels: {0:1, 64:45, 80:54, 96:210, 128:945, 144:770}.
- Chain-kernel dims (MEASURED per chain length, dim Hom measured by
  block, Hom additive over blocks): contact 2 (= V371 kernel),
  ell=1 (2 links) **6**, ell=2 (3 links) **22**, ell=3 (4 links)
  singlet-sector **12**; I_4 = 12/2 = **6** = V100's published 6
  (booked prediction output; non-singlet pair sectors DECLARED
  UNMEASURED at ell=3).
- I_3 = 2 singlet-channel Hom copies = 2 x V100's I_3 = 1.
- V371 contact reproduction: 1-link block sum == T1 full-fiber
  Lanczos kernel 2; 2-link block pattern 2 x {r=1,45,210 blocks}
  sums to 6 == V371's measured kernel; block-built span == V371
  fiber-kernel span, min overlap singular value 1.0000000000.
- Glued-kernel anchor from MEASURED fusion data: Hom 6 x m_45^2
  (90/45 = 2) = 24 == V100/V101 glued kernel, gram ranks 24/24,
  V101 top degeneracy 6.
- Deformation spectra (kernel basis, factored channel form):
  contact -64,-64 (dev 0 = V370 anchor); 2-link exactly
  {-96,-96,-64,-64,0,0} (fiber set reproduced); 3-link
  {-96 x8, -64 x12, 0 x2} - algebraic clean, booked as OUTPUT;
  4-link singlet sector 0 (trivial pair block, declared scope).
  Fiber-route d_k set reproduced to 7.8e-13; fiber JL leak 9.2e-15.
- String Casimir per link: sum_a L_a L_a spectrum flat at 64.000000
  (spread 0.0).
- V368 torus kinematics verbatim (seed 0x51): free-step min gap
  2.017083e-05 reproduced to 0.0e+00; worst K = Gamma (gap_order
  starts at index 0).
- P1 real-space bond propagation: unitarity 2.2e-16, Bloch defect
  1.1e-15 vs SC.bond_product (the V374 pair-hopping machinery,
  gated here).
- Second solver route (materialized link-pair generators), singlet
  channel: 2-link 2 == block 2, 3-link 2 == block 2.

## Machinery lessons (kept for the audit)
1. A reused einsum label in the fiber embedding ("jka,abk->jab")
   silently merges the third link with a block-basis index: the
   3-link deformation spectrum came out NON-algebraic-clean
   (-322.75, -98.09, ...). Detector: the operator is a sum of
   generators of a compact algebra, so every eigenvalue must be
   algebraic-clean; garbage values mean an embedding bug, never
   "interesting physics".
2. Channel dicts keyed by raw eigenvalue floats: census .get(64.0)
   misses 63.999999999999815. Round keys at construction.
3. The deformation operator is block-diagonal over pair channels:
   any Gram over a mixed-sector basis must be built per sector
   (different W dimensions).
4. Unoptimized einsum dispatch on a (945,945) double contraction
   stalls >1 h per tensor; write the two-stage tensordot GEMM form.
5. Sort multisets (np.sort) before comparing concatenated
   sector-blocked spectra to booked ascending lists - a
   sector-ordered spectrum "fails" with a max dev equal to the
   band span.

## What this licenses
The full-lattice pair + straight-geodesic-string operator is now
measured per chain length with certified fusion blocks, both solver
routes agreeing and V371/V368/V100/V101 anchors reproduced from
measured data. V374 (lattice-sum threshold) runs on this machine;
V375 (dressed no-ladder contact) cross-checks these spectra.
