# Claim-to-evidence map

| ID | Claim and scope | Evidence and boundary |
|---|---|---|
| GM-1 | For every fixed n≥3, d≥2 and p>1, the cosine frame has the stated second-order expansion, uniformly on bounded tangent directions. | Written bulk/interface differentiation and locally integrable singularity argument in paper.tex. No finite script proves regularity. |
| GM-2 | The signed-cycle normal form is positive semidefinite with exactly the rotational kernel; strict orbit-local minimality and a local distance-squared bound follow. | Antiperiodic Fourier spectrum and compact-orbit argument in paper.tex. No explicit radius. In d=2 the stated minimum of angular/normal gaps is a lower bound, not necessarily the exact gap. |
| GM-3 | The four-coordinate p=2 specialization agrees with the exact Q(sqrt(2)) identity. | code/check_local_algebra.py checks algebra only, including a deliberate matrix mutation. |
| GM-4 | The supplied four-coordinate endpoint records and sampled evaluator comparisons are internally consistent within stated tolerances. | code/check_results.py; all 200 endpoints, 22 unsuccessful optimizer flags retained. Floating point, not certified optimization. |
| GM-5 | Planar comparison and independent-block special cases supply supporting partial results, not the global rank-three/four theorem. | Written special-case arguments in paper.tex; known radial law and intrinsic-volume antecedents explicitly attributed. |
| GM-6 | Independence is a strict nonglobal local minimum for four-coordinate second moments. | Section 6, Proposition 6.1: written density expansion and positive integral coefficient. Numerical paths are supporting diagnostics, not the proof. |
| GM-7 | Independent-noise floor yields an explicit interior exclusion region, F(C) ≥ lambda_min(C) F(I4), for four-coordinate second moments. | Section 9, Lemma 9.1: elementary Gaussian interval comparison. This is a supporting comparison, not a new general principle or complete global cover. |

The ordinary Hessian is twice the displayed Taylor quadratic coefficient. Common row signs and permutations are discrete symmetries, not extra infinitesimal directions. The range p≤1 is not claimed. The all-order theorem is analytically assessed; the numerical search implementation remains restricted to four rows.
