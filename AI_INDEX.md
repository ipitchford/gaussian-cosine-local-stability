# AI research index

## Entry and claim ceiling

Title: Local stability of cosine configurations for Gaussian minima. Version 1.0.0; scholarly creator Anonymous; candidate tag v1.0.0-candidate. Citation identity is in [CITATION.cff](CITATION.cff).

Read [README.md](README.md), [CLAIMS.md](CLAIMS.md), [paper.tex](paper.tex), and [LIMITATIONS.md](LIMITATIONS.md). The central result is **local** stability for all n≥3 and p>1. Do not summarize it as resolution of the global Gaussian-minimum conjecture. Numerical tools remain four-coordinate diagnostics, not all-order proof checkers.

## Proof dependencies and evidence

GM-1 uses the known cosine radial law, an explicit bulk/interface second variation and integrability for p>1. GM-2 uses the antiperiodic signed-cycle spectrum and an orbit-normal argument. GM-3 is the exact four-coordinate algebraic specialization, checked by [check_local_algebra.py](code/check_local_algebra.py). GM-4 is replayed by [check_results.py](code/check_results.py); [test_guards.py](code/test_guards.py) checks semantic rejection, including NaN and missing batches. See [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md) for actual observed coverage and [ENVIRONMENT.md](ENVIRONMENT.md) for runtime identities.

## Novelty and reusable research directions

Supporting claims GM-6 and GM-7 in [CLAIMS.md](CLAIMS.md) cover the nonglobal local minimum at independence (paper Section 6) and the independent-noise interior exclusion (Section 9). Neither closes the remaining global geometric domain.

[CITATION_AUDIT.md](CITATION_AUDIT.md) separates known identities from the local extension. [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md) records substantive repairs without treating review as proof. Potential worthwhile next steps are a genuine global rank-three/rank-four inequality, or an effective local radius strong enough to support a global cover. More optimizer seeds do not discharge either obligation. Extending to p≤1 requires a new singularity argument, not direct substitution.

## Trust and reuse

[ASSURANCE.md](ASSURANCE.md), [PROVENANCE.md](PROVENANCE.md), and [LICENSES.md](LICENSES.md) identify review/replay scope and component rights. [reproduce.py](reproduce.py) writes only to a new output directory and retains normal/optimized logs. No check count should be promoted into analytic or formal assurance.
