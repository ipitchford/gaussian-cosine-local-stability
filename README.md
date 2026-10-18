# Local stability of cosine configurations for Gaussian minima

Version 1.0.0 · candidate tag `v1.0.0-candidate` · DOI [10.5281/zenodo.23453252](https://doi.org/10.5281/zenodo.23453252). Repository identity: [ipitchford/gaussian-cosine-local-stability](https://github.com/ipitchford/gaussian-cosine-local-stability). These are release identifiers, not an assertion of external peer review.

This candidate develops a local theorem for Gaussian minima in every order n ≥ 3 and every moment p > 1. For the equally spaced planar cosine frame, it computes the full second-order coefficient, including directions that increase covariance rank. Its kernel consists precisely of common ambient rotations. Consequently the frame is a strict local minimum modulo its symmetries, with a positive distance-squared lower bound in a sufficiently small neighbourhood.

The neighbourhood is existential and depends on fixed n, ambient dimension and p. No numerical radius or uniform estimate as p approaches 1 or n increases is certified. **The global Gaussian-minimum conjecture, including the four-coordinate quadratic problem, remains unresolved by this work.**

Read [paper.tex](paper.tex) for the full argument, [CLAIMS.md](CLAIMS.md) for the evidence map, and [LIMITATIONS.md](LIMITATIONS.md) before reusing a claim. The known cosine law is attributed to Kunisky; the intrinsic-volume identity is attributed to Long. The new local argument is not presented as a new proof of those antecedents. [CITATION_AUDIT.md](CITATION_AUDIT.md) records the bounded source comparison.

## What the code contributes

The standard-library checker verifies exact four-coordinate quadratic-form identities in Q(sqrt(2)). Numerical checks compare two evaluators on sampled four-row frames, reevaluate all 200 supplied search endpoints, check local paths and run fail-closed negative controls. These checks do not prove the all-order analytic theorem, its differentiability or Taylor remainder. Neither implementation supplies interval enclosures or a universal evaluator.

Run the checks in fresh copies without altering saved evidence:

```sh
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python3 -B reproduce.py /path/to/new-output-directory
```

This executes the audit and semantic tests in ordinary and optimized Python. See [ENVIRONMENT.md](ENVIRONMENT.md) and [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). The output directory must not exist. Optimizer trajectories are not rerun by this command; unsuccessful historical exits remain visible.

Scholarly creator: Anonymous. This is an unrefereed research candidate, not a formally verified or externally peer-reviewed result. Original prose/data are CC0; original code is MIT, with component boundaries in [LICENSES.md](LICENSES.md).
