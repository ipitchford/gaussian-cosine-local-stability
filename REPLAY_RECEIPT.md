# Observed repair replay

On 18 October 2026 both ordinary and optimized Python executions passed the local algebra checker, numerical endpoint audit and seven-method semantic regression suite. The 200 saved endpoints were preserved byte-for-byte. Dimension three has 100 successful flags; dimension four has 78 successful and 22 unsuccessful flags. Maximum stored-versus-intrinsic discrepancies were 4.551914400963142e-15 and 3.4772185131259903e-12 respectively.

The original seeded sequence of 100 random frames in dimension three followed by 100 in dimension four was checked with explicit finiteness before discrepancy reductions. Eight local perturbation paths and four Sobol sanity comparisons were also replayed. These finite checks do not prove the analytic all-n/all-p theorem, global optimality or numerical conditioning throughout the domain.

The exact original quadratic algebra and deliberate matrix mutation pass with rational Q(sqrt(2)) arithmetic. Semantic tests reject nonfinite evaluator outputs, invalid input shapes, missing batches, duplicate starts and nonfinite saved records. Difficult valid frames test explicit elementary rejection and finite intrinsic evaluation. This is domain disclosure, not universal numerical coverage.

Run `python3 -B reproduce.py NEW_OUTPUT_DIRECTORY` to produce new machine-readable receipts and complete logs. A final clean-archive replay must be distinguished from the development repair replay described here. The original optimizer trajectories were not rerun.
