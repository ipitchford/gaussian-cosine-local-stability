# Reproduction environment

Fresh repair replay used Python 3.11.14, NumPy 2.4.6, SciPy 1.17.1 and SymPy 1.14.0 on macOS arm64 with one BLAS/OpenMP thread. requirements.txt pins these tested numerical-library versions. Original endpoint generation used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0. Endpoints were reevaluated, not regenerated: this is cross-version validation, not trajectory repeatability.

The exact local-algebra script needs only Python's standard library. Other scripts need the pinned dependencies. The public runner copies sources/data into two new directories before running ordinary and -O modes; its receipt includes versions, subprocess statuses and hashes. No network access is used by the checks themselves.
