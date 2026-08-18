# Repository Notes

- Install development dependencies with `Rscript -e "install.packages(c('deSolve', 'bench', 'testthat'))"` if they are missing.
- Run the test suite with `Rscript -e "testthat::test_local()"`.
- Run package checks with `R CMD build .` followed by `R CMD check --no-manual agepi_0.0.0.9000.tar.gz` when practical.
- Run the performance benchmarks with `Rscript development/benchmarks/run_benchmarks.R`.
- Run the profiler with `Rscript development/benchmarks/profile_benchmarks.R`.
- Optimisation work must preserve numerical behaviour unless explicitly authorised otherwise.
- Future performance changes must report benchmark comparisons against the recorded baseline before and after the change.
