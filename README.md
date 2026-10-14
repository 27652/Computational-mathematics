 Computational Mathematics: RPN-CG

C++ implementation of the [Riemannian Proximal Newton-CG method](https://arxiv.org/abs/2405.08365) in ROPTLIB, with the original MATLAB implementation as a reference. The current focus is sparse PCA on the Stiefel manifold.

## Repository overview

- `ROPTLIB/`: C++ library, including `Solvers/RPNCG.cpp` and the test entry points in `test/mainRPNCG.cpp` and `test/testRPNCG.cpp`.
- `RPN-CG/`: upstream MATLAB reference implementation (Git submodule), with SPCA and CM experiments.
- `ROPTLIB/RPNCG - latex/`: implementation notes, requirements, and LaTeX report.
- `paper/`: reference paper.
- `com.md`: detailed completion assessment and remaining validation work.

## Brief roadmap

- [x] Implement the initial C++ RPN-CG solver and sparse PCA test entry point.
- [x] Document class relationships, solver logic, parameters, and memory management.
- [ ] **Check MATLAB/C++ agreement:** share input data and initial points, align retractions and stopping criteria, and investigate differences in inner iterations.
- [ ] **Run fair benchmarks:** compare C++ RPN-CG with a C++ first-order baseline (such as IARPG) and the MATLAB second-order implementation using matched settings and repeated runs.
- [ ] **Publish results and polish:** add timing, accuracy, and operation-count tables to the LaTeX report; record hardware and software settings; clean up solver configuration and code style.

The implementation and technical notes are in place, but numerical equivalence and comparative performance are not yet established. Existing single-run results should not be treated as a speed ranking.
