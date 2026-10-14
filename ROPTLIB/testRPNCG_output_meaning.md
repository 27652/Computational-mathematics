# Meaning of the `testRPNCG` Output

Run the program with:

```bash
./build/bin/testRPNCG
```

This executes a sparse PCA optimization example using the RPNCG solver. The
output is an optimization progress log, not an error report.

## 1. Problem being solved

The test creates a random data matrix `B` of size `50 x 800` and solves the
following problem:

\[
\min_{X^T X = I_5}
F(X) = -\|BX\|_F^2 + 0.8\|X\|_1,
\qquad X \in \mathbb{R}^{800 \times 5}.
\]

The five columns of `X` are orthonormal component directions.

- The term `-||BX||_F^2` encourages the components to capture variation in the
  data.
- The term `0.8||X||_1` encourages sparse components.
- Smaller objective values are better.
- Negative objective values are normal.

The printed `f` is the total objective, including the sparsity penalty.

## 2. Initial configuration

Typical initial output is:

```text
seed:2
GENERAL PARAMETERS:
Tolerance     :          1e-06[YES],
OutputGap     :              1[YES],
Max_Iteration :           1000[YES],
Verbose       :     ITERRESULT[YES],
```

| Output | Meaning |
|---|---|
| `seed:2` | Fixed random seed used to generate the data and initial point. |
| `Tolerance:1e-06` | Stopping threshold: `10^-6`. |
| `Max_Iteration:1000` | Maximum number of outer iterations. |
| `OutputGap:1` | Print every iteration. |
| `Verbose:ITERRESULT` | Print iteration-level diagnostics. |
| `[YES]` | The parameter passes the library's validity check. |

In this implementation, the custom `RPNCG::Run()` stops when

```text
min(nv, nu) <= Tolerance
```

or when `Max_Iteration` is reached.

## 3. Meaning of the iteration output

Each iteration normally prints three lines. For example:

```text
iter:22, ssn_iter:2, inner_it:1, status:4, t:3.276112e-02, tau:1.000000e+02
iter:22, f:-3.736990e+01, nv:3.602351e-02, nu:7.788951e-02, nonzero:952, zero:3048
iter:22, f:-3.742774e+01, nv:3.602351e-02, nu:7.788951e-02, btiter:0, time:0.055626, innertol:1.000000e-11
```

The second line reports the objective before the step. The third line reports
the objective after the step.

| Field | Meaning |
|---|---|
| `iter` | Outer optimization iteration. |
| `ssn_iter` | Semismooth Newton iteration count used to compute the proximal direction `v`. |
| `inner_it` | Truncated conjugate-gradient iteration count for the Newton correction. |
| `status` | Reason the inner tCG calculation stopped. See the table below. |
| `t` | Adaptive proximal step-size parameter. |
| `tau` | Quadratic regularization parameter in the inner model; normally 100 in this test. |
| `f` | Total objective `F(X)`. |
| `nv` | Frobenius norm `||v||_F` of the proximal direction. |
| `nu` | Frobenius norm `||u||_F` of the final search direction. |
| `btiter` | Number of step halvings performed during backtracking. |
| `time` | Cumulative elapsed solver time in seconds. |
| `innertol` | Accuracy threshold for the inner direction calculation. |

For example, an objective changing from `-37.36990` to `-37.42774` means that
the optimization step improved the objective.

## 4. Meaning of `status`

These are inner-solver exit codes, not overall success/failure codes.

| Status | Name | Meaning |
|---:|---|---|
| `-1` | `RPNCG_MAXITER` | Inner CG iteration limit reached. |
| `0` | `RPNCG_EARLY1` | Initial model-decrease check triggered an early return. |
| `1` | `RPNCG_EARLY2` | No eligible correction entries, or the initial curvature/model check triggered an early return. |
| `2` | `RPNCG_NEGATIVE_CURVATURE` | Insufficient positive curvature was detected. |
| `3` | `RPNCG_EARLY3` | A candidate correction failed a model safeguard. |
| `4` | `RPNCG_LINEAR` | The CG residual met the linear-regime stopping threshold. |
| `5` | `RPNCG_SUPERLINEAR` | The CG residual met the superlinear-regime stopping threshold. |

Status `5` is an inner-solver result; by itself it does not prove overall
convergence or that the global optimum was found.

## 5. Meaning of `nonzero` and `zero`

The labels `nonzero` and `zero` can be misleading. They are **not literal
counts of nonzero and zero entries in `X`**.

The implementation divides the `800 x 5 = 4000` entries into two index sets:

- `nonzero` is the size of the set selected for the Newton correction.
- `zero` is the size of its complement.
- They always add up to 4000.

Therefore, output such as

```text
nonzero:0, zero:4000
```

does not mean that `X` is the zero matrix. It means that no entries currently
qualify for the Newton-correction subset.

## 6. Interpreting the final result

A typical final line is similar to:

```text
iter:147, f:-3.850512e+01, nv:8.693090e-08, nu:1.513399e-06, ...
```

This means:

- Final objective: approximately `-38.50512`.
- The stopping measure is

  ```text
  min(8.693090e-08, 1.513399e-06) = 8.693090e-08
  ```

- Since `8.693090e-08 < 1e-06`, the tolerance condition is satisfied.
- The solver stopped because of the tolerance, not because it reached the
  maximum of 1000 iterations.

The objective typically improves from approximately `85.24844` at the first
iteration to approximately `-38.50512` at the end. Since sparse PCA on the
Stiefel manifold is nonconvex, satisfying the stopping condition does not
guarantee that the global optimum has been found.

The program currently does not print a separate final `converged` message; it
simply ends after the last iteration.
