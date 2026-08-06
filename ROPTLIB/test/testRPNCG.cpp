#include "test/TestRPNCGSPCA.h"
#include "Manifolds/Stiefel.h"
#include "Problems/StieSPCA.h"
#include "Solvers/RPNCG.h"

namespace ROPTLIB{


void testRPNCGSPCA(void)
{
    integer n = 800;
    integer m = 50;
    integer p = 5;
    realdp lambda = 0.8;

    Vector B(m, n);
    B.RandGaussian();

    // 和 RPN-CG MATLAB 实验一致：columns zero-mean + normalize
    realdp *Bptr = B.ObtainWriteEntireData();
    for (integer i = 0; i < n; i++)
    {
        realdp mean = 0;
        for (integer j = 0; j < m; j++)
            mean += Bptr[j + i * m];
        mean /= m;

        for (integer j = 0; j < m; j++)
            Bptr[j + i * m] -= mean;

        realdp normcol = 0;
        for (integer j = 0; j < m; j++)
            normcol += Bptr[j + i * m] * Bptr[j + i * m];
        normcol = std::sqrt(normcol);

        for (integer j = 0; j < m; j++)
            Bptr[j + i * m] /= normcol;
    }

    Stiefel Domain(n, p);

    // RPN-CG 的 MATLAB R(x,eta) 是 polar-type retraction；
    // ROPTLIB Stiefel 有 POLAR retraction，可用 Set3/Set4。
    Domain.ChooseParamsSet4();

    Variable X0 = Domain.RandominManifold();

    integer lengthW = 1;
    StieSPCA Prob(B, lambda, n, m, p, lengthW);
    Prob.SetDomain(&Domain);

    RPNCG solver(&Prob, &X0);
    solver.Max_Iteration = 1000;
    solver.Tolerance = 1e-6;
    solver.Verbose = ITERRESULT;
    solver.CheckParams();
    solver.Run();
}
}