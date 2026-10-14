#include <vector>   
#include <limits>   
#include "RPNCG.h"
#include "Element.h"
#include<Solvers.h>
#include "Problem.h"
#include "Problems/StieSPCA.h"
#include "f2c.h"
//#include "Solvers/RPNCG.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>






ROPTLIB::RPNCG::~RPNCG() = default;

/*Define the namespace*/
namespace ROPTLIB{

// MODIFIED: helper for MATLAB/Fortran column-major linear indexing.
static inline integer cmidx(integer row, integer col, integer nrows)
{
    return row + col * nrows;
}

// MODIFIED: gather v(ids) using 0-based linear indices. This replaces the old
// std::vector<realdp> gather(), because tCG/Run need a ROPTLIB Vector.
static Vector GatherVector(const Vector &v, const std::vector<integer> &ids)
{
    Vector out(static_cast<integer>(ids.size()));
    out.SetToZeros();

    const realdp *vptr = v.ObtainReadData();
    realdp *outptr = out.ObtainWritePartialData();

    for (integer k = 0; k < static_cast<integer>(ids.size()); ++k)
        outptr[k] = vptr[ids[k]];

    return out;
}

static Vector Symmetric(const Vector &a)
{
    return (a + a.GetTranspose()) * static_cast<realdp>(0.5);
}

static Vector AstarOp(const Vector &lambda, const Variable &x)
{
    return (x * Symmetric(lambda)) * static_cast<realdp>(-1.0);
}

static Vector Aop(const Vector &eta, const Variable &x)
{
    return (eta.GetTranspose() * x + x.GetTranspose() * eta) * static_cast<realdp>(-0.5);
}

static Vector MakeScalarWeight(realdp value)
{
    Vector weight(1);
    weight.ObtainWriteEntireData()[0] = value;
    return weight;
}

static Vector ProxAt(const Vector &B, realdp t, const Problem *prob)
{
    Vector result(B);
    Vector weight = MakeScalarWeight(static_cast<realdp>(1.0) / t);
    prob->ProxW(B, weight, &result);
    return result;
}

static Vector JacobianAt(const Vector &B, const Vector &eta, realdp t,
                         const Problem *prob)
{
    Vector result(eta);
    Vector weight = MakeScalarWeight(static_cast<realdp>(1.0) / t);
    prob->CalJW(B, eta, weight, &result);
    return result;
}

// MODIFIED: directly builds MATLAB Bbar = reshape(x * tmpp, n*m, d)(id1, :).
static Vector BuildBbar(const Variable &x, const std::vector<integer> &id1);





void RPNCG::UpdateData()
    {
        //Mani->EucGradToGrad(x1, egf, &gf1);
    };

    void RPNCG::PrintInfo()
    {
        SolversNSM::PrintInfo();
    };

    void RPNCG::PrintFinalInfo()
    {
        SolversNSM::PrintFinalInfo();
    };

void RPNCG::SetDefaultParams()
    {
        SolversNSM::SetDefaultParams();
        tau = 1e2;
        rho1 = 1e-3;
        rho2 = 0.5;
        w1 = 1.1;
        w2 = 0.9;
        gamma = 0.01;
        vartheta = 0.01;
        theta = 0.5;
        kappa = 0.1;
        max_innit = 200;
    };


    RPNCG::RPNCG(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};




	void RPNCG::Run(void)
	{
		SolversNSM::Run();
		realdp err = std::numeric_limits<realdp>::infinity();


    f1 = Prob->f(x1) + Prob->g(x1);
    nf++;
    f2=f1;

    Vector egf(x1.Getlength());
    Prob->EucGrad(x1, &egf);
    ng++;

    realdp t = 1.0;
    if (const StieSPCA *spca = dynamic_cast<const StieSPCA *>(Prob))
        t = 1.0 / spca->L;
    realdp t0 = t;

    integer flag = 0;
    RPNCGTCGStatus status = RPNCG_MAXITER; // MODIFIED: initialize before first status check

    const integer n = x1.Getrow();
    const integer r =x1.Getcol();

    Vector Dinitial(r,r);
    Dinitial.SetToZeros();

    realdp nv = 1;
    realdp nu = 1;

    //inner tolence
    realdp innertol =
        std::max(static_cast<realdp>(1e-13),
                 std::min(static_cast<realdp>(1e-11),
                          static_cast<realdp>(1e-3) * std::sqrt(Tolerance) * t * t));


//init iter and ndir
	iter = 0;
    ndir0 = 0;
    ndir1 = 1;


//old var
    Variable xold;
    Vector uold;
    realdp fold = f1;
    realdp nuold = 0;
    realdp nvold = 0;



	while (err > Tolerance && iter < Max_Iteration){

        //innertol
        innertol = std::min(std::max(static_cast<realdp>(1e-30), nv * nv * static_cast<realdp>(1e-8)), innertol);

        if (status >= 4)
            innertol = std::min(std::max(static_cast<realdp>(1e-30), nv * nv * nv), innertol);





      integer ssn_iter = 0;
//find dir - compute search direction using semi-smooth Newton-CG method
        Vector v = FindDir(
            x1, egf, t, &Dinitial, innertol, &ssn_iter
        );
        // nv = norm(v) - Frobenius norm of the search direction
        nv = v.Fnorm();
        const Vector Blambda = x1 * Dinitial;
        const Vector xPlusV = x1 + v;

        const realdp *xData = x1.ObtainReadData();
        const realdp *xvData = xPlusV.ObtainReadData();

        std::vector<integer> id1; // MODIFIED: use std::vector index set for tCG
        std::vector<integer> id2; // MODIFIED: use std::vector index set for tCG

        id1.reserve(x1.Getlength());
        id2.reserve(x1.Getlength());

        for (integer k = 0; k < x1.Getlength(); ++k)
        {
            if (xvData[k] != 0 && std::abs(xData[k]) >= nv)
                id1.push_back(k);
            else
                id2.push_back(k);
        }

        const integer num1 = static_cast<integer>(id1.size()); // MODIFIED: avoid undefined ToIntegerSize
        const integer num2 = static_cast<integer>(id2.size()); // MODIFIED: avoid undefined ToIntegerSize

        max_innit = std::max<integer>(1, std::min<integer>(200,
            static_cast<integer>(std::lround(static_cast<realdp>(1.2) * num1))));

        Vector v_bar = GatherVector(v, id1); // MODIFIED: MATLAB v_bar = v(id1)
        Vector v_hat = GatherVector(v, id2); // MODIFIED: MATLAB v_hat = v(id2)

        const Vector tmp11 = x1.GetTranspose() * Blambda;



        // Initialize u = v + w for final search direction
        integer tcg_iter = 0;
        Vector u(n, r);
        u.SetToZeros();

        // tCG - solve tridiagonal CG subproblem for w
        Vector w = tCG(
            x1, egf, v, nv, id1, id2,
            Blambda, tmp11, t, &tcg_iter, &status);

        // Assemble u = v + w (MATLAB: u(id2) = v_hat; u(id1) = v_bar + w)
        // Since w is only defined on id1 subspace, u = v_bar + w at id1, v_hat at id2

        // First compute v_bar + w on id1 subspace
        Vector v_bar_plus_w = v_bar + w;  // MODIFIED: Element-wise addition on id1 subspace

        // MODIFIED: u(id2) = v_hat, using 0-based linear indices.
        for (integer i = 0; i < num2; ++i)
            u[id2[i]] = v_hat[i];

        // MODIFIED: u(id1) = v_bar + w, using 0-based linear indices.
        for (integer i = 0; i < num1; ++i)
            u[id1[i]] = v_bar_plus_w[i];

        nu = u.Fnorm(); // MODIFIED: MATLAB nu = norm(u, 'fro')



//update t
        if ((static_cast<realdp>(4.0) + static_cast<realdp>(1.0) / t) * nu < nv || status == 0)
            t = std::max(t0, w2 * t);
        else if (status <= 4)
            t = w1 * t;

        if (Verbose >= ITERRESULT && ((iter + 1) % OutputGap == 0))
        {
            printf("iter:%d, ssn_iter:%d, inner_it:%d, status:%d, t:%e, tau:%e\n",
                   iter + 1, ssn_iter, tcg_iter, status, t, tau);
            printf("iter:%d, f:%e, nv:%e, nu:%e, nonzero:%d, zero:%d\n",
                   iter + 1, f1, nv, nu, num1, num2); // MODIFIED: use defined counters
        }

        integer btiter = 0;
//>5
        if (status >= 5 || flag == 1)
        {
            flag = flag + 1;

            if (flag == 1)
            {
                Mani->Retraction(x1, u, &x2);
                nR++;
                f2 = Prob->f(x2) + Prob->g(x2);
                nf++;

                xold = x1;
                fold = f1;
                uold = u;
                nuold = nu;
                nvold = nv;
            }
            else
            {
                flag = 0;

                Mani->Retraction(x1, u, &x2);
                nR++;
                f2 = Prob->f(x2) + Prob->g(x2);
                nf++;

                if (f2 > fold - rho1 * nvold * nvold)
                {
                    realdp alpha = 1;
                    Vector au = uold*alpha;
                    Mani->Retraction(xold, au, &x2);
                    nR++;
                    f2 = Prob->f(x2) + Prob->g(x2);
                    nf++;

                    btiter = 0;
                    while (f2 > fold - rho1 * alpha * nuold * nuold && btiter < 3)
                    {
                        alpha *= rho2;
                        au = uold* alpha;
                        Mani->Retraction(xold, au, &x2);
                        nR++;
                        f2 = Prob->f(x2) + Prob->g(x2);
                        nf++;
                        btiter++;
                    }
                }
            }
        }
        else
        {
            flag = 0;

            realdp alpha = 1;
            Vector au = u*alpha;
            Mani->Retraction(x1, au, &x2);
            nR++;
            f2 = Prob->f(x2) + Prob->g(x2);
            nf++;

            btiter = 0;
            const realdp normDsquared = nu * nu;

            while (f2 > f1 - rho1 * alpha * normDsquared && btiter < 3)
            {
                alpha *= rho2;
                au = u* alpha;
                Mani->Retraction(x1, au, &x2);
                nR++;
                f2 = Prob->f(x2) + Prob->g(x2);
                nf++;
                btiter++;
            }
        }
//update 
        Vector gf2(n, r);
        Prob->EucGrad(x2, &gf2);
        ng++;

        err = std::min(nv, nu);
        iter++;

        ndir1 = nv;
        if (ndir0 == 0)
            ndir0 = ndir1;
//err report
        if (std::isnan(f2) || std::isinf(f2))
        {
            printf("RPNCG stopped: new function value is NaN or Inf.\n");
            break;
        }

        if (Verbose >= ITERRESULT)
        {
            if (iter % OutputGap == 0)
            {
                printf("iter:%d, f:%e, nv:%e, nu:%e, btiter:%d, time:%f, innertol:%e\n",
                       iter, f2, nv, nu, btiter,
                       static_cast<realdp>(getTickCount() - starttime) / CLK_PS,
                       innertol);
            }

        }

        x1 = x2;
        f1 = f2;
        egf = gf2; // MODIFIED: update gradient for next iteration

    }

    ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;


    }

    // ApplyBB - computes the Hessian-vector product for the RPN-CG algorithm
    // BB(y) = Hessian_term + quadratic_term + linear_term
    // where Hessian_term comes from the optimization problem's Hessian
    // quadratic_term comes from the low-rank matrix structure
    // linear_term comes from the current iterate
    Vector RPNCG::ApplyBB(const Vector &y,
           const Variable &x,
           const Problem *Prob,
           const Vector &Blambda,
           const Vector &tmp11) const
    { // MODIFIED: define the class member declared in RPNCG.h
    // MATLAB:
    // tmp  = y' * Blambda;
    // tmp1 = -2 * A' * (A * y);
    // tmp2 = y * tmp11;
    // tmp3 = x * ((tmp + tmp') / 2);
    // out  = tmp1 + tmp2 + tmp3;
    Vector tmp = y.GetTranspose()* Blambda;

    Vector tmp1(y.Getrow(), y.Getcol()); // MODIFIED: Hessian result has same shape as y
    Prob->EucHessianEta(x, y, &tmp1);

    Vector tmp2 = y * tmp11;
    Vector tmp3 = x * ((tmp + tmp.GetTranspose()) * static_cast<realdp>(0.5)); // MODIFIED: use non-mutating transpose

    return tmp1 + tmp2 + tmp3;
    }

// FindDir - semi-smooth Newton-CG method for finding search direction
    // Solves: argmin_p <p, egf> + 0.5 <p, Weight p> + h(R_x(p))
    // where h is the proximal operator and R_x is the retraction
    Vector RPNCG::FindDir(
    const Variable &x,
    const Vector &egf,
    realdp t,
    Vector *Dinitial,
    realdp innertol,
    integer *inneriter
) const {
    if (t <= 0 || Dinitial == nullptr || inneriter == nullptr)
        throw std::invalid_argument("RPNCG::FindDir: invalid step size or output pointer");

    const integer r = x.Getcol();
    const realdp lambda0 = static_cast<realdp>(0.2);
    const realdp nu = static_cast<realdp>(0.99);
    const realdp cg_tol = static_cast<realdp>(0.1);
    const realdp eta1 = static_cast<realdp>(0.2);
    const realdp eta2 = static_cast<realdp>(0.75);
    const realdp gamma1 = static_cast<realdp>(3.0);
    const realdp gamma2 = static_cast<realdp>(5.0);
    const realdp beta = static_cast<realdp>(0.1);
    const integer max_newton = 1000;
    const integer max_cg = 30;

    Vector z(*Dinitial);
    auto makeB = [&](const Vector &multiplier) -> Vector {
        return x - t * (egf - AstarOp(multiplier, x));
    };
    auto residual = [&](const Vector &multiplier, const Vector &B) -> Vector {
        Vector prox = ProxAt(B, t, Prob);
        return Aop(prox - x, x);
    };
    auto jacobian = [&](const Vector &B, const Vector &direction) -> Vector {
        return t * Aop(JacobianAt(B, AstarOp(direction, x), t, Prob), x);
    };

    Vector B = makeB(z);
    Vector F = residual(z, B);
    realdp normF = F.Fnorm();
    const integer q = 5; // nonmonotone window used in driver_rpncg.m
    std::vector<realdp> recent(static_cast<size_t>(q), normF);
    recent[static_cast<size_t>(q - 1)] = normF;
    realdp regularization = lambda0;
    integer iterations = 0;

    while (normF * normF > innertol && iterations < max_newton)
    {
        const realdp mu = regularization * std::max(std::min(normF, static_cast<realdp>(0.1)), static_cast<realdp>(1e-11));
        Vector d(z);
        d.SetToZeros();
        Vector cg_r(-1.0 * F);
        Vector cg_p(cg_r);
        realdp rr = cg_r.DotProduct(cg_r);
        const realdp cg_rhs_scale = regularization * normF;
        integer cg_iterations = 0;

        // The MATLAB stopping test is relative to ||d||, with d initialized
        // at zero.  Its first threshold is therefore zero; perform the first
        // CG step explicitly before applying that test.
        if (rr > 0)
        {
            Vector Ad = jacobian(B, cg_p) + mu * cg_p;
            const realdp denom = cg_p.DotProduct(Ad);
            if (denom > 0 && std::isfinite(denom))
            {
                const realdp alpha = rr / denom;
                d = d + alpha * cg_p;
                cg_r = cg_r - alpha * Ad;
                const realdp rr_new = cg_r.DotProduct(cg_r);
                if (rr > std::numeric_limits<realdp>::min())
                    cg_p = cg_r + (rr_new / rr) * cg_p;
                rr = rr_new;
                ++cg_iterations;
            }
        }

        while (std::sqrt(rr) > cg_tol * std::min(cg_rhs_scale * d.Fnorm(), static_cast<realdp>(1.0)) && cg_iterations < max_cg)
        {
            Vector Ad = jacobian(B, cg_p) + mu * cg_p;
            const realdp denom = cg_p.DotProduct(Ad);
            if (!(denom > 0) || !std::isfinite(denom))
                break;
            const realdp alpha = rr / denom;
            d = d + alpha * cg_p;
            cg_r = cg_r - alpha * Ad;
            const realdp rr_new = cg_r.DotProduct(cg_r);
            if (rr <= std::numeric_limits<realdp>::min())
                break;
            cg_p = cg_r + (rr_new / rr) * cg_p;
            rr = rr_new;
            ++cg_iterations;
        }

        Vector ztrial = z + d;
        Vector Btrial = makeB(ztrial);
        Vector Ftrial = residual(ztrial, Btrial);
        const realdp normTrial = Ftrial.Fnorm();

        if (normTrial < nu * *std::max_element(recent.begin(), recent.end()))
        {
            z = ztrial;
            B = Btrial;
            F = Ftrial;
            normF = normTrial;
            recent[static_cast<size_t>(iterations % q)] = normF;
        }
        else
        {
            const realdp dnorm2 = d.DotProduct(d);
            realdp rho = -Ftrial.DotProduct(d) / std::max(dnorm2, static_cast<realdp>(1e-30));
            if (rho >= eta1 && normTrial > 0)
            {
                Vector v = z - (Ftrial.DotProduct(z - ztrial) / (normTrial * normTrial)) * Ftrial;
                Vector Bv = makeB(v);
                Vector Fv = residual(v, Bv);
                if (Fv.Fnorm() <= normF)
                {
                    z = v; B = Bv; F = Fv; normF = F.Fnorm();
                }
                else
                {
                    z = z - beta * F;
                    B = makeB(z); F = residual(z, B); normF = F.Fnorm();
                }
            }
            if (rho >= eta2)
                regularization = std::max(regularization / 4, static_cast<realdp>(1e-5));
            else if (rho >= eta1)
                regularization *= (1 + gamma1) / 2;
            else
                regularization *= (gamma1 + gamma2) / 2;
        }
        ++iterations;
    }

    *Dinitial = z;
    *inneriter = iterations;
    return ProxAt(B, t, Prob) - x;
}





// tCG - Tridiagonal Conjugate Gradient for solving the RPN-CG subproblem
    // Solves: min l_x^T w + 0.5 w^T BB_1 w s.t. Bx_bar^T w = 0
    // Uses early stopping conditions and adaptive convergence criteria
    Vector RPNCG::tCG(
    const Variable &x,
    const Vector &egf,
    const Vector &v,
    realdp nv,
    const std::vector<integer> &id1, // MODIFIED: match RPNCG.h and Run()
    const std::vector<integer> &id2, // MODIFIED: match RPNCG.h and Run()
    const Vector &Blambda,
    const Vector &tmp11,
    realdp t,
    integer *inner_it,
    RPNCGTCGStatus *status
) const {
    // Initialize status to max iteration
    *status = RPNCG_MAXITER;
    *inner_it = 0;

    const integer n = v.Getrow();
    const integer r = v.Getcol();


    const integer id1_size = static_cast<integer>(id1.size()); // MODIFIED: std::vector size
    const integer id2_size = static_cast<integer>(id2.size()); // MODIFIED: std::vector size

    // With no active entries the constrained Newton correction has no
    // variables.  Avoid constructing a 0-by-d BLAS operand.
    if (id1_size == 0)
    {
        *status = RPNCG_EARLY2;
        return Vector(0);
    }

    // MODIFIED: MATLAB v_bar = v(id1), v_hat = v(id2).
    Vector v_bar = GatherVector(v, id1);
    Vector v_hat = GatherVector(v, id2);

    // z = BB(v)
    Vector z = ApplyBB(v, x, Prob, Blambda, tmp11);

    // tau_nv_hat_2 = tau * sum(v_hat.^2)
    realdp tau_nv_hat_2 = tau * v_hat.Fnorm() * v_hat.Fnorm(); // MODIFIED: v_hat is already gathered

    // v_Bv = v(:)' * z(:)
    realdp v_Bv = v.DotProduct(z);





    //inline



    // tol1 = v_Bv + tau_nv_hat_2
    realdp tol1 = v_Bv + tau_nv_hat_2;

    // tol2 = trace(gf' * v) + 0.5 * tol1
    // Since egf and v are vectors, trace(egf^T * v) = egf^T * v = dot product
    realdp trace_egfv = egf.DotProduct(v);  // Computes dot product
    // G_x(v)-G_x(0), including the nonsmooth part.  Omitting these terms
    // changes the first tCG acceptance test for sparse objectives.
    Variable xPlusV(x);
    xPlusV = x + v;
    realdp tol2 = trace_egfv + 0.5 * tol1
        + Prob->g(xPlusV) - Prob->g(x);





//cheack



    // Early stopping check 1: G_x(v) > G_x(0)
    if (tol2 > 0) {
        *status = RPNCG_EARLY1;
        Vector result(id1_size);
        result.SetToZeros();
        return result;
    }

    // Early stopping check 2
    realdp nv_squared = nv * nv;
    if (tol1 < gamma * nv_squared) {
        *status = RPNCG_EARLY2;
        Vector result(id1_size);
        result.SetToZeros();
        return result;
    }







//buildbbar




// Build Bbar = Bvtmp(id1, :)
Vector Bbar = BuildBbar(x, id1);

// D = Bbar' * Bbar
Vector D = Bbar.GetTranspose() * Bbar;

// Dinv = pinv(D)
Vector Dinv = D.PseudoInverseMatrix();




    // Projection function P_x(eta) = eta - Bbar * (Dinv * (Bbar' * eta))
    // projects eta onto the null space of Bbar
    auto Px = [&](const Vector &eta) -> Vector {
        Vector Bt_eta = Bbar.GetTranspose() * eta;
        Vector Dinv_Bt_eta = Dinv * Bt_eta;
        Vector B_Dinv_Bt_eta = Bbar * Dinv_Bt_eta;
        return eta - B_Dinv_Bt_eta;
    };





    // l_x = -1/t * v_bar + z(id1) - linear term for CG subproblem
    Vector z_slice = GatherVector(z, id1); // MODIFIED: MATLAB z(id1), not z(1:id1_size)
    Vector l_x = (-1.0 / t) * v_bar + z_slice; // MODIFIED: use gathered v_bar directly

    // Initialize CG variables for tridiagonal CG
    Vector w(id1_size);
    w.SetToZeros();

    Vector r_cg = Px(l_x);  // initial residual
    Vector o = r_cg * (-1.0);  // search direction
    realdp delta = r_cg.Fnorm() * r_cg.Fnorm();  // ||r||^2
    realdp r_r = delta;
    realdp norm_r0 = r_cg.Fnorm();  // initial norm

    Vector tt = z;






//loop



    // CG iterations - main loop of tridiagonal CG algorithm
    for (integer j = 0; j < max_innit; j++) {



        // Build VV: zero matrix with ||    VV(id1) = o    ||- trial direction
        Vector VV(n, r);
        VV.SetToZeros();
        realdp *VVptr = VV.ObtainWritePartialData();
        const realdp *optr = o.ObtainReadData();
        for (integer i = 0; i < id1_size; i++) {
            VVptr[id1[i]] = optr[i];
        }



        // ppp = BB(VV) - apply BB operator to trial direction
        Vector ppp = ApplyBB(VV, x, Prob, Blambda, tmp11);



        // p = ppp(id1) - restrict to id1 subspace
        Vector p(id1_size);
        p.SetToZeros();
        const realdp *pppptr = ppp.ObtainReadData();
        realdp *pptr = p.ObtainWritePartialData();
        for (integer i = 0; i < id1_size; i++) {
            pptr[i] = pppptr[id1[i]];
        }

        // o_q = o(:)' * p(:) - curvature condition
        realdp o_q = o.DotProduct(p);

        // Negative curvature check - if insufficient positive curvature, stop
        if (o_q <= vartheta * delta) {
            *status = RPNCG_NEGATIVE_CURVATURE;
            break;
        }

        // Compute step size alpha - using Fletcher-Reeves formula
        realdp alpha = r_r / o_q;

        // Update w - new iterate (manual loop required since we're updating specific elements)
        Vector w_old = w;

        w=w+alpha*o;

        // Update residual: r = r + alpha * q, where q = P_x(p) - projected search direction
        Vector q = Px(p);
        r_cg = r_cg + alpha * q;

        // Build d for error checking - d = v_bar + w at id1, v_hat at id2
        // Using Element.h operations for cleaner code
        Vector d_vec(n, r);
        d_vec.SetToZeros();


// d(id1) = v_bar + w
for (integer k = 0; k < id1_size; ++k)
{
    d_vec[id1[k]] = v_bar[k]+w[k];   // id1 should be 0-based
}

// d(id2) = v_hat
for (integer k = 0; k < id2_size; ++k)
{
    d_vec[id2[k]] = v_hat[k];  // id2 should be 0-based
}






        // Update tt - track accumulated quadratic form
        tt = tt + alpha * ppp;
        // MODIFIED: MATLAB reuses variable name t for tt here; C++ keeps scalar t unchanged.

        // err_1 = d(:)' * tt(:) + tau_nv_hat_2 - energy condition 1
        realdp err_1 = d_vec.DotProduct(tt) + tau_nv_hat_2;


        //x_d = x.main+d;

        // err_2 = gf(:)' * d(:) + 0.5 * err_1 - energy condition 2 (gradient condition)
        Variable xPlusD(x);
        xPlusD = x + d_vec;
        realdp err_2 = egf.DotProduct(d_vec) + 0.5 * err_1
            + Prob->g(xPlusD) - Prob->g(x);


        //err1 and 2 finished

        

        // Early stopping check 3 - if energy conditions violated, backtrack
        realdp d_norm_squared = d_vec.Fnorm() * d_vec.Fnorm();
    
        if (err_1 < gamma * d_norm_squared || err_2 > 0) {
            w = w_old;
            *status = RPNCG_EARLY3;
            break;
        }

        // Update residual norm - monitor convergence
        realdp rold_rold = r_r;
        r_r = r_cg.Fnorm() * r_cg.Fnorm();
        realdp norm_r = r_cg.Fnorm();






        // Stopping criterion check - adaptive tolerance based on progress
        // MODIFIED: MATLAB uses min(norm_r0^theta, kappa).
        realdp threshold = std::min(std::pow(norm_r0, theta), kappa);



        // Check convergence - linear vs superlinear
        if (norm_r <= norm_r0 * threshold) {
            if (kappa < std::pow(norm_r0, theta)) {
                *status = RPNCG_LINEAR;      // kappa stopping condition
            } else {
                *status = RPNCG_SUPERLINEAR; // theta stopping condition
            }
            break;
        }

        // Update search direction
        realdp beta = r_r / rold_rold;
        o = (-1.0) * r_cg + beta * o; // MODIFIED: r_cg and o are already id1-sized vectors
        delta = r_r + beta * beta * delta;  // update direction norm

        *inner_it = j + 1;
    }

    // If max iteration reached, return zero solution
    if (*status == RPNCG_MAXITER) {
        Vector result(id1_size);
        result.SetToZeros();
        return result;
    }

    return w;  // converged solution from CG
}


//func














static Vector BuildBbar(const Variable &x,
                        const std::vector<integer> &id1)
{ // MODIFIED: implementation of MATLAB tmpp/Bvtmp/Bbar block without explicitly forming 3D tmpp
    const integer n = x.Getrow();
    const integer m = x.Getcol();
    const integer id1_size = static_cast<integer>(id1.size());
    const integer d = m * (m + 1) / 2;

    const realdp invsqrt2 = 1.0 / std::sqrt(2.0);

    Vector Bbar(id1_size, d);
    Bbar.SetToZeros();

    const realdp *xptr = x.ObtainReadData();
    realdp *bptr = Bbar.ObtainWritePartialData();

    integer basis = 0;

    // Diagonal symmetric basis:
    // S_i has S(i,i) = 1
    for (integer i = 0; i < m; ++i)
    {
        for (integer k = 0; k < id1_size; ++k)
        {
            integer lin = id1[k];      // 0-based linear index in n-by-m matrix
            integer row = lin % n;
            integer col = lin / n;

            if (col == i)
            {
                bptr[cmidx(k, basis, id1_size)] =
                    xptr[cmidx(row, i, n)];
            }
        }

        basis++;
    }

    // Off-diagonal symmetric basis:
    // S_ij = (e_i e_j^T + e_j e_i^T) / sqrt(2)
    for (integer i = 0; i < m; ++i)
    {
        for (integer j = i + 1; j < m; ++j)
        {
            for (integer k = 0; k < id1_size; ++k)
            {
                integer lin = id1[k];
                integer row = lin % n;
                integer col = lin / n;

                realdp value = 0.0;

                if (col == i)
                {
                    value = xptr[cmidx(row, j, n)] * invsqrt2;
                }
                else if (col == j)
                {
                    value = xptr[cmidx(row, i, n)] * invsqrt2;
                }

                bptr[cmidx(k, basis, id1_size)] = value;
            }

            basis++;
        }
    }

    return Bbar;
}


}; /*end of ROPTLIB namespace*/
