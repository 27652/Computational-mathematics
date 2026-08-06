#ifndef RPNCG_H
#define RPNCG_H

#include "Element.h"
#include "Problem.h"
#include "Solvers/SolversNSM.h"
#include "Others/def.h"
#include "f2c.h"

namespace ROPTLIB {

enum RPNCGTCGStatus {
    RPNCG_MAXITER = -1,
    RPNCG_EARLY1 = 0,
    RPNCG_EARLY2 = 1,
    RPNCG_NEGATIVE_CURVATURE = 2,
    RPNCG_EARLY3 = 3,
    RPNCG_LINEAR = 4,
    RPNCG_SUPERLINEAR = 5
};


struct RPNCGParams
{
    realdp vartheta;
    realdp gamma;
    realdp rho1;
    realdp rho2;
    realdp w1;
    realdp w2;
    realdp theta;
    realdp kappa;
    integer max_innit;

    RPNCGParams()
        : vartheta(0.01),
          gamma(static_cast<realdp>(0.01)),
          rho1(static_cast<realdp>(0.001)),
          rho2(static_cast<realdp>(0.5)),
          w1(static_cast<realdp>(1.1)),
          w2(static_cast<realdp>(0.9)),
          theta(static_cast<realdp>(0.5)),
          kappa(static_cast<realdp>(0.1)),
          max_innit(200)
    {
    }
};


////////mainnnnnnnnnn///
class RPNCG : public SolversNSM {
public:
    RPNCG(const Problem *prob, const Variable * initialx);
    virtual ~RPNCG();

    virtual void SetDefaultParams();
    virtual void Run();

    realdp gamma;
    realdp rho1;
    realdp rho2;
    realdp w1;
    realdp w2;
    realdp theta;
    realdp kappa;
    realdp vartheta;
    realdp tau;
    integer max_innit=200;
    integer max_inner_iter;// param 

protected:// func
    virtual void UpdateData();
    virtual void PrintInfo();
    virtual void PrintFinalInfo();

private://tcg finddir func
    Vector FindDir(
        const Variable &x,
        const Vector &egf,
        realdp t,
        Vector *Dinitial,
        realdp innertol,
        integer *inneriter
    ) const;

    Vector ApplyBB(
        const Vector &y,
        const Variable &x,
        const Problem *Prob,
        const Vector &Blambda,
        const Vector &tmp11
    ) const;

    Vector tCG(
        const Variable &x,
        const Vector &egf,
        const Vector &v,
        realdp nv,
const std::vector<integer> &id1,
const std::vector<integer> &id2,
        const Vector &Blambda,
        const Vector &tmp11,
        realdp t,
        integer *inner_it,
        RPNCGTCGStatus *status
    ) const;
};

};

#endif
