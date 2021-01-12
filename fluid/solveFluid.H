if (pimples.nCorrPimple() <= 1)
{
    #include "rhoEqn.H"
}

#include "UEqn.H"
if (c.average().value() != 1)
{
    #include "cEqn.H"
}

#include "EEqn.H"

// --- Pressure corrector loop
while (pimple.correct())
{
    #include "pEqn.H"
}
if (pimples.pimpleTurbCorr(i))
{
    turbulence.correct();
}
if (pimples.finalPimpleIter())
{
    rho = thermo.rho();
}
