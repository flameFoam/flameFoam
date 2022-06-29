if (pimples.firstPimpleIter())
{
    #include "rhoEqn.H"
}

if (pimple.models())
{
    fvModels.correct();
}

#include "UEqn.H"

tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

mvConvection = tmp<fv::convectionScheme<scalar>>
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,c_h)")
    )
);

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
    thermophysicalTransport.correct();
}
if (pimples.finalPimpleIter())
{
    rho = thermo.rho();
}
