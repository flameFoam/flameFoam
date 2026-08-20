/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

\*---------------------------------------------------------------------------*/

#include "autoIgnition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoIgnition, 0);
    defineRunTimeSelectionTable(autoIgnition, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoIgnition::autoIgnition
(
    const reactionRate& reactRate
)
:
    reactionRate_(reactRate),
    combustionProperties_(reactRate.combModel().coeffs()),
    mesh_(reactRate.mesh()),
    combModel_(reactRate.combModel()),
    debug_(reactRate.debugSwitch()),
    tau_
    (
        IOobject
        (
            "tau",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    Info<< "flameFoam autoIgnition object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoIgnition::~autoIgnition()
{}


// ************************************************************************* //
