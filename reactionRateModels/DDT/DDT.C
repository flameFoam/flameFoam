/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

\*---------------------------------------------------------------------------*/

#include "DDT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "zero.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{
    defineTypeNameAndDebug(DDT, 0);
    addToRunTimeSelectionTable
    (
        reactionRate,
        DDT,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateModels::DDT::DDT
(
    const dictionary& dict,
    const combustionModel& combModel
)
:
    reactionRate(combModel),
    wrinklingCorrelation_
    (
        wrinklingFactor::New(dict, *this)
    ),
    autoIgnition_
    (
        autoIgnition::New(dict, *this)
    ),
    c_(combModel_.thermo().Y("c")),
    rho_(combModel_.rho()),
    tIgn_
    (
        "tIgn",
        dimTime,
        dict.lookupOrDefault<scalar>("tIgn", 0.15e-3)
    )
{
    appendInfo("Reaction rate model: DDT");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateModels::DDT::~DDT()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionRateModels::DDT::correct()
{
    if (debug_)
    {
        Info<< "\tDDT correct:" << endl;
        Info<< "\t\tInitial min/avg/max cSource: "
            << min(cSource_).value() << " "
            << average(cSource_).value() << " "
            << max(cSource_).value() << endl;
    }

    wrinklingCorrelation_->correct();
    autoIgnition_->correct();

    cSource_ =
        rhoU()*wrinklingCorrelation_->burningVelocity()*mag(fvc::grad(c_))
      + rho_*(1 - c_)*max(Zero, autoIgnition_->tau() - 1)/tIgn_;

    if (debug_)
    {
        Info<< "\t\tObtained min/avg/max cSource: "
            << min(cSource_).value() << " "
            << average(cSource_).value() << " "
            << max(cSource_).value() << endl;
        Info<< "\t\tDDT correct finished" << endl;
    }
}


// ************************************************************************* //
