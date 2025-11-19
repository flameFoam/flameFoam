/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2025 Lithuanian Energy Institute

 -------------------------------------------------------------------------------
License
    This file is part of flameFoam, derivative work of OpenFOAM.

    flameFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    flameFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    <http://www.gnu.org/licenses/> for more details.

Disclaimer
    flameFoam is not approved or endorsed by neither the OpenFOAM Foundation
    Limited nor OpenCFD Limited.

\*---------------------------------------------------------------------------*/

#include "TFC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{
    defineTypeNameAndDebug(TFC, 0);
    addToRunTimeSelectionTable
    (
        reactionRate,
        TFC,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateModels::TFC::TFC
(
    const dictionary& dict,
    const combustionModel& combModel
)
:
    reactionRate(combModel),
    turbulentCorrelation_(
        turbulentBurningVelocity::New
        (
            dict,
            *this
        )
    )

{
    appendInfo("Reaction rate model: TFC");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateModels::TFC::~TFC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionRateModels::TFC::correct
(
)
{
    if (debug_)
    {
        Info << "\tTFC correct:" << endl;
        Info << "\t\tInitial min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
    }

    this->correctUnburntProperties();

    turbulentCorrelation_->correct();
    cSource_ = rhoU()*max(turbulentCorrelation_->burningVelocity(), turbulentCorrelation_->getLaminarBurningVelocity())*mag(fvc::grad(combModel_.thermo().Y("c")));

    if (debug_)
    {
        Info << "\t\tObtained min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
        Info << "\t\tTFC correct finished" << endl;
    }
}

char const *Foam::reactionRateModels::TFC::getInfo()
{
    infoString_.append(turbulentCorrelation_().getInfo());
    turbulentCorrelation_().clearInfo();
    return infoString_.c_str();
}


// ************************************************************************* //
