/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2024 Lithuanian Energy Institute

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

#include "FSD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{
    defineTypeNameAndDebug(FSD, 0);
    addToRunTimeSelectionTable
    (
        reactionRate,
        FSD,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateModels::FSD::FSD
(
    const word modelType,
    const dictionary& dict,
    const combustionModel& combModel
)
:
    reactionRate(modelType, dict, combModel),
    wrinklingCorrelation_(
        wrinklingFactor::New
        (
            *this,
            dict
        )
    )
{
    appendInfo("Reaction rate model: FSD");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateModels::FSD::~FSD()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::reactionRateModels::FSD::correct
(
)
{
    if (debug_)
    {
        Info << "\tFSD correct:" << endl;
        Info << "\t\tInitial min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
    }

    wrinklingCorrelation_->correct();

    cSource_ = rhoU()*wrinklingCorrelation_->burningVelocity()*mag(fvc::grad(combModel_.thermo().Y("c")));

    if (debug_)
    {
        Info << "\t\tObtained min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
        Info << "\t\tFSD correct finished" << endl;
    }
}

char const *Foam::reactionRateModels::FSD::getInfo()
{
    infoString_.append(wrinklingCorrelation_().getInfo());
    wrinklingCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
