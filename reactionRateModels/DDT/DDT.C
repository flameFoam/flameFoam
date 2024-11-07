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

#include "DDT.H"
#include "addToRunTimeSelectionTable.H"

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
    ),
    autoIgnition_(
        autoIgnition::New
        (
            *this,
            dict
        )

    ),
    c_(combModel_.thermo().Y("c")),
    rho_(combModel_.rho()),
    tIgn_(dimensionedScalar("tIgn", dimTime, 0.15E-3))
{
    appendInfo("Reaction rate model: DDT");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateModels::DDT::~DDT()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::reactionRateModels::DDT::correct
(
)
{
    if (debug_)
    {
        Info << "\tDDT correct:" << endl;
        Info << "\t\tInitial min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
    }

    wrinklingCorrelation_->correct();
    autoIgnition_->correct();

    cSource_ = rhoU()*wrinklingCorrelation_->burningVelocity()*mag(fvc::grad(c_))
             + rho_*(1-c_)*max(Zero, autoIgnition_->tau()-1)/tIgn_;

    if (debug_)
    {
        Info << "\t\tObtained min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
        Info << "\t\tDDT correct finished" << endl;
    }
}

char const *Foam::reactionRateModels::DDT::getInfo()
{
    infoString_.append(wrinklingCorrelation_().getInfo());
    wrinklingCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
