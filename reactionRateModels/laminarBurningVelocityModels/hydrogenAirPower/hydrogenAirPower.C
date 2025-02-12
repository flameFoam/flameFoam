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

#include "hydrogenAirPower.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{
    defineTypeNameAndDebug(hydrogenAirPower, 0);
    addToRunTimeSelectionTable
    (
        laminarBurningVelocity,
        hydrogenAirPower,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::hydrogenAirPower::hydrogenAirPower
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    laminarBurningVelocity(modelType, reactRate, dict),
    X_H2_0_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2_0")),
    sLaminar0_
    (
        dimensionedScalar
        (
            dimVelocity,
          - 488.9*pow(X_H2_0_, 4)
          + 285.0*pow(X_H2_0_, 3)
          - 21.92*pow(X_H2_0_, 2)
          + 1.351*X_H2_0_
         - 0.04
        )
    ),
    pRef_(dimensionedScalar(dimPressure, 101300)),
    TRef_(dimensionedScalar(dimTemperature, 298)),
    p_(mesh_.lookupObject<volScalarField>("p"))
{
    appendInfo("\tLBV estimation method: power law correlation");
    appendInfo("\tObtained S_L_0: " << sLaminar0_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::hydrogenAirPower::~hydrogenAirPower()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarBurningVelocityModels::hydrogenAirPower::correct
()
{
    if (debug_)
    {
        Info << "\t\t\tPower law correct:" << endl;
        Info << "\t\t\t\tInitial average S_L: "  << average(sLaminar_).value() << endl;
    }

    sLaminar_ = sLaminar0_*pow(reactionRate_.TU()/TRef_, 1.75)*pow(p_/pRef_, -0.2);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tPower law correct finished" << endl;
    }
}

// ************************************************************************* //
