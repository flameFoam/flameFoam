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

#include "LBVPower.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{
    defineTypeNameAndDebug(LBVPower, 0);
    addToRunTimeSelectionTable
    (
        laminarBurningVelocity,
        LBVPower,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::LBVPower::LBVPower
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    laminarBurningVelocity(modelType, reactRate, dict),
    X_H2_0_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2_0")),
    X_H2O_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2O")),

    powerDil_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("powerDil")),
    powerT_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("powerT")),
    powerP_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("powerP")),

    ER_(0.705*X_H2_0_/(0.295*(1-X_H2_0_-X_H2O_))),
    sLaminar0_(dimensionedScalar(dimVelocity, 1.44*ER_*ER_+1.07*ER_-0.29)),
    pRef_(dimensionedScalar(dimPressure, 100000)),
    TRef_(dimensionedScalar(dimTemperature, 298)),
    p_(mesh_.lookupObject<volScalarField>("p"))
{
    appendInfo("\tLBV estimation method: LBVPower correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::LBVPower::~LBVPower()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarBurningVelocityModels::LBVPower::correct
()
{
    if (debug_)
    {
        Info << "\t\t\tLBVPower correct:" << endl;
        Info << "\t\t\t\tInitial average S_L: "  << average(sLaminar_).value() << endl;
    }

    sLaminar_ = sLaminar0_ * pow(1 - X_H2O_, powerDil_) * pow(reactionRate_.TU() / TRef_, powerT_) * pow(p_ / pRef_, powerP_);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tLBVPower correct finished" << endl;
    }
}


