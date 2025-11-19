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
    const dictionary& dict,
    const reactionRate& reactRate
):
    laminarBurningVelocity(reactRate),
    X_H2_0_("X_H2_0", dimless, combustionProperties_),
    X_H2O_("X_H2O", dimless, combustionProperties_),

    powerDil_("powerDil", dimless, dict),
    powerT_("powerT", dimless, dict),
    powerP_("powerP", dimless, dict),

    a2_("a2", dimVelocity, dict),
    a1_("a1", dimVelocity, dict),
    a0_("a0", dimVelocity, dict),

    ER_(0.705*X_H2_0_/(0.295*(1-X_H2_0_-X_H2O_))),
    sLaminar0_((a2_*ER_*ER_ + a1_*ER_ + a0_) * pow(1 - X_H2O_, powerDil_)),
    pRef_("pRef", dimPressure, dict),
    TRef_("TRef", dimTemperature, dict)
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

    const fvMesh& mesh(reactionRate_.mesh());
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    sLaminar_ = sLaminar0_ * pow(reactionRate_.TU() / TRef_, powerT_) * pow(p / pRef_, powerP_);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tLBVPower correct finished" << endl;
    }
}


