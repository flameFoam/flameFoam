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
    const dictionary& dict,
    const reactionRate& reactRate
):
    laminarBurningVelocity(reactRate),
    X_H2_0_("X_H2_0", dimless, combustionProperties_),
    c4_(dimVelocity, -488.9),
    c3_(dimVelocity, 285.0),
    c2_(dimVelocity, -21.92),
    c1_(dimVelocity, 1.351),
    c0_(dimVelocity, -0.04),
    sLaminar0_((c4_*pow(X_H2_0_, 4)+c3_*pow(X_H2_0_, 3)+c2_*pow(X_H2_0_, 2)+c1_*X_H2_0_+c0_)),
    pRef_(dimensionedScalar(dimPressure, 101300)),
    TRef_(dimensionedScalar(dimTemperature, 298))
{
    reactionRate_.appendInfo("\tLBV estimation method: power law correlation");
    OStringStream os;
    os << "Obtained SL_0: " << sLaminar0_;
    reactionRate_.appendInfo(os.str());
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

    const fvMesh& mesh(reactionRate_.mesh());
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    sLaminar_ = sLaminar0_*pow(reactionRate_.TU()/TRef_, 1.75)*pow(p/pRef_, -0.2);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tPower law correct finished" << endl;
    }
}

// ************************************************************************* //
