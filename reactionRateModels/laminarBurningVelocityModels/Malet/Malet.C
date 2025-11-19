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

#include "Malet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{
    defineTypeNameAndDebug(Malet, 0);
    addToRunTimeSelectionTable
    (
        laminarBurningVelocity,
        Malet,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::Malet::Malet
(
    const dictionary& dict,
    const reactionRate& reactRate
):
    laminarBurningVelocity(reactRate),
    X_H2_0_("X_H2_0", dimless, combustionProperties_),
    X_H2O_("X_H2O", dimless, combustionProperties_),
    a2_(dimVelocity, 1.44),
    a1_(dimVelocity, 1.07),
    a0_(dimVelocity, -0.29),
    ER_(0.705*X_H2_0_/(0.295*(1-X_H2_0_-X_H2O_))),
    sLaminar0_((a2_*ER_*ER_ + a1_*ER_ + a0_) * pow(1-X_H2O_,4)),
    pRef_(dimPressure, 100000),
    TRef_(dimTemperature, 298)
{
    appendInfo("\tLBV estimation method: Malet correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::Malet::~Malet()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarBurningVelocityModels::Malet::correct
()
{
    if (debug_)
    {
        Info << "\t\t\tMalet correct:" << endl;
        Info << "\t\t\t\tInitial average S_L: "  << average(sLaminar_).value() << endl;
    }

    const fvMesh& mesh(reactionRate_.mesh());
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    sLaminar_ = sLaminar0_*pow(reactionRate_.TU()/TRef_,2.2)*pow(p/pRef_,-0.5);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tMalet correct finished" << endl;
    }

}


// ************************************************************************* //
