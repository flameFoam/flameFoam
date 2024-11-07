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

#include "turbulentBurningVelocity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turbulentBurningVelocity, 0);
    defineRunTimeSelectionTable(turbulentBurningVelocity, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocity::turbulentBurningVelocity
(
    const word& modelType,
    const reactionRate& reactRate,
    const dictionary& dict
)
:
    reactionRate_(reactRate),
    coeffDict_(dict),
    mesh_(reactionRate_.mesh()),
    combModel_(reactionRate_.combModel()),
    debug_(dict.lookupOrDefault("debug", false)),
    debugFields_(dict.lookupOrDefault("debugFields", false)),
    sTurbulent_
    (
        IOobject
        (
            "TBV",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            debugFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("TBV", dimVelocity, Zero)
    ),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            reactRate,
            combModel_.coeffs()
        )
    )
{
    Info << "flameFoam turbulentBurningVelocity object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocity::~turbulentBurningVelocity()
{}


// NEAIŠKU AR REIKIA, GAL PRAVERS
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::turbulentBurningVelocity::saneEpsilon()
{
    return max(combModel_.turbulence().epsilon(), dimensionedScalar(dimVelocity*dimAcceleration, SMALL));
}

const Foam::volScalarField& Foam::turbulentBurningVelocity::getLaminarBurningVelocity()
{
    return laminarCorrelation_().burningVelocity();
}



// bool Foam::turbulentBurningVelocity::read(const dictionary& dict)
// {
//     dict.lookup("fuel") >> fuel_;
//
//     return true;
// }


// ************************************************************************* //
