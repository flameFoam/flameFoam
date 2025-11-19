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
    const reactionRate& reactRate
)
:
    reactionRate_(reactRate),
    combModel_(reactionRate_.combModel()),
    combustionProperties_(combModel_.coeffs()),
    debug_(combustionProperties_.lookupOrDefault("debug", false)),
    debugFields_(combustionProperties_.lookupOrDefault("debugFields", false)),
    sTurbulent_
    (
        IOobject
        (
            "TBV",
            reactionRate_.mesh().time().name(),
            reactionRate_.mesh(),
            IOobject::NO_READ,
            debugFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        reactionRate_.mesh(),
        dimensionedScalar(dimVelocity, 0)
    ),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            combustionProperties_.subDict("reactionRate"),
            reactRate
        )
    )
{
    Info << "flameFoam turbulentBurningVelocity object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocity::~turbulentBurningVelocity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::turbulentBurningVelocity::saneEpsilon()
{
    return max(combModel_.turbulence().epsilon(), dimensionedScalar(dimVelocity*dimAcceleration, SMALL));
}

const Foam::volScalarField& Foam::turbulentBurningVelocity::getLaminarBurningVelocity()
{
    return laminarCorrelation_().burningVelocity();
}


// ************************************************************************* //
