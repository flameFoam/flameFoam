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

#include "wrinklingFactor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wrinklingFactor, 0);
    defineRunTimeSelectionTable(wrinklingFactor, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wrinklingFactor::wrinklingFactor
(
    const reactionRate& reactRate
)
:
    reactionRate_(reactRate),
    combustionProperties_(reactionRate_.combModel().coeffs()),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            combustionProperties_.subDict("reactionRate"),
            reactRate
        )
    ),
    sTurbulent_
    (
        IOobject
        (
            "TBV",
            reactionRate_.mesh().time().name(),
            reactionRate_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        reactionRate_.mesh(),
        dimensionedScalar(dimVelocity, 0)
    ),
    debug_(combustionProperties_.lookupOrDefault("debug", false))
{
    Info << "flameFoam wrinklingFactor object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wrinklingFactor::~wrinklingFactor()
{}

// ************************************************************************* //
