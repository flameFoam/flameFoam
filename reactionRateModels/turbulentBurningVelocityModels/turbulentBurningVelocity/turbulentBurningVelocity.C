/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
    sTurbulent_
    (
        IOobject
        (
            "TBV",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TBV", dimVelocity, Zero)
    ),
    debug_(dict.lookupOrDefault("debug", false))
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



// bool Foam::turbulentBurningVelocity::read(const dictionary& dict)
// {
//     dict.lookup("fuel") >> fuel_;
//
//     return true;
// }


// ************************************************************************* //
