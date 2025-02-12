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



#include "autoIgnition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoIgnition, 0);
    defineRunTimeSelectionTable
    (
        autoIgnition,
        dictionary
    );

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoIgnition::autoIgnition
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
    debug_(coeffDict_.lookupOrDefault("debug", false)),
    tau_
    (
        IOobject
        (
            "tau",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("tau", dimTime, 0.0)
    )
{
    Info << "flameFoam autoIgnition object initialized" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoIgnition::~autoIgnition()
{}

// ************************************************************************* // 