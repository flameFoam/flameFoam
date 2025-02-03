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
#include "addToRunTimeSelectionTable.H"

//************************************* Static Data Members **************************************//

namespace Foam
{
    defineTypeNameAndDebug(autoIgnition, 0);
    defineRunTimeSelectionTable(autoIgnition, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoIgnition::autoIgnition
(
    const DDT& ddt,
    const dictionary& dict
)
:
    Tref_(dict.lookupOrDefault<scalar>("Tref", 298.15)),
    pref_(dict.lookupOrDefault<scalar>("pref", 101325.0)),
    Ea_(dict.lookupOrDefault<scalar>("Ea", 15000.0)),
    A_(dict.lookupOrDefault<scalar>("A", 0.1))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoIgnition::~autoIgnition()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::autoIgnition::tau() const
{
    return 0.0;
}

// ************************************************************************* // 