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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::turbulentBurningVelocity> Foam::turbulentBurningVelocity::New
(
    const reactionRate& reactRate,
    const dictionary& dict
)
{
    word turbulentBurningVelocityType
    (
        dict.lookup("turbulentBurningVelocity")
    );

    Info<< "Selecting turbulent correlation "
        << turbulentBurningVelocityType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(turbulentBurningVelocityType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown turbulent burning velocity correlation "
            << turbulentBurningVelocityType << endl << endl
            << "Valid turbulent burning velocity correlations are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    const label tempOpen = turbulentBurningVelocityType.find('<');

    const word className = turbulentBurningVelocityType(0, tempOpen);

    return autoPtr<turbulentBurningVelocity>
        (cstrIter()(className, reactRate, dict));
}


// ************************************************************************* //
