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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::turbulentBurningVelocity> Foam::turbulentBurningVelocity::New
(
    const dictionary& reactionRateProperties,
    const reactionRate& reactRate
)
{
    const dictionary& turbulentBurningVelocityProperties =
        reactionRateProperties.subDict("turbulentBurningVelocity");

    const word modelType(turbulentBurningVelocityProperties.lookup("model"));

    Info<< "Selecting turbulent correlation "
        << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(turbulentBurningVelocityProperties)
            << "Unknown turbulent burning velocity correlation "
            << modelType << nl << nl
            << "Valid turbulent burning velocity correlations are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<turbulentBurningVelocity>
        (cstrIter()(turbulentBurningVelocityProperties.optionalSubDict(modelType), reactRate));
}


// ************************************************************************* //
