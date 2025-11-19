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

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laminarBurningVelocity> Foam::laminarBurningVelocity::New
(
    const dictionary& reactRateProperties,
    const reactionRate& reactRate
)
{
    const dictionary& laminarBurningVelocityDict =
        reactRateProperties.subDict("laminarBurningVelocity");

    const word modelType(laminarBurningVelocityDict.lookup("model"));

    Info<< "Selecting laminar correlation "
        << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(reactRateProperties)
            << "Unknown laminar burning velocity correlation "
            << modelType << nl << nl
            << "Valid laminar burning velocity correlations are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<laminarBurningVelocity>
        (cstrIter()(laminarBurningVelocityDict.optionalSubDict(modelType), reactRate));
}


// ************************************************************************* //
