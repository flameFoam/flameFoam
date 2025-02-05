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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::autoIgnition> 
Foam::autoIgnition::New
(
    const reactionRate& reactRate,
    const dictionary& dict
)
{
    word autoIgnitionType
    (
        dict.lookup("autoIgnition")
    );

    Info<< "Selecting auto-ignition model "
        << autoIgnitionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(autoIgnitionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown auto-ignition model "
            << autoIgnitionType << endl << endl
            << "Valid auto-ignition models are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    const label tempOpen = autoIgnitionType.find('<');

    const word className = autoIgnitionType(0, tempOpen);

    return autoPtr<autoIgnition>
        (cstrIter()(className, reactRate, dict));
}


// ************************************************************************* //
