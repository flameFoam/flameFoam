/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laminarBurningVelocity> Foam::laminarBurningVelocity::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel,
    const reactionRate& reactRate
)
{
    word laminarBurningVelocityType
    (
        dict.lookup("laminarBurningVelocity")
    );

    Info<< "Selecting laminar correlation "
        << laminarBurningVelocityType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(laminarBurningVelocityType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown laminar burning velocity correlation "
            << laminarBurningVelocityType << endl << endl
            << "Valid laminar burning velocity correlations are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    const label tempOpen = laminarBurningVelocityType.find('<');

    const word className = laminarBurningVelocityType(0, tempOpen);

    return autoPtr<laminarBurningVelocity>
        (cstrIter()(className, dict, mesh, combModel, reactRate));
}


// ************************************************************************* //
