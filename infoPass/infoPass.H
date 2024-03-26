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

Class
    Foam::infoPass

Description
    Helper class to log flameFoam model data

SourceFiles
    infoPass.C

\*---------------------------------------------------------------------------*/

#ifndef infoPass_H
#define infoPass_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class infoPass Declaration
\*---------------------------------------------------------------------------*/

class infoPass
{
protected:

    // Protected data

        string infoString_;

public:

    //- Runtime type information
    TypeName("infoPass");


    // Constructors

        //- Construct from components
        infoPass();

        //- Disallow default bitwise copy construction
        infoPass(const infoPass&) = delete;


    //- Destructor
    virtual ~infoPass();


    // Member Functions

        //- Return infoString_
        virtual char const *getInfo();

        //- Append new string to infoString_
        void appendInfo(string newData);

        //- Clear infoString_
        void clearInfo();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const infoPass&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

