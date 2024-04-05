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

Class
    Foam::turbulentBurningVelocityModels::Bradley

Description
    Bradley correlation of turbulent burning velocity

SourceFiles
    Bradley.C

\*---------------------------------------------------------------------------*/

#ifndef Bradley_H
#define Bradley_H

#include "turbulentBurningVelocity.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Bradley Declaration
\*---------------------------------------------------------------------------*/

class Bradley
:
    public turbulentBurningVelocity
{
    // Private Data

        //- Molecular Lewis number of fresh mixture
        scalar Le_;

public:

    //- Runtime type information
    TypeName("Bradley");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        Bradley
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        Bradley(const Bradley&) = delete;


    // Destructor

        virtual ~Bradley();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Bradley&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End turbulentBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
