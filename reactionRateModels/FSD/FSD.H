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
    Foam::reactionRateModels::FSD

Description
    FSD model for combustion progress variable equation closure

SourceFiles
    FSD.C

\*---------------------------------------------------------------------------*/

#ifndef FSD_H
#define FSD_H

#include "reactionRate.H"
#include "wrinklingFactor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{

/*---------------------------------------------------------------------------*\
                         Class FSD Declaration
\*---------------------------------------------------------------------------*/

class FSD
:
    public reactionRate
{
    // Private Data

        //- Auto pointer to reaction rate model
        autoPtr<wrinklingFactor> wrinklingCorrelation_;

        //- Proportionality constant for time scale in the FSD Eq.
        // scalar C_;

        //- Proportionality constant for sub-grid strain
        // scalar alpha_;


public:

    //- Runtime type information
    TypeName("FSD");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        FSD
        (
            const word modelType,
            const dictionary& dictCoeffs,
            const fvMesh& mesh,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        FSD(const FSD&) = delete;


    // Destructor

        virtual ~FSD();


    // Member Functions

        //- Correct combustion source
        virtual void correct();

    // Access functions



    // IO

        //- Update properties from given dictionary
        // virtual bool read(const dictionary& dictProperties);

        //- Return infoString_
        virtual char const *getInfo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FSD&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reactionRateFlameAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
