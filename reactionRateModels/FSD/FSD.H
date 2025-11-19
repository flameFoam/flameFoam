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

Class
    Foam::reactionRateModels::FSD

Description
    FSD model for combustion progress variable equation closure

    References:
    \verbatim
        Gubba, S. R., Ibrahim, S., Malalasekera, W., & Masri, A. R. (2007).
        LES modelling of propagating turbulent premixed flames using
        a dynamic flame surface density model.
    \endverbatim

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

public:

    //- Runtime type information
    TypeName("FSD");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        FSD
        (
            const dictionary& dict,
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
