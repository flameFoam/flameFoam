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
    Foam::wrinklingFactorModels::Charlette

Description
    Charlette correlation of flame wrinkling factor

    References:
    \verbatim
        Charlette, F., Meneveau, C., & Veynante, D. (2002).
        A power-law flame wrinkling model for LES of premixed
        turbulent combustion Part I: non-dynamic formulation and
        initial tests. Combustion and Flame, 131(1-2), 159-180.
    \endverbatim

SourceFiles
    Charlette.C

\*---------------------------------------------------------------------------*/

#ifndef Charlette_H
#define Charlette_H

#include "wrinklingFactor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{

/*---------------------------------------------------------------------------*\
                         Class Charlette Declaration
\*---------------------------------------------------------------------------*/

class Charlette
:
    public wrinklingFactor
{
    // Private Data

        const scalar c2_;

        const scalar Ck_;

        const scalar Ck_mult1_;

        const scalar Ck_mult2_;

        const scalar n43_;

        const scalar pi43_;

        const scalar beta_;

public:

    //- Runtime type information
    TypeName("Charlette");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        Charlette
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        Charlette(const Charlette&) = delete;


    // Destructor

        virtual ~Charlette();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Charlette&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End wrinklingFactorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
