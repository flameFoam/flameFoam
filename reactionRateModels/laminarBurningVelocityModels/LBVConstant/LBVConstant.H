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
    Foam::laminarBurningVelocityModels::LBVConstant

Description
    Constant value for laminar burning velocity

SourceFiles
    LBVConstant.C

\*---------------------------------------------------------------------------*/

#ifndef LBVConstant_H
#define LBVConstant_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class LBVConstant Declaration
\*---------------------------------------------------------------------------*/

class LBVConstant
:
    public laminarBurningVelocity
{

    // Private data

        scalar sLaminar0_;

public:

    //- Runtime type information
    TypeName("LBVConstant");


    // Constructors

        //- Construct from dictionary and reaction rate
        LBVConstant
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dictCoeffs
        );

        //- Disallow default bitwise copy construction
        LBVConstant(const LBVConstant&) = delete;


    // Destructor

        virtual ~LBVConstant();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LBVConstant&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
