/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2024 Lithuanian Energy Institute

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
    Foam::laminarBurningVelocityModels::Malet

Description
    Malet correlation of laminar burning velocity

SourceFiles
    Malet.C

\*---------------------------------------------------------------------------*/

#ifndef Malet_H
#define Malet_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Malet Declaration
\*---------------------------------------------------------------------------*/

class Malet
:
    public laminarBurningVelocity
{

    // Private data

        //- Fraction of steam
        const scalar X_H2_0_;
        const scalar X_H2O_;
        const scalar ER_;
        const dimensionedScalar sLaminar0_;
        const dimensionedScalar pRef_;
        const dimensionedScalar TRef_;
        const volScalarField& p_;

public:

    //- Runtime type information
    TypeName("Malet");


    // Constructors

        //- Construct from dictionary and reaction rate
        Malet
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dictCoeffs
        );

        //- Disallow default bitwise copy construction
        Malet(const Malet&) = delete;


    // Destructor

        virtual ~Malet();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Malet&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
