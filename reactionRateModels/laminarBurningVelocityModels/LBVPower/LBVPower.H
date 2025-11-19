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
    Foam::laminarBurningVelocityModels::LBVPower

Description
    LBWPower correlation of laminar burning velocity

SourceFiles
    LBVPower.C

\*---------------------------------------------------------------------------*/

#ifndef LBVPower_H
#define LBVPower_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class LBVPower Declaration
\*---------------------------------------------------------------------------*/

class LBVPower
:
    public laminarBurningVelocity
{

    // Private data

        //- Initial molar fractions of hydrogen and steam
        const dimensionedScalar X_H2_0_;
        const dimensionedScalar X_H2O_;

        //- Exponents for dilution, temperature and pressure
        const dimensionedScalar powerDil_;
        const dimensionedScalar powerT_;
        const dimensionedScalar powerP_;

        //- Quadratic polynomial coefficients for sLaminar0
        const dimensionedScalar a2_;
        const dimensionedScalar a1_;
        const dimensionedScalar a0_;

        //- Equivalence ratio
        const dimensionedScalar ER_;

        //- Laminar burning velocity at pRef and TRef
        const dimensionedScalar sLaminar0_;

        //- Reference pressure and temperature
        const dimensionedScalar pRef_;
        const dimensionedScalar TRef_;

public:

    //- Runtime type information
    TypeName("LBVPower");


    // Constructors

        //- Construct from dictionary and reaction rate
        LBVPower
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        LBVPower(const LBVPower&) = delete;


    // Destructor

        virtual ~LBVPower();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LBVPower&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
