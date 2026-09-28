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
    Foam::laminarBurningVelocityModels::hydrogenAirPower

Description
    hydrogen-air power law based correlation for laminar burning velocity

SourceFiles
    hydrogenAirPower.C

\*---------------------------------------------------------------------------*/

#ifndef hydrogenAirPower_H
#define hydrogenAirPower_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class hydrogenAirPower Declaration
\*---------------------------------------------------------------------------*/

class hydrogenAirPower
:
    public laminarBurningVelocity
{

    // Private data

        const dimensionedScalar X_H2_0_;

        //- Polynomial coefficients for sLaminar0
        const dimensionedScalar c4_;
        const dimensionedScalar c3_;
        const dimensionedScalar c2_;
        const dimensionedScalar c1_;
        const dimensionedScalar c0_;

        //- Laminar burning velocity at pRef and TRef
        const dimensionedScalar sLaminar0_;


        //- Reference pressure and temperature
        const dimensionedScalar pRef_;
        const dimensionedScalar TRef_;

public:

    //- Runtime type information
    TypeName("hydrogenAirPower");


    // Constructors

        //- Construct from dictionary and reaction rate
        hydrogenAirPower
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        hydrogenAirPower(const hydrogenAirPower&) = delete;


    // Destructor

        virtual ~hydrogenAirPower();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hydrogenAirPower&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
