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
    Foam::turbulentBurningVelocityModels::Zimont

Description
    Zimont correlation of turbulent burning velocity

    References:
    \verbatim
        Zimont, V. L. (2000).
        Gas premixed combustion at high turbulence.
        Turbulent flame closure combustion model.
        Experimental thermal and fluid science, 21(1-3), 179-186.
    \endverbatim

SourceFiles
    Zimont.C

\*---------------------------------------------------------------------------*/

#ifndef Zimont_H
#define Zimont_H

#include "turbulentBurningVelocity.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Zimont Declaration
\*---------------------------------------------------------------------------*/

class Zimont
:
    public turbulentBurningVelocity
{
    // Private Data

        //- Molecular Lewis number of fresh mixture
        const dimensionedScalar ZimontA_;

        const dimensionedScalar Le_;

        const dimensionedScalar ACalpha_;

public:

    //- Runtime type information
    TypeName("Zimont");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        Zimont
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        Zimont(const Zimont&) = delete;


    // Destructor

        virtual ~Zimont();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Zimont&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End turbulentBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
