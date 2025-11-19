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
    Foam::turbulentBurningVelocityModels::Bray

Description
    Bray correlation of turbulent burning velocity

    References:
    \verbatim
        Bray, K. N. C. (1990).
        Studies of the turbulent burning velocity.
        Proceedings of the Royal Society of London. Series A:
        Mathematical and Physical Sciences, 431(1882), 315-335.
    \endverbatim

SourceFiles
    Bray.C

\*---------------------------------------------------------------------------*/

#ifndef Bray_H
#define Bray_H

#include "turbulentBurningVelocity.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Bray Declaration
\*---------------------------------------------------------------------------*/

class Bray
:
    public turbulentBurningVelocity
{
    // Private Data
    const scalar c1_;

public:

    //- Runtime type information
    TypeName("Bray");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        Bray
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        Bray(const Bray&) = delete;


    // Destructor

        virtual ~Bray();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Bray&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End turbulentBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
