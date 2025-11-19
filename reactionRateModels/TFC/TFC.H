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
    Foam::reactionRateModels::TFC

Description
    TFC model for combustion progress variable equation closure

    References:
    \verbatim
        Zimont, V. L. (2000).
        Gas premixed combustion at high turbulence.
        Turbulent flame closure combustion model.
        Experimental thermal and fluid science, 21(1-3), 179-186.
    \endverbatim

SourceFiles
    TFC.C

\*---------------------------------------------------------------------------*/

#ifndef TFC_H
#define TFC_H

#include "reactionRate.H"
#include "turbulentBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{

/*---------------------------------------------------------------------------*\
                         Class TFC Declaration
\*---------------------------------------------------------------------------*/

class TFC
:
    public reactionRate
{
    // Private Data

        //- Auto pointer to reaction rate model
        autoPtr<turbulentBurningVelocity> turbulentCorrelation_;


public:

    //- Runtime type information
    TypeName("TFC");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        TFC
        (
            const dictionary& dict,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        TFC(const TFC&) = delete;


    // Destructor

        virtual ~TFC();


    // Member Functions

        //- Correct combustion source
        virtual void correct();

    // IO

        //- Return infoString_
        virtual char const *getInfo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TFC&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reactionRateFlameAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
