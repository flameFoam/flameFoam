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
    Foam::reactionRateModels::ETFC

Description
    ETFC model for combustion progress variable equation closure

    References:
    \verbatim
        Lipatnikov, A. N., & Chomiak, J. (2002).
        Turbulent flame speed and thickness: phenomenology,
        evaluation, and application in multi-dimensional simulations.
        Progress in energy and combustion science, 28(1), 1-74.
    \endverbatim

SourceFiles
    ETFC.C

\*---------------------------------------------------------------------------*/

#ifndef ETFC_H
#define ETFC_H

#include "reactionRate.H"
#include "turbulentBurningVelocity.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{

/*---------------------------------------------------------------------------*\
                         Class ETFC Declaration
\*---------------------------------------------------------------------------*/

class ETFC
:
    public reactionRate
{
    // Private Data

        //- Auto pointers to reaction rate models
        autoPtr<turbulentBurningVelocity> turbulentCorrelation_;

        autoPtr<laminarBurningVelocity> laminarCorrelation_;

	volScalarField Dt_inf_;

	volScalarField DEffByRho_;

	volScalarField TauByT_;

	volScalarField expFactor_;

	volScalarField cLam_;

protected:

    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;

            const dimensionedScalar Le_;

public:

    //- Runtime type information
    TypeName("ETFC");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        ETFC
        (
            const dictionary& dict,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        ETFC(const ETFC&) = delete;


    // Destructor

        virtual ~ETFC();


    // Member Functions

        //- Correct combustion source
        virtual void correct();

    // IO

        //- Return infoString_
        virtual char const *getInfo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ETFC&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reactionRateFlameAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
