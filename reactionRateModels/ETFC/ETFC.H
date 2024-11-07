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
    Foam::reactionRateModels::ETFC

Description
    ETFC model for combustion progress variable equation closure

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

	const volScalarField& c_;

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

            dimensionedScalar alpha_u_;

            dimensionedScalar Le_;

public:

    //- Runtime type information
    TypeName("ETFC");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        ETFC
        (
            const word modelType,
            const dictionary& dictCoeffs,
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

        //- Update properties from given dictionary
        // virtual bool read(const dictionary& dictProperties);

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
