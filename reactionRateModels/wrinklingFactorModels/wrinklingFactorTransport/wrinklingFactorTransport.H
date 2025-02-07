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
    Foam::wrinklingFactorModels::wrinklingFactorTransport

Description
    wrinklingFactorTransport correlation of flame wrinkling factor

SourceFiles
    wrinklingFactorTransport.C

\*---------------------------------------------------------------------------*/

#ifndef wrinklingFactorTransport_H
#define wrinklingFactorTransport_H

#include "wrinklingFactor.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{

/*---------------------------------------------------------------------------*\
                         Class wrinklingFactorTransport Declaration
\*---------------------------------------------------------------------------*/

class wrinklingFactorTransport
:
    public wrinklingFactor
{
    // Private Data

        autoPtr<laminarBurningVelocity> laminarCorrelation_;

        //subgrid scale flame wrinkling factor
        volScalarField Xi_;

        //- Minimum laminar flame speed
        dimensionedScalar SuMin_;
        
        //- Maximum laminar flame speed
        dimensionedScalar SuMax_;
        
        //- Xi model coefficient
        scalar XiCoef_;
        
        //- Xi shape coefficient
        scalar XiShapeCoef_;
        
        //- Turbulence coefficient
        scalar uPrimeCoef_;
        
        //- Strain rate coefficient
        dimensionedScalar sigmaExt_;

        //- Laminar burning velocity field
        volScalarField Su;
        
        //- Progress variable field (usually from composition)
        const volScalarField* b_;


        //- Debug flag
        const Switch debug_;

public:

    //- Runtime type information
    TypeName("wrinklingFactorTransport");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        wrinklingFactorTransport
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wrinklingFactorTransport(const wrinklingFactorTransport&) = delete;


    // Destructor

        virtual ~wrinklingFactorTransport();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wrinklingFactorTransport&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End wrinklingFactorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
