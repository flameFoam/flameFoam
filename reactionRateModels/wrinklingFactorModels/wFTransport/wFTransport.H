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
    Foam::wrinklingFactorModels::wFTransport

Description
    transport equation for flame wrinkling factor

SourceFiles
    wFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef wFTransport_H
#define wFTransport_H

#include "reactionRate.H"
#include "wrinklingFactor.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{

/*---------------------------------------------------------------------------*\
                         Class wFTransport Declaration
\*---------------------------------------------------------------------------*/

class wFTransport
:
    public wrinklingFactor
{
    // Private Data

        autoPtr<laminarBurningVelocity> laminarCorrelation_;

        //subgrid scale flame wrinkling factor
        volScalarField Xi_;

        //- Le
        dimensionedScalar Le_;

        //- ReT
        dimensionedScalar ReT_;

        const volScalarField& rho_;
        const surfaceScalarField& phi_;

        //- Pressure
        const volScalarField& p_;

        //- Pressure
        const dimensionedScalar p0_;

        //- Debug flag
        const Switch debug_;


protected:

    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;
public:

    //- Runtime type information
    TypeName("wFTransport");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        wFTransport
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wFTransport(const wFTransport&) = delete;


    // Destructor

        virtual ~wFTransport();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wFTransport&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End wrinklingFactorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
