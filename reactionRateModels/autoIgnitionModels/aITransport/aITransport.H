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
    Foam::wrinklingFactorModels::aITransport

Description
    transport equation for autoignition delay

SourceFiles
    aITransport.C

\*---------------------------------------------------------------------------*/

#ifndef aITransport_H
#define aITransport_H

#include "reactionRate.H"
#include "autoIgnition.H"
#include "UniformTable2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace autoIgnitionModels
{

/*---------------------------------------------------------------------------*\
                         Class aITransport Declaration
\*---------------------------------------------------------------------------*/

class aITransport
:
    public autoIgnition
{
    // Private Data
        const volScalarField& rho_;
        const surfaceScalarField& phi_;

        //- Pressure
        const volScalarField& p_;

        volScalarField ADT_;

        //- Debug flag
        const Switch debug_;

protected:

    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;
public:

    //- Runtime type information
    TypeName("aITransport");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        aITransport
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        aITransport(const aITransport&) = delete;


    // Destructor

        virtual ~aITransport();


    // Member Functions

        //- Correct TBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const aITransport&) = delete;
        //- Return infoString_
        virtual char const *getInfo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End autoIgnitionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
