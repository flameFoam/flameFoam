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
    Foam::laminarBurningVelocityModels::ANN

Description
    ANN correlation of laminar burning velocity

SourceFiles
    ANN.C

\*---------------------------------------------------------------------------*/

#ifndef ANN_H
#define ANN_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class ANN Declaration
\*---------------------------------------------------------------------------*/

class ANN
:
    public laminarBurningVelocity
{

    // Private data

        //- Fraction of steam
        const scalar X_H2_0_;
        const scalar X_H2O_;
        const scalar ER_;

        const volScalarField& p_;
        PtrList<volScalarField> par_;
        const scalarListList W0_;
        const scalarList B0_;
        const scalarListList W1_;
        const scalarList B1_;
        const scalarListList W2_;
        const scalarList B2_;
        const scalarListList W3_;
        const scalarList B3_;
        const scalarListList W4_;
        const dimensionedScalar B4_;

        // Input layer
        PtrList<volScalarField> L0_;
        List<string> L0_names_;
        PtrList<volScalarField> Y0out_;
        List<string> Y0out_names_;

        // Layer 1
        PtrList<volScalarField> L1_;
        List<string> L1_names_;
        PtrList<volScalarField> Y1out_;
        List<string> Y1out_names_;

        // Layer 2
        PtrList<volScalarField> L2_;
        List<string> L2_names_;
        PtrList<volScalarField> Y2out_;
        List<string> Y2out_names_;

        // Layer 3
        PtrList<volScalarField> L3_;
        List<string> L3_names_;
        PtrList<volScalarField> Y3out_;
        List<string> Y3out_names_;

        // Layer 4
        volScalarField L4_;

public:

    //- Runtime type information
    TypeName("ANN");


    // Constructors

        //- Construct from dictionary and reaction rate
        ANN
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dictCoeffs
        );

        //- Disallow default bitwise copy construction
        ANN(const ANN&) = delete;


    // Destructor

        virtual ~ANN();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ANN&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
