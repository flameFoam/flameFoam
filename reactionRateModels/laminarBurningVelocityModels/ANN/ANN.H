/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
        const dimensionedScalar sLaminar0_;
        const dimensionedScalar pRef_;
        const dimensionedScalar TRef_;
        const volScalarField& p_;
        const volScalarField& T_;
        
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

        //- Construct from dictionary, mesh and combustion model
        ANN
        (
            const word modelType,
            const dictionary& dictCoeffs,
            const fvMesh& mesh,
            const combustionModel& combModel
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
