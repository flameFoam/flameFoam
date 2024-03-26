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
    Foam::laminarBurningVelocityModels::Malet

Description
    Malet correlation of laminar burning velocity

SourceFiles
    Malet.C

\*---------------------------------------------------------------------------*/

#ifndef Malet_H
#define Malet_H

#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Malet Declaration
\*---------------------------------------------------------------------------*/

class Malet
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

public:

    //- Runtime type information
    TypeName("Malet");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        Malet
        (
            const word modelType,
            const dictionary& dictCoeffs,
            const fvMesh& mesh,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        Malet(const Malet&) = delete;


    // Destructor

        virtual ~Malet();


    // Member Functions

        //- Correct LBV
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Malet&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarBurningVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
