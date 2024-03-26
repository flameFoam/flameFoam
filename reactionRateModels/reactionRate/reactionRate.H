/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::reactionRate

Description
    Abstract class for reaction rate model

SourceFiles
    reactionRate.C
    reactionRateNew.C

\*---------------------------------------------------------------------------*/

#ifndef reactionRate_H
#define reactionRate_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "combustionModel.H"
#include "fvmSup.H"
#include "infoPass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reactionRate Declaration
\*---------------------------------------------------------------------------*/

class reactionRate
:
    public infoPass
{

protected:

    // Protected data

        //- Dictionary
        const dictionary& coeffDict_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combustion model owner
        const combustionModel& combModel_;

        //- Mixture
        const basicSpecieMixture& mixture_;

        //- Fuel consumption rate per unit of flame area
        volScalarField cSource_;

        const volScalarField& p_;

        const dimensionedScalar molarH2_;

        const scalar X_H2_0_;
        const scalar Y_H2_0_;
        const scalar Y_H2_99_;

        //- User supplied combustion heat
        const scalar H0_;

        //- Effective combustion heat
        const dimensionedScalar HEff_;

        const label yIndex_;

        const dimensionedScalar WU_;
        const dimensionedScalar p0_;
        const dimensionedScalar rho0_;

        bool debug_;


public:

    //- Runtime type information
    TypeName("reactionRate");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reactionRate,
            dictionary,
            (
                const word modelType,
                const dictionary& dict,
                const fvMesh& mesh,
                const combustionModel& combModel
            ),
            (modelType, dict, mesh, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRate
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        reactionRate(const reactionRate&) = delete;


    // Selector

        static autoPtr<reactionRate> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~reactionRate();


    // Member Functions

        //- Correct progress variable source
        virtual void correct() = 0;


        //- Access functions
        Foam::tmp<Foam::volScalarField> TU() const;
        Foam::tmp<Foam::volScalarField> rhoU() const;
        Foam::tmp<Foam::volScalarField> muU() const;

        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        inline const combustionModel& combModel() const
        {
            return combModel_;
        }

        //- Return progress variable source field
        virtual tmp<volScalarField::Internal> R(const label speciei) const;

        //- Return progress variable source matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Return heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        // NEAIŠKU AR REIKIA, GAL PRAVERS
        //- Update from dictionary
        // virtual bool read(const dictionary& dictProperties);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactionRate&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
