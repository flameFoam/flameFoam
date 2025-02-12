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

        //- Combustion model owner
        const combustionModel& combModel_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Fuel consumption rate per unit of flame area (field)
        volScalarField cSource_;

        //- Other fields
        const volScalarField& p_;
        const volScalarField& T_;

        //- Unburnt component index
        const label yIndex_;

        //- Unburnt temperature correction switch
        const Switch Tu_;

        //- Debug switches
        const Switch debug_;
        const Switch debugFields_;

        //- Initial values (fields)
        const dimensionedScalar p0_;
        const dimensionedScalar rho0_;

        //- Molar mass of hydrogen
        const dimensionedScalar WH2_;

        //- Molas mass of unburned mixture
        const dimensionedScalar WU_;

        //- Hydrogen concentration selections
        const scalar X_H2_0_;
        const scalar Y_H2_99_;
        const scalar Y_H2_0_;

        //- User supplied combustion heat
        const scalar H0_;

        //- Effective combustion heat
        const dimensionedScalar HEff_;




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
                const combustionModel& combModel
            ),
            (modelType, dict, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRate
        (
            const word& modelType,
            const dictionary& dict,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        reactionRate(const reactionRate&) = delete;


    // Selector

        static autoPtr<reactionRate> New
        (
            const dictionary& dict,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~reactionRate();


    // Member Functions

        //- Correct progress variable source
        virtual void correct() = 0;

        // Sanitized epsilon value for use in divisions
        Foam::tmp<Foam::volScalarField> saneEpsilon() const;

        //- Unburned mixture property functions
        Foam::tmp<Foam::volScalarField> TU() const;
        Foam::tmp<Foam::volScalarField> rhoU() const;
        Foam::tmp<Foam::volScalarField> muU() const;
        Foam::tmp<Foam::volScalarField> nuU() const;

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
