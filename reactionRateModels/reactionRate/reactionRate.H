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

        //- Combustion model owner
        const combustionModel& combModel_;

        //- Dictionary of global mixture parameters
        const dictionary& combustionProperties_;

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
        const dimensionedScalar X_H2_0_;
        const dimensionedScalar Y_H2_99_;
        const dimensionedScalar Y_H2_0_;

        //- User supplied combustion heat
        const dimensionedScalar H0_;

        //- Effective combustion heat
        const dimensionedScalar HEff_;

        //- Unburnt mixture properties
        volScalarField rhoU_;
        volScalarField TU_;
        volScalarField muU_;

        //- Unburnt thermal diffusivity and correction switch
        volScalarField alphaU_;
        const Switch calculateAlphaU_;



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
                const dictionary& dict,
                const combustionModel& combModel
            ),
            (dict, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRate
        (
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        reactionRate(const reactionRate&) = delete;


    // Selector

        static autoPtr<reactionRate> New
        (
            const dictionary& combustionProperties,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~reactionRate();


    // Member Functions

        //- Correct progress variable source
        virtual void correct() = 0;

        void correctUnburntProperties();

        //- Access functions
        inline const volScalarField& TU() const
        {
            return TU_;
        }

        inline const volScalarField& rhoU() const
        {
            return rhoU_;
        }

        inline const volScalarField& muU() const
        {
            return muU_;
        }

        inline const volScalarField& alphaU() const
        {
            return alphaU_;
        }

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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactionRate&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
