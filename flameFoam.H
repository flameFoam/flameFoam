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
    Foam::combustionModels::flameFoam

Description
    Collection of combustion models for turbulent premixed combustion of
    hydrogen/air mixtures for OpenFOAM-13.

    Models provided for laminar burning velocity estimation:
        Constant, user defined value
        Malet correlation for lean mixtures [1]
        Power law correlation with user defined coefficients
        Deep neural network model based on experimental data [2]

    Models provided for turbulent burning velocity estimation:
        Turbulent Flame speed Closure (TFC) [3] and Extented TFC [4] models
        for RAS with Bradley [5], Zimont [6] and Bray [7] correlations
        Flame Surface Density (FSD) [8] model for LES with Charlette [9]
        correlation

    References:

SourceFiles
flameFoam.C

\*---------------------------------------------------------------------------*/

#ifndef flameFoam_H
#define flameFoam_H

#include "combustionModel.H"
#include "reactionRate.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                           Class flameFoam Declaration
\*---------------------------------------------------------------------------*/

class flameFoam
:
    public combustionModel
{
    // Private Data

        //- Auto pointer to reaction rate model
        autoPtr<reactionRate> reactionRate_;


        //- progress variable specie index
        label cIndex_;

        OFstream runInfo_;

        bool debug_;


public:

    //- Runtime type information
    TypeName("flameFoam");


    // Constructors

        //- Construct from components
        flameFoam
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        flameFoam(const flameFoam&) = delete;


    //- Destructor
    virtual ~flameFoam();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Specie consumption rate field
        virtual tmp<volScalarField::Internal> R(const label speciei) const;

        //- Specie consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();

        void outputSubInfo();

   // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const flameFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
