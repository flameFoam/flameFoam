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
    Foam::turbulentBurningVelocity

Description
    Abstract class for turbulent burning velocity correlations

SourceFiles
    turbulentBurningVelocity.C
    turbulentBurningVelocityNew.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentBurningVelocity_H
#define turbulentBurningVelocity_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "combustionModel.H"
#include "fvmSup.H"
#include "infoPass.H"
#include "reactionRate.H"
#include "laminarBurningVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class turbulentBurningVelocity Declaration
\*---------------------------------------------------------------------------*/

class turbulentBurningVelocity
:
    public infoPass
{

protected:

    // Protected data

        const reactionRate& reactionRate_;

        //- Dictionary
        dictionary coeffDict_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combustion model owner
        const combustionModel& combModel_;

        //- Debug switches
        const Switch debug_;
        const Switch debugFields_;

        //- Turbulent burning velocity
        volScalarField sTurbulent_;

        //- Laminar corrrelation
        autoPtr<laminarBurningVelocity> laminarCorrelation_;


public:

    //- Runtime type information
    TypeName("turbulentBurningVelocity");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            turbulentBurningVelocity,
            dictionary,
            (
                const word modelType,
                const reactionRate& reactRate,
                const dictionary& dict
            ),
            (modelType, reactRate, dict)
        );


    // Constructors

        //- Construct from components
        turbulentBurningVelocity
        (
            const word& modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        turbulentBurningVelocity(const turbulentBurningVelocity&) = delete;


    // Selector

        static autoPtr<turbulentBurningVelocity> New
        (
            const reactionRate& reactRate,
            const dictionary& dict
        );


    // Destructor

        virtual ~turbulentBurningVelocity();


    // Member Functions

        // Sanitized epsilon value for use in divisions
        Foam::tmp<Foam::volScalarField> saneEpsilon();

        //- Access functions

            //- Return turbulent burning velocity
            const volScalarField& burningVelocity() const
            {
                return sTurbulent_;
            }

            //- Return laminar burning velocity
            const volScalarField& getLaminarBurningVelocity();

        //- Correct progress variable source
        // virtual void correct() = 0;

        virtual void correct() = 0;

        // NEAIŠKU AR REIKIA, GAL PRAVERS
        //- Update from dictionary
        // virtual bool read(const dictionary& dictProperties);


    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const turbulentBurningVelocity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
