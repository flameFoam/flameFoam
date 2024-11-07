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
    Foam::wrinklingFactor

Description
    Abstract class for wrinklingFactor correlations

SourceFiles
    wrinklingFactor.C
    wrinklingFactorNew.C

\*---------------------------------------------------------------------------*/

#ifndef wrinklingFactor_H
#define wrinklingFactor_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "combustionModel.H"
#include "fvmSup.H"
#include "infoPass.H"
#include "reactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wrinklingFactor Declaration
\*---------------------------------------------------------------------------*/

class wrinklingFactor
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

        //- Wrinkling factor
        volScalarField sTurbulent_;

        bool debug_;


public:

    //- Runtime type information
    TypeName("wrinklingFactor");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wrinklingFactor,
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
        wrinklingFactor
        (
            const word& modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wrinklingFactor(const wrinklingFactor&) = delete;


    // Selector

        static autoPtr<wrinklingFactor> New
        (
            const reactionRate& reactRate,
            const dictionary& dict
        );


    // Destructor

        virtual ~wrinklingFactor();


    // Member Functions

        //- Access functions

            //- Return wrinkling factor
            const volScalarField& burningVelocity() const
            {
                return sTurbulent_;
            }


        //- Correct progress variable source
        // virtual void correct() = 0;

        virtual void correct() = 0;

    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const wrinklingFactor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
