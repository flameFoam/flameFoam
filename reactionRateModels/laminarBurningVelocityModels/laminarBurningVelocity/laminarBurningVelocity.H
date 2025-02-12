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
    Foam::laminarBurningVelocity

Description
    Abstract class for laminar burning velocity correlations

SourceFiles
    laminarBurningVelocity.C
    laminarBurningVelocityNew.C

\*---------------------------------------------------------------------------*/

#ifndef laminarBurningVelocity_H
#define laminarBurningVelocity_H

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
                    Class laminarBurningVelocity Declaration
\*---------------------------------------------------------------------------*/

class laminarBurningVelocity
:
    public infoPass
{

protected:

    // Protected data

        // Reaction rate
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

        //- laminar burning velocity
        volScalarField sLaminar_;


public:

    //- Runtime type information
    TypeName("laminarBurningVelocity");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            laminarBurningVelocity,
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
        laminarBurningVelocity
        (
            const word& modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        laminarBurningVelocity(const laminarBurningVelocity&) = delete;


    // Selector

        static autoPtr<laminarBurningVelocity> New
        (
            const reactionRate& reactRate,
            const dictionary& dict
        );


    // Destructor

        virtual ~laminarBurningVelocity();


    // Member Functions

        //- Access functions

            //- Return laminar burning velocity
            const volScalarField& burningVelocity() const
            {
                return sLaminar_;
            }


        //- Correct LBV
        virtual void correct() = 0;

    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const laminarBurningVelocity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
