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
#include "laminarBurningVelocity.H"

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

        //- Full combustionProperties dictionary
        const dictionary& combustionProperties_;

        autoPtr<laminarBurningVelocity> laminarCorrelation_;

        //- Wrinkling factor
        volScalarField sTurbulent_;

        const Switch debug_;


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
                const dictionary& dict,
                const reactionRate& reactRate
            ),
            (dict, reactRate)
        );


    // Constructors

        //- Construct from components
        wrinklingFactor
        (
            const reactionRate& reactRate
        );

        //- Disallow default bitwise copy construction
        wrinklingFactor(const wrinklingFactor&) = delete;


    // Selector

        static autoPtr<wrinklingFactor> New
        (
            const dictionary& reactionRateProperties,
            const reactionRate& reactRate
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
