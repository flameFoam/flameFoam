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

        //- Dictionary
        dictionary coeffDict_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combustion model owner
        const combustionModel& combModel_;

        //- laminar burning velocity
        volScalarField sLaminar_;

        bool debug_;


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
                const dictionary& dict,
                const fvMesh& mesh,
                const combustionModel& combModel
            ),
            (modelType, dict, mesh, combModel)
        );


    // Constructors

        //- Construct from components
        laminarBurningVelocity
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        laminarBurningVelocity(const laminarBurningVelocity&) = delete;


    // Selector

        static autoPtr<laminarBurningVelocity> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
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

        // virtual void correct(
        //     const volScalarField& k,
        //     const volScalarField& epsilon,
        //     const volScalarField& rho,
        //     const volScalarField& mu
        // ) = 0;

        // NEAIŠKU AR REIKIA, GAL PRAVERS
        //- Update from dictionary
        // virtual bool read(const dictionary& dictProperties);


    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const laminarBurningVelocity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
