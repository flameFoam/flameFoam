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
    Foam::autoIgnition

Description
    Auto-ignition delay time model based on temperature and pressure.

SourceFiles
    autoIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef autoIgnition_H
#define autoIgnition_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "combustionModel.H"
#include "fvmSup.H"
#include "infoPass.H"
#include "reactionRate.H"

namespace Foam
{

namespace reactionRateModels
{
/*---------------------------------------------------------------------------*\
                           Class autoIgnition Declaration
\*---------------------------------------------------------------------------*/

class autoIgnition
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

    bool debug_;

    //- Auto-ignition delay time field
    volScalarField tau_;

public:
    //- Runtime type information
    TypeName("autoIgnition");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        autoIgnition,
        dictionary,
        (
            const word modelType,
            const reactionRate& reactRate,
            const dictionary& dict
        ),
        (modelType, reactRate, dict)
    );

    // Constructors
    autoIgnition
    (
        const word& modelType,
        const reactionRate& reactRate,
        const dictionary& dict
    );

    //- Disallow default bitwise copy construction
    autoIgnition(const autoIgnition&) = delete;

    // Selector
    static autoPtr<autoIgnition> New
    (
        const reactionRate& reactRate,
        const dictionary& dict
    );

    // Destructor
    virtual ~autoIgnition();

    // Member Functions
    const volScalarField& tau() const
    {
        return tau_;
    }

    virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const autoIgnition&) = delete;
};
} // End namespace reactionRateModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 