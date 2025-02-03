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

#include "volFields.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "combustionModel.H"

namespace Foam
{
namespace reactionRateModels
{

class DDT;  // Forward declaration

/*---------------------------------------------------------------------------*\
                           Class autoIgnition Declaration
\*---------------------------------------------------------------------------*/

class autoIgnition
{
    // Private data
    scalar Tref_;
    scalar pref_;
    scalar Ea_;
    scalar A_;

public:
    TypeName("autoIgnition");

    // Constructors

        //- Construct from dictionary
        autoIgnition
        (
            const DDT& ddt,
            const dictionary& dict
        );

    //- Destructor
    virtual ~autoIgnition();

    // Member Functions

        //- Return the auto-ignition delay time
        scalar tau() const;
        virtual void correct() {}

    static autoPtr<autoIgnition> New
    (
        const DDT& ddt,
        const dictionary& dict
    )
    {
        return autoPtr<autoIgnition>(new autoIgnition(ddt, dict));
    }
};

} // End namespace reactionRateModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 