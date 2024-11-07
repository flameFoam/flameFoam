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
    Foam::infoPass

Description
    Helper class to log flameFoam model data

SourceFiles
    infoPass.C

\*---------------------------------------------------------------------------*/

#ifndef infoPass_H
#define infoPass_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class infoPass Declaration
\*---------------------------------------------------------------------------*/

class infoPass
{
protected:

    // Protected data

        string infoString_;

public:

    //- Runtime type information
    TypeName("infoPass");


    // Constructors

        //- Construct from components
        infoPass();

        //- Disallow default bitwise copy construction
        infoPass(const infoPass&) = delete;


    //- Destructor
    virtual ~infoPass();


    // Member Functions

        //- Return infoString_
        virtual char const *getInfo();

        //- Append new string to infoString_
        void appendInfo(string newData);

        //- Clear infoString_
        void clearInfo();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const infoPass&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

