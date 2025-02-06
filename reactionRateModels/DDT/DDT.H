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
    Foam::reactionRateModels::DDT

Description
    DDT model for combustion progress variable equation closure

    Original implementation:
    Karanam, A., Verma, V. (2024). Numerical Modeling of Hydrogen Flame
    Acceleration and DDT Using Generalized Turbulent Flame Closure Approach.
    Combustion Science and Technology, 1–32.
    https://doi.org/10.1080/00102202.2024.2408643

    References:
    Ettner, F., Vollmer, K.G., Sattelmayer T.



SourceFiles
    DDT.C

\*---------------------------------------------------------------------------*/

#ifndef DDT_H
#define DDT_H

#include "reactionRate.H"
#include "wrinklingFactor.H"
#include "autoIgnition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{

/*---------------------------------------------------------------------------*\
                         Class DDT Declaration
\*---------------------------------------------------------------------------*/

class DDT
:
    public reactionRate
{
    // Private Data

        //- Auto pointer to reaction rate model
        autoPtr<wrinklingFactor> wrinklingCorrelation_;

        //- Auto pointer to autoignition delay model
        autoPtr<autoIgnition> autoIgnition_;

        const volScalarField& c_;

        const volScalarField& rho_;

        const dimensionedScalar tIgn_;

public:

    //- Runtime type information
    TypeName("DDT");


    // Constructors

        //- Construct from dictionary, mesh and combustion model
        DDT
        (
            const word modelType,
            const dictionary& dictCoeffs,
            const combustionModel& combModel
        );

        //- Disallow default bitwise copy construction
        DDT(const DDT&) = delete;


    // Destructor

        virtual ~DDT();


    // Member Functions

        //- Correct combustion source
        virtual void correct();

    // Access functions



    // IO

        //- Update properties from given dictionary
        // virtual bool read(const dictionary& dictProperties);

        //- Return infoString_
        virtual char const *getInfo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DDT&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reactionRateFlameAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
