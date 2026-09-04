/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

Class
    Foam::reactionRateModels::DDT

Description
    DDT model: flame wrinkling plus auto-ignition source.

    Original implementation:
    Karanam, A., Verma, V. (2024). Numerical Modeling of Hydrogen Flame
    Acceleration and DDT Using Generalized Turbulent Flame Closure Approach.
    Combustion Science and Technology.

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

        autoPtr<wrinklingFactor> wrinklingCorrelation_;

        autoPtr<autoIgnition> autoIgnition_;

        const volScalarField& c_;

        const volScalarField& rho_;

        const dimensionedScalar tIgn_;


public:

    TypeName("DDT");


    // Constructors

        DDT
        (
            const dictionary& dict,
            const combustionModel& combModel
        );

        DDT(const DDT&) = delete;


    //- Destructor
    virtual ~DDT();


    // Member Functions

        virtual void correct();

        void operator=(const DDT&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRateModels
} // End namespace Foam

#endif

// ************************************************************************* //
