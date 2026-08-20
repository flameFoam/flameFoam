/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

Class
    Foam::autoIgnitionModels::aITransport

Description
    Transport equation for auto-ignition progress, with bilinear ADT lookup.

SourceFiles
    aITransport.C

\*---------------------------------------------------------------------------*/

#ifndef aITransport_H
#define aITransport_H

#include "autoIgnition.H"
#include "scalarTable2D.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace autoIgnitionModels
{

/*---------------------------------------------------------------------------*\
                       Class aITransport Declaration
\*---------------------------------------------------------------------------*/

class aITransport
:
    public autoIgnition
{
    // Private Data

        const volScalarField& rho_;
        const surfaceScalarField& phi_;
        const volScalarField& p_;

        volScalarField ADT_;

        scalarTable2D adtTable_;

        dimensionedScalar Sct_;


public:

    TypeName("aITransport");


    // Constructors

        aITransport
        (
            const dictionary& dict,
            const reactionRate& reactRate
        );

        aITransport(const aITransport&) = delete;


    //- Destructor
    virtual ~aITransport();


    // Member Functions

        virtual void correct();

        void operator=(const aITransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace autoIgnitionModels
} // End namespace Foam

#endif

// ************************************************************************* //
