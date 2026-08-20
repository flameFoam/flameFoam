/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

Class
    Foam::autoIgnition

Description
    Abstract auto-ignition delay-time model.

SourceFiles
    autoIgnition.C
    autoIgnitionNew.C

\*---------------------------------------------------------------------------*/

#ifndef autoIgnition_H
#define autoIgnition_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "reactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class autoIgnition Declaration
\*---------------------------------------------------------------------------*/

class autoIgnition
{
protected:

    // Protected data

        const reactionRate& reactionRate_;

        const dictionary& combustionProperties_;

        const fvMesh& mesh_;

        const combustionModel& combModel_;

        const Switch debug_;

        volScalarField tau_;


public:

    TypeName("autoIgnition");

    declareRunTimeSelectionTable
    (
        autoPtr,
        autoIgnition,
        dictionary,
        (
            const dictionary& dict,
            const reactionRate& reactRate
        ),
        (dict, reactRate)
    );


    // Constructors

        autoIgnition(const reactionRate& reactRate);

        autoIgnition(const autoIgnition&) = delete;


    // Selector

        static autoPtr<autoIgnition> New
        (
            const dictionary& reactRateProperties,
            const reactionRate& reactRate
        );


    //- Destructor
    virtual ~autoIgnition();


    // Member Functions

        const volScalarField& tau() const
        {
            return tau_;
        }

        virtual void correct() = 0;

        void operator=(const autoIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
