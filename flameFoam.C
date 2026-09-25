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

\*---------------------------------------------------------------------------*/

#include "flameFoam.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(flameFoam, 0);
    addToRunTimeSelectionTable(combustionModel, flameFoam, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::flameFoam::flameFoam
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, combustionProperties),
    reactionRate_(
        reactionRate::New
        (
            this->coeffs(),
            *this
        )
    ),
    cIndex_(thermo.specieIndex(thermo.Y("c"))),
    runInfo_("flameFoam." + this->mesh().name() + ".combustionInfo"),
    debug_(this->coeffs().lookupOrDefault("debug", false)),
    debugFields_(this->coeffs().lookupOrDefault("debugFields", false))
{
    #include "../version.H"

    logCritical("flameFoam combustion model selected");
    logCritical(string("flameFoam library version: ") + flameFoamVersion);
    logCritical
    (
        string("Mesh size: ")
      + Foam::name(returnReduce(this->mesh().cells().size(), sumOp<label>()))
    );

    logCritical("Average initial values:");
    logCritical(string("\tp_rgh: ") + Foam::name(average(db().lookupObject<volScalarField>("p_rgh")).value()));
    logCritical(string("\tT: ") + Foam::name(average(thermo.T()).value()));
    logCritical(string("\trho: ") + Foam::name(average(thermo.rho()).value()));
    logCritical(string("\tmu: ") + Foam::name(average(thermo.mu()).value()));
    logCritical(string("\tc: ") + Foam::name(average(thermo.Y("c")).value()));
    logCritical
    (
        string("debug: ") + Switch(debug_).asText()
      + ", debugFields: " + Switch(debugFields_).asText()
    );

    outputSubInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::flameFoam::~flameFoam()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::flameFoam::correct()
{
    if (debug_)
    {
        Info << "flameFoam correct: " << endl;
    }
    reactionRate_->correct();
    if (debug_)
    {
        Info << "\tflameFoam correct finished" << endl;
    }
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::flameFoam::R(const label speciei) const
{
    if (debug_)
    {
        Info << "flameFoam R(" << speciei << "): " << endl;
    }
    if (speciei == cIndex_)
    {
        return reactionRate_->R(cIndex_);
    }
    else
    {
        return
        volScalarField::Internal::New
        (
            typedName("R_" + this->thermo().Y()[speciei].name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::flameFoam::R(volScalarField& Y) const
{
    if (debug_)
    {
        Info << "flameFoam R(" << Y.name() << "): " << endl;
    }
    if (this->thermo().specieIndex(Y) == cIndex_)
    {
        return reactionRate_->R(Y);
    }
    else
    {
        return tmp<fvScalarMatrix>(new fvScalarMatrix(Y, dimMass/dimTime));
    }
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::flameFoam::Qdot() const
{
    if (debug_)
    {
        Info << "flameFoam Qdot: " << endl;
    }
    return reactionRate_->Qdot();
}

void Foam::combustionModels::flameFoam::logCritical(const string& msg)
{
    Info<< msg << endl;
    runInfo_ << msg.c_str() << endl;
}


void Foam::combustionModels::flameFoam::outputSubInfo()
{
    const char* info = reactionRate_().getInfo();
    if (info && info[0] != '\0')
    {
        Info<< info << endl;
        runInfo_ << info << endl;
        reactionRate_().clearInfo();
    }
}

bool Foam::combustionModels::flameFoam::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
