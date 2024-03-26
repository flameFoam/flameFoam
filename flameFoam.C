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
            this->mesh(),
            *this
        )
    ),
    mixture_(dynamic_cast<const basicSpecieMixture&>(this->thermo().composition())),
    cIndex_(mixture_.index(mixture_.Y("c"))),
    runInfo_("flameFoam." + this->mesh().name() + ".combustionInfo"),
    debug_(this->coeffs().lookupOrDefault("debug", false))
{
    // Log flameFoam combustion model information
    runInfo_ << "flameFoam combustion model selected" << endl;

    #include "../version.H"
    runInfo_ << "flameFoam library version: " << flameFoamVersion << endl;

    runInfo_ << "\nMesh size: " << returnReduce(this->mesh().cells().size(), sumOp<label>()) << endl;

    runInfo_ << "\nAverage initial values:" << endl;
    runInfo_ << "\tp_rgh: " << average(db().lookupObject<volScalarField>("p_rgh")).value() << endl;
    runInfo_ << "\tT: "     << average(thermo.T()).value() << endl;
    runInfo_ << "\trho: "   << average(thermo.rho()).value() << endl;
    runInfo_ << "\tmu: "    << average(thermo.mu()).value() << endl;
    runInfo_ << "\tc: "     << average(mixture_.Y("c")).value() << endl;

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
            typedName("R_" + mixture_.Y()[speciei].name()),
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
    if (mixture_.index(Y) == cIndex_)
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

void Foam::combustionModels::flameFoam::outputSubInfo()
{
    if (strcmp(reactionRate_().getInfo(), "") != 0)
    {
        runInfo_ << reactionRate_().getInfo() << endl;
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
