/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "ETFC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateModels
{
    defineTypeNameAndDebug(ETFC, 0);
    addToRunTimeSelectionTable
    (
        reactionRate,
        ETFC,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateModels::ETFC::ETFC
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    reactionRate(modelType, dict, mesh, combModel),
    turbulentCorrelation_(
        turbulentBurningVelocity::New
        (
            *this,
            dict
        )
    ),
    c_(combModel_.thermo().Y("c")),
    Dt_inf_
    (
    	IOobject
        (
            "Dt_inf",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Dt_inf", dimKinematicViscosity, Zero)
    ),
    TauByT_
    (
    	IOobject
        (
            "TauByT",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TauByT", dimless, Zero)
    ),
    expFactor_
    (
    	IOobject
        (
            "expFactor",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("expFactor", dimless, Zero)
    ),
    cLam_
    (
    	IOobject
        (
            "cLam",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cLam", dimDensity/dimTime, Zero)
    ),
    Sct_("Sct", dimless, 0),
    alpha_u_("alpha_u", dimKinematicViscosity, this->coeffDict_),
    Le_("Le", dimless, this->coeffDict_)
{
    IOdictionary thermophysicalTransportDict
    (
        IOobject
        (
            "thermophysicalTransport",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    Sct_ = thermophysicalTransportDict.lookup<scalar>("Sct");

    appendInfo("Reaction rate model: ETFC");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateModels::ETFC::~ETFC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionRateModels::ETFC::correct
(
)
{
    if (debug_)
    {
        Info << "\tETFC correct:" << endl;
        Info << "\t\tInitial min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
    }

    turbulentCorrelation_->correct();

    Dt_inf_ = combModel_.turbulence().nut()/Sct_; // TODO: check against dev2-efix-ZimontLe-no0-fix2/

    TauByT_ = max(1.5* Dt_inf_/(combModel_.turbulence().k()*mesh_.time()), SMALL);  // TODO: check against dev2-efix-ZimontLe-no0-fix2/

    expFactor_ = 1 - exp(-1/TauByT_);

    cLam_ = 0.25*pow(turbulentCorrelation_->getLaminarBurningVelocity(), 2)*
        rhoU()*max(c_-SMALL*mesh_.time().deltaT().value(),0.0)*(1-c_)/            // TODO: check if deltaT or absolute value should be used
        (alpha_u_/Le_+Dt_inf_*expFactor_);

    cSource_ =
        rhoU()*
        turbulentCorrelation_->burningVelocity()*
        mag(fvc::grad(c_))*
        pow((max(1-expFactor_*TauByT_, 0.0)),0.5)  + cLam_;

    if (debug_)
    {
        Info << "\t\tObtained min/avg/max cSource: " << min(cSource_).value() << " " << average(cSource_).value() << " " << max(cSource_).value() << endl;
        Info << "\t\tETFC correct finished" << endl;
    }
}

char const *Foam::reactionRateModels::ETFC::getInfo()
{
    infoString_.append(turbulentCorrelation_().getInfo());
    turbulentCorrelation_().clearInfo();
    return infoString_.c_str();
}
