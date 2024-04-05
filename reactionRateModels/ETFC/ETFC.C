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

    Sct_("Sct", dimless, this->coeffDict_),
    alpha_u_("alpha_u", dimViscosity, this->coeffDict_),
    Le_("Le", dimless, this->coeffDict_)
{
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

    volScalarField Dt_inf = combModel_.turbulence().nut()/Sct_;

    volScalarField TauByT = max(1.5* Dt_inf/(combModel_.turbulence().k()*mesh_.time()), SMALL);
    
    volScalarField expFactor = 1 - exp(-1/TauByT);
    
    volScalarField c = combModel_.thermo().composition().Y("c");
    
    volScalarField dFactor = rhoU()*(alpha_u_/Le_+Dt_inf*expFactor);
    
    volScalarField cLam = 0.25*pow(turbulentCorrelation_->laminarCorrelation().burningVelocity(), 2)*
        rhoU()*max(c-SMALL*mesh_.time().deltaT().value(),0.0)*(1-c)/
        (alpha_u_/Le_+Dt_inf*expFactor);
  
    cSource_ =
        rhoU()*
        turbulentCorrelation_->burningVelocity()*
        mag(fvc::grad(c))*
        pow((max(1-expFactor*TauByT, 0.0)),0.5)  + cLam;

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
