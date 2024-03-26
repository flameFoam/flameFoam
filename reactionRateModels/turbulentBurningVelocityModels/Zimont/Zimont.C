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

#include "Zimont.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{
    defineTypeNameAndDebug(Zimont, 0);
    addToRunTimeSelectionTable
    (
        turbulentBurningVelocity,
        Zimont,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Zimont::Zimont
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    turbulentBurningVelocity(modelType, reactRate, dict),
    ZimontA_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("ZimontA")),
    alpha_u_(dimensionedScalar(dimViscosity, dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("alpha_u"))),
    Le_(dict.optionalSubDict(modelType + "Coeffs").lookupOrDefault<scalar>("Le", 1.0)),
    ACalpha_(ZimontA_*Foam::pow(0.37, 0.25)*Foam::pow(alpha_u_, -0.25)*Foam::pow(Le_, -0.3)),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            combModel_.coeffs(),
            this->mesh_,
            combModel_
        )
    )
{
    appendInfo("\tTBV estimation method: Zimont correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Zimont::~Zimont()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbulentBurningVelocityModels::Zimont::correct()
{
    if (debug_)
    {
        Info << "\t\tZimont correct:" << endl;
        Info << "\t\t\tInitial average S_T: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();
    sTurbulent_ = max(
        laminarCorrelation_->burningVelocity(),
        ACalpha_*pow(2.0/3*combModel_.turbulence().k(), 0.75)*pow(saneEpsilon(), -0.25)*pow(laminarCorrelation_->burningVelocity(), 0.5));

    if (debug_)
    {
        Info << "\t\t\tObtained average S_T: "  << average(sTurbulent_).value() << endl;
        Info << "\t\t\tZimont correct finished" << endl;
    }
}

char const *Foam::turbulentBurningVelocityModels::Zimont::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}
