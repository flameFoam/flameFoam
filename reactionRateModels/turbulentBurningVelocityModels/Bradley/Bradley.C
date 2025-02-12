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

#include "Bradley.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{
    defineTypeNameAndDebug(Bradley, 0);
    addToRunTimeSelectionTable
    (
        turbulentBurningVelocity,
        Bradley,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bradley::Bradley
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    turbulentBurningVelocity(modelType, reactRate, dict),
    Le_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("Le"))
{
    appendInfo("\tTBV estimation method: Bradley correlation");
    appendInfo("\t\tLe: " + name(Le_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bradley::~Bradley()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbulentBurningVelocityModels::Bradley::correct()
{
    if (debug_)
    {
        Info << "\t\tBradley correct:" << endl;
        Info << "\t\t\tInitial average S_T: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();
    sTurbulent_ =
        1.37179015019233  //0.88*0.157^(-0.3)*(2/3)^(0.275)
      * pow(laminarCorrelation_->burningVelocity(), 0.6)
      * pow(reactionRate_.saneEpsilon() * reactionRate_.nuU(), -0.15)
      * pow(Le_, -0.3)
      * pow(combModel_.turbulence().k(), 0.5);

    if (debug_)
    {
        Info << "\t\t\tObtained average S_T: "  << average(sTurbulent_).value() << endl;
        Info << "\t\t\tBradley correct finished" << endl;
    }

}

char const *Foam::turbulentBurningVelocityModels::Bradley::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
