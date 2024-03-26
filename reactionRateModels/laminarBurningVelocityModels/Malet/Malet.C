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

#include "Malet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{
    defineTypeNameAndDebug(Malet, 0);
    addToRunTimeSelectionTable
    (
        laminarBurningVelocity,
        Malet,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::Malet::Malet
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
):
    laminarBurningVelocity(modelType, dict, mesh, combModel),
    X_H2_0_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2_0")),
    X_H2O_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2O")),
    ER_(0.705*X_H2_0_/(0.295*(1-X_H2_0_-X_H2O_))),
    sLaminar0_(dimensionedScalar(dimVelocity, 1.44*ER_*ER_+1.07*ER_-0.29)),
    pRef_(dimensionedScalar(dimPressure, 100000)),
    TRef_(dimensionedScalar(dimTemperature, 298)),
    p_(mesh.lookupObject<volScalarField>("p")),
    T_(mesh.lookupObject<volScalarField>("T"))
{
    appendInfo("\tLBV estimation method: Malet correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::Malet::~Malet()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarBurningVelocityModels::Malet::correct
()
{
    if (debug_)
    {
        Info << "\t\t\tMalet correct:" << endl;
        Info << "\t\t\t\tInitial average S_L: "  << average(sLaminar_).value() << endl;
    }

    sLaminar_ = sLaminar0_*pow(1-X_H2O_,4)*pow(T_/TRef_,2.2)*pow(p_/pRef_,-0.5);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tMalet correct finished" << endl;
    }


    // dimensionedScalar omega0
    // (
    //     "omega0",
    //     dimensionSet(1, -2, -1, 0, 0, 0, 0),
    //     correlation_.omega0()
    // );
    //
    // dimensionedScalar sigmaExt
    // (
    //     "sigmaExt",
    //     dimensionSet(0, 0, -1, 0, 0, 0, 0),
    //     correlation_.sigmaExt()
    // );
    //
    // dimensionedScalar omegaMin
    // (
    //     "omegaMin",
    //     omega0.dimensions(),
    //     1e-4
    // );
    //
    // dimensionedScalar kMin
    // (
    //     "kMin",
    //     sqr(dimVelocity),
    //     small
    // );
    //
    // const compressibleMomentumTransportModel& turbulence =
    //     combModel_.turbulence();
    //
    // // Total strain
    // const volScalarField sigmaTotal
    // (
    //     sigma + alpha_*turbulence.epsilon()/(turbulence.k() + kMin)
    // );
    //
    // const volScalarField omegaInf(correlation_.omega0Sigma(sigmaTotal));
    //
    // dimensionedScalar sigma0("sigma0", sigma.dimensions(), 0.0);
    //
    // const volScalarField tau(C_*mag(sigmaTotal));
    //
    // volScalarField Rc
    // (
    //     (tau*omegaInf*(omega0 - omegaInf) + sqr(omegaMin)*sigmaExt)
    //    /(sqr(omega0 - omegaInf) + sqr(omegaMin))
    // );
    //
    // const volScalarField& rho = combModel_.rho();
    // const tmp<surfaceScalarField> tphi = combModel_.phi();
    // const surfaceScalarField& phi = tphi();
    //
    // solve
    // (
    //      fvm::ddt(rho, omega_)
    //    + fvm::div(phi, omega_)
    //   ==
    //      rho*Rc*omega0
    //    - fvm::SuSp(rho*(tau + Rc), omega_)
    // );
    //
    // omega_.min(omega0);
    // omega_.max(0.0);


}


// bool  Foam::laminarBurningVelocityModels::Malet::read
// (
//     const dictionary& dict
// )
// {
//     if (reactionRateFlameArea::read(dict))
//     {
//         coeffDict_ = dict.optionalSubDict(typeName + "Coeffs");
//         coeffDict_.lookup("C") >> C_;
//         coeffDict_.lookup("alpha") >> alpha_;
//         correlation_.read
//         (
//             coeffDict_.subDict(fuel_)
//         );
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }

// ************************************************************************* //
