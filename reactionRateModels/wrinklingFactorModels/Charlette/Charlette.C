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

#include "Charlette.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcLaplacian.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{
    defineTypeNameAndDebug(Charlette, 0);
    addToRunTimeSelectionTable
    (
        wrinklingFactor,
        Charlette,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::Charlette::Charlette
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    wrinklingFactor(modelType, reactRate, dict),
    U_(mesh_.lookupObject<volVectorField>("U")),
    delta_(mesh_.objectRegistry::lookupObject<volScalarField>("delta")),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            combModel_.coeffs(),
            mesh_,
            combModel_
        )
    ),
    c2_(0.5),
    Ck_(1.5),
    Ck_mult1_(0.245454545454545), // 27/110
    Ck_mult2_(0.648567745274438), // 4*sqrt(27/110)*18/55
    n43_(4.0/3.0),
    pi43_(Foam::pow(Foam::constant::mathematical::pi, n43_)),
    beta_(0.5)
{
    appendInfo("\tWrinkling factor estimation method: Charlette correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::Charlette::~Charlette()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::wrinklingFactorModels::Charlette::correct()
{
    if (debug_)
    {
        Info << "\t\tCharlette correct:" << endl;
        Info << "\t\t\tInitial average TBV: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();

    volScalarField ud_ = c2_*pow3(delta_)*mag(fvc::curl(fvc::laplacian(U_)));

    volScalarField lf_ = 4*reactionRate_.muU()/(laminarCorrelation_->burningVelocity() * reactionRate_.rhoU());

    volScalarField udByLBV_ = ud_/laminarCorrelation_->burningVelocity();

    volScalarField deltaBylf_ = delta_/lf_;

    volScalarField Red_ = 4*deltaBylf_*udByLBV_+SMALL;

    volScalarField fu_ = Ck_mult2_*Foam::pow(Ck_, 1.5)*pow(udByLBV_, 2);
    volScalarField fd_ = pow(Ck_mult1_*Ck_*pi43_*max(0.0, pow(deltaBylf_, n43_) - 1), 0.5);
    volScalarField fRe_ = pow(0.163636363636364*exp(-1.5*Ck_*pi43_/Red_), 0.5)*pow(Red_, 0.5);

    volScalarField d_ = 0.6 + 0.2*exp(-0.1*udByLBV_)-0.2*exp(-0.01*deltaBylf_);

    volScalarField gamma_ =
    pow(
        pow(
            pow(
                pow(fu_+SMALL, -d_)
                +
                pow(fd_+SMALL, -d_),
                -1/d_),
            -1.4)
        +
        pow(fRe_+SMALL, -1.4),
        -0.714285714285714
    );

    volScalarField xi_ = pow(1 + min(deltaBylf_, gamma_*udByLBV_), beta_);

    sTurbulent_ = xi_*laminarCorrelation_->burningVelocity();

    if (debug_)
    {
        Info << "\t\t\tObtained average TBV: "  << average(sTurbulent_).value() << endl;
        Info << "\t\t\tCharlette correct finished" << endl;
    }

}

char const *Foam::wrinklingFactorModels::Charlette::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
