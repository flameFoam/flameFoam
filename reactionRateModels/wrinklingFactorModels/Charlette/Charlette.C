/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2025 Lithuanian Energy Institute

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
    const dictionary& dict,
    const reactionRate& reactRate
):
    wrinklingFactor(reactRate),
    c2_(2.0),
    Ck_(1.5),
    Ck_mult1_(27.0/110.0),
    Ck_mult2_(4*pow(Ck_mult1_, 0.5)*18.0/55.0),
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

    const fvMesh& mesh_(reactionRate_.mesh());

    const volVectorField& U_ = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& delta_ = mesh_.objectRegistry::lookupObject<volScalarField>("delta");

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
