volScalarField& T_u = (Tu[i]) ? TunFluid[i] : thermo.T();
Info << "Min p value: " << min(p).value();
if (Tu[i])
{
    T_u = W_u*p/(rho0*constant::physicoChemical::R*pow(p/p0,1.0/thermo.gamma()));
}

if (ignite[i])
{
    run_info << "Ignition activated:\n" << "\tIgnition radius:\t" << Foam::sqrt(ignR2[i]) << endl;

    forAll (mesh.C(), celli)
    {
        const scalar r2 = pow(mesh.C()[celli].x()-ignX[i], 2)+pow(mesh.C()[celli].y()-ignY[i], 2)+pow(mesh.C()[celli].z()-ignZ[i], 2);
        if (r2 <= ignR2[i])
        {
            Pout << "Ignition in cell number " << celli << endl;
            if (model[i] == "ETFC")
            {
                c[celli] = 0.5;
            }
            else
            {
                c[celli] = 1;
            }
        }
    }

    Pout << "Updating c boundary" << endl;
    c.correctBoundaryConditions();

    ignite[i] = false;
}

if (CorrelationL[i] == "ANN")
{
    #include "ANN.H"
}
else
{
    const dimensionedScalar S_L0("S_L0", dimVelocity, LBV_param);
    S_L = S_L0;
    if (CorrelationL[i] == "Malet")
    {
        S_L = S_L*pow(1-X_H2O[i],4)*pow(T_u/Tref,2.2)*pow(p/pref,-0.5);
    }
}

if (model[i] != "ETFC")
{
    dFactor = turbulence.muEff()/ScT[i];
}

if (turbulenceType[i] == "RAS")
{
    #include "RANS.H"
}
else if (turbulenceType[i] == "LES")
{
    #include "LES.H"
}
else
{
    S_T = S_L;
    cSource = rho0*S_T*mag(fvc::grad(c));
}

fvScalarMatrix cEqn
(
    fvm::ddt(rho, c) + mvConvection->fvmDiv(phi, c)
  ==
    fvm::laplacian(dFactor, c)
  + cSource
);

cEqn.relax();
cEqn.solve();
c = c*max(c - 1E-9, scalar(0)) / (max(c - 1E-9, scalar(0)) + VSMALL);