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
        S_L = S_L*pow(1-X_H2O[i],4)*pow(T/Tref,2.2)*pow(p/pref,-0.5);
    }
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
    fvm::laplacian(turbulence.muEff()/ScT[i], c)
  + cSource
);

cEqn.relax();
cEqn.solve();
