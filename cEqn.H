if (CorrelationL[i] == "none")
{
    S_L = S_L0;
}
else if (CorrelationL[i] == "Malet")
{
    S_L = S_L0*pow(1-X_H2O[i],4)*pow(T/Tref,2.2)*pow(p/pref,-0.5);
}

if (max(turbulence.k()).value() > 0)
{
    volScalarField sanitized_epsilon = max(dimensionedScalar(dimVelocity*dimAcceleration, 1E-10), turbulence.epsilon());
    if (Correlation[i] == "Zimont")
    {
        S_T = max(S_L, AC[i]*pow(2.0/3*turbulence.k(), 0.75)*pow(sanitized_epsilon, -0.25)*pow(S_L, 0.5)*alpha_u);
    }
    else if (Correlation[i] == "Bradley")
    {
        S_T = max(S_L, 0.88*pow(0.157*Foam::pow(2.0/3, 0.75)*pow(sanitized_epsilon*mu_u/rho0, 0.5)/pow(S_L, 2)*Le[i], -0.3)*Foam::pow(2.0/3, 0.5)*pow(turbulence.k(), 0.5));
    }
    else if (Correlation[i] == "Bray")
    {
        S_T = max(S_L, 0.875*pow( 0.157*2/3/pow(S_L, 2)*pow(pow(pow(3/2,-1), 0.5)/sanitized_epsilon/mu_u*rho0 ,-0.5) ,-0.392)*pow(2.0/3*turbulence.k(), 0.5));
    }
}
else
{
    S_T = S_L;
}

cSource = rho0*S_T*mag(fvc::grad(c));

if (ignite[i])
{
    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            fluidRegions[i],
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    run_info << "Ignition activated:\n" << endl;

    const scalar ignX = readScalar(combustionProperties.lookup("ignX"));
    const scalar ignY = readScalar(combustionProperties.lookup("ignY"));
    const scalar ignZ = readScalar(combustionProperties.lookup("ignZ"));
    const scalar ignR = readScalar(combustionProperties.lookup("ignitionRadius"));
    run_info << "\tIgnition radius:\t" << ignR << endl;

    forAll (mesh.C(), celli)
    {
        const scalar X2 = (mesh.C()[celli].x()-ignX)*(mesh.C()[celli].x()-ignX);
        const scalar Y2 = (mesh.C()[celli].y()-ignY)*(mesh.C()[celli].y()-ignY);
        const scalar Z2 = (mesh.C()[celli].z()-ignZ)*(mesh.C()[celli].z()-ignZ);
        const scalar r2 = X2+Y2+Z2;


        if (r2 <= ignR*ignR)
        {
            Pout << "Ignition in cell number " << celli << endl;
            c[celli] = 1;

        }
    }

    Pout << "Correcting boundary, ";
    cSource.correctBoundaryConditions();
    Pout << "S_c done, ";
    c.correctBoundaryConditions();
    Pout << "c done." << endl;

    ignite[i] = false;
}

fvScalarMatrix cEqn
(
    fvm::ddt(rho, c)
  + fvm::div(phi, c)
  ==
     fvm::laplacian(turbulence.muEff()/ScT[i], c)
  +  cSource
);

cEqn.relax();
cEqn.solve();
