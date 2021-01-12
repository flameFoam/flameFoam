const dimensionedScalar molarH2 = dimensionedScalar("molarH2",dimensionSet(1,0,0,0,-1,0,0), 0.002);
const dimensionedScalar Tref("Tref", dimTemperature, 298);
const dimensionedScalar pref("pref", dimPressure, 100000);

Info<< "Reading combustion properties\n" << endl;
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

ScT[i] = readScalar(combustionProperties.lookup("ScT"));
X_H2_0[i] = readScalar(combustionProperties.lookup("X_H2_0"));
run_info << "Initial hydrogen concentration: " << X_H2_0[i] << endl;

Y_H2_0[i] = dimensionedScalar("Y_H2_0", molarH2*X_H2_0[i]*thermoFluid[i].p().average()/(rhoFluid[i].average()*constant::physicoChemical::R*thermoFluid[i].T().average()));
Y_H2_99[i] = 0;
run_info << "Initial hydrogen mass fraction: " << Y_H2_0[i] << endl;

X_H2O[i] = readScalar(combustionProperties.lookup("X_H2O"));
run_info << "Initial steam concentration: " << X_H2O[i] << endl;

if ((X_H2_0_c>1) || (X_H2O_c>1))
{
    WarningInFunction
    << "Hydrogen and/or steam concentrations not present in thermophysicalProperties." << nl << endl;
    run_info << "Hydrogen and/or steam concentrations not present in thermophysicalProperties. " << endl;

}
else
{
    if ((X_H2_0_c!=X_H2_0[i]) || (X_H2O_c!=X_H2O[i]))
    {
        FatalErrorInFunction
        << "Hyrogen and/or steam concentrations are not equal in combustionProperties and thermophysicalProperties." <<  endl
        << exit(FatalError);
    }
}

dimensionedScalar S_L0(combustionProperties.lookup("S_L0"));
S_L0_value[i] = S_L0.value();
CorrelationL[i] = "none";
if (S_L0_value[i] == 0)
{
    CorrelationL[i] = "Malet";
    const scalar ER = 0.705*X_H2_0[i]/(0.295*(1-X_H2_0[i]-X_H2O[i]));
    S_L0 = dimensionedScalar("S_L0", dimVelocity,1.44*ER*ER+1.07*ER-0.29);
    Info << "Malet's correlation used for the referential laminar velocity. Calculated equivalence ratio: " << ER << ", calculated S_L0: " << S_L0 << endl;
    run_info << "Malet's correlation used, calculated equivalence ratio: " << ER << endl;
    S_L0_value[i] = S_L0.value();
}
Tref_value[i] = Tref.value();
pref_value[i] = pref.value();
run_info << "Reference laminar flame speed: " << S_L0 << endl;

Le[i] = readScalar(combustionProperties.lookup("Le"));
run_info << "Le: " << Le[i] << endl;

Switch ignition(combustionProperties.lookup("ignite"));
ignite[i] = ignition;

S_LFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "S_L",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        S_L0*pow(1-X_H2O[i],4)*pow(thermoFluid[i].T()/Tref,2.2)*pow(thermoFluid[i].p()/pref,-0.5)
    )
);

S_TFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "S_T",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        S_LFluid[i]
    )
);

cSourceFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "cSource",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho0*S_TFluid[i]*mag(fvc::grad(cFluid[i]))
    )
);

Correlation[i] = combustionProperties.lookup("correlation");
Info << "Correlation: " << Correlation[i] << endl;
run_info << "Correlation: " << Correlation[i] << endl;

dimensionedScalar alpha_u = dimensionedScalar("alpha_u", dimensionSet(0,-0.5,0.25,0,0,0,0), 0);
if (Correlation[i] == "Zimont")
{
    const scalar ZimontA = readScalar(combustionProperties.lookup("ZimontA"));
    const dimensionedScalar alpha_read(combustionProperties.lookup("alpha"));
    run_info << "\tZimont coefficient: " << ZimontA << endl;
    run_info << "\tUnburned alpha: " << alpha_read << endl;
    AC[i] = ZimontA*Foam::pow(0.37, 0.25);
    alpha_u = pow(alpha_read, -0.25);
}
alpha_u_value[i] = alpha_u.value();

Switch visc_switch(combustionProperties.lookup("visciousDis"));
visciousDis[i] = visc_switch;
run_info << "Viscious dissipation: " << visciousDis[i] << endl;
