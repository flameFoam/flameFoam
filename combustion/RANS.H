volScalarField sanitized_epsilon = max(turbulence.epsilon(), dimensionedScalar(dimVelocity*dimAcceleration, SMALL));
if (Correlation[i] == "Zimont")
{
    S_T = ACalpha*pow(2.0/3*turbulence.k(), 0.75)*pow(sanitized_epsilon, -0.25)*pow(S_L, 0.5)*Foam::pow(Le[i], -0.3);
}
else if (Correlation[i] == "Bradley")
{
    S_T = 1.37179015019233  //0.88*0.157^(-0.3)*(2/3)^(0.275)
    *pow(S_L, 0.6)
    *pow(sanitized_epsilon*mu_u/rho0, -0.15)
    *Foam::pow(Le[i], -0.3)
    *pow(turbulence.k(), 0.5);
}
else if (Correlation[i] == "Bray")
{
    S_T = 0.875
    *pow( 0.157*2/3/pow(S_L, 2)
    *pow(pow(pow(3/2,-1), 0.5)/sanitized_epsilon/mu_u*rho0 ,-0.5) ,-0.392)
    *pow(2.0/3*turbulence.k(), 0.5);
}

if (model[i] != "ETFC")
{
    S_T = max(S_T,S_L);
}

cSource = rho0*S_T*mag(fvc::grad(c));

if (model[i] == "ETFC")
{
    volScalarField Dt_inf = 0.09*pow(turbulence.k(),2)/(sanitized_epsilon*ScT[i]);
    volScalarField TauByT = max(0.135*turbulence.k()/(ScT[i]*sanitized_epsilon*mesh.time()), SMALL);
    volScalarField expFactor = 1 - exp(-1/TauByT);

    dFactor = rho*(alpha_u/Le[i]+Dt_inf*expFactor);
    cLam = 0.25*pow(S_L, 2)*rho0*max(c-SMALL*1E-3,0.0)*(1-c)/(alpha_u/Le[i]+Dt_inf*expFactor);
    cSource = cSource*pow((max(1-expFactor*TauByT, 0.0)),0.5) + cLam;
}

if (wallQuenching[i])
{
    Pe_Q = 1.9/ER[i] * pow(p/dimensionedScalar(dimPressure, 3.0e5), 0.26*min(1.0, pow(ER[i], -2)));
    Pe = wallDistance * S_L * Le[i] / alpha_u;
    D_Q = max((Pe-Pe_Q)/(Pe_max-Pe_Q), SMALL);
    quenchAlpha = quenchBeta / max(1.0 - 0.6 * D_Q, SMALL);
    I_QR = max(1.0 - 2.0 / (1+exp(min(D_Q * quenchAlpha, scalar(50.0) ))), 0.0);
    cSource = cSource * I_QR;
}

if (quenching[i])
{
    volScalarField sigma = 0.26*Foam::log
    (
        0.37*Foam::pow
        (
            Foam::sqrt
            (
                max(turbulence.k(), k_SMALL)
            )
            /
            Foam::pow
            (
                thermo.nu()*sanitized_epsilon,
                0.25
            ),
            3
        )
    );
    sigma.max(SMALL);

    G = 0.5*Foam::erfc
    (
        -
        (
            Foam::log
            (
                15.0*thermo.nu()*Foam::pow(gcr,2)
                /
                sanitized_epsilon
            )
            +
            0.5*sigma
        )
        /
        Foam::sqrt(2.0*sigma)
    );

    cSource = G*cSource;
}