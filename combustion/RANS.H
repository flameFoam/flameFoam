volScalarField sanitized_epsilon = max(turbulence.epsilon(), dimensionedScalar(dimVelocity*dimAcceleration, SMALL));
if (Correlation[i] == "Zimont")
{
    S_T = max(S_L, ACalpha*pow(2.0/3*turbulence.k(), 0.75)*pow(sanitized_epsilon, -0.25)*pow(S_L, 0.5));
}
else if (Correlation[i] == "Bradley")
{
    S_T = max(
        S_L,
        1.37179015019233* //0.88*0.157^(-0.3)*(2/3)^(0.275)
        pow(S_L, 0.6)*pow(sanitized_epsilon*mu_u/rho0, -0.15)*Foam::pow(Le[i], -0.3)*pow(turbulence.k(), 0.5)
    );
}
else if (Correlation[i] == "Bray")
{
    S_T = max(S_L, 0.875*pow( 0.157*2/3/pow(S_L, 2)*pow(pow(pow(3/2,-1), 0.5)/sanitized_epsilon/mu_u*rho0 ,-0.5) ,-0.392)*pow(2.0/3*turbulence.k(), 0.5));
}

cSource = rho0*S_T*mag(fvc::grad(c));

if (model[i] == "ETFC")
{
    volScalarField Dt_inf = 0.09*pow(turbulence.k(),2)/(sanitized_epsilon*ScT[i]);
    volScalarField TauByT = max(0.135*turbulence.k()/(ScT[i]*sanitized_epsilon*mesh.time()), SMALL);
    volScalarField expFactor = 1 - exp(-1/TauByT);

    cSource = cSource*pow((1-(expFactor*TauByT)),0.5) + 0.25*pow(S_L, 2)*rho0*c*(1-c)/(alpha_u+Dt_inf*expFactor);
}
