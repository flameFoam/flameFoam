const volScalarField& delta = mesh.objectRegistry::lookupObject<volScalarField>("delta");
scalar c2 = 0.5;

volScalarField ud = c2*pow3(delta)*mag(fvc::curl(fvc::laplacian(U)));

volScalarField lf = 4*mu_u/(S_L*rho0);

volScalarField udByS_L = ud/S_L;
volScalarField deltaBylf = delta/lf;

if (Correlation[i] == "PitschDuchamp" )
{
    scalar b1 = 2;
    scalar b3 = 1;
    scalar Scd = 0.5;
    scalar b3b1= pow(b3,2)/pow(b1,2);
    volScalarField udDa = S_L*deltaBylf;
    Xi = 1 + udByS_L*b3*sqrt(udDa/(ud*Scd+b3b1*udDa));
}
else if (Correlation[i] == "Charlette")
{
    scalar Ck = 1.5;
    scalar Ck_mult1 = 0.245454545454545; // 27/110
    scalar Ck_mult2 = 0.648567745274438; // 4*sqrt(27/110)*18/55

    scalar n43 = 4.0/3.0;
    scalar pi43 = Foam::pow(Foam::constant::mathematical::pi, n43);

    volScalarField Red = 4*deltaBylf*udByS_L+SMALL;

    scalar beta = 0.5;

    volScalarField fu = Ck_mult2*Foam::pow(Ck, 1.5)*pow(udByS_L, 2);
    volScalarField fd = pow(Ck_mult1*Ck*pi43*max(0.0, pow(deltaBylf, n43) - 1), 0.5);
    volScalarField fRe = pow(0.163636363636364*exp(-1.5*Ck*pi43/Red), 0.5)*pow(Red, 0.5);

    volScalarField d = 0.6 + 0.2*exp(-0.1*udByS_L)-0.2*exp(-0.01*deltaBylf);

    volScalarField gamma =
    pow(
        pow(
            pow(
                pow(fu+SMALL, -d)
                +
                pow(fd+SMALL, -d),
                -1/d),
            -1.4)
        +
        pow(fRe+SMALL, -1.4),
        -0.714285714285714
    );

    Xi = pow(1 + min(deltaBylf, gamma*udByS_L), beta);
}
else
{
    FatalErrorInFunction
    << "Unknown turbulent flame speed model." <<  endl
    << exit(FatalError);
}

S_T = S_L*Xi;

cSource = rho0*S_T*mag(fvc::grad(c));
