{

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
            ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
            : -dpdt
        )
      - fvm::laplacian(turbulence.alphaEff(), he)
        ==
        rho*(U&g)
      + cSource*(Y_H2_0[i]-Y_H2_99[i])*H0*(min(c, scalar(1))-c.oldTime())/(c-c.oldTime()+1E-10) // renormalized to max c = 1
      + fvOptions(rho, he)
    );

    if (visciousDis[i])
    {
        volTensorField gradU = fvc::grad(U);
        volTensorField tau = turbulence.muEff() * (gradU + gradU.T());

        EEqn = EEqn + fvc::div(tau & U);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    c.min(1);
    Info << "Average c: " << c.average().value() << endl;

    thermo.composition().Y(0) = c;
    Ydefault = 1-c;

    fvOptions.correct(he);

    thermo.correct();

}
