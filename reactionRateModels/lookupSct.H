#ifndef lookupSct_H
#define lookupSct_H

#include "IOdictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Read Sct from the already-registered thermophysicalTransport dictionary.
// Combustion is constructed before the transport model, so this must be
// called from correct(), not from the ETFC/aITransport constructors.
inline dimensionedScalar lookupSctFromRegistry(const fvMesh& mesh)
{
    const objectRegistry* obr = &mesh;

    if (!obr->foundObject<IOdictionary>("thermophysicalTransport"))
    {
        obr = &mesh.time();
    }

    if (!obr->foundObject<IOdictionary>("thermophysicalTransport"))
    {
        FatalErrorInFunction
            << "thermophysicalTransport is not in the object registry. "
            << "Sct is read from the existing transport object "
            << "(no extra IOdictionary). Call this after transport construction."
            << exit(FatalError);
    }

    const IOdictionary& ttm =
        obr->lookupObject<IOdictionary>("thermophysicalTransport");

    const dictionary& coeffs =
        ttm.found("RAS") ? ttm.subDict("RAS")
      : ttm.found("LES") ? ttm.subDict("LES")
      : ttm;

    return dimensionedScalar("Sct", dimless, coeffs.lookup<scalar>("Sct"));
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
