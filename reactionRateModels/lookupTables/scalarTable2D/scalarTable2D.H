/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

Class
    Foam::scalarTable2D

Description
    Sorted 2-D lookup table with bilinear interpolation.
    Loaded from constant/ADT directory (one pressure per .ADT file).

SourceFiles
    scalarTable2D.C

\*---------------------------------------------------------------------------*/

#ifndef scalarTable2D_H
#define scalarTable2D_H

#include "scalarField.H"
#include "fileName.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarTable2D Declaration
\*---------------------------------------------------------------------------*/

class scalarTable2D
{
    // Private data

        //- Sorted x samples (pressure)
        scalarField x_;

        //- y samples (temperature) for each x
        List<scalarField> y_;

        //- Table values for each x
        List<scalarField> values_;

        //- Missing-data sentinel
        scalar missingValue_;


    // Private member functions

        //- Last index i such that samples[i] <= s; 0 if s is below range
        static label findLower(const scalarField& samples, const scalar s);

        //- 1-D linear interpolation with end clamping
        static scalar interpolate1D
        (
            const scalarField& samples,
            const scalarField& vals,
            const scalar s
        );


public:

    // Constructors

        scalarTable2D();

        explicit scalarTable2D
        (
            const fileName& directory,
            const scalar missingValue = 1e9
        );


    // Member Functions

        //- Load *.ADT files from directory (portable, no dirent)
        void loadAdtDirectory(const fileName& directory);

        //- Bilinear interpolation in (x, y) = (p, T)
        scalar interpolate(const scalar x, const scalar y) const;

        label nX() const
        {
            return x_.size();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
