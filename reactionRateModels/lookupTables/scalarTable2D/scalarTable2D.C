/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

\*---------------------------------------------------------------------------*/

#include "scalarTable2D.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "OSspecific.H"
#include "ListOps.H"
#include "DynamicList.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarTable2D::scalarTable2D()
:
    x_(),
    y_(),
    values_(),
    missingValue_(1e9)
{}


Foam::scalarTable2D::scalarTable2D
(
    const fileName& directory,
    const scalar missingValue
)
:
    x_(),
    y_(),
    values_(),
    missingValue_(missingValue)
{
    loadAdtDirectory(directory);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::scalarTable2D::findLower
(
    const scalarField& samples,
    const scalar s
)
{
    if (samples.empty())
    {
        return 0;
    }

    label lo = 0;
    label hi = samples.size() - 1;

    if (s <= samples[0])
    {
        return 0;
    }
    if (s >= samples[hi])
    {
        return hi;
    }

    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;
        if (samples[mid] <= s)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return lo;
}


Foam::scalar Foam::scalarTable2D::interpolate1D
(
    const scalarField& samples,
    const scalarField& vals,
    const scalar s
)
{
    if (samples.empty())
    {
        return 0;
    }
    if (samples.size() == 1)
    {
        return vals[0];
    }

    const label i = findLower(samples, s);
    if (i == samples.size() - 1)
    {
        return vals[i];
    }

    const scalar x0 = samples[i];
    const scalar x1 = samples[i + 1];
    const scalar w = (s - x0)/max(x1 - x0, VSMALL);
    return (1 - w)*vals[i] + w*vals[i + 1];
}


void Foam::scalarTable2D::loadAdtDirectory(const fileName& directory)
{
    const fileNameList files(readDir(directory, fileType::file));

    DynamicList<scalar> pBuf;
    DynamicList<scalarField> yTmp;
    DynamicList<scalarField> vTmp;

    forAll(files, fi)
    {
        const fileName& f = files[fi];
        if (f.ext() != "ADT")
        {
            continue;
        }

        const word stem(f.lessExt());
        scalar pVal = 0;
        {
            IStringStream iss(stem);
            iss >> pVal;
            if (iss.bad())
            {
                WarningInFunction
                    << "Skipping ADT file with non-numeric name: " << f
                    << endl;
                continue;
            }
        }

        IFstream is(directory/f);
        if (!is.good())
        {
            WarningInFunction
                << "Cannot open " << directory/f << endl;
            continue;
        }

        string line;
        for (label header = 0; header < 4 && is.good(); ++header)
        {
            is.getLine(line);
        }

        DynamicList<scalar> Tbuf;
        DynamicList<scalar> vbuf;

        while (is.good())
        {
            is.getLine(line);
            if (line.empty())
            {
                continue;
            }

            IStringStream iss(line);
            scalar T = 0;
            scalar tau = 0;
            iss >> T >> tau;
            if (iss.bad() || tau < 0)
            {
                continue;
            }

            Tbuf.append(T);
            vbuf.append(tau);
        }

        if (Tbuf.empty())
        {
            continue;
        }

        labelList Torder(Tbuf.size());
        sortedOrder(Tbuf, Torder);

        scalarField Tsorted(Tbuf.size());
        scalarField vSorted(Tbuf.size());
        forAll(Torder, i)
        {
            Tsorted[i] = Tbuf[Torder[i]];
            vSorted[i] = vbuf[Torder[i]];
        }

        pBuf.append(pVal);
        yTmp.append(Tsorted);
        vTmp.append(vSorted);
    }

    if (pBuf.empty())
    {
        FatalErrorInFunction
            << "No valid .ADT files in " << directory
            << exit(FatalError);
    }

    labelList pOrder(pBuf.size());
    sortedOrder(pBuf, pOrder);

    x_.setSize(pBuf.size());
    y_.setSize(pBuf.size());
    values_.setSize(pBuf.size());

    forAll(pOrder, i)
    {
        x_[i] = pBuf[pOrder[i]];
        y_[i] = yTmp[pOrder[i]];
        values_[i] = vTmp[pOrder[i]];
    }

    Info<< "scalarTable2D: loaded " << x_.size()
        << " pressure samples from " << directory << endl;
}


Foam::scalar Foam::scalarTable2D::interpolate
(
    const scalar x,
    const scalar y
) const
{
    if (x_.empty())
    {
        return missingValue_;
    }

    const label i = findLower(x_, x);

    const scalar v0 = interpolate1D(y_[i], values_[i], y);

    if (i == x_.size() - 1 || x <= x_[i] + VSMALL)
    {
        return v0;
    }

    const scalar v1 = interpolate1D(y_[i + 1], values_[i + 1], y);
    const scalar w = (x - x_[i])/max(x_[i + 1] - x_[i], VSMALL);

    return (1 - w)*v0 + w*v1;
}


// ************************************************************************* //
