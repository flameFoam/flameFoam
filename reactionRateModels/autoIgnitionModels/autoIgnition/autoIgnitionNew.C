/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

\*---------------------------------------------------------------------------*/

#include "autoIgnition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::autoIgnition> Foam::autoIgnition::New
(
    const dictionary& reactRateProperties,
    const reactionRate& reactRate
)
{
    const dictionary& autoIgnitionDict =
        reactRateProperties.subDict("autoIgnition");

    const word modelType(autoIgnitionDict.lookup("model"));

    Info<< "Selecting auto-ignition model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(reactRateProperties)
            << "Unknown auto-ignition model "
            << modelType << nl << nl
            << "Valid auto-ignition models are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<autoIgnition>
    (
        cstrIter()(autoIgnitionDict.optionalSubDict(modelType), reactRate)
    );
}


// ************************************************************************* //
