// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "settlingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(settlingModel, 0);
    defineRunTimeSelectionTable(settlingModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModel::settlingModel
(
    const dictionary& sedimentDict
)
:
    sedimentDict_(sedimentDict)
{
        Info << "settlingModel object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModel::~settlingModel()
{}

// * * * * * * * * * * * * * * * *  Selector   * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::settlingModel> Foam::settlingModel::New
(
    const dictionary& UfallDict
)
{
    word settlingModelType
    (
        UfallDict.get<word>("settlingModel")
    );

    info << "Selecting settlingModel"
        << settlingModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(settlingModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "settlingModel::New : " << endl
                << "    unknown settlingModelType type "
                << settlingModelType
                << ", constructor not in hash table" << endl << endl
                << "    Valid settlingModel types are : " << endl;
        Info << dictionaryConstructorTablePtr_->sortedToc()
             << abort(FatalError);
    }

    return cstrIter()(ppDict);
}
    */
