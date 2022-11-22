// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SettlingModel.H"

namespace Foam
{
    defineTypeNameAndDebug(SettlingModel, 0);
    defineRunTimeSelectionTable(SettlingModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SettlingModel::SettlingModel
(
    const dictionary& sedimentDict
)
:
    sedimentDict_(sedimentDict)
{
        Info << "settlingModel object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SettlingModel::~SettlingModel()
{}

// * * * * * * * * * * * * * * * *  Selector   * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::SettlingModel> Foam::SettlingModel::New
(
    const dictionary& sedimentDict
)
{
    word settlingModelType
    (
        sedimentDict.get<word>("settlingModel")
    );

    Info << "Selecting settlingModel"
        << settlingModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(settlingModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "SettlingModel::New : " << endl
                << "    unknown settlingModelType type "
                << settlingModelType
                << ", constructor not in hash table" << endl << endl
                << "    Valid settlingModel types are : " << endl;
        Info << dictionaryConstructorTablePtr_->sortedToc()
             << abort(FatalError);
    }

    return cstrIter()(sedimentDict);
}
