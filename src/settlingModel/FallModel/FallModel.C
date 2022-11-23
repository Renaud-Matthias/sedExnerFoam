#include "FallModel.H"


namespace Foam
{
    defineTypeNameAndDebug(FallModel, 0);
    defineRunTimeSelectionTable(FallModel, dictionary);
}

// Constructor

Foam::FallModel::FallModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// Destructor

Foam::FallModel::~FallModel()
{}

// Selector

Foam::autoPtr<Foam::FallModel> Foam::FallModel::New
(
    const dictionary& dict
)
{
    word fallModelType(dict.get<word>("fallModel"));

    Info<< "Selecting fallModel "
        << fallModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(fallModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "FallModel::New(const dictionary&) : " << endl
            << "    unknown fallModelType type "
            << fallModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid fallModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<FallModel>(cstrIter()(dict));
}
