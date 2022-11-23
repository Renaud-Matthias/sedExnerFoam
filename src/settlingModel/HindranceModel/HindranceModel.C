#include "HindranceModel.H"


namespace Foam
{
    defineTypeNameAndDebug(HindranceModel, 0);
    defineRunTimeSelectionTable(HindranceModel, dictionary);
}

// Constructor

Foam::HindranceModel::HindranceModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// Destructor

Foam::HindranceModel::~HindranceModel()
{}

// Selector

Foam::autoPtr<Foam::HindranceModel> Foam::HindranceModel::New
(
    const dictionary& dict
)
{
    word hindranceModelType(dict.get<word>("hindranceModel"));

    Info<< "Selecting hindranceModel "
        << hindranceModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(hindranceModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "HindranceModel::New(const dictionary&) : " << endl
            << "    unknown hindranceModelType type "
            << hindranceModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid hindranceModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<HindranceModel>(cstrIter()(dict));
}
