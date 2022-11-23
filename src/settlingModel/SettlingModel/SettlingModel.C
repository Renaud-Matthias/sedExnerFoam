// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SettlingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SettlingModel::SettlingModel
(
    const dictionary& sedimentDict
)
:
    sedimentDict_(sedimentDict)
{
        Info << "Initialization of SettlingModel" << endl;
        Info << "Initialization of FallModel" << endl;
        FallModel_ = FallModel::New (sedimentDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SettlingModel::~SettlingModel()
{}

