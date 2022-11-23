#include "Deigaard.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Deigaard, 0);
    addToRunTimeSelectionTable(FallModel, Deigaard, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Deigaard::Deigaard(const dictionary& dict)
:
    FallModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Deigaard::~Deigaard()
{}

