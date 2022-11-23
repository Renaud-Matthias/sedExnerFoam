// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef SettlingModel_H
#define SettlingModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "HindranceModel.H"
#include "FallModel.H"

namespace Foam
{

class SettlingModel
{
protected:

    const dictionary& sedimentDict_;

    autoPtr<FallModel> FallModel_;

    //autoPtr<HindranceModel> HindranceModel_;
    
public:

    //- Constructor
    SettlingModel
    (
        const dictionary& sedimentDict
    );


    //- Destructor
    virtual ~SettlingModel();

};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
