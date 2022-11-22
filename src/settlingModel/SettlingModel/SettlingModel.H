// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef SettlingModel_H
#define SettlingModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
//#include "HindranceModel.H"

namespace Foam
{

class SettlingModel
{
protected:

    const dictionary& sedimentDict_;
    
public:

    //- Runtime type information
    TypeName("SettlingModel");


    // Declare runtime construction

    declareRunTimeSelectionTable
    (
        autoPtr,
        SettlingModel,
        dictionary,
        (
            const dictionary& sedimentDict
        ),
        (sedimentDict)
    );

    //- Constructor
    SettlingModel
    (
        const dictionary& sedimentDict
    );


    //- Destructor
    virtual ~SettlingModel();

    // Selectors
    
        static autoPtr<SettlingModel> New
        (
            const dictionary& sedimentDict
        );
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
