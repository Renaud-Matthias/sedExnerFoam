// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef settlingModel_H
#define settlingModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class settlingModel
{
protected:

    const dictionary& sedimentDict_;
    
public:

    //- Runtime type information
    TypeName("settlingModel");


    // Declare runtime construction

    declareRunTimeSelectionTable
        (
            autoPtr,
	    settlingModel,
            dictionary,
            (const dictionary& sedimentDict),
            (sedimentDict)
        );

    //- Constructor
    settlingModel
    (
        const dictionary& sedimentDict
    );


    //- Destructor
    virtual ~settlingModel();

    // Selectors
    
        static autoPtr<settlingModel> New
        (
            const dictionary& UfallDict
        );
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
