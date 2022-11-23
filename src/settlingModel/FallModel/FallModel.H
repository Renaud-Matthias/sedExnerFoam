#ifndef FallModel_H
#define FallModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class FallModel
{
protected:
    const dictionary& dict_;
    
    public:
    
    //- Runtime type information
    TypeName("fallModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FallModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
    
    // Constructor
    FallModel(const dictionary& dict);

    // Selectors
    static autoPtr<FallModel> New
    (
        const dictionary& dict
    );

    //- Destructor
    virtual ~FallModel();
       
};
} // End namespace Foam

#endif
