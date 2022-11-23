#ifndef HindranceModel_H
#define HindranceModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class HindranceModel
{
protected:
    const dictionary& dict_;
    
    public:
    
    //- Runtime type information
    TypeName("hindranceModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HindranceModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
    
    // Constructor
    HindranceModel(const dictionary& dict);

    // Selectors
    static autoPtr<HindranceModel> New
    (
        const dictionary& dict
    );

    //- Destructor
    virtual ~HindranceModel();
       
};
} // End namespace Foam

#endif
