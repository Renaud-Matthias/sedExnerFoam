#ifndef Deigaard_H
#define Deigaard_H

#include "FallModel.H"

namespace Foam
{

class Deigaard
:
    public FallModel
{

public:

    //- Runtime type information
    TypeName("Deigaard");


    // Constructors

        //- Construct from components
    Deigaard(const dictionary& dict);


    //- Destructor
    virtual ~Deigaard();

};

}

#endif
