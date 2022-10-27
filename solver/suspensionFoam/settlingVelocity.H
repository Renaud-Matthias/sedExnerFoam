Info<< "Compute sediment settling velocity Ufall\n" << endl;

/*
law to compute settling velocity

float functionTest( float a, float b )
{
    float result;
    if (a > b)
        result = a;
    else
        result = b;
    return result;
}

Info<< functionTest(7.9,3.2) << endl

*/
Info<< "TEST OF FUNCTION\n" << endl;



dimensionedVector Ufall("Ufall", dimVelocity, vector (0, 0, -0.68));
