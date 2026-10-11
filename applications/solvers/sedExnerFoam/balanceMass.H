// ============================================================
// Sediment mass conservation check
// ============================================================

// phiqbProj is already integrated over each edge:
//   phiqbProj > 0 : flux leaving the domain
//   phiqbProj < 0 : flux entering the domain

scalar totalBedloadIn  = 0.0;
scalar totalBedloadOut = 0.0;

// Only boundary edges
for
(
    label edgei = meshProj.nInternalEdges();
    edgei < meshProj.nEdges();
    ++edgei
)
{
    const scalar flux = phiqbProj[edgei];

    if (flux > 0)
    {
        totalBedloadOut += flux;
    }
    else
    {
        totalBedloadIn -= flux;
    }
}

// Net flux leaving the domain
const scalar netBedloadFlux =
    totalBedloadOut - totalBedloadIn;


// ============================================================
// Volume variation computed by Exner
// ============================================================

// dVsed has dimensions m3/s.
// Integrate over the projected face areas.
const scalar dVsedRate =
    gSum
    (
        dVsed.primitiveField() * SfProj
    );


// ============================================================
// Conservation error
// ============================================================

// Exner convention:
//   positive boundary flux = sediment leaving the domain
//   therefore dVsed/dt = - net boundary flux

const scalar conservationError =
    dVsedRate + netBedloadFlux;


// ============================================================
// Output
// ============================================================

Info<< "----------------------------------------" << endl;
Info<< "Sediment mass conservation:" << endl;

Info<< "  Bedload IN  = "
    << totalBedloadIn
    << " m3/s" << endl;

Info<< "  Bedload OUT = "
    << totalBedloadOut
    << " m3/s" << endl;

Info<< "  Net flux    = "
    << netBedloadFlux
    << " m3/s" << endl;

Info<< "  dVsed/dt    = "
    << dVsedRate
    << " m3/s" << endl;

Info<< "  Error       = "
    << conservationError
    << " m3/s" << endl;

if (CsPtr.valid())
{
    const volScalarField& Cs = CsPtr();

    const scalar Vsed =
        fvc::domainIntegrate(Cs).value();

    const scalar CsMin =
        gMin(Cs.primitiveField());

    const scalar CsMaxValue =
        gMax(Cs.primitiveField());

    Info<< "----------------------------------------" << endl;
    Info<< "Suspended sediment mass:" << endl;
    Info<< "  Integral(Cs dV) = "
        << Vsed << " m3" << endl;
    Info<< "  Min(Cs)         = "
        << CsMin << endl;
    Info<< "  Max(Cs)         = "
        << CsMaxValue << endl;
    Info<< "----------------------------------------" << endl;
}
scalar dVsedTotal = gSum(dVsed.primitiveField());

Info<< "dVsed      = " << dVsedTotal << " m3/s" << endl;

Info<< "----------------------------------------" << endl;
