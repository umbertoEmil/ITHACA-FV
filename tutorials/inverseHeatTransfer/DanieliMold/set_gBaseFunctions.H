volScalarField& T = _T();
fvMesh& mesh = _mesh();
Info << "Radial Basis Functions are used." << endl;
// The center of each function is the projection of each thermocouple
// on the boundary hotSide
if (thermocouplesCellID.size() == 0)
{
    readThermocouples();
}

gBaseFunctions.resize(thermocouplesNum);
gWeights.resize(thermocouplesNum);
//for (label funcI = 0; funcI < thermocouplesNum; funcI++)
forAll(gBaseFunctions, funcI)
{
    gBaseFunctions[funcI].resize(T.boundaryField()[hotSide_ind].size());
    scalar thermocoupleX = thermocouplesPos[funcI].x();
    scalar thermocoupleZ = thermocouplesPos[funcI].z();
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
        scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
        scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - thermocoupleX) +
                                   (faceZ - thermocoupleZ) * (faceZ - thermocoupleZ));

        if (faceX > - castWidth / 2 && faceX < castWidth / 2)
        {
            if (faceZ > meniscusHeight)
            {
                //Fast decay of the basis
                gBaseFunctions[funcI][faceI] = 0;
                //gBaseFunctions[funcI][faceI] = Foam::exp(- (10 *
                //                               shapeParameter
                //                               * radius * radius));
            }
            else
            {
                gBaseFunctions[funcI][faceI] = Foam::exp(- (shapeParameter *
                                               shapeParameter
                                               * radius * radius));
            }
        }
        else
        {
            gBaseFunctions[funcI][faceI] = 0;
        }
    }
}
