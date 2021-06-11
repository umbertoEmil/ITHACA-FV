Info << "Setting BC\n" << endl;
fvMesh& mesh = _mesh();
volScalarField& T = _T();
// Setting BC at the cold side
coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
label coldSideSize = T.boundaryField()[coldSide_ind].size();
Tf.resize(coldSideSize);
set_valueFraction();
forAll(Tf, faceI)
{
    scalar faceZ =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
    Tf[faceI] = 298 - faceZ / 1.2 * 0.4; //Cooling mold
    //Tf[faceI] = 298 - faceZ / 1.2 * 8.7; //Steady mold
    //Tf[faceI] = 0;
}
// Setting BC at hotSide
hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
label hotSideSize = T.boundaryField()[hotSide_ind].size();
gTrue.resize(hotSideSize);
forAll(gTrue, faceI)
{
    scalar faceX =
        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
    scalar faceZ =
        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
    gTrue[faceI] = - 100000 * Foam::sqrt(faceX * faceX + faceZ * faceZ);
    //gTrue[faceI] = 0;
}
forAll(mesh.boundaryMesh(), patchI)
{
    if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
    {
        ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
    }
    else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
    {
        ITHACAutilities::assignBC(T, patchI, - gTrue / k);
    }
    else
    {
        ITHACAutilities::assignBC(T, patchI, homogeneousBC);
    }
}
