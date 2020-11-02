fvMesh& mesh = _mesh();
volScalarField& T(_T());
// Setting BC at the cold side
coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
label coldSideSize = T.boundaryField()[coldSide_ind].size();
Tf.resize(coldSideSize);
refGrad.resize(coldSideSize);
valueFraction.resize(coldSideSize);
forAll(Tf, faceI)
{
    scalar faceZ =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
    scalar faceY =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].y();
    scalar faceX =
        mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].x();

    Tf[faceI] = (Tf_delta * (1.2 - faceZ) + Tf_0);
}
