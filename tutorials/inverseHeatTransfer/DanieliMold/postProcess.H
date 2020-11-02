Info << "Computing errors" << endl;
fvMesh& mesh = _mesh();
gFullField.resize(0);
ITHACAstream::read_fields(gFullField, heatFluxFieldName,
                          folder);
label Nsolutions = gFullField.size();

if(syntheticData)
{
    volScalarField gTrueField = list2Field(gTrue, innerField);
    Eigen::MatrixXd heatFluxL2norm;
    heatFluxL2norm.resize(Nsolutions, 1);
    Eigen::MatrixXd heatFluxLinfNorm = heatFluxL2norm;
    forAll(gFullField, solutionI)
    {
        volScalarField gDiffField = gFullField[solutionI] - gTrueField;
        ITHACAstream::exportSolution(gDiffField,
                                     std::to_string(solutionI + 1), folder,
                                     "gDiffField");
        ITHACAstream::exportSolution(gTrueField,
                                     std::to_string(solutionI + 1), folder,
                                     "gTrueField");
        heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh, gDiffField,
                                    "hotSide") / ITHACAutilities::L2normOnPatch(mesh, gTrueField, "hotSide");
        heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh, gDiffField,
                                      "hotSide") / ITHACAutilities::LinfNormOnPatch(mesh, gTrueField, "hotSide");
    }
    ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                               folder);
    ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                               folder);
}
else
{
    Eigen::VectorXd totalHeatFluxVec(gFullField.size());
    forAll(gFullField, solutionI)
    {
        //Access the mesh information for the boundary
        label patchID = mesh.boundaryMesh().findPatchID("hotSide");
        const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
        //List of cells close to a boundary
        const labelUList& faceCells = cPatch.faceCells();
        scalar totalArea = 1.687; //m2
        scalar totalHeatFlux = 0;
        forAll(cPatch, faceI)
        {
            //id of the owner cell having the face
            label faceOwner = faceCells[faceI] ;
            scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
            totalHeatFlux += gFullField[solutionI][faceOwner] * faceArea;
        }
        totalHeatFlux /= totalArea;
        Info << "Heat flux = " << totalHeatFlux << " W/m2" << endl;
        Info << "HotSide area = " << totalArea << " m2" << endl;
        totalHeatFluxVec(solutionI) = totalHeatFlux; 
    }
    ITHACAstream::exportMatrix(totalHeatFluxVec, "totalHeatFlux", "eigen",
                               folder);

}
Info << "Postprocessing done\n" << endl;
