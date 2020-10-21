Info << "Computing errors" << endl;
gFullField.resize(0);
ITHACAstream::read_fields(gFullField, heatFluxFieldName,
                          folder);
volScalarField gTrueField = list2Field(gTrue, innerField);
label Nsolutions = gFullField.size();
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
    scalar EPS = 1e-6;
    volScalarField relativeErrorField(gDiffField);

    for (label i = 0; i < relativeErrorField.internalField().size(); i++)
    {
        if (std::abs(gTrueField.ref()[i]) < EPS)
        {
            relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
        }
        else
        {
            relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / std::abs(
                                              gTrueField.ref()[i]);
        }
    }

    ITHACAstream::exportSolution(relativeErrorField,
                                 std::to_string(solutionI + 1), folder,
                                 "relativeErrorField");
    heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh,
                                relativeErrorField,
                                "hotSide");
    heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                  relativeErrorField,
                                  "hotSide");
}
ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                           folder);
ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                           folder);
