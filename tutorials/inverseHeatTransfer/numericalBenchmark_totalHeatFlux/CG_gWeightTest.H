Info << endl;
Info << "*********************************************************" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/CG_gWeightTest/";
volScalarField gTrueField = example_CG.list2Field(example_CG.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;

int Ntests = 300;
Eigen::VectorXd gWeightsVector = expList(1e-10, 1e2, Ntests);
ITHACAstream::exportMatrix(gWeightsVector, "gWeightsVector", "eigen",
                           outputFolder);

for (int i = 0; i < Ntests; i++)
{
    Info << "Test " << i + 1 << endl;
    example_CG.gIntegralWeight = gWeightsVector(i);

    if (example_CG.conjugateGradient())
    {
        Info << "CG converged" << endl;
        volScalarField heatFluxField = example_CG.list2Field(
                                           example_CG.gList[example_CG.gList.size() - 1]);
        ITHACAstream::exportSolution(heatFluxField,
                                     std::to_string(i + 1), outputFolder,
                                     "g_CG");
        Info << "************************************" << endl;
        Info << endl << endl;
    }
    else
    {
        Info << "CG did not converged" << endl;
        Info << "************************************" << endl;
        Info << endl << endl;
    }
}

example_CG.postProcess(outputFolder, "g_CG");
Info << "*********************************************************" << endl;
Info << endl;
