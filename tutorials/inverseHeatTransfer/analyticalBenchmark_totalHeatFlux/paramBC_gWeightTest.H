Info << endl;
Info << "*********************************************************" << endl;
Info << endl;
word outputFolder = "./ITHACAoutput/paramBC_gWeightTest/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");

int Ntests = 150;
Info << "debug 1" << endl;
Eigen::VectorXd gWeightsVector = expList(1e-5, 1e3, Ntests);
ITHACAstream::exportMatrix(gWeightsVector, "gWeightsVector", "eigen",
                           outputFolder);
example_paramBC.set_gParametrized("rbf", 0.1);
example_paramBC.parameterizedBCoffline();


for (int i = 0; i < Ntests; i++)
{
    Info << "Test " << i + 1 << endl;
    example_paramBC.gIntegralWeight = gWeightsVector(i);
    example_paramBC.parameterizedBC("fullPivLU");
    volScalarField gParametrizedField = example_paramBC.list2Field(
                                            example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i + 1), outputFolder,
                                 "g_paramBC");
}

example_paramBC.postProcess(outputFolder, "g_paramBC");
Info << "*********************************************************" << endl;
Info << endl;
