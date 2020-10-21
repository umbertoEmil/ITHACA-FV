Info << endl;
Info << "*********************************************************" << endl;
Info << "Testing parameterized BC inverse solver with NOISY data" <<
     endl;
word outputFolder = "./ITHACAoutput/CGandParamErrorTest/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
List<List<scalar>> heatFluxWeights;
Eigen::VectorXd residualNorms;
scalar innerField = 1.0;
example_paramBC.set_gParametrized("rbf", rbfShapePar);
example_paramBC.parameterizedBCoffline();

Info << "Introducing error in the measurements" << endl;
Info << "Performing " << Ntests << " tests." << endl;
Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;


Eigen::VectorXd varianceVect = example_paramBC.Tmeas * noiseLevel;
varianceVect.cwiseProduct(varianceVect);
example_CG.Jtol = example_CG.Tmeas.size() * varianceVect.maxCoeff() *
                  varianceVect.maxCoeff() / 2;
std::cout << "Jtol = " << example_CG.Jtol << std::endl;
example_CG.saveSolInLists = 1;

for (label i = 0; i < Ntests; i++)
{
    Info << "Test " << i << endl;
    example_paramBC.addNoise(noiseLevel);
    example_CG.addNoise(noiseLevel);
    //TSVD solution
    example_paramBC.parameterizedBC("TSVD", 3);
    volScalarField gParametrizedField = example_paramBC.list2Field(
                                            example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i * 3 + 1),
                                 outputFolder,
                                 "g");
    //LU solution
    example_paramBC.parameterizedBC("fullPivLU");
    gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(i * 3 + 2),
                                 outputFolder,
                                 "g");
    //CG solution
    example_CG.conjugateGradient();
    volScalarField heatFluxField = example_CG.list2Field(
                                       example_CG.gList[example_CG.gList.size() - 1]);
    ITHACAstream::exportSolution(heatFluxField,
                                 std::to_string(i * 3 + 3), outputFolder,
                                 "g");
    example_paramBC.Tmeas = TmeasOrig;
    example_CG.Tmeas = TmeasOrig;
}

example_paramBC.postProcess(outputFolder, "g", innerField);
Info << "*********************************************************" << endl;
Info << endl;
