Info << endl;
Info << "*********************************************************" << endl;
Info << "Testing inverse solver with NOISY data" << endl;
Info << "Performing " << Ntests << " tests." << endl;
Info << endl;
word outputFolder;
outputFolder = "./ITHACAoutput/ParamBCnoiseLevelTest/";

Eigen::VectorXd noiseLevelVec = Eigen::VectorXd::LinSpaced(16, 0.01, 0.5);
ITHACAstream::exportMatrix(noiseLevelVec, "noiseLevelVec", "eigen",
                           outputFolder);

example_paramBC.set_gParametrized("rbf", RBFshapePar);
example_paramBC.parameterizedBCoffline();
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;



for (label NLi = 0; NLi < noiseLevelVec.size() ; NLi++)
{
    noiseLevel = noiseLevelVec(NLi);

    for (label i = 0; i < Ntests; i++)
    {
        Info << "Test " << i << endl;
        example_paramBC.addNoise(noiseLevel);
        example_paramBC.parameterizedBC("TSVD", 5);
        volScalarField gParametrizedField = example_paramBC.list2Field(
                                                example_paramBC.g);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(NLi * Ntests + i + 1),
                                     outputFolder,
                                     "g");
        example_paramBC.Tmeas = TmeasOrig;
    }
}

example_paramBC.postProcess(outputFolder, "g");
Info << "*********************************************************" << endl;
Info << endl;
