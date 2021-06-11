Info << endl;
Info << "*********************************************************" << endl;
Info << "Testing parameterized BC TSVD with NOISY data" <<
     endl;
word outputFolder = "./ITHACAoutput/parameterizedBCnoiseTest_TSVD/";
volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
List<List<scalar>> heatFluxWeights;
scalar innerField = 1.0;
example_paramBC.set_gParametrized("rbf", RBFshapePar);
example_paramBC.parameterizedBCoffline();

Info << "Introducing error in the measurements" << endl;
Info << "Performing " << Ntests << " tests." << endl;
Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;

Eigen::VectorXi TSVDtruc(25);
TSVDtruc << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
         20, 22, 24, 26, 28, 30;

ITHACAstream::exportMatrix(TSVDtruc, "TSVDtrucVector", "eigen",
                           outputFolder);

std::cout << "debug : Tmeas orig = \n" << example_paramBC.Tmeas.transpose() <<
          std::endl;

int Ntrunc = TSVDtruc.size();

for (label i = 0; i < Ntests; i++)
{
    Info << "Test " << i << endl;
    example_paramBC.addNoise(noiseLevel);
    std::cout << "debug : Tmeas = \n" << example_paramBC.Tmeas.transpose() <<
              std::endl;
    List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
    List<scalar> solutionNorms;

    for (int truncI = 0; truncI < Ntrunc; truncI++)
    {
        example_paramBC.parameterizedBC("TSVD", TSVDtruc(truncI));
        volScalarField gParametrizedField = example_paramBC.list2Field(
                                                example_paramBC.g);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(i * Ntrunc + truncI + 1),
                                     outputFolder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example_paramBC.T,
                                     std::to_string(i * Ntrunc + truncI + 1),
                                     outputFolder,
                                     "T");
    }

    example_paramBC.Tmeas = TmeasOrig;
}

example_paramBC.postProcess(outputFolder, "gParametrized", innerField);
Info << "*********************************************************" << endl;
Info << endl;
