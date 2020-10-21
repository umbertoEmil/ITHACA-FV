Info << endl;
Info << "*********************************************************" << endl;
Info << "Testing CG inverse solver with NOISY data" << endl;
Info << "Performing " << Ntests << " tests." << endl;
Info << endl;
word outputFolder;

if (CGnoiseLevelTest)
{
    outputFolder = "./ITHACAoutput/CGnoiseLevelTest/";
}
else if (ParamBCnoiseLevelTest)
{
    outputFolder = "./ITHACAoutput/ParamBCnoiseLevelTest/";
    example_paramBC.set_gParametrized("rbf", RBFshapePar);
    example_paramBC.parameterizedBCoffline();
}

volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");
example_CG.saveSolInLists = 1;
Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;
double gIntegral_measOrig = example_paramBC.gIntegral_meas;

Eigen::VectorXd noiseLevelVec(10);
noiseLevelVec << .1, .2, .3, .4, .5, .6, .7, .8, .9, 1;


for (label NLi = 0; NLi < noiseLevelVec.size() ; NLi++)
{
    double noiseLevel = noiseLevelVec(NLi);

    for (label i = 0; i < Ntests; i++)
    {
        Info << "Test " << i << endl;
        example_paramBC.addNoise(noiseLevelTC, noiseLevel);
        example_CG.Tmeas = example_paramBC.Tmeas;
        example_CG.gIntegral_meas = example_paramBC.gIntegral_meas;
        Info << "Stopping for J < " << example_CG.Jtol << endl;

        if (CGnoiseLevelTest)
        {
            if (example_CG.conjugateGradient())
            {
                Info << "CG converged" << endl;
                volScalarField heatFluxField = example_CG.list2Field(
                                                   example_CG.gList[example_CG.gList.size() - 1]);
                ITHACAstream::exportSolution(heatFluxField,
                                             std::to_string(NLi * Ntests + i + 1), outputFolder,
                                             "g");
                Info << "debug = " << NLi* noiseLevelVec.size() + i + 1 << endl;
                Info << "************************************" << endl;
                Info << endl << endl;
            }
            else
            {
                Info << "CG did not converged" << endl;
                Info << "************************************" << endl;
                i--;
                Info << endl << endl;
            }
        }
        else if (ParamBCnoiseLevelTest)
        {
            example_paramBC.parameterizedBC("TSVD", TSVDtrunc);
            volScalarField gParametrizedField = example_paramBC.list2Field(
                                                    example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(NLi * Ntests + i + 1),
                                         outputFolder,
                                         "g");
        }

        example_paramBC.Tmeas = TmeasOrig;
        example_CG.Tmeas = TmeasOrig;
        example_paramBC.gIntegral_meas = gIntegral_measOrig;
        example_CG.gIntegral_meas = gIntegral_measOrig;
    }
}

example_paramBC.postProcess(outputFolder, "g");
Info << "*********************************************************" << endl;
Info << endl;
