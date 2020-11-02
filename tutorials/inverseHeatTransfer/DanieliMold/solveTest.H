Info << "debug: solving test" << endl;
volScalarField& T = _T();
ITHACAstream::exportSolution(T, "2", "./ITHACAoutput/true/", "test");
simpleControl& simple = _simple();

ITHACAutilities::assignIF(T, homogeneousBC);
Info << "DT = " << DT << endl;

while (simple.loop())
{
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::laplacian(DT, T)
        );
        TEqn.solve();
    }
}

Info << "Exporting" << endl;
ITHACAstream::exportSolution(T, "1", "./ITHACAoutput/true/", "test");
//Reinitializing
restart();
