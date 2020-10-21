//Reinitializing
//restart();
setBC();

volScalarField& T = _T();
simpleControl& simple = _simple();
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
ITHACAstream::exportSolution(T, "1", "./ITHACAoutput/true/", T.name());
