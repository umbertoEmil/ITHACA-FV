volScalarField& T = _T();
fvMesh& mesh = _mesh();
Foam::Time& runTime = _runTime();
//ITHACAutilities::assignIF(T, homogeneousBC);
simpleControl simple(mesh);
Info << "DT = " << DT << endl;
assignDirectBC();
word outputFolder = "./ITHACAoutput/true/";

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
ITHACAstream::exportSolution(T, "1", outputFolder, T.name());
volScalarField gTrueField = list2Field(g);
ITHACAstream::exportSolution(gTrueField,
                             "1", outputFolder,
                             "gTrue");


//Reinitializing runTime
instantList Times = runTime.times();
runTime.setTime(Times[1], 0);
