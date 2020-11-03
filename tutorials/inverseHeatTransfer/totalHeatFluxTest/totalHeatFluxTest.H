#ifndef totalHeatFluxTest_H
#define totalHeatFluxTest_H



class totalHeatFluxTest_CG: public inverseLaplacianProblemTotalHeatMeasure_CG
{
    public:
        explicit totalHeatFluxTest_CG(int argc, char* argv[])
            :
            inverseLaplacianProblemTotalHeatMeasure_CG(argc, argv),
            T(_T()),
            lambda(_lambda()),
            deltaT(_deltaT()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }
        volScalarField& T;
        volScalarField& lambda;
        volScalarField& deltaT;
        fvMesh& mesh;
        Time& runTime;

        double g_0;
        double g_Z;
        double g_X;
        double Tf_0;
        double Tf_delta;
        double d;

        //Boudary Conditions
        List<scalar> heatFlux_hotSide;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include "postProcess.H"
        }

	void set_Tf()
        {
#include "setTf.H"
	}

        void solveTrue()
        {
#include"solveTrue.H"
        }

};

class totalHeatFluxTest_paramBC: public inverseLaplacianProblemTotalHeatMeasure_paramBC
{
    public:
        explicit totalHeatFluxTest_paramBC(int argc, char* argv[])
            :
            inverseLaplacianProblemTotalHeatMeasure_paramBC(argc, argv),
            T(_T()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
        }
        volScalarField& T;
        fvMesh& mesh;
        Time& runTime;

        double g_0;
        double g_Z;
        double g_X;
        double Tf_0;
        double Tf_delta;
        double d;

        //Boudary Conditions
        List<scalar> heatFlux_hotSide;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include "postProcess.H"
        }

	void set_Tf()
        {
#include "setTf.H"
	}

        void solveTrue()
        {
#include"solveTrue.H"
        }

};
#endif
