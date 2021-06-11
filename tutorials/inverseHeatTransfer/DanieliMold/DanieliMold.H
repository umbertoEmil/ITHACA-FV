#ifndef DanieliMold_H
#define DanieliMold_H


namespace ITHACAutilities
{
Eigen::MatrixXd colWiseMean(Eigen::MatrixXd matIn)
{
    Eigen::MatrixXd matOut;
    matOut.resize(matIn.cols(), 1);

    for (int col = 0; col < matIn.cols(); col++)
    {
        matOut(col) = 0;

        for (int row = 0; row < matIn.rows(); row++)
        {
            matOut(col) += matIn(row, col);
        }

        matOut(col) /= matIn.rows();
    }

    return matOut;
}
}

class DanieliMold_CG: public inverseLaplacianProblem_CG
{
    public:
        explicit DanieliMold_CG(int argc, char* argv[])
            :
            inverseLaplacianProblem_CG(argc, argv)
        {
            fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> lambdaFullField;
        PtrList<volScalarField> deltaTfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;

        scalar meniscusHeight = 1.11; //[m] from bottom of the mold
        scalar castWidth;

        List<int> thermocouplesPlaneCellID;
        const std::vector<label>* surfaceCellList;

        void getThermocouplesPlane()
        {
            Info << "getThermocouplesPlane not implemented, exiting" << endl;
            exit (120);
        }

        void readTCfromFile(word TCfile)
        {
#include"readTCfromFile.H"
        }

        Eigen::MatrixXd danieliTCtoMyNumeration(Eigen::MatrixXd Tdanieli)
        {
#include"danieliTCtoMyNumeration.H"
        }

        void postProcess(word folder, word heatFluxFieldName, bool syntheticData = 1,
                         scalar innerField = 0.0)
        {
#include"postProcess.H"
        }

        void setBC()
        {
#include"setBC.H"
        }

        void solveTrue()
        {
#include"solveTrue.H"
        }

        void solveTest()
        {
#include"solveTest.H"
        }
};

class DanieliMold_paramBC: public inverseLaplacianProblem_paramBC
{
    public:
        explicit DanieliMold_paramBC(int argc, char* argv[])
            :
            inverseLaplacianProblem_paramBC(argc, argv)
        {
            fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
        }

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;

        scalar meniscusHeight = 1.11; //[m] from bottom of the mold
        scalar castWidth;

        List<int> thermocouplesPlaneCellID;
        const std::vector<label>* surfaceCellList;

        void getThermocouplesPlane()
        {
            Info << "getThermocouplesPlane not implemented, exiting" << endl;
            exit (120);
        }

        void set_gBaseFunctions()
        {
#include "set_gBaseFunctions.H"
        }


        void readTCfromFile(word TCfile)
        {
#include"readTCfromFile.H"
        }

        Eigen::MatrixXd danieliTCtoMyNumeration(Eigen::MatrixXd Tdanieli)
        {
#include"danieliTCtoMyNumeration.H"
        }

        void postProcess(word folder, word heatFluxFieldName, bool syntheticData = 1,
                         scalar innerField = 0.0)
        {
#include"postProcess.H"
        }

        void setBC()
        {
#include"setBC.H"
        }

        void solveTrue()
        {
#include"solveTrue.H"
        }

        void solveTest()
        {
#include"solveTest.H"
        }
};

#endif
