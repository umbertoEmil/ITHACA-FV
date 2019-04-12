#ifndef DanieliMold_H
#define DanieliMold_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "triSurfaceMesh.H"
#include <vector>
#include "intersectionSurface.H"


class DanieliMold: public inverseLaplacianProblem
{
    public:
        explicit DanieliMold(int argc, char* argv[])
            :
            inverseLaplacianProblem(argc, argv)
        {
	    fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }
        //volScalarField& T;
        //volScalarField& lambda;
        //volScalarField& deltaT;
        //fvMesh& mesh;
        //Time& runTime;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> lambdaFullField;
        PtrList<volScalarField> deltaTfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;

        List<int> thermocouplesPlaneCellID;
        const std::vector<label>* surfaceCellList;

        void getThermocouplesPlane()
        {
            Info << "getThermocouplesPlane not implemented, exiting" << endl;
	    exit (120);
            //fvMesh& mesh = _mesh();
	    //volScalarField& T = _T();
	    //intersectionSurface thermocouplesSurf("thermocouplesSurface", mesh); 
            //ITHACAutilities::assignIF(T, homogeneousBC);
	    //volScalarField intersectionSurfField(T);
	    //thermocouplesSurf.createImmersedBody(intersectionSurfField);
            //ITHACAstream::exportSolution(intersectionSurfField, "1", "./ITHACAoutput/debug/",
            //                             "intersectionSurfField");
	    //surfaceCellList = thermocouplesSurf.getSurfaceCellList();

	    //ITHACAutilities::assignIF(intersectionSurfField, homogeneousBC);
	    //for(unsigned int cellI=0; cellI<surfaceCellList->size(); cellI++)
	    //{
	    //    intersectionSurfField.ref()[surfaceCellList->operator[](cellI)] = 1;
	    //}
            //ITHACAstream::exportSolution(intersectionSurfField, "1", "./ITHACAoutput/debug/",
            //                             "intersectionSurfField2");


        }

	void set_gBaseFunctions(word type, scalar shapeParameter = 1.0)
        {
            volScalarField& T = _T();
            fvMesh& mesh = _mesh();
            Info << "Radial Basis Functions are used." << endl;
            // The center of each function is the projection of each thermocouple
            // on the boundary hotSide
            if (thermocouplesCellID.size() == 0)
            {
                readThermocouples();
            }

            gBaseFunctions.resize(thermocouplesNum);
            gWeights.resize(thermocouplesNum);
            //for (label funcI = 0; funcI < thermocouplesNum; funcI++)
            forAll(gBaseFunctions, funcI)
            {
                gBaseFunctions[funcI].resize(T.boundaryField()[hotSide_ind].size());
		scalar thermocoupleX = thermocouplesPos[funcI].x();
                scalar thermocoupleZ = thermocouplesPos[funcI].z();
                forAll (T.boundaryField()[hotSide_ind], faceI)
                {
                    scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                    scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                    scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - thermocoupleX) +
                                               (faceZ - thermocoupleZ) * (faceZ - thermocoupleZ));
                    if(faceX > - castWidth / 2 && faceX < castWidth / 2)
                    {
                        gBaseFunctions[funcI][faceI] = Foam::exp(- (shapeParameter *
                                                       shapeParameter
                                                       * radius * radius));
                    }
                    else
                    {
                        gBaseFunctions[funcI][faceI] = 0;
                    }
                }
            }
	}


        List<int> EigenMatrix2List ( Eigen::MatrixXi matrix )
	{
	    if(matrix.cols() == 1)
	    {
	        List<int> list(matrix.rows());
	        for(label i = 0; i < matrix.rows(); i++)
	        {
	            list[i] = matrix(i,0);
	        }
	        return list;
	    }
	    else
	    {
	        Info << "I still must code this function, exiting" << endl;
		exit(11);
	    }
	}

        //void readThermocouples()
        //{
        //    if (!thermocouplesRead)
        //    {
	//        word fileName = "./thermocouplesCellsID"; 
	//        if(ITHACAutilities::check_file(fileName + "_mat.txt"))
	//	{
	//	    Info << "Reading thermocouples cells from file" << endl;
	//	    thermocouplesCellID = EigenMatrix2List(ITHACAstream::readMatrix(fileName + "_mat.txt").cast <int> ());
	//	}
	//	else
	//	{
        //            Info << "Defining positions of thermocouples" << endl;
        //            fvMesh& mesh = _mesh();
        //            volScalarField& T = _T();
	//            thermocouplesCellID.resize(thermocouplesPos.size());
	//	    forAll(thermocouplesPos, tcI)
	//	    {
	//	        thermocouplesCellID[tcI] = mesh.findCell(thermocouplesPos[tcI]);
	//	    }
	//	    
	//	    volScalarField thermocouplesField(T);
	//	    ITHACAutilities::assignIF(thermocouplesField, homogeneousBC);
        //            forAll(thermocouplesCellID, tcI)
        //            {
        //                thermocouplesField.ref()[thermocouplesCellID[tcI]] = 1;
        //            }
        //            ITHACAstream::exportSolution(thermocouplesField, "1", "./ITHACAoutput/debug/",
        //                                         "thermocouplesField,");

        //            Eigen::MatrixXi thermocouplesCellID_eigen = List2EigenMatrix(thermocouplesCellID);
        //            ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
        //                         "eigen", "./");

        //            ////Define thermocouples plane
        //            //defineThermocouplesPlane();
        //            //Info << "Exporting debug field" << endl;
        //            //ITHACAutilities::assignIF(T, homogeneousBC);
        //            //forAll(thermocouplesPlaneCellID, cellI)
        //            //{
        //            //    T.ref()[thermocouplesPlaneCellID[cellI]] = 1000;
        //            //}
        //            //ITHACAstream::exportSolution(T, "1", "./ITHACAoutput/debug/",
        //            //                             "thermocouplesSurface");
	//	}
	//	thermocouplesRead = 1;
        //    }
        //    else
        //    {
        //        WarningInFunction << "readThermocouples function called twice." << endl;
        //        WarningInFunction << "I am not doing the second reading." << endl;
        //    }
        //}

	Eigen::MatrixXd colWiseMean(Eigen::MatrixXd matIn)
	{
	    Eigen::MatrixXd matOut;
	    matOut.resize(matIn.cols(),1);
	    for(int col = 0; col < matIn.cols(); col++)
	    {
		matOut(col) = 0;
	        for(int row = 0; row < matIn.rows(); row++)
		{
		    matOut(col) += matIn(row,col);
		}
		matOut(col) /= matIn.rows();
	    }
	    return matOut;
	}

	void readTCfromFile(word TCfile)
	{
	    // The TC file must have each TC on a column
	    // This function will mean the temperatures on the rows
	    Eigen::MatrixXd Tdanieli;
	    Tdanieli = cnpy::load(Tdanieli, TCfile, "colMajor");
	    Info << "Thermocouples temperature matrix shape = " << Tdanieli.rows() << " x " << Tdanieli.cols() << endl;
	    //std::cout << "Tdanieli = " << Tdanieli << std::endl;
	    Tdanieli = Tdanieli.array() + 273.15; // C to K
	    Tdanieli = colWiseMean(Tdanieli);
	    //std::cout << "Tdanieli = " << Tdanieli << std::endl;
	    label TCnumber = thermocouplesCellID.size();
	    std::string message = "The file " + TCfile + " must contain " + std::to_string(TCnumber) + " elements. It has " + std::to_string(Tdanieli.size());
	    M_Assert(Tdanieli.size() == TCnumber,
	    message.c_str());
	    Tmeas = danieliTCtoMyNumeration(Tdanieli);
	}

	Eigen::MatrixXd danieliTCtoMyNumeration(Eigen::MatrixXd Tdanieli)
	{
	    Eigen::MatrixXd Tout = Tdanieli;
	    label rows = 10;
	    label cols = 9;
	    M_Assert(Tdanieli.size() == rows * cols,
	    "Reprogram the function danieliTCtoMyNumeration with the right TC numbers.");
	    label i = 0;
	    for(int col = 0; col < cols; col++)
	    {
	        for(int row = 0; row < rows; row++)
		{
		    //Info << i << " = " << (rows - 1 - row) * cols + col << endl;
		    Tout(i) = Tdanieli((rows - 1 - row) * cols + col);
		    i++;
		}
	    }
	    return Tout;
	}

        volScalarField analyticalErrorField (volScalarField& field1,
                                             volScalarField& field2)
        {
            volScalarField fieldDiff(field1 - field2);
            return fieldDiff;
        }

        void computeTdiff()
        {
            TdiffField = TredField;
            forAll(TdiffField, fieldI)
            {
                TdiffField[fieldI] = TredField[fieldI] - TfullField[0];
            }
        }

        int isInPlane(double cx, double cz)
        {
            return (cx >= interpolationPlane.minX -
                    interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz >= interpolationPlane.minZ - interpolationPlane.thermocoupleCellDim[2] / 4 &&
                    cx <= interpolationPlane.maxX + interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz <= interpolationPlane.maxZ + interpolationPlane.thermocoupleCellDim[2] / 4
                   );
        }

        void postProcess(word folder, word heatFluxFieldName, bool syntheticData = 1, scalar innerField = 0.0)
        {
            Info << "Computing errors" << endl;
            fvMesh& mesh = _mesh();
            gFullField.resize(0);
            ITHACAstream::read_fields(mesh, gFullField, heatFluxFieldName,
                                      folder);
            label Nsolutions = gFullField.size();
            
	    if(syntheticData)
	    {
	        volScalarField gTrueField = list2Field(gTrue, innerField);
                Eigen::MatrixXd heatFluxL2norm;
                heatFluxL2norm.resize(Nsolutions, 1);
                Eigen::MatrixXd heatFluxLinfNorm = heatFluxL2norm;
                forAll(gFullField, solutionI)
                {
                    volScalarField gDiffField = gFullField[solutionI] - gTrueField;
                    ITHACAstream::exportSolution(gDiffField,
                                                 std::to_string(solutionI + 1), folder,
                                                 "gDiffField");
                    ITHACAstream::exportSolution(gTrueField,
                                                 std::to_string(solutionI + 1), folder,
                                                 "gTrueField");
                    heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh, gDiffField,
                                                "hotSide") / ITHACAutilities::L2normOnPatch(mesh, gTrueField, "hotSide");
                    heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh, gDiffField,
                                                  "hotSide") / ITHACAutilities::LinfNormOnPatch(mesh, gTrueField, "hotSide");
                }
                ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                                           folder);
                ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                                           folder);
	    }
	    else
	    {
	        forAll(gFullField, solutionI)
                {
	            //Access the mesh information for the boundary
                    label patchID = mesh.boundaryMesh().findPatchID("hotSide");
                    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
                    //List of cells close to a boundary
                    const labelUList& faceCells = cPatch.faceCells();
		    scalar totalArea = 1.687; //m2
		    scalar totalHeatFlux = 0;
                    forAll(cPatch, faceI)
                    {
                        //id of the owner cell having the face
                        label faceOwner = faceCells[faceI] ;
                        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];
                        totalHeatFlux += gFullField[solutionI][faceOwner] * faceArea;
                    }
		    totalHeatFlux /= totalArea;
		    Info << "Heat flux = " << totalHeatFlux << " W/m2" << endl;
		    Info << "HotSide area = " << totalArea << " m2" << endl;
		}

	    }
            Info << "Postprocessing done\n" << endl;
	}

	void setBC()
	{
	    Info << "Setting BC\n" << endl;
            fvMesh& mesh = _mesh();
	    volScalarField& T = _T();
            // Setting BC at the cold side
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            label coldSideSize = T.boundaryField()[coldSide_ind].size();
            Tf.resize(coldSideSize);
            set_valueFraction();
            forAll(Tf, faceI)
            {
                scalar faceZ =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
                //Tf[faceI] = 323 - faceZ / 1.2 * 0.4; //Cooling mold
                Tf[faceI] = 323 - faceZ / 1.2 * 8.7; //Steady mold
		//Tf[faceI] = 0;
            }
            // Setting BC at hotSide
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
	    gTrue.resize(hotSideSize);
            forAll(gTrue, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                scalar faceZ =
                    mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                gTrue[faceI] = - 100000 * Foam::sqrt(faceX * faceX + faceZ * faceZ);
		//gTrue[faceI] = 0;
            }
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(T, patchI, - gTrue / k);
                }
                else
                {
                    ITHACAutilities::assignBC(T, patchI, homogeneousBC);
                }
            }
	}

        void solveTrue()
        {
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
        }

        void solveTest()
        {
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
        }
};

#endif
