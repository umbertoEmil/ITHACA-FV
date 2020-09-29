#ifndef EnKFtest_1DheatTransfer_H
#define EnKFtest_1DheatTransfer_H
#include "simpleControl.H"
#include "fvOptions.H"


class EnKFtest_1DheatTransfer: public laplacianProblem
{
    public:
        explicit EnKFtest_1DheatTransfer(int argc, char* argv[], int Nseeds_)
            :
            laplacianProblem(argc, argv)
        {
            fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
            simpleControl& simple = _simple();
            Time& runTime = _runTime();
#include "createFvOptions.H"
            left_ind = mesh.boundaryMesh().findPatchID("left");
            Nseeds = Nseeds_;
            Tensemble = PtrList<volScalarField>(Nseeds);
            volScalarField& T(_T());
            stateSize = T.size();
            priorSamples = Eigen::MatrixXd::Zero(stateSize, Nseeds);
            posteriorSamples = priorSamples;

            startTime = runTime.startTime().value();
            deltaTime = runTime.deltaTValue();
            endTime = runTime.endTime().value();
            Ntimes = (endTime - startTime) / deltaTime;
            Info << "startTime = " << startTime << endl;
            Info << "endTime = " << endTime << endl;
            Info << "deltaTime = " << deltaTime << endl;
            Info << "Ntimes = " << Ntimes << endl;
            timeVector.resize(Ntimes);
            probe_mean.resize(Nprobes, Ntimes);
            probe_true = probe_mean;
            probe_MaxConfidence = probe_mean;
            probe_minConfidence = probe_mean;
        }

        autoPtr<simpleControl> _simple;
        autoPtr<fv::options> _fvOptions;
        double k;
        double rho;
        double Cp;
        label left_ind;
        Eigen::MatrixXd measurements;
        int NtimeStepsBetweenSamps = 3;
        int Nseeds;
        PtrList<volScalarField> Tensemble;
        int stateSize;
        
        std::shared_ptr<muq::Modeling::Gaussian> priorDensity;
        std::shared_ptr<muq::Modeling::Gaussian> modelErrorDensity;
        std::shared_ptr<muq::Modeling::Gaussian> measNoiseDensity;
        double meas_cov;

        Eigen::MatrixXd priorSamples;
        Eigen::MatrixXd posteriorSamples;

        double startTime; 
        double deltaTime; 
        double endTime;   
        int Ntimes;
        
        Eigen::VectorXd timeVector;
        Foam::vector probePosition = Foam::vector(0.5, 0.1, 0.005);
        int Nprobes = 1;
        Eigen::MatrixXd probe_true;
        Eigen::MatrixXd probe_mean;
        Eigen::MatrixXd probe_MaxConfidence;
        Eigen::MatrixXd probe_minConfidence;

        int samplingTimeI = 0;
        int sampleI = 0;
        int sampleFlag = 0;

        double updateBC(bool reconstruction = 0)
        {
            Time& runTime = _runTime();
            scalar time = runTime.value();
            double BC = time * 5 + 2;
            if(reconstruction)
            {
                BC = time * 5 + 0;
            }
            return BC;
        }

        void resetSamplingCounters()
        {
            samplingTimeI = 0;
            sampleI = 0;
            sampleFlag = 0;
        }

        Eigen::VectorXd observe(volScalarField field)
        {
            fvMesh& mesh = _mesh();
            Time& runTime = _runTime();
            IOdictionary measurementsDict
            (
                IOobject
                (
                    "measurementsDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            
            List<vector> measurementPoints(measurementsDict.lookup("positions"));
            Eigen::VectorXd measures(measurementPoints.size());

            forAll(measurementPoints, pntI)
            {
                //TODO
                //add check if I put two measurements on the same cell
                measures(pntI) = field[mesh.findCell(measurementPoints[pntI])];
            }
            return measures;
        }

        Eigen::MatrixXd solve(volScalarField& T, word outputFolder, word outputFieldName, bool reconstruction = 0)
        {
            Eigen::MatrixXd obsMat;
            fvMesh& mesh = _mesh();
            Foam::Time& runTime = _runTime();
            simpleControl& simple = _simple();
	    fv::options& fvOptions(_fvOptions());

	    dimensionedScalar diffusivity("diffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0), k / (rho * Cp));
	    Info << "Thermal diffusivity = " << diffusivity << " m2/s"<< endl;

            int timeI = 0;
            while (runTime.loop())
            {
                scalar time = runTime.value();
                if(!reconstruction)
                {
                    sampleFlag++;
                    timeVector(timeI) = time;
                    timeI++;
                }

                ITHACAutilities::assignBC(T, mesh.boundaryMesh().findPatchID("left"), updateBC(reconstruction));
            
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::ddt(T) - fvm::laplacian(diffusivity, T)
                    );
            
                    fvOptions.constrain(TEqn);
                    TEqn.solve();
                    fvOptions.correct(T);
                }

                if(reconstruction)
                {
                    //Adding model error
                    forAll(T.internalField(), cellI)
                    {
                        T.ref()[cellI] += modelErrorDensity->Sample()(cellI);
                    }
                }
                else
                {
                    probe_true.col(timeI - 1) = sampleField(T, probePosition);
                }
                ITHACAstream::exportSolution(T, runTime.timeName(),
                                 outputFolder,
                                 outputFieldName);
                if(sampleFlag == NtimeStepsBetweenSamps)
                {
                    Info << "Sampling at time = " << runTime.timeName() << nl << endl;
                    obsMat.conservativeResize(observe(T).size(), obsMat.cols()+1);
                    obsMat.col(sampleI) = observe(T);
                    if(!reconstruction)
                    {
                        sampleFlag = 0;
                        sampleI++;
                    }
                }
            
                runTime.write();
            }

            if(!reconstruction)
            {
                ITHACAstream::exportVector(timeVector, "timeVector", "eigen", outputFolder);
                ITHACAstream::exportMatrix(probe_true, "probe_true", "eigen", outputFolder);
            }
            return obsMat;
        }

	void solveDirect()
	{
            Info << "\n****************************************\n" << endl;
            Info << "\nPerforming true solution\n" << endl;

	    word outputFolder = "./ITHACAoutput/direct/";
            restart();

            volScalarField& T(_T());
            measurements = solve(T, outputFolder, "Tdirect");

            ITHACAstream::exportMatrix(measurements, "trueMeasurements", "eigen", outputFolder);
            std::cout << "Number of samples in time = " << measurements.cols() << std::endl;
            std::cout << "Number of sample in space = " << measurements.rows() << std::endl;
            resetSamplingCounters();
            Info << "\nEND true solution\n" << endl;
            Info << "\n****************************************\n" << endl;
	}

        void reconstruction()
        {
            restart();
            Time& runTime = _runTime();
            resetSamplingCounters();
            std::vector<Eigen::MatrixXd> observationTensor;
            word outputFolder = "ITHACAoutput/reconstuction";

            Info << "\n****************************************\n" << endl;
            Info << "\nStarting reconstruction\n" << endl;
            int samplesCounter = 0;
            for(int timeI = 0; timeI < Ntimes; timeI++)
            {

                double initialTime = timeI * deltaTime + startTime;
                Info << "\n\nTime " << initialTime + deltaTime << endl; 

                Eigen::MatrixXd forecastOutput = forecastStep(initialTime, timeI, initialTime + deltaTime); 

                if(forecastOutput.cols() > 0)
                {
                    samplesCounter++;
                    Info << "debug : YES sampling step" << endl;
                    //This is a sampling step
                    observationTensor.push_back(forecastOutput);
                    sampleFlag = 0;
                    Eigen::VectorXd meas = measurements.col(samplesCounter - 1);

                    //Kalman filter
                    posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(Tensemble, meas, Eigen::MatrixXd::Identity(meas.size(), meas.size()) * meas_cov, forecastOutput);

                    for(int i = 0; i < Nseeds; i++)
                    {
                        restart();
                        volScalarField& T = _T();
                        Eigen::VectorXd internalField = posteriorSamples.col(i);
                        assignIF(T, internalField);
                        Tensemble.set(i, T);
                    }
                }
                else
                {
                    Info << "debug : NOT a sampling step" << endl;
                }

                Eigen::MatrixXd ensambleProbe(Nprobes, Nseeds);
                for(int i = 0; i < Nseeds; i++)
                {
                    ensambleProbe.col(i) = sampleField(Tensemble[i], probePosition);
                }
                probe_mean.col(timeI) = ensambleProbe.rowwise().mean();
                
                probe_MaxConfidence.col(timeI) = ITHACAmuq::muq2ithaca::quantile(ensambleProbe, 0.95);
                probe_minConfidence.col(timeI) = ITHACAmuq::muq2ithaca::quantile(ensambleProbe, 0.05);
            }
            ITHACAstream::exportMatrix(probe_mean, "probe_mean", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe_MaxConfidence, "probe_MaxConfidence", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe_minConfidence, "probe_minConfidence", "eigen", outputFolder);
            Info << "\n****************************************\n" << endl;

        }


        Eigen::MatrixXd forecastStep(double startTime_, int startIndex_, double endTime_)
        {
            word outputFolder = "ITHACAoutput/forwardSamples";
            bool modelError = 1;
            Eigen::MatrixXd observationMat;
            sampleFlag++;
            for(int i = 0; i < Nseeds; i++)
            {
                word fieldName = "Tsample" + std::to_string(i);
                volScalarField T(_T());
                if(startIndex_ == 0)
                {
                    restart();
                    T = _T();
                    Eigen::VectorXd internalField = priorSamples.col(i);
                    assignIF(T, internalField);
                }
                else
                {
                    T = Tensemble.release(i);
                }
                resetRunTime(startTime_, startIndex_, endTime_);


                if(observationMat.cols() == 0)
                {
                    observationMat = solve(T, outputFolder, fieldName, modelError);
                }
                else
                {
                    observationMat.conservativeResize(observationMat.rows(), observationMat.cols()+1);
                    Eigen::MatrixXd observationForCheck = solve(T, outputFolder, fieldName, modelError);
                    M_Assert(observationForCheck.cols() == 1, "The forecast step passed trough more than a sampling step");
                    observationMat.col(i) = observationForCheck;
                }

                Tensemble.set(i, T);
            }
            return observationMat;
        }

        void assignIF(volScalarField& field_, Eigen::VectorXd internalField_)
        {
            for(int i = 0; i < internalField_.size(); i++)
            {
                field_.ref()[i] = internalField_(i);
            }
        }

        Eigen::VectorXd sampleField(volScalarField field_, Foam::vector probeLocation_)
        {
            Foam::fvMesh& mesh = _mesh();
            Eigen::VectorXd output(1);
            output(0) = field_[mesh.findCell(probeLocation_)];

            return output;
        }

        void restart()
        {
            Time& runTime = _runTime();
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 0);
            _simple.clear();
        
            _T.clear();
        
            Foam::fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
                 (
                     new simpleControl
                     (
                         mesh
                     )
                 );
        
            _T = autoPtr<volScalarField>
                 (
                     new volScalarField
                     (
                         IOobject
                         (
                             "T",
                             runTime.timeName(),
                             mesh,
                             IOobject::MUST_READ,
                             IOobject::AUTO_WRITE
                         ),
                         mesh
                     )
                 );
        }

        void resetRunTime(double startTime_, int startIndex_, double endTime_)
        {
            Time& runTime = _runTime();
            instantList Times = runTime.times();
            runTime.setTime(startTime_, startIndex_);
            runTime.setEndTime(endTime_);

            _simple.clear();
            Foam::fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
                 (
                     new simpleControl
                     (
                         mesh
                     )
                 );
        }

        void priorSetup(double mean, double cov)
        {
            volScalarField& T(_T());
            int stateSize = T.size();
            Eigen::VectorXd prior_mu = Eigen::MatrixXd::Ones(stateSize,1) * mean;
            Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(stateSize, stateSize) * cov;
            priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu, prior_cov);
        }

        void priorSampling() 
        {
            for(int i = 0; i < Nseeds; i++)
            {
                priorSamples.col(i) = priorDensity->Sample();
            }
            posteriorSamples = priorSamples;
        }


        void modelErrorSetup(double mean, double cov)
        {
            volScalarField& T(_T());
            int stateSize = T.size();
            Eigen::VectorXd modelError_mu = Eigen::MatrixXd::Ones(stateSize,1) * mean;
            Eigen::MatrixXd modelError_cov = Eigen::MatrixXd::Identity(stateSize, stateSize) * cov;
            modelErrorDensity = std::make_shared<muq::Modeling::Gaussian>(modelError_mu, modelError_cov);
        }

        void measNoiseSetup(double mean, double cov)
        {
            meas_cov = cov;
            int obsSize = measurements.rows();
            M_Assert(obsSize > 0, "Read measurements before setting up the noise");
            Eigen::VectorXd measNoise_mu = Eigen::MatrixXd::Ones(obsSize,1) * mean;
            Eigen::MatrixXd measNoise_cov = Eigen::MatrixXd::Identity(obsSize, obsSize) * cov;
            measNoiseDensity = std::make_shared<muq::Modeling::Gaussian>(measNoise_mu, measNoise_cov);
        }
};

#endif
