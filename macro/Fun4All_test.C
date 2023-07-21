// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLTEST_C
#define MACRO_FUN4ALLTEST_C

#include </sphenix/user/shuhangli/etacalib/analysis/pi0ClusterAna/src/pi0ClusterAna.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <caloreco/RawClusterBuilderGraph.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/RawTowerCalibration.h>
#include <g4centrality/PHG4CentralityReco.h>
#include </sphenix/user/shuhangli/etacalib/source/GetCaloInfo.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libGetCaloInfo.so)
R__LOAD_LIBRARY(libpi0ClusterAna.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libqa_modules.so)
R__LOAD_LIBRARY(libg4centrality.so)



void Fun4All_test(
	const int nEvents = 1,
	const string &trackFile = "dst_tracks.list",
	const string &clusterFile = "dst_calo_cluster.list",
	const string &truthFile = "dst_truth.list",
	const string &globalFile = "dst_bbc_g4hit.list",
	const string &outputFile = "cemcCluster_Out.root")
{
	// this convenience library knows all our i/o objects so you don't
	// have to figure out what is in each dst type
	gSystem->Load("libg4dst.so");

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(1); // set it to 1 if you want event printouts
	PHG4CentralityReco *cent = new PHG4CentralityReco();
	cent->Verbosity(0);
	cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
	se->registerSubsystem(cent);

	Fun4AllInputManager *inCluster = new Fun4AllDstInputManager("DSTClusters");
	std::cout << "Adding file list " << clusterFile << std::endl;
	inCluster->AddListFile(clusterFile, 1);
	se->registerInputManager(inCluster);

	Fun4AllInputManager *truthCalo = new Fun4AllDstInputManager("DSTCaloTruth");
	std::cout << "Adding file list " << truthFile << std::endl;
	truthCalo->AddListFile(truthFile, 1);
	se->registerInputManager(truthCalo);

	Fun4AllInputManager *Global = new Fun4AllDstInputManager("DSTGlobal");
	std::cout << "Adding file list " << globalFile << std::endl;
	Global->AddListFile(globalFile, 1);
	se->registerInputManager(Global);

	Fun4AllInputManager *inTrack = new Fun4AllDstInputManager("DSTTracks");
	std::cout << "Adding file list " << trackFile << std::endl;
	inTrack->AddListFile(trackFile, 1);
	se->registerInputManager(inTrack);

	/*
	RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
	ClusterBuilder->Detector("CEMC");
	ClusterBuilder->Verbosity(1);
	ClusterBuilder->set_threshold_energy(0.030);  // This threshold should be the same as in CEMCprof_Thresh**.root file below
	std::string emc_prof = getenv("CALIBRATIONROOT");
	emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
	ClusterBuilder->LoadProfile(emc_prof);
	//    ClusterBuilder->set_UseTowerInfo(1); // to use towerinfo objects rather than old RawTower
	se->registerSubsystem(ClusterBuilder);

	RawClusterPositionCorrection *clusterCorrection = new RawClusterPositionCorrection("CEMC");
	clusterCorrection->Get_eclus_CalibrationParameters().ReadFromFile("CEMC_RECALIB", "xml", 0, 0,
			//raw location
			string(getenv("CALIBRATIONROOT")) + string("/CEMC/PositionRecalibration_EMCal_9deg_tilt/"));
	clusterCorrection->Get_ecore_CalibrationParameters().ReadFromFile("CEMC_ECORE_RECALIB", "xml", 0, 0,
			//raw location
			string(getenv("CALIBRATIONROOT")) + string("/CEMC/PositionRecalibration_EMCal_9deg_tilt/"));

	clusterCorrection->Verbosity(1);
	se->registerSubsystem(clusterCorrection);
	*/
	GetCaloInfo *getcalo = new GetCaloInfo("dummy", outputFile);
	se->registerSubsystem(getcalo);

	se->run(nEvents);
	se->End();

	delete se;
	cout << "Analysis Completed" << endl;

	gSystem->Exit(0);
	
	
	
}

#endif // MACRO_FUN4ALLTEST_C
