// $Id: THypGeAnalysis.cxx 754 2011-05-18 11:04:52Z adamczew $
//-----------------------------------------------------------------------
//			 The GSI Online Offline Object Oriented (Go4) Project
//				 Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//										 Planckstr. 1, 64291 Darmstadt, Germany
// Contact:						http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "THypGeAnalysis.h"
#include "defines.h"						// defines (~globals) are in this file, add this line to every file
#include <stdlib.h>
#include "Riostream.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TSystem.h"

extern "C" {
	 #include "s_filhe_swap.h"
	 #include "s_bufhe_swap.h"
	 #include "f_ut_utime.h"
}

#include "TGo4Fitter.h"
#include "TGo4FitterEnvelope.h"
#include "TGo4AnalysisStep.h"
#include "TGo4StepFactory.h"
#include "THypGeParameter.h"
#include "THypGeUnpackEvent.h"
#include "THypGeAnlEvent.h"
#include "TGo4Version.h"

#include "Go4EventServer.h"

#include "THypGeSpectrumAnalyser.h"




//***********************************************************
THypGeAnalysis::THypGeAnalysis() :
	 TGo4Analysis(),
	 fMbsEvent(0),
	 fRawEvent(0),
	 fCalEvent(0)
{
	cout << "Wrong constructor THypGeAnalysis()!" << endl;
}

//***********************************************************
// this constructor is called by go4analysis executable
THypGeAnalysis::THypGeAnalysis(int argc, char** argv) :
	 TGo4Analysis(argc, argv),
	 fMbsEvent(0),
	 fRawEvent(0),
	 fCalEvent(0),
	 fEvents(0),
	 fLastEvent(0)
{
	 if (!TGo4Version::CheckVersion(__GO4BUILDVERSION__)) {
			cout << "****	Go4 version mismatch" << endl;
			exit(-1);
	 }
	 cout << "Argc = " << argc << endl;
	 cout << "**** THypGeAnalysis: Create " << argv[0] << endl;
		// get parameter for analysis from user input
	 MWDm = 200;			// M of MWD
	 MAl = 100;				// L of MA
	 NoS = 100;				// Number of smoothings of mean and WA filter
	 Width = 3;				// Width of mean filter
	 sigmaGaus = 11;			// sigma of gaussian shaper
	 sigmaBil = 1500;			// sigma of gaussian shaper
	 tau = 5383;			// tau of deconvolution
	 EnableMA = 0;			// Switch for second moving average filter
	 SmoothingMethod = 4;	// Choose Smoothing Filter: 0 = Off, 1 = Mean, 2 = WA, 3 = Gaus, 4 = Bil
	 EnableBaselineCorrection = 1; 	//Switch baseline correction on or off
	 
	 if (argc>1)
	 {
		MWDm = atoi(argv[1]);
		cout << "MWDm\t"<<atoi(argv[1]) << "   " <<MWDm << endl;
	}
	if (argc>2)
	{
		MAl = atoi(argv[2]);
		cout << "MAl\t" << argv[2] << endl;
	}
	if (argc>3)
	{
		NoS = atoi(argv[3]);
		cout << "Number of Smoothing\t" << argv[3] << endl;
	}
	if (argc>4)
	{
		Width = atoi(argv[4]);
		cout << "Width\t" << argv[4] << endl;
	}
	if (argc>5)
	{
		sigmaGaus = atoi(argv[5]);
		cout << "sigmaGaus\t" << argv[5] << endl;
	}
	if (argc>6)
	{
		sigmaBil = atoi(argv[6]);
		cout << "sigmaBil\t" << argv[6] << endl;
	}
	if (argc>7)
	{
		tau = atoi(argv[7]);
		cout << "tau\t" << argv[7] << endl;
	}
	if (argc>8)
	{
		EnableMA = atoi(argv[8]);
		cout << "Enable MA\t" << argv[8] << endl;
	}
	if (argc>9)
	{
		SmoothingMethod = atoi(argv[9]);
		cout << "Enable Smoothing\t" << argv[9] << endl;
	}
	if (argc>10)
	{
		EnableBaselineCorrection = atoi(argv[10]);
		cout << "EnableBaselineCorrection\t" << argv[10] << endl;
	}
		//TGo4MbsFileParameter* input = new TGo4MbsFileParameter(userinput);
	 TString kind, input, out1, out2;

// Create step 1 Unpack.
	 TGo4StepFactory* factory1 = new TGo4StepFactory("UnpackFactory");
	 factory1->DefEventProcessor("UnpackProc", "THypGeUnpackProc");// object name, class name
	 factory1->DefOutputEvent("UnpackEvent", "THypGeUnpackEvent"); // object name, class name
	 TGo4AnalysisStep* step1 = new TGo4AnalysisStep("Unpack",factory1,0,0,0);
	 step1->SetErrorStopEnabled(kTRUE);
	 AddAnalysisStep(step1);
// These settings will be overwritten by setup.C
	 step1->SetSourceEnabled(kTRUE);
	 step1->SetStoreEnabled(kFALSE);
	 step1->SetProcessEnabled(kTRUE);

// Create step 2 Analysis.
	 TGo4StepFactory* factory2 = new TGo4StepFactory("AnalysisFactory");
	 factory2->DefInputEvent("UnpackEvent", "THypGeUnpackEvent"); // object name, class name
	 factory2->DefEventProcessor("AnlProc", "THypGeAnlProc"); // object name, class name
	 factory2->DefOutputEvent("AnlEvent", "THypGeAnlEvent"); // object name, class name
	 TGo4AnalysisStep* step2		= new TGo4AnalysisStep("Analysis",factory2,0,0,0);
	 step2->SetErrorStopEnabled(kTRUE);
	 AddAnalysisStep(step2);
// These settings will be overwritten by setup.C
	 step2->SetSourceEnabled(kFALSE);
	 step2->SetStoreEnabled(kFALSE);
	 step2->SetProcessEnabled(kTRUE);
	
	// Parameters of the analysis
	fPar = new THypGeParameter("HypGeParameter");
	 cout << "ParAdded " << AddParameter(fPar)  << endl;
	 fPar->SetParameters( MWDm, MAl,NoS, Width ,sigmaGaus,sigmaBil, tau, EnableMA, SmoothingMethod, EnableBaselineCorrection);

	
	char chis[100], chead[100];
	for(Int_t i=0;i<FADC_CHAN;i++)
		{
			snprintf(chis,15,"Trace%02d",i+1);	
			snprintf(chead,63,"Trace channel %2d; time [#mus];Amplitude [a.u.]",i+1);
			fhTrace[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			fhTrace[i]->GetXaxis()->CenterTitle();
			fhTrace[i]->GetYaxis()->CenterTitle();
			
			AddHistogram(fhTrace[i],"V1724");
		}
		
	//create histogram for energy spectrum
	fhEnergySpectrum = new TH1D("Energy","Energy",4000,0,4000);
	AddHistogram(fhEnergySpectrum,"V1724/Energyspectrum");

		//risetime histos
	fhRisetime1090 = new TH1D("Risetime1090","Risetime1090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	AddHistogram(fhRisetime1090,"V1724/Risetime1090");
	fhRisetime3090 = new TH1D("Risetime3090","Risetime3090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	AddHistogram(fhRisetime3090,"V1724/Risetime3090");
	fhEnergyRise1090Corr = new TH2D("EnergyRise1090Corr","Enegy-Risetime1090-Correlation;Rt;E",100,0,1000,2000,0,2000);
	AddHistogram(fhEnergyRise1090Corr,"V1724/EnergyRise1090Corr");
	fhEnergyRise3090Corr = new TH2D("EnergyRise3090Corr","Enegy-Risetime3090-Correlation;Rt;E",100,0,1000,2000,0,2000);
	AddHistogram(fhEnergyRise3090Corr,"V1724/EnergyRise3090Corr");
	
	fhEnergyTimeSinceLastPulse = new TH2D("EnergyTimeSinceLastPulse","Energy- Time since last Pulse - Correlation;E ;Time since last pulse [#mus]",20000,0,2000,200,0,200);
	AddHistogram(fhEnergyTimeSinceLastPulse,"V1724/EnergyTimeSinceLastPulse");
	
	cout << "All global histograms created"<< endl;
	
	
	//SetDynListInterval(1000);					// set the auto save interval. the value is the time between 2 automatic saves in seconds
	//SetDynListInterval (100);					// change this value to get more updates on the drawing of histograms. the value is the number of events processed before a TTree::Draw() is called. only in GUI mode
}
//***********************************************************
THypGeAnalysis::~THypGeAnalysis()
{
	cout << "**** THypGeAnalysis: Delete" << endl;
}
//***********************************************************

//-----------------------------------------------------------
Int_t THypGeAnalysis::UserPreLoop()
{
	cout << "**** THypGeAnalysis: PreLoop" << endl;
	Print(); // printout the step settings
	cout << "**************************************" << endl;
	 // we update the pointers to the current event structures here:
	 fMbsEvent = dynamic_cast<TGo4MbsEvent*>		(GetInputEvent("Unpack"));	 // of step "Unpack"
	 fRawEvent = dynamic_cast<THypGeUnpackEvent*> (GetOutputEvent("Unpack"));
	 fCalEvent = dynamic_cast<THypGeAnlEvent*>		(GetOutputEvent("Analysis"));
	 fEvents=0;
	 fLastEvent=0;

	 // create histogram for UserEventFunc
	 // At this point, the histogram has been restored
	 // from auto-save file if any.
	 fPar->SetParameters( MWDm, MAl,NoS, Width, sigmaGaus, sigmaBil, tau, EnableMA, SmoothingMethod, EnableBaselineCorrection);
	 fPar->PrintParameters();
	 return 0;
}
//-----------------------------------------------------------
Int_t THypGeAnalysis::UserPostLoop()	// fitting can be done here
{
	cout << "**** THypGeAnalysis: PostLoop" << endl;
	cout << "Last event: " << fLastEvent << " Total events: " << fEvents << endl;
	if(fMbsEvent)
		{
			// we can check some properties of last event here:
			//fMbsEvent->PrintEvent(); // header and data content

			// fileheader structure:
			s_filhe* fileheader=fMbsEvent->GetMbsSourceHeader();
			if(fileheader) {
					 cout <<"\nInput file was: "<<fileheader->filhe_file << endl;
					 cout <<"Tapelabel:\t" << fileheader->filhe_label<<endl;
					 cout <<"UserName:\t" << fileheader->filhe_user<<endl;
					 cout <<"RunID:\t" << fileheader->filhe_run<<endl;
					 cout <<"\tExplanation: "<<fileheader->filhe_exp <<endl;
					 cout <<"\tComments: "<<endl;
					 Int_t numlines=fileheader->filhe_lines;
					 for(Int_t i=0; i<numlines;++i)
							 cout<<"\t\t"<<fileheader->s_strings[i].string << endl;
				 }

			// mbs buffer header structure:
			s_bufhe* bufheader=fMbsEvent->GetMbsBufferHeader();
			if(bufheader) {
				 char sbuf[1000];
				 f_ut_utime(bufheader->l_time[0], bufheader->l_time[1], sbuf);
				 cout <<"Last Buffer:"<<endl;
				 cout <<"\tNumber: " << bufheader->l_buf << endl;
				 cout <<"\tTime: " << sbuf << endl;
			 }


		}

	fMbsEvent = 0; // reset to avoid invalid pointer if analysis is changed in between
	fRawEvent = 0;
	fCalEvent = 0;
	fEvents=0;

	// ana of energyspectrum

	//  (12.02.2014) fitting is done afterwards
	
	
	//THypGeSpectrumAnalyser* SpecAna = new THypGeSpectrumAnalyser(fhEnergySpectrum,"co60",10);			// 2 is placeholder
		//SpecAna->SetLinearCalibration(1);
		////SpecAna->SetTxtFileOutputName(TxtFilename);
		//SpecAna->SetTxtFileOutputName("test.txt");
		////SpecAna->SetRootFileOutputName(RootFilename);
		//SpecAna->SetRootFileOutputName("test.root");
		//SpecAna->SetGaussianFitting();
		//SpecAna->EnableGo4Mode();
		//SpecAna->SetAnalysisObject(this);
		//SpecAna->AnalyseSpectrum();
		//fhEnergySpectrum->Draw();
	return 0;
}

//-----------------------------------------------------------
Int_t THypGeAnalysis::UserEventFunc()
{
//// This function is called once for each event.
	
	 return 0;
}
