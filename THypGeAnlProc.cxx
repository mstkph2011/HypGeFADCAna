// $Id: THypGeAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "THypGeAnlProc.h"

#include "defines.h"						// defines (~globals) are in this file, add this line to every file
#include <algorithm>


#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"

#include "TGo4WinCond.h"
#include "TGo4Analysis.h"

#include "THypGeAnlEvent.h"
#include "THypGeUnpackEvent.h"
#include "THypGeParameter.h"

//-----------------------------------------------------------
THypGeAnlProc::THypGeAnlProc() :
   TGo4EventProcessor()
{
}
//-----------------------------------------------------------
THypGeAnlProc::THypGeAnlProc(const char* name) :
   TGo4EventProcessor(name)
{
	EventCounter = 0;
	
	char chis[100],chead[100];
	for (Int_t i =0;i < FADC_CHAN; i++)
	{
		//get trace histograms
		snprintf(chis,63,"V1724/Trace%02d",i+1);
		fhTrace[i] = (TH1D*) GetHistogram(chis);
		
		// create histograms for smoothed traces
		snprintf(chis,15,"Trace_smoothed_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after smooting; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_Smoothed[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_Smoothed[i]->GetXaxis()->CenterTitle();
		fhTrace_Smoothed[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_Smoothed[i],"V1724");
		
		//create histograms for baseline corrected traces
		sprintf(chis,"Trace_bc_%02d",i+1);  
		sprintf(chead,"Trace channel %2d after smoothing and baseline correction; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_BaseCorr[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_BaseCorr[i]->GetXaxis()->CenterTitle();
		fhTrace_BaseCorr[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_BaseCorr[i],"V1724");
		
		snprintf(chis,15,"Trace_deconv_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after deconvolution; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_deconv[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_deconv[i]->GetXaxis()->CenterTitle();
		fhTrace_deconv[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_deconv[i],"V1724");
		
		snprintf(chis,15,"Trace_MWD_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after MWD filter; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_MWD[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MWD[i]->GetXaxis()->CenterTitle();
		fhTrace_MWD[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MWD[i],"V1724");
		
		snprintf(chis,15,"Trace_MA_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after MA filter; time [#mus];Amplitude [a.u.] ",i+1);
		fhTrace_MA[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MA[i]->GetXaxis()->CenterTitle();
		fhTrace_MA[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MA[i],"V1724");
	}
	if (!fhTrace[0])
		cout << "fhTrace[0] not found"<< endl;
		
		// get histo for energy spectrum 
	fhEnergySpectrum = (TH1D*) GetHistogram("V1724/Energyspectrum/Energy");
	cout << fhEnergySpectrum << endl;
		// get histo for risetime and correlations
	fhRisetime1090 =  (TH1D*) GetHistogram("V1724/Risetime1090/Risetime1090");		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	fhRisetime3090 =  (TH1D*) GetHistogram("V1724/Risetime3090/Risetime3090");		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	fhEnergyRise1090Corr = (TH2D*) GetHistogram("V1724/EnergyRise1090Corr/EnergyRise1090Corr");
	fhEnergyRise3090Corr = (TH2D*) GetHistogram("V1724/EnergyRise3090Corr/EnergyRise3090Corr");
	
	fhEnergyTimeSinceLastPulse= (TH2D*) GetHistogram("V1724/EnergyTimeSinceLastPulse/EnergyTimeSinceLastPulse");
	
	fhAmplBaselinegradient= new TH1D ("Baselinegrad","Gradient of the start of the baseline in amplitude signal",1000,0,10);
	AddHistogram(fhAmplBaselinegradient,"V1724/BaseGrad");
		
	fhEnergySpectrum_withCut = (TH1D*) GetHistogram("V1724/Energyspectrum/Energy_withCut");
	
		// get parameters
	fHypPar = (THypGeParameter*)  GetParameter("HypGeParameter");
	cout << "fHypPar" <<fHypPar << endl;
		//cout << "\tMWDm " << fHypPar->GetMWDm() << endl;
		// real analysis object
	fMWDAna = new THypGeMWD(TRACE_LENGTH);
	fMWDAna->ConnectHistograms(fhEnergySpectrum_withCut);
	cout << "**** THypGeAnlProc: Create" << endl;
}
//-----------------------------------------------------------
THypGeAnlProc::~THypGeAnlProc()						// 25.3.14 something here gives an error when shuting down the analysis
{
	//delete[] fhTrace_Smoothed;
	//delete[] fhTrace_BaseCorr;
	//delete[] fhTrace_deconv;
	//delete[] fhTrace_MWD;
	//delete[] fhTrace_MA;
	
	//delete fhAmplBaselinegradient;
	//delete fMWDAna;
   cout << "**** THypGeAnlProc: Delete" << endl;

}
//-----------------------------------------------------------
Bool_t THypGeAnlProc::BuildEvent(TGo4EventElement* dest)
{

	Bool_t isValid=kFALSE; // validity of output event

   THypGeUnpackEvent* inp_evt  = (THypGeUnpackEvent*) GetInputEvent();
   THypGeAnlEvent* out_evt = (THypGeAnlEvent*) dest;

   // see comments in UnpackProc
   if((inp_evt==0) || !inp_evt->IsValid()){ // input invalid
	  out_evt->SetValid(isValid); // invalid
	  return isValid; // must be same is for SetValid
   }
   isValid=kTRUE;

	//add ana code here
	
	fMWDAna->SetUseMWD(1);
	if (fHypPar->GetParametersChanged())
		fMWDAna->SetParameters(fHypPar->GetMWDm(),fHypPar->GetMAl(),fHypPar->GetNoOfSmoothing(),fHypPar->GetWidth() ,fHypPar->GetSigmaGaus(),fHypPar->GetSigmaBil(),fHypPar->GetTau(), fHypPar->GetEnableMA(),fHypPar->GetSmoothingMethod(),fHypPar->GetEnableBaselineCorrection());

	if (fMWDAna->FullAnalysis(fhTrace[0],fhTrace_Smoothed[0],fhTrace_BaseCorr[0], fhTrace_deconv[0],fhTrace_MWD[0],fhEnergySpectrum,fhRisetime1090,fhRisetime3090,fhTrace_MA[0],(TH2D*) fhEnergyRise1090Corr,(TH2D*) fhEnergyRise3090Corr, fhEnergyTimeSinceLastPulse) != -1)				// some error here
		fhAmplBaselinegradient->Fill((fhTrace_deconv[0]->GetBinContent(1)-fhTrace_deconv[0]->GetBinContent(301))/300);
	
	//this shows number of real events
	
	//EventCounter++;
	//cout << "EventCounter " << EventCounter << endl;
	
	// see comments in UnpackProc
	out_evt->SetValid(isValid);
	return isValid;
}
