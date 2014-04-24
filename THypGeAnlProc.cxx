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
		snprintf(chis,63,"Traces/Trace%02d",i+1);
		fhTrace[i] = (TH1D*) GetHistogram(chis);
		
		// create histograms for smoothed traces
		snprintf(chis,15,"Trace_smoothed_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after smooting; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_Smoothed[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_Smoothed[i]->GetXaxis()->CenterTitle();
		fhTrace_Smoothed[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_Smoothed[i],"Traces");
		
		//create histograms for baseline corrected traces
		sprintf(chis,"Trace_bc_%02d",i+1);  
		sprintf(chead,"Trace channel %2d after smoothing and baseline correction; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_BaseCorr[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_BaseCorr[i]->GetXaxis()->CenterTitle();
		fhTrace_BaseCorr[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_BaseCorr[i],"Traces");
		
		snprintf(chis,15,"Trace_deconv_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after deconvolution; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_deconv[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_deconv[i]->GetXaxis()->CenterTitle();
		fhTrace_deconv[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_deconv[i],"Traces");
		
		snprintf(chis,15,"Trace_MWD_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after MWD filter; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_MWD[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MWD[i]->GetXaxis()->CenterTitle();
		fhTrace_MWD[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MWD[i],"Traces");
		
		snprintf(chis,15,"Trace_MA_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after MA filter; time [#mus];Amplitude [a.u.] ",i+1);
		fhTrace_MA[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MA[i]->GetXaxis()->CenterTitle();
		fhTrace_MA[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MA[i],"Traces");
			
		snprintf(chis,15,"Trace_Direct_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after Direct filter; time [#mus];Amplitude [a.u.] ",i+1);
		fhTrace_Direct[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_Direct[i]->GetXaxis()->CenterTitle();
		fhTrace_Direct[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_Direct[i],"Traces");
	}
	
	//create histogram for energy spectrum
	fhEnergySpectrum = new TH1D("Energy","Energy",4000,0,4000);
		AddHistogram(fhEnergySpectrum,"Energyspectrum");
	//create histogram for energy spectrum with a cut in the pile up time
	fhEnergySpectrum_withCut = new TH1D("Energy_withCut","Energy_withCut",4000,0,4000);
		AddHistogram(fhEnergySpectrum_withCut,"Energyspectrum");
	//risetime histos
	fhRisetime1090 = new TH1D("Risetime1090","Risetime1090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
		AddHistogram(fhRisetime1090,"Risetime1090");
	fhRisetime3090 = new TH1D("Risetime3090","Risetime3090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
		AddHistogram(fhRisetime3090,"Risetime3090");
	fhEnergyRise1090Corr = new TH2D("EnergyRise1090Corr","Energy-Risetime1090-Correlation;Rt;ADC Value [a.u.]",100,0,1000,2000,0,2000);
		AddHistogram(fhEnergyRise1090Corr,"EnergyRise1090Corr");
	
	fhEnergyRise3090Corr = new TH2D("EnergyRise3090Corr","Energy-Risetime3090-Correlation;Rt;ADC Value [a.u.]",100,0,1000,2000,0,2000);
		AddHistogram(fhEnergyRise3090Corr,"EnergyRise3090Corr");
	
	
	// histogram to see the correlation of energy and the time between two pulses to examine the effect of pile up
	fhEnergyTimeSinceLastPulse = new TH2D("EnergyTimeSinceLastPulse","Energy- Time since last Pulse - Correlation;Time since last pulse [#mus];ADC Value [a.u.]",200,0,200,20000,0,2000);
		AddHistogram(fhEnergyTimeSinceLastPulse,"EnergyTimeSinceLastPulse");
	for (int i = 0; i < 20; i++)
	{
		snprintf(chis,100,"EnergyTimeSinceLastPulse_withCut_%02d",i+1);  
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation with cut on ADC value of previous signal (%d - %d);Time since last pulse [#mus];ADC Value [a.u.]",i*100,(i+1)*100-1);
		fhEnergyTimeSinceLastPulse_withCut[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
			AddHistogram(fhEnergyTimeSinceLastPulse_withCut[i],"EnergyTimeSinceLastPulse");
		//cout << fhEnergyTimeSinceLastPulse_withCut[i] << endl;
	}
	fhEnergyTimeSinceLastPulseCorr = new TH2D("EnergyTimeSinceLastPulseCorrected","Energy- Time since last Pulse - Correlation with correction;Time since last pulse [#mus];ADC Value [a.u.]",200,0,200,20000,0,2000);
		AddHistogram(fhEnergyTimeSinceLastPulseCorr,"EnergyTimeSinceLastPulseCorrected");
	for (int i = 0; i < 20; i++)
	{
		snprintf(chis,100,"EnergyTimeSinceLastPulseCorrected_withCut_%02d",i+1);  
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation with correction and cut on ADC value of previous signal (%d - %d);Time since last pulse [#mus];ADC Value [a.u.]",i*100,(i+1)*100-1);
		fhEnergyTimeSinceLastPulseCorr_withCut[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
		AddHistogram(fhEnergyTimeSinceLastPulseCorr_withCut[i],"EnergyTimeSinceLastPulseCorrected");
		//cout << fhEnergyTimeSinceLastPulse_withCut[i] << endl;
	}
		
	fhAmplBaselinegradient= new TH1D ("Baselinegrad","Gradient of the start of the baseline in amplitude signal",1000,0,10);
		AddHistogram(fhAmplBaselinegradient,"BaseGrad");
		
	
	
		// get parameters
	fHypPar = (THypGeParameter*)  GetParameter("HypGeParameter");
	
		//cout << "\tMWDm " << fHypPar->GetMWDm() << endl;
		// real analysis object
	fMWDAna = new THypGeMWD(TRACE_LENGTH);
		fMWDAna->ConnectTraceHistograms(fhTrace, fhTrace_Smoothed, fhTrace_BaseCorr, fhTrace_deconv, fhTrace_MWD, fhTrace_MA, fhTrace_Direct);
		fMWDAna->Connect1DEnergySpectraHistograms(fhEnergySpectrum,fhEnergySpectrum_withCut);
		fMWDAna->Connect1DRisetimeHistograms(fhRisetime1090, fhRisetime3090);
		fMWDAna->Connect2DEnergyRisetimeHistograms(fhEnergyRise1090Corr, fhEnergyRise3090Corr);
		fMWDAna->Connect2DEnergyTimeSinceLastPulseHistograms(fhEnergyTimeSinceLastPulse, fhEnergyTimeSinceLastPulseCorr, fhEnergyTimeSinceLastPulse_withCut, fhEnergyTimeSinceLastPulseCorr_withCut, 20);
	
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

	//if (fMWDAna->FullAnalysis(fhTrace[0],fhTrace_Smoothed[0],fhTrace_BaseCorr[0], fhTrace_deconv[0],fhTrace_MWD[0],fhEnergySpectrum,fhRisetime1090,fhRisetime3090,fhTrace_MA[0],(TH2D*) fhEnergyRise1090Corr,(TH2D*) fhEnergyRise3090Corr, fhEnergyTimeSinceLastPulse) != -1)				// some error here
	if (fMWDAna->FullAnalysis() != -1)				// some error here
		fhAmplBaselinegradient->Fill((fhTrace_deconv[0]->GetBinContent(1)-fhTrace_deconv[0]->GetBinContent(301))/300);
	
	//this shows number of real events
	
	//EventCounter++;
	//cout << "EventCounter " << EventCounter << endl;
	
	// see comments in UnpackProc
	out_evt->SetValid(isValid);
	return isValid;
}
