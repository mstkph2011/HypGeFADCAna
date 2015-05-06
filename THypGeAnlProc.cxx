// $Id: THypGeAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
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

#define FADC_CHAN 1
// value from defines.h overwritten !!!!!

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
		snprintf(chis,63,"Traces/Trace_%02d",i+1);
		fhTrace[i] = (TH1D*) GetHistogram(chis);
		
		// create histograms for smoothed traces
		snprintf(chis,30,"Trace_Deri_%02d",i+1);
		snprintf(chead,63,"Derivative of trace channel %2d ; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_Deri[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_Deri[i]->GetXaxis()->CenterTitle();
		fhTrace_Deri[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_Deri[i],"Traces");

		// create histograms for smoothed traces
		snprintf(chis,30,"Trace_smoothed_%02d",i+1);
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
		
		snprintf(chis,30,"Trace_deconv_%02d",i+1);
		snprintf(chead,63,"Trace channel %2d after deconvolution; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_deconv[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_deconv[i]->GetXaxis()->CenterTitle();
		fhTrace_deconv[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_deconv[i],"Traces");
		
		snprintf(chis,30,"Trace_MWD_%02d",i+1);
		snprintf(chead,63,"Trace channel %2d after MWD filter; time [#mus];Amplitude [a.u.]",i+1);
		fhTrace_MWD[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MWD[i]->GetXaxis()->CenterTitle();
		fhTrace_MWD[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MWD[i],"Traces");
		
		snprintf(chis,30,"Trace_MA_%02d",i+1);
		snprintf(chead,63,"Trace channel %2d after MA filter; time [#mus];Amplitude [a.u.] ",i+1);
		fhTrace_MA[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_MA[i]->GetXaxis()->CenterTitle();
		fhTrace_MA[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_MA[i],"Traces");
			
		snprintf(chis,30,"Trace_Direct_%02d",i+1);
		snprintf(chead,63,"Trace channel %2d after Direct filter; time [#mus];Amplitude [a.u.] ",i+1);
		fhTrace_Direct[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
		fhTrace_Direct[i]->GetXaxis()->CenterTitle();
		fhTrace_Direct[i]->GetYaxis()->CenterTitle();
			AddHistogram(fhTrace_Direct[i],"Traces");

			//create histogram for energy spectrum from MA
		snprintf(chis,30,"Energy_%02d",i+1);
		snprintf(chead,63,"Energy spectrum from MA signal channel %2d; ADC channel [a.u.];Counts [a.u.] ",i+1);
		fhEnergySpectrum[i] = new TH1D(chis,chead,16000,0,8000);
			AddHistogram(fhEnergySpectrum[i],"Energyspectrum");
		//create histogram for energy spectrum
		snprintf(chis,30,"EnergyCorr_%02d",i+1);
		snprintf(chead,63,"Energy spectrum from MA signal channel with risetime correction %2d; ADC channel [a.u.];Counts [a.u.] ",i+1);
		fhEnergySpectrumCorr[i] = new TH1D(chis,chead,16000,0,8000);
			AddHistogram(fhEnergySpectrumCorr[i],"Energyspectrum");
		//create histograms for energy spectrum filtering on first or other pulses in the trace
		snprintf(chis,30,"EnergyFirstPulse_%02d",i+1);
		snprintf(chead,63,"Energy spectrum from MA signal of first pulse in trace channel %2d; ADC channel [a.u.];Counts [a.u.] ",i+1);
		fhEnergySpectrumFirstPulse[i] = new TH1D(chis,chead,16000,0,8000);
			AddHistogram(fhEnergySpectrumFirstPulse[i],"Energyspectrum");
		snprintf(chis,30,"EnergyOtherPulses_%02d",i+1);
		snprintf(chead,63,"Energy spectrum from MA signal of all but first pulses in trace channel %2d; ADC channel [a.u.];Counts [a.u.] ",i+1);
		fhEnergySpectrumOtherPulses[i] = new TH1D(chis,chead,16000,0,8000);
			AddHistogram(fhEnergySpectrumOtherPulses[i],"Energyspectrum");
		//create histogram for energy spectrum with a cut in the pile up time
		snprintf(chis,30,"EnergyWithCut_%02d",i+1);
		snprintf(chead,63,"Energy spectrum with cut channel %2d; ADC channel [a.u.];Counts [a.u.]",i+1);
		fhEnergySpectrum_withCut[i] = new TH1D(chis,chead,4000,0,4000);
			AddHistogram(fhEnergySpectrum_withCut[i],"Energyspectrum");

		//risetime histos
		snprintf(chis,30,"Risetime1090_%02d",i+1);
		snprintf(chead,63,"Risetime1090 %2d; Risetime 1090 [ns];Counts [a.u.]",i+1);
		fhRisetime1090[i] = new TH1D(chis,chead,100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRisetime1090[i],"Risetime1090");
		snprintf(chis,30,"Risetime1090FirstPulse_%02d",i+1);
		snprintf(chead,63,"Risetime1090 of first pulse in trace %2d; Risetime 1090 [ns];Counts [a.u.]",i+1);
		fhRisetime1090Co1332FirstPulse[i] = new TH1D(chis,chead,100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRisetime1090Co1332FirstPulse[i],"Risetime1090Co1332Only");
		snprintf(chis,30,"Risetime1090OtherPulses_%02d",i+1);
		snprintf(chead,63,"Risetime1090 of all but first pulse in trace %2d; Risetime 1090 [ns];Counts [a.u.]",i+1);
		fhRisetime1090Co1332OtherPulses[i] = new TH1D(chis,chead,100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRisetime1090Co1332OtherPulses[i],"Risetime1090Co1332Only");
		snprintf(chis,30,"Risetime3090_%02d",i+1);
		snprintf(chead,63,"Risetime3090 %2d; Risetime 1090 [ns];Counts [a.u.]",i+1);
		fhRisetime3090[i] = new TH1D(chis,chead,100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRisetime3090[i],"Risetime3090");
		snprintf(chis,30,"Risetime1090Co1332Only_%02d",i+1);
		snprintf(chead,63,"Risetime1090Co1332Only %2d; Risetime 1090 [ns];Counts [a.u.]",i+1);
		fhRisetime1090Co1332Only[i] = new TH1D(chis,chead,100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRisetime1090Co1332Only[i],"Risetime1090Co1332Only");



		// risetime, energy correlation histograms
		snprintf(chis,30,"EnergyRt1090_%02d",i+1);
		snprintf(chead,63,"Energy-Risetime1090-Correlation channel %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRt1090[i] = new TH2D(chis,chead,100,0,1000,4000,0,4000);
			AddHistogram(fhEnergyRt1090[i],"EnergyRt1090");
		snprintf(chis,30,"EnergyRt1090CorrectionRt_%02d",i+1);
		snprintf(chead,63,"Energy-Risetime1090-Correlation with Rt correction channel %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRt1090CorrectionRt[i] = new TH2D(chis,chead,100,0,1000,4000,0,4000);
			AddHistogram(fhEnergyRt1090CorrectionRt[i],"EnergyRt1090");
		snprintf(chis,30,"EnergyRt1090Co1332Only_%02d",i+1);
		snprintf(chead,63,"Energy-Risetime1090-Correlation of Co 1332 keV line of channel %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRt1090Co1332Only[i] = new TH2D(chis,chead,100,0,1000,4000,0,4000);
			AddHistogram(fhEnergyRt1090Co1332Only[i],"EnergyRt1090");


		// histograms using the Rratio
		snprintf(chis,30,"RratioCo1332Only_%02d",i+1);
		snprintf(chead,63,"RratioCo1332Only %2d; Rratio;Counts [a.u.]",i+1);
		fhRratioCo1332Only[i] = new TH1D(chis,chead,50,0,1);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
			AddHistogram(fhRratioCo1332Only[i],"RratioHistograms");
		snprintf(chis,30,"EnergyRratio_%02d",i+1);
		snprintf(chead,63,"Energy-Rratio-Correlation channel %02d;Rratio;ADC Value [a.u.]",i+1);
		fhEnergyRratio[i] = new TH2D(chis,chead,50,0,1,4000,0,4000);
			AddHistogram(fhEnergyRratio[i],"RratioHistograms");
		snprintf(chis,30,"EnergyRratioCorrectionRratio_%02d",i+1);
		snprintf(chead,63,"Energy-Rratio-Correlation with Rratio correction channel %02d;Rratio;ADC Value [a.u.]",i+1);
		fhEnergyRratioCorrectionR[i] = new TH2D(chis,chead,50,0,1,4000,0,4000);
			AddHistogram(fhEnergyRratioCorrectionR[i],"RratioHistograms");
		snprintf(chis,30,"EnergyRratioCo1332Only_%02d",i+1);
		snprintf(chead,63,"Energy-Risetime1090-Correlation of Co 1332 keV line of channel %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRratioCo1332Only[i] = new TH2D(chis,chead,50,0,1,4000,0,4000);
			AddHistogram(fhEnergyRratioCo1332Only[i],"RratioHistograms");


		snprintf(chis,30,"Rt1030Rt1090Co1332Only_%02d",i+1);
		snprintf(chead,100,"Risetime1030-Risetime1090-Correlation of Co 1332 keV line of channel %02d;Risetime 1030 [ns];Risetime 1090 [ns]",i+1);
		fhRt1030Rt1090Co1332Only[i] = new TH2D(chis,chead,100,0,1000,100,0,1000);
			AddHistogram(fhRt1030Rt1090Co1332Only[i],"RtRtCorrelation");
		snprintf(chis,30,"Rt1030Rt80100Co1332Only_%02d",i+1);
		snprintf(chead,100,"Risetime1030-Risetime80100-Correlation of Co 1332 keV line of channel %02d;Risetime 1030 [ns];Risetime 80100 [ns]",i+1);
		fhRt1030Rt80100Co1332Only[i] = new TH2D(chis,chead,100,0,1000,100,0,1000);
			AddHistogram(fhRt1030Rt80100Co1332Only[i],"RtRtCorrelation");


			// rt, energy correlation histogram with ballistic deficit correction
		snprintf(chis,40,"EnergyRt1090_ballistic_%02d",i+1);
		snprintf(chead,100,"Energy-Risetime1090-Correlation channel without ballistic %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRt1090_ballistic[i] = new TH2D(chis,chead,1001,0,1000,4000,0,4000);
			AddHistogram(fhEnergyRt1090_ballistic[i],"EnergyRt1090");

		snprintf(chis,30,"EnergyRt3090_%02d",i+1);
		snprintf(chead,100,"Energy-Risetime3090-Correlation channel %02d;Risetime [ns];ADC Value [a.u.]",i+1);
		fhEnergyRt3090[i] = new TH2D(chis,chead,1001,0,1000,4000,0,4000);
			AddHistogram(fhEnergyRt3090[i],"EnergyRt3090");


		// histogram to see the correlation of energy and the time between two pulses to examine the effect of pile up
		snprintf(chis,30,"EnergyTimeSinceLastPulse_%02d",i+1);
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation channel %02d;Time since last pulse [#mus];ADC Value [a.u.]",i+1);
		fhEnergyTimeSinceLastPulse[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
			AddHistogram(fhEnergyTimeSinceLastPulse[i],"EnergyTimeSinceLastPulse");
		snprintf(chis,50,"EnergyTimeSinceLastPulseCorrected_%02d",i+1);
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation with correction channel %02d;Time since last pulse [#mus];ADC Value [a.u.]",i+1);
		fhEnergyTimeSinceLastPulseCorr[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
			AddHistogram(fhEnergyTimeSinceLastPulseCorr[i],"EnergyTimeSinceLastPulseCorrected");


		snprintf(chis,50,"TraceDeriMaximumPosition_%02d",i+1);
		snprintf(chead,100,"TraceDeriMaximumPosition %02d;TraceDeriMaximumPosition [#mus];counts [a.u.]",i+1);
		fhTraceDeriMaximumPosition[i]  =new TH1D (chis,chead,100,0,100);
			AddHistogram(fhTraceDeriMaximumPosition[i],"TraceDeriMaxPos");


	}
	for (int i = 0; i < 20; i++)
	{
		snprintf(chis,100,"EnergyTimeSinceLastPulse_withCut_%02d",i+1);  
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation with cut on ADC value of previous signal (%d - %d);Time since last pulse [#mus];ADC Value [a.u.]",i*100,(i+1)*100-1);
		fhEnergyTimeSinceLastPulse_withCut[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
			AddHistogram(fhEnergyTimeSinceLastPulse_withCut[i],"EnergyTimeSinceLastPulse");
		//cout << fhEnergyTimeSinceLastPulse_withCut[i] << endl;
	}

	for (int i = 0; i < 20; i++)
	{
		snprintf(chis,100,"EnergyTimeSinceLastPulseCorrected_withCut_%02d",i+1);  
		snprintf(chead,100,"Energy- Time since last Pulse - Correlation with correction and cut on ADC value of previous signal (%d - %d);Time since last pulse [#mus];ADC Value [a.u.]",i*100,(i+1)*100-1);
		fhEnergyTimeSinceLastPulseCorr_withCut[i] = new TH2D(chis,chead,200,0,200,20000,0,2000);
		AddHistogram(fhEnergyTimeSinceLastPulseCorr_withCut[i],"EnergyTimeSinceLastPulseCorrected");
		//cout << fhEnergyTimeSinceLastPulse_withCut[i] << endl;
	}
		

		fhTrace_deri1 = new TH1D ("deri1","deri1",TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			AddHistogram(fhTrace_deri1,"Test");
		fhTrace_deri2 = new TH1D ("deri2","deri2",TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			AddHistogram(fhTrace_deri2,"Test");
		fhTrace_deri3 = new TH1D ("deri3","deri3",TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			AddHistogram(fhTrace_deri3,"Test");
		fhTrace_deri4 = new TH1D ("deri4","deri4",TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			AddHistogram(fhTrace_deri4,"Test");

	for (Int_t i =0;i < 1000; i++)
	{
		snprintf(chis,40,"OutPutTrace_%i",i);
		snprintf(chead,100,"Trace %i;Samples [a.u.];ADC Value [a.u.]",i);
		hTraceOutput[i] = new TH1D(chis,chead,TRACE_LENGTH,0,TRACE_LENGTH * TIME_RESOLUTION_FACTOR);
			AddHistogram(hTraceOutput[i],"OutputTraces");
	}

		// get parameters
	fHypPar = (THypGeParameter*)  GetParameter("HypGeParameter");
	


		//cout << "\tMWDm " << fHypPar->GetMWDm() << endl;
		// real analysis object
	fMWDAna = new THypGeMWD(TRACE_LENGTH,FADC_CHAN);
		fMWDAna->ConnectTraceHistograms(fhTrace, fhTrace_Deri, fhTrace_Smoothed, fhTrace_BaseCorr, fhTrace_deconv, fhTrace_MWD, fhTrace_MA, fhTrace_Direct);
		fMWDAna->Connect1DEnergySpectraHistograms(fhEnergySpectrum, fhEnergySpectrumCorr, fhEnergySpectrum_withCut);
		fMWDAna->Connect1DRisetimeHistograms(fhRisetime1090, fhRisetime3090,fhRisetime1090Co1332Only);
		fMWDAna->Connect2DEnergyRisetimeHistograms(fhEnergyRt1090,fhEnergyRt1090CorrectionRt, fhEnergyRt1090Co1332Only, fhEnergyRt3090,fhEnergyRt1090_ballistic);
		fMWDAna->ConnectRratioHistograms(fhRratioCo1332Only, fhEnergyRratio, fhEnergyRratioCorrectionR, fhEnergyRratioCo1332Only);
		fMWDAna->Connect2DEnergyTimeSinceLastPulseHistograms(fhEnergyTimeSinceLastPulse, fhEnergyTimeSinceLastPulseCorr, fhEnergyTimeSinceLastPulse_withCut, fhEnergyTimeSinceLastPulseCorr_withCut, 20);
		fMWDAna->ConnectTestHistograms(fhTrace_deri1,fhTrace_deri2,fhTrace_deri3,fhTrace_deri4);
		fMWDAna->ConnectOutputTraces(hTraceOutput);
		fMWDAna->ConnectPulseFilteredHistograms(fhEnergySpectrumFirstPulse, fhEnergySpectrumOtherPulses, fhRisetime1090Co1332FirstPulse, fhRisetime1090Co1332OtherPulses);
		fMWDAna->ConnectRtCorrelationHistograms(fhRt1030Rt1090Co1332Only,fhRt1030Rt80100Co1332Only);
		fMWDAna->ConnectTraceDeriMaximumHistograms(fhTraceDeriMaximumPosition);

		fMWDAna->SetParameters(fHypPar->GetMWDm(),fHypPar->GetMAl(),fHypPar->GetNoOfSmoothing(),fHypPar->GetWidth() ,fHypPar->GetSigmaGaus(),fHypPar->GetSigmaBil(),fHypPar->GetTau(), fHypPar->GetEnableMA(),fHypPar->GetSmoothingMethod(),fHypPar->GetEnableBaselineCorrection());
		if (fHypPar->GetSecondAnalysisRound())
		{
			fMWDAna->IsSecondRun();
			fMWDAna->SetSecondRunParametersFileName(fHypPar->GetParameterFileName());
		}
		fMWDAna->Init();
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
	
	//delete fhTraceDeriMaximumPosition;
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
	//cout << fHypPar->GetMAl() << endl;
	fMWDAna->SetUseMWD(1);
	//if (fHypPar->GetParametersChanged())
	//{
		fMWDAna->SetParameters(fHypPar->GetMWDm(),fHypPar->GetMAl(),fHypPar->GetNoOfSmoothing(),fHypPar->GetWidth() ,fHypPar->GetSigmaGaus(),fHypPar->GetSigmaBil(),fHypPar->GetTau(), fHypPar->GetEnableMA(),fHypPar->GetSmoothingMethod(),fHypPar->GetEnableBaselineCorrection());
		//fHypPar->SetParametersChanged(0);
	//}
	fMWDAna->FullAnalysis();
	//if (fMWDAna->FullAnalysis() != -1)				// some error here
		//fhTraceDeriMaximumPosition->Fill((fhTrace_deconv[0]->GetBinContent(1)-fhTrace_deconv[0]->GetBinContent(301))/300);
	
	//this shows number of real events
	
	//EventCounter++;
	//cout << "EventCounter " << EventCounter << endl;
	
	// see comments in UnpackProc
	out_evt->SetValid(isValid);
	return isValid;
}
