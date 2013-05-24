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
		snprintf(chead,63,"Trace channel %2d after smooting",i+1);
		fhTrace_Smoothed[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH);
		AddHistogram(fhTrace_Smoothed[i],"V1724/Trace_smoothed");
		
		//create histograms for baseline corrected traces
		snprintf(chis,15,"Trace_bc_%02d",i+1);  
		snprintf(chead,63,"Trace channel %2d after smoothing and baseline correction",i+1);
		fhTrace_BaseCorr[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH);
		AddHistogram(fhTrace_BaseCorr[i],"V1724/Trace_baselineCorrected");
	}
	if (!fhTrace[0])
		cout << "fhTrace[0] not found"<< endl;
		
	//add elements of this ana step here
		
		//MWD histos
	fhBaseline = new TH1D ("Baseline","Baseline",10000,0,10000);		//Kai edit
	AddHistogram(fhBaseline,"V1724/Baseline");
	fhAmplitude = new TH1D ("Amplitude","A",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fhAmplitude,"V1724/Ampl");
	fhMWD = new TH1D ("MWD","MWD",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fhMWD,"V1724/MWD");
	fhMWD_MA = new TH1D("MWDMA","MWDMA",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fhMWD_MA,"V1724/MWDMA");

		// get histo for energy spectrum
	fhEnergySpectrum = (TH1D*) GetHistogram("V1724/Energyspectrum/Energy");
	cout << "1090 " << fhEnergySpectrum << endl;
		// get risetime histos
	fhRisetime1090 = (TH1D*) GetHistogram("V1724/Rt1090/Rt1090");
	fhRisetime3090 = (TH1D*) GetHistogram("V1724/Rt3090/Rt3090");
	
		//get energy risetime correlation histos
		
	fhEnergyRise1090Corr = (TH2D*) GetHistogram("V1724/ERt1090Corr/ERt1090Corr");
	cout << "1090 " << fhEnergyRise1090Corr << endl;
	fhEnergyRise3090Corr = (TH2D*) GetHistogram("V1724/ERt3090Corr/ERt3090Corr");
		
		// correlation spectrum
		// get parameters
	fParam1 = (THypGeParameter*)  GetParameter("HypGeParameter");
		// real analysis object
	fMWDAna = new THypGeMWD(TRACE_LENGTH);
	
	//oldMax is used to prevent the analysis of a trace if the trace does not change (as in the actual (bad) data 2.4.14
	oldMax=0;
	cout << "**** THypGeAnlProc: Create" << endl;

}
//-----------------------------------------------------------
THypGeAnlProc::~THypGeAnlProc()
{
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
	
	//for(Int_t i = 1; i< TRACE_LENGTH;i++)
	//	fhTrace_BaseCorr[0]->SetBinContent(i,fhTrace[0]->GetBinContent(i));
	
	
	//cout << "M " << fParam1->M << endl;
	//cout << fhTrace[0]->GetMaximum() << endl;
	
	//if (oldMax !=fhTrace[0]->GetMaximum())	// change from here 21.5.13; maybe adapt histograms
	//{
		//oldMax = fhTrace[0]->GetMaximum();
		//fMWDAna->Test(fhTrace[0],fhTrace_Smoothed[0],fhTrace_BaseCorr[0], fhAmplitude, fhMWD, fhEnergySpectrum, fhRisetime1090,fhRisetime3090,fhMWD_MA); 
		fMWDAna->FullAnalysis(fhTrace[0],fhTrace_Smoothed[0],fhTrace_BaseCorr[0], fhAmplitude, fhMWD, fhEnergySpectrum, fhRisetime1090,fhRisetime3090,fhMWD_MA,(TH2D*) fhEnergyRise1090Corr,(TH2D*) fhEnergyRise3090Corr); 
	
		//this shows number of real events
		//EventCounter++;
		//cout << "EventCounter " << EventCounter << endl;
	//}			// to here
	

	
	EventCounter++;
	//cout << "EventCounter " << EventCounter << endl;
	
	// see comments in UnpackProc
	out_evt->SetValid(isValid);
	return isValid;
}
