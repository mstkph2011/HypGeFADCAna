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
#ifdef EXA_CODE
   ,
   fSum1(0), fSum2(0), fSum3(0),
   fParam1(0), fWinCon(0)
#endif
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
					//snprintf(chis,15,"Trace_smoothed_%02d",i+1);  
					//snprintf(chead,63,"Trace channel %2d after smooting",i+1);
					//fhTrace_Smoothed[i] = new TH1D (chis,chead,TRACE_LENGTH,0,TRACE_LENGTH);
					//AddHistogram(fhTrace_Smoothed[i],"V1724/Trace_smoothed");
		
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
					//fhBaseline = new TH1D ("Baseline","Baseline",10000,0,10000);		//Kai edit
					//AddHistogram(fhBaseline,"V1724/Baseline");
					//fhAmplitude = new TH1D ("Amplitudeo","A",TRACE_LENGTH,0,TRACE_LENGTH);
					//AddHistogram(fhAmplitude,"V1724/Ampl");
	fhMWD = new TH1D ("MWDo","MWD",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fhMWD,"V1724/MWD");
	fhMWD_MA = new TH1D("MWDMAo","MWDMA",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fhMWD_MA,"V1724/MWDMA");

		//risetime histos
	
	fhRisetime1090 = new TH1D("Risetime1090","Risetime1090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	AddHistogram(fhRisetime1090,"V1724/Risetime1090");
	fhRisetime3090 = new TH1D("Risetime3090","Risetime3090",100,0,1000);		// risetime is multiplied by 10 (10 ns/bin)	--> 1 bin covers 10 ns intervall
	AddHistogram(fhRisetime3090,"V1724/Risetime3090");

		// get histo for energy spectrum
	fhEnergySpectrum = (TH1D*) GetHistogram("V1724/Energyspectrum/Energy");

				//fRisetime = new TH1D ("Risetime", "Risetime of Trace 1", 20000,0,20000);			//Marcell_edt
				//AddHistogram(fRisetime,"V1724/Risetime");
				//fBaseline = new TH1D ("Baseline","Baseline",10000,0,10000);		//Kai edit
				//AddHistogram(fBaseline,"V1724/Baseline");
				//fMovAvg = new TH1D ("MovAvg","Moving Average of Trace 1",TRACE_LENGTH,0,TRACE_LENGTH);
				//AddHistogram(fMovAvg,"V1724/MovAvg");
	MovAvgArr = new int[TRACE_LENGTH];
	fSmoothedTrace = new TH1D ("SmoothedTrace","SmoothedTrace",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fSmoothedTrace,"V1724/SmoothedTrace");
	fAmplitude = new TH1D ("Amplitude","A",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fAmplitude,"V1724/Ampl");
	fMWD = new TH1D ("MWD","MWD",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fMWD,"V1724/MWD");
				//fEnergy = new TH1D ("Energyo","Energy",20000,0,2000);
				//AddHistogram(fEnergy,"V1724/Energy");
				//fRisetime1090 = new TH1D("Risetime1090o","Risetime1090",10000,0,1000);
				//AddHistogram(fRisetime1090,"V1724/Risetime1090");
				//fRisetime3090 = new TH1D("Risetime3090o","Risetime3090",10000,0,1000);
				//AddHistogram(fRisetime3090,"V1724/Risetime3090");
	fMWDMA = new TH1D("MWDMA","MWDMA",TRACE_LENGTH,0,TRACE_LENGTH);
	AddHistogram(fMWDMA,"V1724/MWDMA");
	fEnergyRise1090Corr = new TH2D("EnergyRise1090Corr","Enegy-Risetime1090-Correlation;Rt;E",100,0,1000,2000,0,2000);
	AddHistogram(fEnergyRise1090Corr,"V1724/EnergyRise1090Corr");
	fEnergyRise3090Corr = new TH2D("EnergyRise3090Corr","Enegy-Risetime3090-Correlation;Rt;E",100,0,1000,2000,0,2000);
	AddHistogram(fEnergyRise3090Corr,"V1724/EnergyRise3090Corr");
				//fEnergyRiseCorr = MakeTH2('D',"ERC","ERC",2000,0,2000,100,0,1000);
	
	
	fhAmplBaselinegradient= new TH1D ("Baselinegrad","Gradient of the start of the baseline in amplitude signal",1000,0,10);
	AddHistogram(fhAmplBaselinegradient,"V1724/BaseGrad");
		
		// get parameters
	fHypPar = (THypGeParameter*)  GetParameter("HypGeParameter");
	cout << "fHypPar" <<fHypPar << endl;
		//cout << "\tMWDm " << fHypPar->GetMWDm() << endl;
		// real analysis object
	fMWDAna = new THypGeMWD(TRACE_LENGTH);
	
	//fMWDAna->FullAnalysis(fhTrace[0],fSmoothedTrace,fhTrace_BaseCorr[0], fAmplitude,fMWD,fhEnergySpectrum,fhRisetime1090,fhRisetime3090,fMWDMA,(TH2D*) fEnergyRise1090Corr,(TH2D*) fEnergyRise3090Corr);
	/*
	Int_t NoOfSmoothing = 100;
	
	
	S = new TMatrixD (0,TRACE_LENGTH-1,0,TRACE_LENGTH-1);
	T = new TMatrixD(0,TRACE_LENGTH-1,0,TRACE_LENGTH-1);
	TMatrixD *N = new TMatrixD(0,TRACE_LENGTH-1,0,TRACE_LENGTH-1);
	w = new TVectorD(TRACE_LENGTH);
	v = new TVectorD(TRACE_LENGTH);
	
	for(Int_t i=0;i<TRACE_LENGTH;i++)
	{
		(*S)(i,i)=0.5;
		if (i)
			(*S)(i,i-1)=0.25;
		if (i != TRACE_LENGTH-1)	
			(*S)(i,i+1)=0.25;
	}
	(*S)(0,0)=0.5;
	(*S)(0,1)=0.5;
	(*S)(TRACE_LENGTH-1,TRACE_LENGTH-2)=0.5;
	(*S)(TRACE_LENGTH-1,TRACE_LENGTH-1)=0.5;
	(*T) = (*S);
	Double_t result;
	for(Int_t i=1;i<NoOfSmoothing;i++)
	{
		cout << "Loaded " << i*100/NoOfSmoothing << "%" << endl;
		for(Int_t k=0;k<TRACE_LENGTH;k++)
		{
			//cout << "k " << k << endl;
			for(Int_t j=max(0,k-1-i);j<=min(k+1+i,TRACE_LENGTH-1);j++)
			{
			//	cout << "j " << j << endl;
					result = 0;
					for(Int_t l=max(0,k-1-i);l<=min(k+1+i,TRACE_LENGTH-1);l++)
					{
						//cout << "l " << l << endl;
						result += (*S)(k,l) * (*T)(l,j);
					}
				(*N)(k,j) = result;
				//cout << "\t\tresult" << result << endl;
			}
		}	
		(*S) = (*N);
		
	}
	cout << (*S)(0,0) << "\t" << (*S)(0,1) << endl;
	cout << (*N)(0,0) << "\t" << (*N)(0,1) << endl;
	//S->Print();
	fMWDAna->GetS(S);
	*/
	//oldMax is used to prevent the analysis of a trace if the trace does not change (as in the actual (bad) data 2.4.14
	oldMax=0;
	cout << "**** THypGeAnlProc: Create" << endl;
#ifdef EXA_CODE
   //// init user analysis objects:
   
////////////////////////////////////////////////////
// uncomment following lines if you want to use external macro to set parameter values:
//   // if unpack was enabled, parameters have been printed already
//   if(TGo4Analysis::Instance()->GetAnalysisStep("Unpack")->IsProcessEnabled())
//	    gROOT->ProcessLine(".x setparam.C(0)");
//   else gROOT->ProcessLine(".x setparam.C(1)");
///////////////////////////////////////////////////

   if(fParam1->fbHisto){
      // this one is created in THypGeAnalysis, because it is used in both steps
      fWinCon = (TGo4WinCond *) GetAnalysisCondition("wincon1");
      if (fWinCon) fWinCon->PrintCondition(true);
      fSum1     = MakeTH1('I', "Sum1", "Sum over 8 channels", 5000, 1., 5001.);
      fSum2     = MakeTH1('I', "Sum2", "Sum over 8 channels shift 1", 5000, 1., 5001.);
      fSum3     = MakeTH1('I', "Sum3", "Sum over 8 channels shift 2", 5000, 1., 5001.);
   }
#endif
}
//-----------------------------------------------------------
THypGeAnlProc::~THypGeAnlProc()
{
   cout << "**** THypGeAnlProc: Delete" << endl;
#ifdef EXA_CODE
   if(fParam1->fbHisto){
      if (fWinCon) fWinCon->PrintCondition(true);
   }
#endif
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

#ifdef EXA_CODE
   Int_t cnt(0);
   for(Int_t ii=0;ii<HypGe_NUM_CHAN/2;ii++) {
      out_evt->frData[cnt] = (Float_t)inp_evt->fiCrate1[ii];
      cnt++;
   }
   for(Int_t ii=0; ii<HypGe_NUM_CHAN/2; ii++) {
      out_evt->frData[cnt] = (Float_t)inp_evt->fiCrate2[ii];
      cnt++;
   }
   if(fParam1->fbHisto) { // histogramming
      for(Int_t ii=0;ii<HypGe_NUM_CHAN;ii++)
         if(out_evt->frData[ii]) {
            if(fWinCon && fWinCon->Test(out_evt->frData[ii])) fSum1->Fill(out_evt->frData[ii]);
            if (fParam1) fSum2->Fill(out_evt->frData[ii] + fParam1->frP1);
            if (fParam1) fSum3->Fill(out_evt->frData[ii] + fParam1->frP2);
         }
   }
#endif
	
	//add ana code here
	
	//for(Int_t i = 1; i< TRACE_LENGTH;i++)
	//	fhTrace_BaseCorr[0]->SetBinContent(i,fhTrace[0]->GetBinContent(i));
	
	
	//cout << "M " << fParam1->M << endl;
	//cout << fhTrace[0]->GetMaximum() << endl;
	
	//fHypPar->PrintParameters();
	fMWDAna->SetUseMWD(1);
	fMWDAna->SetParameters(fHypPar->GetMWDm(),fHypPar->GetMAl(),fHypPar->GetNoOfSmoothing(),fHypPar->GetWidth() ,fHypPar->GetSigmaGaus(),fHypPar->GetSigmaBil(),fHypPar->GetTau(), fHypPar->GetEnableMA(),fHypPar->GetEnableSmoothing(),fHypPar->GetEnableBaselineCorrection());

	if (fMWDAna->FullAnalysis(fhTrace[0],fSmoothedTrace,fhTrace_BaseCorr[0], fAmplitude,fMWD,fhEnergySpectrum,fhRisetime1090,fhRisetime3090,fMWDMA,(TH2D*) fEnergyRise1090Corr,(TH2D*) fEnergyRise3090Corr) != -1)
		fhAmplBaselinegradient->Fill((fAmplitude->GetBinContent(1)-fAmplitude->GetBinContent(301))/300);
	
	
		//this shows number of real events
		//EventCounter++;
		//cout << "EventCounter " << EventCounter << endl;

	

	
	EventCounter++;
	//cout << "EventCounter " << EventCounter << endl;
	
	// see comments in UnpackProc
	out_evt->SetValid(isValid);
	return isValid;
}
