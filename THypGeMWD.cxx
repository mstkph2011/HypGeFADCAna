//      THypGeMWD.cpp
//      
//      Copyright 2013 Kai Thomas Rittgen <rittgen@lynxsv>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include "THypGeMWD.h"
#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVirtualFFT.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <cassert>

using namespace std;

//default constructor DON'T use!
THypGeMWD::THypGeMWD()
{
	
}

THypGeMWD::THypGeMWD(Int_t TraceLength_ext)            //constructor
{
	TraceLength=TraceLength_ext;
	cout << "Analysis created"<<endl;
	hTraceBuffer =new TH1D();
	M = 200;					// window width for MWD
	L = 100;					// top width for MA
	NoOfSmoothing = 100;
	Width = 3;
	Sigma =11;					// sigma for gaussian smoothing
	Sigma2 = 1500;
	tau = 5383;				// tau of pre amp in samples (1 sample = 10 ns @ 100 MSa/s)
	EnableMA = 0;			// Switch for second moving average filter
	SmoothingMethod = 3;	// Switch for smoothing methods
	EnableBaselineCorrection = 1; 	//Switch baseline correction on or off
	PileUpTimeThreshold = 20;	// in µs, no ext parameter yet (28.03.14)
	
	Aarray = new Double_t[TraceLength+1];
	MWDarray = new Double_t[TraceLength+1];
	GradMWD1array = new Double_t[TraceLength+1];
	GradMWD2array = new Double_t[TraceLength+1];
	Sarray = new Double_t[TraceLength+1];		// array to store the result of the direct filter
	Parray = new Double_t[TraceLength+1];		// array to store an intermediate result of the direct filter
	useMWD= 1;
	GausNorm = new Double_t[TraceLength];
	CalculateGausCoeff();
	
	EnergyRtCorrFunc = new TF1("EnergyRtCorrFunc","pol4",0,2000);
		EnergyRtCorrFunc->SetParameters(1,-9.1237e-05, 7.95504e-07, -3.0778e-09, 4.39764e-12);
		
	EnergyPileUpTimeCorrFunc = new TF1("EnergyPileUpTimeCorrFunc","[0]*(1-[1]/pow(x,[2]))",0,2000);
		EnergyPileUpTimeCorrFunc->SetParameters(1,3.74495e-02,1.95113);
	//CalculateSecondGausCoeff();
}

THypGeMWD::~THypGeMWD()            //destructor
{
	cout << "HypGeMWD deconstructed!" << endl;
}

//Double_t THypGeMWD::FullAnalysis (TH1D* hTrace[0]_ext, TH1D* hSmoothedTrace[0], TH1D* hTrace_bc[0], TH1D* hAmplitude[0],TH1D* hMWD[0], TH1D* hEnergySpectrum, TH1D* hRisetime1090, TH1D* hRisetime3090, TH1D* hMWDMA[0], TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr, TH2D* hEnergyTimeSinceLastPulse)
//{
	//bool UseTimer = 0;
	//if (UseTimer)
	//{ 
		//timer.Reset();
		//timer.Start();
	//}
	////cout << SmoothingMethod << endl;
	//if ( SmoothingMethod)
	//{
		//Smoothing(hTrace[0]_ext, hSmoothedTrace[0]);
		//if (UseTimer) 
		//{
			//timer.Stop();
			//cout << "Smoothing step took " << timer.RealTime() << " seconds" << endl;
		//}
	//}
	//if (UseTimer)
	//{ 
		//timer.Reset();
		//timer.Start();
	//}
	//Baseline(hTrace[0]_ext, hTrace_bc[0]);
	//if (UseTimer) 
	//{	
		//timer.Stop();
		//cout << "Baseline correction step took " << timer.RealTime() << " seconds" << endl;
	//}
	
	
	//if (UseTimer)
	//{ 
		//timer.Reset();
		//timer.Start();
	//}
	//if (MWD(hTrace[0]_ext, hMWD[0], hAmplitude[0]) == -1)
		//return -1;
	
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Devonvolution and MWD step took " << timer.RealTime() << " seconds" << endl;
	//}
	
	//if (EnableMA)
	//{
		//if (UseTimer)
		//{ 
			//timer.Reset();
			//timer.Start();
		//}
		//MA(hTrace[0]_ext, hMWDMA[0], hMWD[0]);
		//if (UseTimer)
		//{	
			//timer.Stop();
			//cout << "MA step took " << timer.RealTime() << " seconds" << endl;
		//}
	//}
	
	//if (UseTimer)
	//{ 
		//timer.Reset();
		//timer.Start();
	//}
	//if ( useMWD)
		//Energyspectrum(hTrace[0]_ext, hEnergySpectrum, hMWD[0]);				//needs MWD
	//else
		//Energyspectrum(hTrace[0]_ext, hEnergySpectrum, hAmplitude[0]);		//needs MWD
	
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Energyspectrum step took " << timer.RealTime() << " seconds" << endl;

		//timer.Reset();
		//timer.Start();
	//}
	//Risetime(hTrace[0]_ext, hRisetime1090, hRisetime3090);		//needs Energyspectrum
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Risetime step took " << timer.RealTime() << " seconds" << endl;

		//timer.Reset();
		//timer.Start();
	//}
	//ERC(hEnergyRise1090Corr, hEnergyRise3090Corr);			//needs Energyspectrum and Risetime
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Energy Rt correlation step took " << timer.RealTime() << " seconds" << endl;

		//timer.Reset();
		//timer.Start();
	//}
	//AnaStep_EnergyTimeSinceLastPulseCorrelation(hEnergyTimeSinceLastPulse);
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Energy TSLPulse correlation step took " << timer.RealTime() << " seconds" << endl;
		
		//timer.Reset();
		//timer.Start();
	//}
	//AnaStep_FillEnergySpectrumWithPileUpCut(PileUpTimeThreshold);
	//if (UseTimer)
	//{	
		//timer.Stop();
		//cout << "Energy Spectrum with cut step took " << timer.RealTime() << " seconds" << endl;
	//}
	//return 0;
//}

Double_t THypGeMWD::FullAnalysis ()
{
	//cout << "entering ana" << endl;
	bool UseTimer = 0;
	if (UseTimer)
	{ 
		timer.Reset();
		timer.Start();
	}
	//cout << SmoothingMethod << endl;
	if ( SmoothingMethod)
	{
		AnaStep_Smoothing();
		if (UseTimer) 
		{
			timer.Stop();
			cout << "Smoothing step took " << timer.RealTime() << " seconds" << endl;
		}
	}
	if (UseTimer)
	{ 
		timer.Reset();
		timer.Start();
	}
	AnaStep_BaselineCorrection();
	if (UseTimer) 
	{	
		timer.Stop();
		cout << "Baseline correction step took " << timer.RealTime() << " seconds" << endl;
	}
	
	
	if (UseTimer)
	{ 
		timer.Reset();
		timer.Start();
	}
	if (AnaStep_DoMovingWindowDeconvolution() == -1)
		return -1;
	
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Devonvolution and MWD step took " << timer.RealTime() << " seconds" << endl;
	}
	
	if (EnableMA)
	{
		if (UseTimer)
		{ 
			timer.Reset();
			timer.Start();
		}
		AnaStep_DoMovingAverageFilter();
		if (UseTimer)
		{	
			timer.Stop();
			cout << "MA step took " << timer.RealTime() << " seconds" << endl;
		}
	}
	
	if (UseTimer)
	{ 
		timer.Reset();
		timer.Start();
	}
	
	AnaStep_FillEnergyspectrum();				//needs MWD
	
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Energyspectrum step took " << timer.RealTime() << " seconds" << endl;

		timer.Reset();
		timer.Start();
	}
	AnaStep_ExtractRisetime();		//needs Energyspectrum
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Risetime step took " << timer.RealTime() << " seconds" << endl;

		timer.Reset();
		timer.Start();
	}
	AnaStep_DoEnergyRisetimeCorrelation();			//needs Energyspectrum and Risetime
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Energy Rt correlation step took " << timer.RealTime() << " seconds" << endl;

		timer.Reset();
		timer.Start();
	}
	AnaStep_EnergyTimeSinceLastPulseCorrelation();
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Energy TSLPulse correlation step took " << timer.RealTime() << " seconds" << endl;
		
		timer.Reset();
		timer.Start();
	}
	AnaStep_FillEnergySpectrumWithPileUpCut(PileUpTimeThreshold);
	if (UseTimer)
	{	
		timer.Stop();
		cout << "Energy Spectrum with cut step took " << timer.RealTime() << " seconds" << endl;
	}
	AnaStep_DoDirectFilter();
	return 0;
}
	
	//Smoothing
	
Int_t THypGeMWD::AnaStep_Smoothing()
{
	
	SetTrace(hTrace[0]);
	
	if( SmoothingMethod == 1)
		DoMeanFilter();
	if( SmoothingMethod == 2)
		DoWeightedAverageFilter();
	if( SmoothingMethod == 3)
		DoGaussianFilter();
	if( SmoothingMethod == 4)	
		DoBilateralFilter();
		

	for (Int_t i =0; i <= hSmoothedTrace[0]->GetNbinsX(); i++)
		hTraceBuffer->SetBinContent(i,hSmoothedTrace[0]->GetBinContent(i));
	return 0;
}
	
	//Baseline
	
Int_t THypGeMWD::AnaStep_BaselineCorrection()
{
	if (hTraceBuffer == 0)
		SetTrace(hTrace[0]);
		
	if (hTraceBuffer->GetMaximum() == 0)
		return -1;
	Double_t sumoffset = 0;
	Int_t BaseLineCorrectionBins = 20;//100;		//building average of first 100 bins
	
		for(Int_t i=1;i<=BaseLineCorrectionBins;i++)
		{
			sumoffset = sumoffset + hTraceBuffer->GetBinContent(i);
		}

		 offset_av=hTraceBuffer->GetMinimum();//sumoffset/BaseLineCorrectionBins;			// changed in Jülich, trace must be positive to avoid a strong negative drift of the deconvoluted signal!, remaining offset is corrected after deconvolution
		
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hTraceBuffer->SetBinContent(i,hTraceBuffer->GetBinContent(i) - offset_av);
			hTrace_bc[0]->SetBinContent(i,hTraceBuffer->GetBinContent(i));
		}
	return 0;
}
	
	//Moving-Window-Deconvolution
	
Int_t THypGeMWD::AnaStep_DoMovingWindowDeconvolution()
{
	//if (hTraceBuffer == 0)
	//SetTrace(hTrace[0]);
	
	Aarray[0] = hTraceBuffer->GetBinContent(1);			// first value of baseline corrected trace
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Aarray[i] = hTraceBuffer->GetBinContent(i) - hTraceBuffer->GetBinContent(i-1) * (1.-(1./tau)) + Aarray[i-1];
	}
	//bc correction of amp signal
	Aarray[0]=0;
	Double_t SumAmpOffset=0;
	Int_t LengthSumOffsetAverage=20;
	for(Int_t i = 1; i <= LengthSumOffsetAverage; i++)
		SumAmpOffset+= Aarray[i];
	SumAmpOffset= SumAmpOffset/LengthSumOffsetAverage;
		
	// try to correct slope in baseline - not working (yet)
							/*	Double_t mCorrA[301];			// mCorrA[0] = sum and later on the mean
								mCorrA[0] = 0;
								for (Int_t i=1; i <= 300; i++)
								{
									mCorrA[i] = (Aarray[1]-Aarray[i+1])/i;				// get gradient of first 300 slope triangles
									mCorrA[0] += mCorrA[i];
								}
								mCorrection = mCorrA[0]/300;					// mean of slope
	*/
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Aarray[i]= Aarray[i]-SumAmpOffset; //+ mCorrection*i;						// bc correction of amplitude signal
		if (i > Int_t(M))
			MWDarray[i] = Aarray[i] - Aarray[i-Int_t(M)];
		else
			MWDarray[i] = Aarray[i]-Aarray[1];				// added to cover edge efects 
		hAmplitude[0]->SetBinContent(i,Aarray[i]);
		hMWD[0]->SetBinContent(i,MWDarray[i]);
	}
	GradMWD1array[0]=0;
	GradMWD2array[0]=0;
	for(Int_t i=1;i<=TraceLength;i++)
	{
		GradMWD1array[i] = hMWD[0]->GetBinContent(i+1) - hMWD[0]->GetBinContent(i);
	}
	for(Int_t i=1;i<=TraceLength;i++)
	{
		GradMWD2array[i] = GradMWD1array[i+1] - GradMWD1array[i];
	}	
	return 0;
}
	
	//Moving-Average
	
Int_t THypGeMWD::AnaStep_DoMovingAverageFilter()
{
	if (hTraceBuffer == 0)
		SetTrace(hTrace[0]);
	
	Double_t  MWDMA[TraceLength+1];
	
	MWDMA[0] = hMWD[0]->GetBinContent(1);
	for(Int_t n=1;n<=TraceLength;n++)
	{
		if (n > Int_t(L))
			MWDMA[n] = MWDMA[n-1] + 1./L * (hMWD[0]->GetBinContent(n) - hMWD[0]->GetBinContent(n-Int_t(L)));
		else
			MWDMA[n] = MWDMA[n-1] + 1./L * (hMWD[0]->GetBinContent(n) - hMWD[0]->GetBinContent(1));
		
		hTraceBuffer->SetBinContent(n,MWDMA[n]);
		hMWDMA[0]->SetBinContent(n,MWDMA[n]);
	}

	return 0;
}

Int_t		THypGeMWD::AnaStep_DoDirectFilter()
{
	// this filter replaces MWD + MA (other version of it, after Jordanov/Knoll), do this step aber smoothing and baseline correction
	
	Double_t d_kl;		// variable to calculate a needed value from the trace
	Sarray[0]= 0;
	Parray[0] = 0;
	
	for(Int_t n = 1;n <= TraceLength;n++)
	{
		d_kl=hTrace_bc[0]->GetBinContent(n);
		if ( n-L > 0)
			d_kl = d_kl - hTrace_bc[0]->GetBinContent(n-L);
		if (n-M >0) 
			d_kl = d_kl - hTrace_bc[0]->GetBinContent(n-M);
		if (n-M-L > 0)
			d_kl += hTrace_bc[0]->GetBinContent(n-M-L);
	
		Parray[n]=Parray[n-1]+d_kl;
		Sarray[n]=Sarray[n-1]+Parray[n]+ tau* d_kl;
		hTrace_Direct[0]->SetBinContent(n,Sarray[n]);
	}
	return 0;
}	
	//PileupCompensation and Energyspectrum
		
Int_t THypGeMWD::AnaStep_FillEnergyspectrum()
{
	
	SetTrace(hMWD[0]);
	if (useMWD)
	{
		EvaluateMWD();
	}
	else
	{
		EvaluateAmplitude();
	}
	
	return 0;
}
	
	//Risetime
	
Int_t THypGeMWD::AnaStep_ExtractRisetime()
{
	//cout << "new Trace" << endl;
	if (hTrace[0] == 0)
		SetTrace(hTrace[0]);
	/*
	risetime1090.clear();
	risetime3090.clear();
	
	if(leftborder.size() == rightborder.size())
	{	
		Double_t Grad1[TraceLength];
		Double_t Grad2[TraceLength];
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Grad1[i] = hTrace[0]_ext->GetBinContent(i+5) - hTrace[0]_ext->GetBinContent(i);
		}
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Grad2[i] = Grad1[i+5] - Grad1[i];
		}
		/*
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hRisetime1090->SetBinContent(i+5,Grad2[i]);
		}
		
		for(Int_t k=0;k<energy.size();k++)
		{
			//cout << "k " << k << endl;
				
			Int_t leftborderrisegrad = FindLocalMaximumBin(Grad2,max(leftborder[k]-20,1),min(leftborder[k]+40,TraceLength));
			Int_t rightborderrisegrad = FindLocalMinimumBin(Grad2,max(rightborder[k]-40,1),min(rightborder[k]+20, TraceLength));
			
			Double_t leftborderrise = FindLocalMinimumBin(hTrace[0]_ext,max(leftborderrisegrad-5,1),min(leftborderrisegrad+5,TraceLength));
			Double_t rightborderrise = FindLocalMaximumBin(hTrace[0]_ext,max(rightborderrisegrad-5,1),min(rightborderrisegrad+5, TraceLength));
			
			//Double_t rightborderrise = FindFirstBinAboveInterpolated(Grad2,0,FindLocalMinimumBin(Grad2,max(rightborder[k]-40,1),min(rightborder[k]+20, TraceLength)),min(FindLocalMinimumBin(Grad2,max(rightborder[k]-40,1),min(rightborder[k]+20, TraceLength))+20,TraceLength));
			//cout << "lfbr" << hTrace[0]_ext->GetBinContent(leftborderrise) << " rfbr "  << hTrace[0]_ext->GetBinContent(rightborderrise) <<endl;
			Double_t mean_low = 0;
			for(Int_t i=leftborderrise-9;i<=leftborderrise;i++)
			{
				mean_low += hTrace[0]_ext->GetBinContent(i);
			}
			mean_low = mean_low/10.;
			//cout << "mv1 " << mean_low << endl;
			Double_t mean_up = 0;
			for(Int_t i=rightborderrise;i<=rightborderrise+9;i++)
			{
				mean_up += hTrace[0]_ext->GetBinContent(i);
			}
			mean_up = mean_up/10.;
			//cout << "mv2 " << mean_up << endl;
			if(mean_up > mean_low)
			{
				Double_t l = (mean_up-mean_low)/10;		//Calculates 10% of height of the rising signal
				//cout << "llllll " << mean_low+9*l << endl;
				//Double_t RiseX1090 = FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+1*l,leftborderrise,rightborderrise);
				//Double_t RiseX1090 = Double_t(FindFirstBinAbove(hRisetime1090,leftborderrise+mean_low+9*l,leftborderrise,rightborderrise)) - Double_t(FindFirstBinAbove(hRisetime1090,leftborderrise+mean_low+l,leftborderrise,rightborderrise));
				//Double_t RiseX3090 = FindFirstBinAbove(hTrace[0]_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAbove(hTrace[0]_ext,mean_low+3*l,leftborderrise,rightborderrise);
				//cout << "bla" <<RiseX3090 << endl;	
				Double_t RiseX1090 = FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+l,leftborderrise,rightborderrise);
				Double_t RiseX3090 = FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace[0]_ext,mean_low+3*l,leftborderrise,rightborderrise);
				//Double_t RiseX1030 = RiseX1090 -  RiseX3090;
				//cout << leftborderrise<< " blubb " << rightborderrise<< " blubb "<<RiseX3090 << endl;
				Double_t RiseT1090 = RiseX1090 * 10;
				Double_t RiseT3090 = RiseX3090 * 10;
				
				//if (RiseX3090 > 1)
				{
					hRisetime1090->Fill(RiseT1090);
					hRisetime3090->Fill(RiseT3090);
			
					risetime1090.push_back(RiseT1090);		//used for Energy-Risetime-Correlation
					risetime3090.push_back(RiseT3090);		//used for Energy-Risetime-Correlation
				}
			}
			
		}
	}
	*/
	/*	
	risetime1090.clear();
	risetime3090.clear();
	Int_t abort = 0;
	Int_t i = 0;
	Double_t rightborder = -9;
	Double_t tresholdRisetime = 50;
	Double_t gradRisetime = 5;
	
	do
	{
		i = rightborder+10;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort=1;
				break;
			}
			i++;
		}
		while(hTrace[0]_ext->GetBinContent(i+1) - hTrace[0]_ext->GetBinContent(i) < gradRisetime && hTrace[0]_ext->GetBinContent(i+2) - hTrace[0]_ext->GetBinContent(i+1) < gradRisetime && hTrace[0]_ext->GetBinContent(i+3) - hTrace[0]_ext->GetBinContent(i+2) < gradRisetime);
		Double_t leftborder = hTrace[0]_ext->GetBinContent(i);
		i = leftborder+10;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort=1;
				break;
			}
			i++;
		}
		while(hTrace[0]_ext->GetBinContent(i+1) - hTrace[0]_ext->GetBinContent(i) >= 0 && hTrace[0]_ext->GetBinContent(i+2) - hTrace[0]_ext->GetBinContent(i+1) >= 0 && hTrace[0]_ext->GetBinContent(i+3) - hTrace[0]_ext->GetBinContent(i+2) >= 0);
		Double_t rightborder = hTrace[0]_ext->GetBinContent(i);
		
		if(hTrace[0]_ext->GetBinContent(rightborder) - hTrace[0]_ext->GetBinContent(leftborder) > tresholdRisetime && rightborder - leftborder < 100)
		{
			Double_t mean_low = 0;
			for(Int_t i=leftborder-29;i<=leftborder;i++)
			{
				mean_low += hTrace[0]_ext->GetBinContent(i);
			}
			mean_low = mean_low/30;
		
			Double_t mean_up = 0;
			for(Int_t i=rightborder;i<=rightborder+29;i++)
			{
				mean_up += hTrace[0]_ext->GetBinContent(i);
			}
			mean_up = mean_up/30;
		
			cout << leftborder << endl;
			cout << rightborder << endl;
		
			Double_t l = (mean_up-mean_low)/10;		//Calculates 10% of height of the rising signal
		
			Double_t RiseX1090 = Double_t(FindFirstBinAbove(hTrace[0]_ext,leftborder+mean_low+9*l,leftborder,rightborder)) - Double_t(FindFirstBinAbove(hTrace[0]_ext,leftborder+mean_low+l,leftborder,rightborder));
			Double_t RiseX3090 = Double_t(FindFirstBinAbove(hTrace[0]_ext,leftborder+mean_low+l,leftborder,rightborder)) - Double_t(FindFirstBinAbove(hTrace[0]_ext,leftborder+mean_low+3*l,leftborder,rightborder));
						
			Double_t RiseT1090 = RiseX1090 * 10;		//Conversion into nanoseconds (FADC 100 MS/s)
			Double_t RiseT3090 = RiseX3090 * 10;
	
			hRisetime1090->Fill(RiseT1090);
			hRisetime3090->Fill(RiseT3090);
		
			risetime1090.push_back(RiseT1090);		//used for Energy-Risetime-Correlation
			risetime3090.push_back(RiseT3090);		//used for Energy-Risetime-Correlation
		}
	}
	while(abort == 0);
		
	/*
	for(Int_t k=0;k<=Int_t(leftborder.size());k++)
	{
		Int_t d = leftborder[k]+20;
		do
		{
			if(d+50 >= TraceLength)
			{
				abort=1;
				break;
			}
			d++;
		}	
		while(hTrace[0]_ext->GetBinContent(d+1) - hTrace[0]_ext->GetBinContent(d) >= 0 && abort==0);
		if(abort==0)
		{
		Double_t rightborder = hTrace[0]_ext->GetBin(d);
		

		}
	}
	*/
	if (hTrace[0] == 0)
	SetTrace(hTrace[0]);
		
	risetime1090.clear();
	risetime3090.clear();

	Double_t threshold_Risetime = 50;		//threshold up where signals were identified as useful for Risetime calculation

	for(Int_t k=0;k<Int_t(leftborder.size());k++)
	{		
		Double_t rightborder = FindLocalMaximumBin(hTrace[0],leftborder[k],leftborder[k]+40);
						
		if(hTrace[0]->GetBinContent(rightborder) > threshold_Risetime+hTrace[0]->GetBinContent(leftborder[k]))
		{
			Double_t l = (hTrace[0]->GetBinContent(rightborder)-hTrace[0]->GetBinContent(leftborder[k]))/10;		//Calculates 10% of height of the rising signal
						
			Double_t RiseX1090 = Double_t(FindFirstBinAbove(hTrace[0],hTrace[0]->GetBinContent(rightborder)-l,leftborder[k],rightborder)) - Double_t(FindFirstBinAbove(hTrace[0],hTrace[0]->GetBinContent(rightborder)-9*l,leftborder[k],rightborder));
			Double_t RiseX3090 = Double_t(FindFirstBinAbove(hTrace[0],hTrace[0]->GetBinContent(rightborder)-l,leftborder[k],rightborder)) - Double_t(FindFirstBinAbove(hTrace[0],hTrace[0]->GetBinContent(rightborder)-7*l,leftborder[k],rightborder));
			
			Double_t RiseT1090 = RiseX1090 * 10;		//Conversion into nanoseconds (FADC 100 MS/s)
			Double_t RiseT3090 = RiseX3090 * 10;
	
			hRisetime1090->Fill(RiseT1090);
			hRisetime3090->Fill(RiseT3090);
			
			risetime1090.push_back(RiseT1090);		//used for Energy-Risetime-Correlation
			risetime3090.push_back(RiseT3090);		//used for Energy-Risetime-Correlation
		}
	}
	return 0;
}
	
	//Energy-Risetime-Correlation
	
Int_t THypGeMWD::AnaStep_DoEnergyRisetimeCorrelation()
{
		
	if(Int_t(energy.size()) == Int_t(risetime1090.size()))
	{
		for(Int_t n=0;n < Int_t(energy.size());n++)
		{
			hEnergyRise1090Corr->Fill(risetime1090[n],energy[n]);
		}
	}
	
	if(Int_t(energy.size()) == Int_t(risetime3090.size()))
	{
		for(Int_t n=0;n < Int_t(energy.size());n++)
		{
			hEnergyRise3090Corr->Fill(risetime3090[n],energy[n]);
		}
	}
	return 0; 
}
	


TH1D* THypGeMWD::GetTrace()
{
	return hTraceBuffer;
}

void THypGeMWD::SetTrace(TH1D* hTrace_ext)
{

	if (hTraceBuffer )
		delete hTraceBuffer;
	hTraceBuffer = (TH1D*) hTrace_ext->Clone("hTraceBuffer");
	
}



Int_t THypGeMWD::FindFirstBinAbove(TH1D* fHisto,Double_t threshold,Int_t low, Int_t high)
{
   //find first bin with content > threshold between low and high
   //if no bins with content > threshold is found the function returns -1.
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (fHisto->GetBinContent(bin) > threshold) 
				return bin;
   }
   return -1;
}

Double_t THypGeMWD::FindFirstBinAboveInterpolated(TH1D* fHisto, Double_t threshold,Int_t low, Int_t high)
{
   //find first bin with content > threshold between low and high
   //after that linear interpolation is done to improve the resolution
   //if no bins with content > threshold is found the function returns -1.
   Double_t InterBin = -1;
   for (Int_t bin=low;bin<=high+1;bin++) 
   {
	   //cout << bin << " ";
		if (fHisto->GetBinContent(bin) > threshold)
		{
			InterBin = bin-1 + (threshold -fHisto->GetBinContent(bin-1))/(fHisto->GetBinContent(bin) -fHisto->GetBinContent(bin-1));
			//cout << "InterBin " << InterBin << endl;
			return InterBin;			
		}
   }
   return InterBin;
}
Double_t THypGeMWD::FindFirstBinAboveInterpolated(Double_t *Array, Double_t threshold,Int_t low, Int_t high)
{
   //find first bin with content > threshold between low and high
   //after that linear interpolation is done to improve the resolution
   //if no bins with content > threshold is found the function returns -1.
   Double_t InterBin = -1;
   for (Int_t bin=low;bin<=high+1;bin++) 
   {
	   //cout << bin << " ";
		if (Array[bin] > threshold)
		{
			InterBin = bin-1 + (threshold -Array[bin-1])/(Array[bin] -Array[bin-1]);
			//cout << "InterBin " << InterBin << endl;
			return InterBin;			
		}
   }
   return InterBin;
}
Double_t THypGeMWD::FindLocalMaximum(TH1D* fHisto,Int_t low, Int_t high)
{
	//Find local maximum between low and high
   Double_t Max = 0;
      
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (fHisto->GetBinContent(bin) > Max) 
				Max = fHisto->GetBinContent(bin);
   }
   return Max;
}
Double_t THypGeMWD::FindLocalMaximum(Double_t *Array,Int_t low, Int_t high)
{
	//Find local maximum between low and high
   Double_t Max = 0;
      
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (Array[bin] > Max) 
				Max = Array[bin];
   }
   return Max;
}

Double_t THypGeMWD::FindLocalMinimum(Double_t *Array,Int_t low, Int_t high)
{
	//Find local minimum between low and high
   Double_t Min = 0;
      
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (Array[bin] < Min) 
				Min = Array[bin];
   }
   return Min;
}

Double_t THypGeMWD::FindLocalMinimum(TH1D* fHisto2,Int_t low2, Int_t high2)
{
	//Find local minimum between low and high
   Double_t Min = 0;
      
   for (Int_t bin=low2;bin<=high2;bin++) 
   {
      if (fHisto2->GetBinContent(bin) < Min) 
				Min = fHisto2->GetBinContent(bin);
   }
   return Min;
}


Int_t THypGeMWD::FindLocalMaximumBin(TH1D* fHisto,Int_t low, Int_t high)
{
	//Find local maximum bin between low and high
   Double_t Max = 0;
   Int_t MaxBin = low;
   
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (fHisto->GetBinContent(bin) > Max)
      {
				Max = fHisto->GetBinContent(bin);
				MaxBin = bin;
			}
   }
   return MaxBin;
}
Int_t THypGeMWD::FindLocalMaximumBin(Double_t *Array,Int_t low, Int_t high)
{
	//Find local maximum bin between low and high
   Double_t Max = 0;
  Int_t MaxBin = low;
   
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (Array[bin] > Max)
      {
				Max = Array[bin];
				MaxBin = bin;
			}
   }
   return MaxBin;
}

Int_t THypGeMWD::FindLocalMinimumBin(Double_t *Array,Int_t low, Int_t high)
{
	//Find local minimum bin between low and high
   Double_t Min = 0;
   Int_t MinBin = low;
   
   for (Int_t bin=low;bin<=high;bin++) 
   {
      if (Array[bin] < Min)
      {
				Min = Array[bin];
				MinBin = bin;
			}
   }
   return MinBin;
}

Int_t THypGeMWD::FindLocalMinimumBin(TH1D* fHisto2,Int_t low2, Int_t high2)				// whatever this function is used for (steinen)
{
	//Find local minimum bin between low and high, starting at the high value
   Double_t Min = 0;
   Int_t MinBin = low2;
   
   for (Int_t bin=high2 ; bin >=low2; --bin) 								// 26.2.14 : correction of loop header (steinen)
   {
      if (fHisto2->GetBinContent(bin) < Min)
      {
				Min = fHisto2->GetBinContent(bin);
				MinBin = bin;
			}
   }
   return MinBin;
}



void THypGeMWD::SetParameters( Int_t M_ext, Int_t L_ext,Int_t NoS_ext, Int_t Width_ext, Int_t Sigma_ext, Int_t SigmaBil_ext, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC)
{
	M = M_ext;					// window width for MWD
	L = L_ext;					// top width for MA
	if (Sigma != Sigma_ext)
	{
		Sigma = Sigma_ext;					// sigma for gaussian smoothing
		CalculateGausCoeff();
	}
	NoOfSmoothing= NoS_ext;
	Width= Width_ext;
	Sigma2 = SigmaBil_ext;
	tau = tau_ext;
	EnableMA = EnaMA;			// Switch for second moving average filter
	SmoothingMethod = EnaSmo;	// Switch smoothing on or off
	EnableBaselineCorrection = EnaBC; 	//Switch baseline correction on or off

}

void THypGeMWD::EvaluateAmplitude()
{
	
	//for MWD without MA
	Double_t threshold_MWD = 50;		//threshold up where signals were identified as useful for MWD
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high;
	Int_t abort1 = 0;
	leftborder.clear();
	energy.clear();


	
	do
	{
		Double_t Sumenergy = 0;
		Double_t Sumnoise_left = 0;
		Double_t Energy_av = 0;
		Double_t Noise_av = 0;
	
		Int_t i = leftborder_high+10;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hAmplitude[0]->GetBinContent(i+1) - hAmplitude[0]->GetBinContent(i) < grad_MWD && hAmplitude[0]->GetBinContent(i+2) - hAmplitude[0]->GetBinContent(i+1) < grad_MWD && hAmplitude[0]->GetBinContent(i+3) - hAmplitude[0]->GetBinContent(i+2) < grad_MWD);
		
		leftborder_low = hAmplitude[0]->GetBin(i);
						
		//leftborder.push_back(leftborder_low);		//used for Risetime (add +10 to leftborder_low at high smoothing rates -> results from the discrepance between smoothed and unsmoothed MWD-shape)

		i = leftborder_low+5;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hAmplitude[0]->GetBinContent(i+1) - hAmplitude[0]->GetBinContent(i) >= grad_MWD && hAmplitude[0]->GetBinContent(i+2) - hAmplitude[0]->GetBinContent(i+1) >= grad_MWD && hAmplitude[0]->GetBinContent(i+3) - hAmplitude[0]->GetBinContent(i+2) >= grad_MWD);
		
		leftborder_high = hAmplitude[0]->GetBin(i);
				
				
		if(hAmplitude[0]->GetBinContent(leftborder_low) < threshold_MWD && hAmplitude[0]->GetBinContent(leftborder_high) > threshold_MWD)
		{
			
			for(i=leftborder_high+5;i<=leftborder_high+15;i++)
			{
				Sumenergy = Sumenergy + hAmplitude[0]->GetBinContent(i);
			}
			
			Energy_av = Sumenergy/11;
			
			Int_t k = 0;
			
			for(i=leftborder_low-10;i<leftborder_low;i++)
			{	
				if(hAmplitude[0]->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_left = Sumnoise_left + hAmplitude[0]->GetBinContent(i);
					k++;
				}
			}
									
			Noise_av = Sumnoise_left/k;
		
			Energy_av = Energy_av - Noise_av;
			if (Energy_av > threshold_MWD)
			{
				hEnergySpectrum->Fill(Energy_av);
				leftborder.push_back(leftborder_low);		//used for Risetime
				rightborder.push_back(leftborder_high);		//used for Risetime 
				energy.push_back(Energy_av);		//used for Energy-Risetime-Correlation 
			}
		}
	}
	while(abort1==0);

}

/*
for(Int_t i=0;i<=TraceLength;i++)
{
	if(hTraceBuffer->GetBinContent(i) > threshold_MWD)
	{
		Double_t leftborder_low = FindLocalMinimumBin(Grad2_MWD,Grad2_MWD->GetBin(i)-30,hTraceBuffer->GetBin(i)+30);
		Double_t leftborder_high = FindLocalMaximumBin(Grad2_MWD,Grad2_MWD->GetBin(i)-30,hTraceBuffer->GetBin(i)+30);

		Double_t leftborder_low_base = FindLocalMinimumBin(Grad2_MWD,leftborder_low-30,leftborder_low);
		Double_t leftborder_high_base = FindLocalMaximumBin(Grad2_MWD,leftborder_high,leftborder_high+30);

		Double_t rightborder_high = FindLocalMaximumBin(Grad2_MWD,leftborder_high,leftborder_high+M);
		Double_t rightborder_low = FindLocalMinimumBin(Grad2_MWD,rightborder_high,rightborder_high+30);
		
		Double_t rightborder_high_base = FindLocalMaximumBin(Grad2_MWD,rightborder_high-30,rightborder_high);
		Double_t rightborder_low_base = FindLocalMinimumBin(Grad2_MWD,rightborder_low,rightborder_low+30);
*/		
		
		
		

void THypGeMWD::EvaluateMWD()
{
	//for MWD
 
	Double_t threshold_MWD = 150;		//threshold value for energy of signals 
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high, rightborder_high;
	/*
	Double_t Grad3[TraceLength];
	Double_t Grad4[TraceLength];
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Grad3[i] = hTraceBuffer->GetBinContent(i+5) - hTraceBuffer->GetBinContent(i);
	}
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Grad4[i] = Grad3[i+5] - Grad3[i];
	}
	
	Int_t treshold = 10;
	border_max.clear();
	border_min.clear();

	Int_t i=1-M/2;
	do
	{
		i=i+M/2;

		do
		{
			i++;
		}
		while(abs(Grad4[i]) < treshold);

		border_max.push_back(FindLocalMaximumBin(Grad4,i,i+M/2));
		border_min.push_back(FindLocalMinimumBin(Grad4,i,i+M/2));
	}
	while(i<=TraceLength);
for(Int_t k=0;k<=border_max.size();k++)
{
	cout << "border_max" << k << " " << border_max[k] << endl;
}
for(Int_t k=0;k<=border_min.size();k++)
{
	cout << "border_min" << k << " " << border_min[k] << endl;
}

if(border_max.size() == border_min.size())
	{
		border_max[border_max.size()+1] = TraceLength;
		border_min[border_min.size()+1] = TraceLength;
		
		for(Int_t k=0;k<=border_max.size();k++)
		{
			
			if(border_max[k] < border_min[k] && border_max[k+1] > border_min[k+1])
			{
				Double_t Sumenergy = 0;
				Double_t Sumnoise = 0;
				for(Int_t l=border_min[k]+10;l<=border_min[k+1]-10;l++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_min[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-20;i<=border_max[k]-10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1]+10;i<=border_max[k+1]+20;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergySpectrum->Fill(Energy);
				}
			}
			if(border_max[k] < border_min[k] && border_max[k+1] < border_min[k+1] && border_max[k+2] > border_min[k+2] && border_max[k+3] > border_min[k+3])
			{
				
				Double_t Sumnoise = 0;
				Int_t p = 0;
				for(i=border_max[k]-20;i<=border_max[k]-10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+3]+10;i<=border_max[k+3]+20;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Double_t Sumenergy = 0;
				for(Int_t i=border_min[k]+10;i<=border_max[k+1]-10;i++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(i);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergySpectrum->Fill(Energy);
				}
				Sumenergy = 0;
				for(Int_t i=border_max[k+2]+10;i<=border_min[k+3]-10;i++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(i);
				}
				Energy = Sumenergy/(border_min[k+3]-border_max[k+2]);
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergySpectrum->Fill(Energy);
				}
			}
			
		}
	}
*/
	/*if(border_max.size() == border_min.size())
	{
		border_max[border_max.size()+1] = TraceLength;
		border_min[border_min.size()+1] = TraceLength;
		Double_t Sumnoise = 0;
		for(Int_t k=0;k<=border_max.size();k++)
		{
			if(border_max[k] < border_min[k] && border_max[k+1] > border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_min[k];l<=border_min[k+1];l++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_min[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-10;i<=border_max[k];i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1];i<=border_max[k+1]+10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergySpectrum->Fill(Energy);
			}
			if(border_max[k] > border_min[k] && border_max[k+1] < border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_max[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_max[k]);
				Int_t p = 0;
				for(i=border_min[k]-10;i<=border_min[k];i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_min[k+1];i<=border_min[k+1]+10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergySpectrum->Fill(Energy);
			}
			if(border_max[k] < border_min[k] && border_max[k+1] < border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_min[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-10;i<=border_max[k];i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_min[k+1];i<=border_min[k+1]+10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergySpectrum->Fill(Energy);
			}
			if(border_max[k] > border_min[k] && border_max[k+1] > border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_min[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTraceBuffer->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_min[k]-10;i<=border_min[k];i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1];i<=border_max[k+1]+10;i++)
				{
					if(hTraceBuffer->GetBinContent(i)<treshold)
					{
						Sumnoise += hTraceBuffer->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergySpectrum->Fill(Energy);
			}
		}
	}
	*/
	/*
	for(Int_t i=0;i<=TraceLength;i++)
	{
		hEnergySpectrum->SetBinContent(i,Grad4[i]);
	}
	*/
	
	Double_t rightborder_low = -9;
	Int_t abort1 = 0;
	leftborder.clear();
	rightborder.clear();
	SignalTime.clear();
	energy.clear();

	Double_t Sumenergy = 0;
	Double_t Sumnoise_left = 0;
	Double_t Sumnoise_right = 0;
	Double_t Energy_av = 0;
	Double_t Noise_av = 0;
		
	for(;;)			// loop is finished internally
	{
		Sumenergy = 0;
		Sumnoise_left = 0;
		Sumnoise_right = 0;
		Energy_av = 0;
		Noise_av = 0;
	
		Int_t i = rightborder_low+10;
		while(hMWD[0]->GetBinContent(i+1) - hMWD[0]->GetBinContent(i) < grad_MWD && hMWD[0]->GetBinContent(i+2) - hMWD[0]->GetBinContent(i+1) < grad_MWD && hMWD[0]->GetBinContent(i+3) - hMWD[0]->GetBinContent(i+2) < grad_MWD)	// check if derivative of trace is smaller than a threshold for 3 continous times, when this fails the start of a signal is found
		{
			if(i+10 >= TraceLength)		// break of the outer loop
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		if (abort1) break;
		
		leftborder_low = hMWD[0]->GetBin(i);			// start of MWD signal
						
		//leftborder.push_back(leftborder_low);		//used for Risetime (add +10 to leftborder_low at high smoothing rates -> results from the discrepance between smoothed and unsmoothed MWD-shape)
		//cout << leftborder_low << endl;
		i = leftborder_low+5;
		while(hMWD[0]->GetBinContent(i+1) - hMWD[0]->GetBinContent(i) >= grad_MWD && hMWD[0]->GetBinContent(i+2) - hMWD[0]->GetBinContent(i+1) >= grad_MWD && hMWD[0]->GetBinContent(i+3) - hMWD[0]->GetBinContent(i+2) >= grad_MWD) // check if derivative of trace is bigger than a threshold for 3 continous times, when this fails the top of a signal is found
		{
			if(i+10 >= TraceLength) // break of the outer loop
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		if (abort1) break;
		
		leftborder_high = hMWD[0]->GetBin(i);			// upper left edge of MWD signal
				
		i = leftborder_high+5;	// little offset just as a speedup
		
		while(hMWD[0]->GetBinContent(i+1) - hMWD[0]->GetBinContent(i) > -grad_MWD && hMWD[0]->GetBinContent(i+2) - hMWD[0]->GetBinContent(i+1) > -grad_MWD && hMWD[0]->GetBinContent(i+3) - hMWD[0]->GetBinContent(i+2) > -grad_MWD) // check if derivative of trace is bigger than a threshold for 3 continous times, when this fails the end of the flat top of a signal is found
		{
			if(i+10 >= TraceLength) // break of the outer loop
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		if (abort1) break;
	
		rightborder_high = hMWD[0]->GetBin(i);		// upper right edge of MWD signal
		
		i = rightborder_high+5;
		while(hMWD[0]->GetBinContent(i+1) - hMWD[0]->GetBinContent(i) <= -grad_MWD && hMWD[0]->GetBinContent(i+2) - hMWD[0]->GetBinContent(i+1) <= -grad_MWD && hMWD[0]->GetBinContent(i+3) - hMWD[0]->GetBinContent(i+2) <= -grad_MWD) // check if derivative of trace is smaller than a threshold for 3 continous times, when this fails the end of the MWD signal is found
		{
			if(i+10 >= TraceLength) // break of the outer loop
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		if (abort1) break;
	
		rightborder_low = hMWD[0]->GetBin(i);		// end of MWD signal
	
		if(hMWD[0]->GetBinContent(leftborder_low) < threshold_MWD && hMWD[0]->GetBinContent(rightborder_low) < threshold_MWD && hMWD[0]->GetBinContent(leftborder_high) > threshold_MWD && hMWD[0]->GetBinContent(rightborder_high) > threshold_MWD)
		{
			
			for(i=leftborder_high+10;i<=rightborder_high-10;i++)			// calculate average value of flat top
			{
				Sumenergy = Sumenergy + hMWD[0]->GetBinContent(i);
			}
			
			Energy_av = Sumenergy/(rightborder_high-leftborder_high-20);
			
			Int_t k = 0;
			
			for(i=leftborder_low-10;i<leftborder_low;i++)							// calculate average value of base before pulse
			{	
				if(hMWD[0]->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_left = Sumnoise_left + hMWD[0]->GetBinContent(i);	
					k++;
				}
			}
			
			for(i=rightborder_low+1;i<=rightborder_low+10;i++)				// calculate average value of base after pulse
			{	
				if(hMWD[0]->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_right = Sumnoise_right + hMWD[0]->GetBinContent(i);
					k++;
				}
			}
			
			Noise_av = (Sumnoise_left + Sumnoise_right)/k;					// average of complete noise (left + right)/NoAllBaselinePoints
			
			Energy_av = Energy_av - Noise_av; //+ M*mCorrection/50;			// added correction for inclined amplitude signal --> most likely not working
			
			
			
			if (Energy_av > threshold_MWD)
			{
				hEnergySpectrum->Fill(Energy_av);
				energy.push_back(Energy_av);		//used for Energy-Risetime-Correlation 
				leftborder.push_back(leftborder_low);		//used for Risetime
				rightborder.push_back(leftborder_high);		//used for Risetime 
				SignalTime.push_back((leftborder_low+leftborder_high)/2/100);	// used for Energy - time before pulse Correlation, time value in µs
				
			}
		}
	}
	
}

void THypGeMWD::SetUseMWD(Bool_t useMWD_ext)
{
	useMWD = useMWD_ext;
}

void THypGeMWD::CalculateGausCoeff()
{
	/**
	* Function to calculate the coefficients of the gaussian smoothing filter. This values are also used for the first gaussian of the bilateral filter.
	*/
	Double_t limit = 0.999;
	
	GausBreakUp = 0;
	Double_t Norm =sqrt(2*TMath::Pi());
	Double_t GaussInt = 0;
	
	do{
		g[GausBreakUp] = (TMath::Exp(-pow((Double_t(GausBreakUp)/(sqrt(2)*Sigma)),2)))/Double_t(Sigma)/Norm;
		if(GausBreakUp==0)
			GaussInt += g[GausBreakUp];
		else
			GaussInt += 2*g[GausBreakUp];
		GausBreakUp++;
		//cout << g[GausBreakUp] <<endl;
	}
	while( (GaussInt < limit ) );
	
	// the norm for every point should be calculated in order to prevent effects at the beginning and the end of the trace
	
	for (Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
	{
		GausNorm[i-1]=0;
		for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)
		{
			if (j <1 )			// condition to fit the algorithm to the beginning of the trace
				continue;
			if (j > TraceLength) 		// condition to fit the algorithm to the end of the trace
				continue;
			GausNorm[i-1] += g[abs(j-i)];
		}
	}
	cout << "coefficients of gausian filter (re)calculated" << endl;
}

Double_t THypGeMWD::EnergyRtCorrection(Double_t Rt, Double_t EnergyUncorr )
{
	Double_t EnergyCorr = EnergyUncorr / EnergyRtCorrFunc->Eval(Rt);
	return EnergyCorr;
}


Double_t	THypGeMWD::EnergyPileUpTimeCorrection(Double_t PileUpTime, Double_t EnergyUncorr)
{
	// energy correction for high rates (close signals)
	Double_t EnergyCorr = EnergyUncorr / EnergyPileUpTimeCorrFunc->Eval(PileUpTime);
	return EnergyCorr;
}


Double_t THypGeMWD::Gaus(Double_t x)				// function used by bilateral filter to calculate the value of den second gausian
{
	Double_t Norm =sqrt(2*TMath::Pi());
	return (TMath::Exp(-pow((x/(sqrt(2)*Sigma2)),2)))/Double_t(Sigma2)/Norm;
}

void THypGeMWD::DoMeanFilter()
{
	//Rectangular filter
	for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Double_t Content = 0;
			if(i <= Width)
			{
				for(Int_t j=1-i;j<=Width;j++)
				{
					Content += hTraceBuffer->GetBinContent(i+j);
				}
				Content = Content/(Width+i);
			}
			else
			{
				if(i >= TraceLength-Width)
				{
					for(Int_t j=-Width;j<=TraceLength-i;j++)
					{
						Content += hTraceBuffer->GetBinContent(i+j);
					}
					Content = Content/(Width+TraceLength-i+1);
				}	
				else
				{
					for(Int_t j=-Width;j<=Width;j++)
					{
						Content += hTraceBuffer->GetBinContent(i+j);
					}
					Content = Content/(2*Width+1);
				}
			}
			hSmoothedTrace[0]->SetBinContent(i,Content);
		}
	}
}
void THypGeMWD::DoWeightedAverageFilter()
{
	// weighted average
	for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)	
		{
			if (i == 1)	
				hSmoothedTrace[0]->SetBinContent(i,Double_t(hTraceBuffer->GetBinContent(i))/2 + Double_t(hTraceBuffer->GetBinContent(i+1))/2);
			else
				if (i == TraceLength)
					hSmoothedTrace[0]->SetBinContent(i,Double_t(hTraceBuffer->GetBinContent(i-1))/2 + Double_t(hTraceBuffer->GetBinContent(i))/2);
				else
					hSmoothedTrace[0]->SetBinContent(i,Double_t(hTraceBuffer->GetBinContent(i-1))/4 + Double_t(hTraceBuffer->GetBinContent(i))/2 + Double_t(hTraceBuffer->GetBinContent(i+1))/4);
		}
	}
}
void THypGeMWD::DoGaussianFilter()
{
	Double_t value;

	for(Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
	{	
		value = 0;
		for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)			// smoothing loop
		{
			if (j <1  || j > TraceLength)			// condition to fit the algorithm to the beginning and the end of the trace
				continue;
			value += g[abs(j-i)] * hTraceBuffer->GetBinContent (j);			// real gaussian smoothing
			// since we don't integrate over the whole gaussian we have to renormate the value, this has to be calcutated only once --> GausNorm[] in CalculateGausCoeff
			//cout << g[j] << endl;
		}
		value = value/GausNorm[i-1];							// renormalization of the value, because not the whole Gaus is included (GausNorm < 1)
		
		//cout << value << "\t" << GausNorm[i-1] << endl;
		hSmoothedTrace[0]->SetBinContent(i,value);
	}
}
void THypGeMWD::DoBilateralFilter()
{
// bilateral
	Double_t value;
	Double_t norm;
	for (Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
	{
		value = 0;
		norm = 0;
		for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)			// smoothing loop
		{
			if (j <1  || j > TraceLength)			// condition to fit the algorithm to the beginning and the end of the trace
				continue;
			value += g[abs(j-i)]* Gaus(hTraceBuffer->GetBinContent (i)-hTraceBuffer->GetBinContent (j)) * hTraceBuffer->GetBinContent (j);				// bilateral (2d gaussian) smoothing
			norm += g[abs(j-i)]* Gaus(hTraceBuffer->GetBinContent (i)-hTraceBuffer->GetBinContent (j));				// here we have to recalculate the norm every time because of the non-linearity of the bilateral filter (norm is not const for position i)
			
		}
		
		value = value/norm;				
		//cout << value << "\t" << norm << endl;
		hSmoothedTrace[0]->SetBinContent(i,value);
	}
}

void THypGeMWD::DoFourierTransformation()
{
	
}
void THypGeMWD::DoBandStopFilter()
{
	
}
void THypGeMWD::DoFourierBackTransformation()
{
	
}

Int_t THypGeMWD::AnaStep_EnergyTimeSinceLastPulseCorrelation()
{
	Double_t PUTime;
	for (Int_t i =1; i <Int_t(SignalTime.size()); i++)
	{
		PUTime= SignalTime[i]- SignalTime[i-1];
		hEnergyTimeSinceLastPulse->Fill(PUTime,energy[i]);
		hEnergyTimeSinceLastPulseCorr->Fill(PUTime,EnergyPileUpTimeCorrection(PUTime, energy[i]));
		for(Int_t j = 0; j < NumberOfPileUpTimeHistograms; j++)
		{
			if (energy[i-1] > j*100 && energy[i-1] < (j+1)*100-1 )
				hEnergyTimeSinceLastPulse_WithCuts[j]->Fill(PUTime,energy[i]);
		}
	}
	return 0;
}

void THypGeMWD::ConnectTraceHistograms(TH1D** hTrace_ext, TH1D** hSmoothedTrace_ext, TH1D** hTrace_bc_ext, TH1D** hAmplitude_ext,TH1D** hMWD_ext,TH1D** hMWDMA_ext,TH1D** hTrace_Direct_ext)
{
	
		hTrace = hTrace_ext; 
		hSmoothedTrace = hSmoothedTrace_ext;
		hTrace_bc = hTrace_bc_ext;
		hAmplitude = hAmplitude_ext;
		hMWD = hMWD_ext;
		hMWDMA = hMWDMA_ext;
		hTrace_Direct = hTrace_Direct_ext;
}
void THypGeMWD::Connect1DEnergySpectraHistograms(TH1D *hEnergySpectrum_ext,TH1D *hEnergySpectrumWithCut_ext)
{
	hEnergySpectrum = hEnergySpectrum_ext;
	hEnergySpectrumWithCut = hEnergySpectrumWithCut_ext;
}
void THypGeMWD::Connect1DRisetimeHistograms(TH1D* hRisetime1090_ext, TH1D* hRisetime3090_ext)
{
	hRisetime1090 = hRisetime1090_ext;
	hRisetime3090 = hRisetime3090_ext;
}
void THypGeMWD::Connect2DEnergyRisetimeHistograms(TH2D* hEnergyRise1090Corr_ext, TH2D* hEnergyRise3090Corr_ext)
{
	hEnergyRise1090Corr = hEnergyRise1090Corr_ext;
	hEnergyRise3090Corr = hEnergyRise3090Corr_ext;
}
void THypGeMWD::Connect2DEnergyTimeSinceLastPulseHistograms(TH2D* hEnergyTimeSinceLastPulse_ext, TH2D* hEnergyTimeSinceLastPulseCorr_ext ,TH2D** hEnergyTimeSinceLastPulse_WithCuts_ext, TH2D** hEnergyTimeSinceLastPulseCorr_WithCuts_ext, Int_t NumberOfCuts)
{
	NumberOfPileUpTimeHistograms = NumberOfCuts;
	hEnergyTimeSinceLastPulse = hEnergyTimeSinceLastPulse_ext;
	hEnergyTimeSinceLastPulseCorr = hEnergyTimeSinceLastPulseCorr_ext;
	hEnergyTimeSinceLastPulse_WithCuts = hEnergyTimeSinceLastPulse_WithCuts_ext;
	hEnergyTimeSinceLastPulseCorr_WithCuts = hEnergyTimeSinceLastPulseCorr_WithCuts_ext;
	//for (int i = 0; i < NumberOfCuts;i++)
	//{
		//cout << hEnergyTimeSinceLastPulse_WithCuts[i] << "\t\t" << hEnergyTimeSinceLastPulse_WithCuts_ext[i] << endl;
	//}
}

Int_t		THypGeMWD::AnaStep_FillEnergySpectrumWithPileUpCut(Double_t CutValue)
{
	for (Int_t i =1; i <Int_t(SignalTime.size()); i++)
	{
		if (SignalTime[i]- SignalTime[i-1] > CutValue )
			hEnergySpectrumWithCut->Fill(energy[i]);
		//cout << i << "\t" <<energy[i]<< "\t" << SignalTime[i]- SignalTime[i-1] << endl;
	}
	return 0;
}