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

#include "TGo4Analysis.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <cassert>

using namespace std;

//default constructor DON'T use!
THypGeMWD::THypGeMWD()
{
	THypGeMWD(16384,1);
}

THypGeMWD::THypGeMWD(Int_t TraceLength_ext,Int_t NumberOfChannels_ext)            //constructor
{
	TraceLength=TraceLength_ext;
	NumberOfChannels = NumberOfChannels_ext;
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

	Aarray = new Double_t* [NumberOfChannels];
	MWDarray = new Double_t* [NumberOfChannels];
	Derivative1array = new Double_t* [NumberOfChannels];
	Derivative2array = new Double_t* [NumberOfChannels];
	Derivative3array = new Double_t* [NumberOfChannels];
	Derivative4array = new Double_t* [NumberOfChannels];
	MWDMAarray = new Double_t* [NumberOfChannels];
	Sarray = new Double_t* [NumberOfChannels];		// array to store the result of the direct filter
	Parray = new Double_t* [NumberOfChannels];		// array to store an intermediate result of the direct filter

	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		Aarray[ChanNumber] = new Double_t[TraceLength+1];
		MWDarray[ChanNumber] = new Double_t[TraceLength+1];
		Derivative1array[ChanNumber] = new Double_t[TraceLength+1];
		Derivative2array[ChanNumber] = new Double_t[TraceLength+1];
		Derivative3array[ChanNumber] = new Double_t[TraceLength+1];
		Derivative4array[ChanNumber] = new Double_t[TraceLength+1];
		MWDMAarray[ChanNumber] = new Double_t[TraceLength+1];
		Sarray[ChanNumber] = new Double_t[TraceLength+1];		// array to store the result of the direct filter
		Parray[ChanNumber] = new Double_t[TraceLength+1];		// array to store an intermediate result of the direct filter
	}
	useMWD= 1;
	GausNorm = new Double_t[TraceLength];
	CalculateGausCoeff();
	
	EnergyRtCorrFunc = new TF1("EnergyRtCorrFunc","pol4",0,2000);
		EnergyRtCorrFunc->SetParameters(1,-9.1237e-05, 7.95504e-07, -3.0778e-09, 4.39764e-12);
		
	EnergyPileUpTimeCorrFunc = new TF1("EnergyPileUpTimeCorrFunc","[0]*(1-[1]/pow(x,[2]))",0,2000);
		EnergyPileUpTimeCorrFunc->SetParameters(1,3.74495e-02,1.95113);
	//CalculateSecondGausCoeff();
	cout << "Analysis created"<<endl;
}

THypGeMWD::~THypGeMWD()            //destructor
{
	cout << "HypGeMWD deconstructed!" << endl;
}

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
	//cout << "finished direct filter"<< endl;
	return 0;
}
	

	//Smoothing
	
Int_t THypGeMWD::AnaStep_Smoothing()
{
	// hSmoothedTrace will be filled be the filter functions
	if( SmoothingMethod == 1)
		DoMeanFilter();
	if( SmoothingMethod == 2)
		DoWeightedAverageFilter();
	if( SmoothingMethod == 3)
		DoGaussianFilter();
	if( SmoothingMethod == 4)	
		DoBilateralFilter();
		
	return 0;
}
	
	//Baseline
	
Int_t THypGeMWD::AnaStep_BaselineCorrection()
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
//		if (hSmoothedTrace[ChanNumber]->GetMaximum() == 0)  	// changed in Jülich, trace must be positive to avoid a strong negative drift of the deconvoluted signal!, remaining offset is corrected after deconvolution
//			return -1;
//		Double_t sumoffset = 0;
//		Int_t BaseLineCorrectionBins = 20;//100;		//building average of first <value> bins
//
//		for(Int_t i=1;i<=BaseLineCorrectionBins;i++)
//		{
//			sumoffset = sumoffset + hSmoothedTrace[ChanNumber]->GetBinContent(i);
//		}
//		offset_av= sumoffset/BaseLineCorrectionBins;

		offset_av=hSmoothedTrace[ChanNumber]->GetMinimum();			// changed in Jülich, trace must be positive to avoid a strong negative drift of the deconvoluted signal!, remaining offset is corrected after deconvolution
		
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hTrace_bc[ChanNumber]->SetBinContent(i,hSmoothedTrace[ChanNumber]->GetBinContent(i) - offset_av);
		}
	}
	return 0;
}
	
	//Moving-Window-Deconvolution
	
Int_t THypGeMWD::AnaStep_DoMovingWindowDeconvolution()
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		Aarray[ChanNumber][0] = hTrace_bc[ChanNumber]->GetBinContent(1);			// first value of baseline corrected trace
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Aarray[ChanNumber][i] = hTrace_bc[ChanNumber]->GetBinContent(i) - hTrace_bc[ChanNumber]->GetBinContent(i-1) * (1.-(1./tau)) + Aarray[ChanNumber][i-1];
		}
		//bc correction of amp signal
		// the idea here is, that after the deconvolution, the lowest point is right at the beginning and this is than corrected to 0
		Aarray[ChanNumber][0]=0;
		Double_t SumAmpOffset=0;
		Int_t LengthSumOffsetAverage=20;
		for(Int_t i = 1; i <= LengthSumOffsetAverage; i++)
			SumAmpOffset+= Aarray[ChanNumber][i];
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
			Aarray[ChanNumber][i]= Aarray[ChanNumber][i]-SumAmpOffset; //+ mCorrection*i;						// bc correction of amplitude signal
			if (i > Int_t(M))
				MWDarray[ChanNumber][i] = Aarray[ChanNumber][i] - Aarray[ChanNumber][i-Int_t(M)];
			else
				MWDarray[ChanNumber][i] = Aarray[ChanNumber][i]-Aarray[ChanNumber][1];				// added to cover edge efects
			hAmplitude[ChanNumber]->SetBinContent(i,Aarray[ChanNumber][i]);
			hMWD[ChanNumber]->SetBinContent(i,MWDarray[ChanNumber][i]);
		}
	}
	return 0;
}
	
	//Moving-Average
	
Int_t THypGeMWD::AnaStep_DoMovingAverageFilter()
{
	//cout << L << endl;
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		MWDMAarray[ChanNumber][0] = hMWD[ChanNumber]->GetBinContent(1);
		for(Int_t n=1;n<=TraceLength;n++)
		{
			if (n > Int_t(L))
				MWDMAarray[ChanNumber][n] = MWDMAarray[ChanNumber][n-1] + 1./L * (hMWD[ChanNumber]->GetBinContent(n) - hMWD[ChanNumber]->GetBinContent(n-Int_t(L)));
			else
				MWDMAarray[ChanNumber][n] = MWDMAarray[ChanNumber][n-1] + 1./L * (hMWD[ChanNumber]->GetBinContent(n) - hMWD[ChanNumber]->GetBinContent(1));

			hMWDMA[ChanNumber]->SetBinContent(n,MWDMAarray[ChanNumber][n]);
		}
	}

	return 0;
}
Int_t THypGeMWD::AnaStep_DoDirectFilter()
{
	// this filter replaces MWD + MA (other version of it, after Jordanov/Knoll), do this step aber smoothing and baseline correction
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		Double_t d_kl;		// variable to calculate a needed value from the trace
		Sarray[ChanNumber][0]= 0;
		Parray[ChanNumber][0] = 0;

		for(Int_t n = 1;n <= TraceLength;n++)
		{
			d_kl=hTrace_bc[ChanNumber]->GetBinContent(n);
			if ( n-L > 0)
				d_kl = d_kl - hTrace_bc[ChanNumber]->GetBinContent(n-L);
			if (n-M >0)
				d_kl = d_kl - hTrace_bc[ChanNumber]->GetBinContent(n-M);
			if (n-M-L > 0)
				d_kl += hTrace_bc[ChanNumber]->GetBinContent(n-M-L);

			Parray[ChanNumber][n]=Parray[ChanNumber][n-1]+d_kl;
			Sarray[ChanNumber][n]=Sarray[ChanNumber][n-1]+Parray[ChanNumber][n]+ tau* d_kl;
			hTrace_Direct[ChanNumber]->SetBinContent(n,Sarray[ChanNumber][n]);
		}
	}
	return 0;
}	
	//PileupCompensation and Energyspectrum
		
Int_t THypGeMWD::AnaStep_FillEnergyspectrum()
{
	if (useMWD)
	{
		EvaluateMWD();
	}
	else
	{
		EvaluateAmplitude();
	}
	if (EnableMA)
		EvaluateMA();
	return 0;
}
	
	//Risetime
	
Int_t THypGeMWD::AnaStep_ExtractRisetime()
{
	//cout << "new Trace" << endl;
		

	Double_t RiseX1090, RiseX3090, RiseT1090, RiseT3090;
	Double_t threshold_Risetime = 50;		//threshold up where signals were identified as useful for Risetime calculation

	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		risetime1090[ChanNumber].clear();
		risetime3090[ChanNumber].clear();
		for(Int_t k=0;k<Int_t(leftborder[ChanNumber].size());k++)
		{
			Double_t RightEdge = FindLocalMaximumBin(hTrace[ChanNumber],leftborder[ChanNumber][k],leftborder[ChanNumber][k]+40);

			if(hTrace[ChanNumber]->GetBinContent(RightEdge) > threshold_Risetime+hTrace[ChanNumber]->GetBinContent(leftborder[ChanNumber][k]))
			{
				Double_t l = (hTrace[ChanNumber]->GetBinContent(RightEdge)-hTrace[ChanNumber]->GetBinContent(leftborder[ChanNumber][k]))/10;		//Calculates 10% of height of the rising signal

				RiseX1090 = Double_t(FindFirstBinAbove(hTrace[ChanNumber],hTrace[ChanNumber]->GetBinContent(RightEdge)-l,leftborder[ChanNumber][k],RightEdge)) - Double_t(FindFirstBinAbove(hTrace[ChanNumber],hTrace[ChanNumber]->GetBinContent(RightEdge)-9*l,leftborder[ChanNumber][k],RightEdge));
				RiseX3090 = Double_t(FindFirstBinAbove(hTrace[ChanNumber],hTrace[ChanNumber]->GetBinContent(RightEdge)-l,leftborder[ChanNumber][k],RightEdge)) - Double_t(FindFirstBinAbove(hTrace[ChanNumber],hTrace[ChanNumber]->GetBinContent(RightEdge)-7*l,leftborder[ChanNumber][k],RightEdge));
				
				RiseT1090 = RiseX1090 * 10;		//Conversion into nanoseconds (FADC 100 MS/s)
				RiseT3090 = RiseX3090 * 10;

				hRisetime1090[ChanNumber]->Fill(RiseT1090);
				hRisetime3090[ChanNumber]->Fill(RiseT3090);
			
				risetime1090[ChanNumber].push_back(RiseT1090);		//used for Energy-Risetime-Correlation
				risetime3090[ChanNumber].push_back(RiseT3090);		//used for Energy-Risetime-Correlation
			}
		}
	}
	return 0;
}
	
	//Energy-Risetime-Correlation
	
Int_t THypGeMWD::AnaStep_DoEnergyRisetimeCorrelation()
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		if(Int_t(energy[ChanNumber].size()) == Int_t(risetime1090[ChanNumber].size()))
		{
			for(Int_t n=0;n < Int_t(energy[ChanNumber].size());n++)
			{
				hEnergyRise1090Corr[ChanNumber]->Fill(risetime1090[ChanNumber][n],energy[ChanNumber][n]);
			}
		}

		if(Int_t(energy[ChanNumber].size()) == Int_t(risetime3090[ChanNumber].size()))
		{
			for(Int_t n=0;n < Int_t(energy[ChanNumber].size());n++)
			{
				hEnergyRise3090Corr[ChanNumber]->Fill(risetime3090[ChanNumber][n],energy[ChanNumber][n]);
			}
		}
	}
	return 0; 
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
	//cout << "new parameters set"<< endl;
}

void THypGeMWD::EvaluateAmplitude()
{
	
	//for MWD without MA
	// 2.5.14 still old do ... while style
	Double_t threshold_MWD = 50;		//threshold up where signals were identified as useful for MWD
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high;
	Int_t abort1 = 0;



	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		CalculateDerivatives(hAmplitude[ChanNumber],ChanNumber);
		leftborder[ChanNumber].clear();
		energy[ChanNumber].clear();
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
			while(hAmplitude[ChanNumber]->GetBinContent(i+1) - hAmplitude[ChanNumber]->GetBinContent(i) < grad_MWD && hAmplitude[ChanNumber]->GetBinContent(i+2) - hAmplitude[ChanNumber]->GetBinContent(i+1) < grad_MWD && hAmplitude[ChanNumber]->GetBinContent(i+3) - hAmplitude[ChanNumber]->GetBinContent(i+2) < grad_MWD);
			
			leftborder_low = hAmplitude[ChanNumber]->GetBin(i);

			//leftborder[ChanNumber].push_back(leftborder_low);		//used for Risetime (add +10 to leftborder_low at high smoothing rates -> results from the discrepance between smoothed and unsmoothed MWD-shape)

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
			while(hAmplitude[ChanNumber]->GetBinContent(i+1) - hAmplitude[ChanNumber]->GetBinContent(i) >= grad_MWD && hAmplitude[ChanNumber]->GetBinContent(i+2) - hAmplitude[ChanNumber]->GetBinContent(i+1) >= grad_MWD && hAmplitude[ChanNumber]->GetBinContent(i+3) - hAmplitude[ChanNumber]->GetBinContent(i+2) >= grad_MWD);
			
			leftborder_high = hAmplitude[ChanNumber]->GetBin(i);


			if(hAmplitude[ChanNumber]->GetBinContent(leftborder_low) < threshold_MWD && hAmplitude[ChanNumber]->GetBinContent(leftborder_high) > threshold_MWD)
			{

				for(i=leftborder_high+5;i<=leftborder_high+15;i++)
				{
					Sumenergy = Sumenergy + hAmplitude[ChanNumber]->GetBinContent(i);
				}

				Energy_av = Sumenergy/11;

				Int_t k = 0;

				for(i=leftborder_low-10;i<leftborder_low;i++)
				{
					if(hAmplitude[ChanNumber]->GetBinContent(i) < threshold_MWD)
					{
						Sumnoise_left = Sumnoise_left + hAmplitude[ChanNumber]->GetBinContent(i);
						k++;
					}
				}

				Noise_av = Sumnoise_left/k;
			
				Energy_av = Energy_av - Noise_av;
				if (Energy_av > threshold_MWD)
				{
					hEnergySpectrum[ChanNumber]->Fill(Energy_av);
					leftborder[ChanNumber].push_back(leftborder_low);		//used for Risetime
					rightborder[ChanNumber].push_back(leftborder_high);		//used for Risetime
					energy[ChanNumber].push_back(Energy_av);		//used for Energy-Risetime-Correlation
				}
			}
		}
		while(abort1==0);
	}

}
void THypGeMWD::EvaluateMWD()
{
	//for MWD

	Double_t threshold_MWD = 150;		//threshold value for energy[ChanNumber] of signals
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high, rightborder_high;
	Double_t rightborder_low = -9;
	Int_t abort1 = 0;

	Double_t Sumenergy = 0;
	Double_t Sumnoise_left = 0;
	Double_t Sumnoise_right = 0;
	Double_t Energy_av = 0;
	Double_t Noise_av = 0;

	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		CalculateDerivatives(hMWD[ChanNumber],ChanNumber);
		leftborder[ChanNumber].clear();
		rightborder[ChanNumber].clear();
		SignalTime[ChanNumber].clear();
		energy[ChanNumber].clear();
		for(;;)			// loop is finished internally
		{
			Sumenergy = 0;
			Sumnoise_left = 0;
			Sumnoise_right = 0;
			Energy_av = 0;
			Noise_av = 0;
		
			Int_t i = rightborder_low+10;
			while(hMWD[ChanNumber]->GetBinContent(i+1) - hMWD[ChanNumber]->GetBinContent(i) < grad_MWD && hMWD[ChanNumber]->GetBinContent(i+2) - hMWD[ChanNumber]->GetBinContent(i+1) < grad_MWD && hMWD[ChanNumber]->GetBinContent(i+3) - hMWD[ChanNumber]->GetBinContent(i+2) < grad_MWD)	// check if derivative of trace is smaller than a threshold for 3 continous times, when this fails the start of a signal is found
			{
				if(i+10 >= TraceLength)		// break of the outer loop
				{
					abort1 = 1;
					break;
				}
				i++;
			}
			if (abort1) break;

			leftborder_low = hMWD[ChanNumber]->GetBin(i);			// start of MWD signal

			//leftborder[ChanNumber].push_back(leftborder_low);		//used for Risetime (add +10 to leftborder_low at high smoothing rates -> results from the discrepance between smoothed and unsmoothed MWD-shape)
			//cout << leftborder_low << endl;
			i = leftborder_low+5;
			while(hMWD[ChanNumber]->GetBinContent(i+1) - hMWD[ChanNumber]->GetBinContent(i) >= grad_MWD && hMWD[ChanNumber]->GetBinContent(i+2) - hMWD[ChanNumber]->GetBinContent(i+1) >= grad_MWD && hMWD[ChanNumber]->GetBinContent(i+3) - hMWD[ChanNumber]->GetBinContent(i+2) >= grad_MWD) // check if derivative of trace is bigger than a threshold for 3 continous times, when this fails the top of a signal is found
			{
				if(i+10 >= TraceLength) // break of the outer loop
				{
					abort1 = 1;
					break;
				}
				i++;
			}
			if (abort1) break;

			leftborder_high = hMWD[ChanNumber]->GetBin(i);			// upper left edge of MWD signal

			i = leftborder_high+5;	// little offset just as a speedup
			
			while(hMWD[ChanNumber]->GetBinContent(i+1) - hMWD[ChanNumber]->GetBinContent(i) > -grad_MWD && hMWD[ChanNumber]->GetBinContent(i+2) - hMWD[ChanNumber]->GetBinContent(i+1) > -grad_MWD && hMWD[ChanNumber]->GetBinContent(i+3) - hMWD[ChanNumber]->GetBinContent(i+2) > -grad_MWD) // check if derivative of trace is bigger than a threshold for 3 continous times, when this fails the end of the flat top of a signal is found
			{
				if(i+10 >= TraceLength) // break of the outer loop
				{
					abort1 = 1;
					break;
				}
				i++;
			}
			if (abort1) break;

			rightborder_high = hMWD[ChanNumber]->GetBin(i);		// upper right edge of MWD signal

			i = rightborder_high+5;
			while(hMWD[ChanNumber]->GetBinContent(i+1) - hMWD[ChanNumber]->GetBinContent(i) <= -grad_MWD && hMWD[ChanNumber]->GetBinContent(i+2) - hMWD[ChanNumber]->GetBinContent(i+1) <= -grad_MWD && hMWD[ChanNumber]->GetBinContent(i+3) - hMWD[ChanNumber]->GetBinContent(i+2) <= -grad_MWD) // check if derivative of trace is smaller than a threshold for 3 continous times, when this fails the end of the MWD signal is found
			{
				if(i+10 >= TraceLength) // break of the outer loop
				{
					abort1 = 1;
					break;
				}
				i++;
			}
			if (abort1) break;

			rightborder_low = hMWD[ChanNumber]->GetBin(i);		// end of MWD signal

			if(hMWD[ChanNumber]->GetBinContent(leftborder_low) < threshold_MWD && hMWD[ChanNumber]->GetBinContent(rightborder_low) < threshold_MWD && hMWD[ChanNumber]->GetBinContent(leftborder_high) > threshold_MWD && hMWD[ChanNumber]->GetBinContent(rightborder_high) > threshold_MWD)
			{

				for(i=leftborder_high+10;i<=rightborder_high-10;i++)			// calculate average value of flat top
				{
					Sumenergy = Sumenergy + hMWD[ChanNumber]->GetBinContent(i);
				}

				Energy_av = Sumenergy/(rightborder_high-leftborder_high-20);

				Int_t k = 0;

				for(i=leftborder_low-10;i<leftborder_low;i++)							// calculate average value of base before pulse
				{
					if(hMWD[ChanNumber]->GetBinContent(i) < threshold_MWD)
					{
						Sumnoise_left = Sumnoise_left + hMWD[ChanNumber]->GetBinContent(i);
						k++;
					}
				}

				for(i=rightborder_low+1;i<=rightborder_low+10;i++)				// calculate average value of base after pulse
				{
					if(hMWD[ChanNumber]->GetBinContent(i) < threshold_MWD)
					{
						Sumnoise_right = Sumnoise_right + hMWD[ChanNumber]->GetBinContent(i);
						k++;
					}
				}
				
				Noise_av = (Sumnoise_left + Sumnoise_right)/k;					// average of complete noise (left + right)/NoAllBaselinePoints

				Energy_av = Energy_av - Noise_av; //+ M*mCorrection/50;			// added correction for inclined amplitude signal --> most likely not working



				if (Energy_av > threshold_MWD)
				{
					hEnergySpectrum[ChanNumber]->Fill(Energy_av);
					energy[ChanNumber].push_back(Energy_av);		//used for Energy-Risetime-Correlation
					leftborder[ChanNumber].push_back(leftborder_low);		//used for Risetime
					rightborder[ChanNumber].push_back(leftborder_high);		//used for Risetime
					SignalTime[ChanNumber].push_back((leftborder_low+leftborder_high)/2/100);	// used for Energy - time before pulse Correlation, time value in µs

				}
			}
		}
	}
}

void THypGeMWD::EvaluateMA()
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
			CalculateDerivatives(hMWDMA[ChanNumber],ChanNumber);
			for(Int_t i=1;i<=TraceLength;i++)
			{
				hTraceDeri1 ->SetBinContent(i, Derivative1array[ChanNumber][i]*100);
				hTraceDeri2 ->SetBinContent(i, Derivative2array[ChanNumber][i]*1000);
				hTraceDeri3 ->SetBinContent(i, Derivative3array[ChanNumber][i]*10000);
				hTraceDeri4 ->SetBinContent(i, Derivative4array[ChanNumber][i]*100000);

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
	// energy[ChanNumber] correction for high rates (close signals)
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
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
		{
			for(Int_t i=1;i<=TraceLength;i++)
			{
				Double_t Content = 0;
				if(i <= Width)
				{
					for(Int_t j=1-i;j<=Width;j++)
					{
						Content += hTrace[ChanNumber]->GetBinContent(i+j);
					}
					Content = Content/(Width+i);
				}
				else
				{
					if(i >= TraceLength-Width)
					{
						for(Int_t j=-Width;j<=TraceLength-i;j++)
						{
							Content += hTrace[ChanNumber]->GetBinContent(i+j);
						}
						Content = Content/(Width+TraceLength-i+1);
					}
					else
					{
						for(Int_t j=-Width;j<=Width;j++)
						{
							Content += hTrace[ChanNumber]->GetBinContent(i+j);
						}
						Content = Content/(2*Width+1);
					}
				}
				hSmoothedTrace[0]->SetBinContent(i,Content);
			}
		}
	}
}
void THypGeMWD::DoWeightedAverageFilter()
{
	// weighted average
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
		{
			for(Int_t i=1;i<=TraceLength;i++)
			{
				if (i == 1)
					hSmoothedTrace[ChanNumber]->SetBinContent(i,Double_t(hTrace[ChanNumber]->GetBinContent(i))/2 + Double_t(hTrace[ChanNumber]->GetBinContent(i+1))/2);
				else
					if (i == TraceLength)
						hSmoothedTrace[ChanNumber]->SetBinContent(i,Double_t(hTrace[ChanNumber]->GetBinContent(i-1))/2 + Double_t(hTrace[ChanNumber]->GetBinContent(i))/2);
					else
						hSmoothedTrace[ChanNumber]->SetBinContent(i,Double_t(hTrace[ChanNumber]->GetBinContent(i-1))/4 + Double_t(hTrace[ChanNumber]->GetBinContent(i))/2 + Double_t(hTrace[ChanNumber]->GetBinContent(i+1))/4);
			}
		}
	}
}
void THypGeMWD::DoGaussianFilter()
{
	Double_t value;
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for(Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
		{
			value = 0;
			for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)			// smoothing loop
			{
				if (j <1  || j > TraceLength)			// condition to fit the algorithm to the beginning and the end of the trace
					continue;
				value += g[abs(j-i)] * hTrace[ChanNumber]->GetBinContent (j);			// real gaussian smoothing
				// since we don't integrate over the whole gaussian we have to renormate the value, this has to be calcutated only once --> GausNorm[] in CalculateGausCoeff
				//cout << g[j] << endl;
			}
			value = value/GausNorm[i-1];							// renormalization of the value, because not the whole Gaus is included (GausNorm < 1)

			//cout << value << "\t" << GausNorm[i-1] << endl;
			hSmoothedTrace[ChanNumber]->SetBinContent(i,value);
		}
	}
}
void THypGeMWD::DoBilateralFilter()
{
// bilateral
	Double_t value;
	Double_t norm;
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for (Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
		{
			value = 0;
			norm = 0;
			for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)			// smoothing loop
			{
				if (j <1  || j > TraceLength)			// condition to fit the algorithm to the beginning and the end of the trace
					continue;
				value += g[abs(j-i)]* Gaus(hTrace[ChanNumber]->GetBinContent (i)-hTrace[ChanNumber]->GetBinContent (j)) * hTrace[ChanNumber]->GetBinContent (j);				// bilateral (2d gaussian) smoothing
				norm += g[abs(j-i)]* Gaus(hTrace[ChanNumber]->GetBinContent (i)-hTrace[ChanNumber]->GetBinContent (j));				// here we have to recalculate the norm every time because of the non-linearity of the bilateral filter (norm is not const for position i)

			}
			
			value = value/norm;
			//cout << value << "\t" << norm << endl;
			hSmoothedTrace[ChanNumber]->SetBinContent(i,value);
		}
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
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for (Int_t i =1; i <Int_t(SignalTime[ChanNumber].size()); i++)
		{
			PUTime= SignalTime[ChanNumber][i]- SignalTime[ChanNumber][i-1];
			hEnergyTimeSinceLastPulse[ChanNumber]->Fill(PUTime,energy[ChanNumber][i]);
			hEnergyTimeSinceLastPulseCorr[ChanNumber]->Fill(PUTime,EnergyPileUpTimeCorrection(PUTime, energy[ChanNumber][i]));
			if (ChanNumber ==1)		// Cut version of analysis is at the moment only done for the first channel
			{
				for(Int_t j = 0; j < NumberOfPileUpTimeHistograms; j++)
				{
					if (energy[ChanNumber][i-1] > j*100 && energy[ChanNumber][i-1] < (j+1)*100-1 )
						hEnergyTimeSinceLastPulse_WithCuts[j]->Fill(PUTime,energy[ChanNumber][i]);
				}
			}
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
void THypGeMWD::Connect1DEnergySpectraHistograms(TH1D **hEnergySpectrum_ext,TH1D **hEnergySpectrumWithCut_ext)
{
	hEnergySpectrum = hEnergySpectrum_ext;
	hEnergySpectrumWithCut = hEnergySpectrumWithCut_ext;
}
void THypGeMWD::Connect1DRisetimeHistograms(TH1D** hRisetime1090_ext, TH1D** hRisetime3090_ext)
{
	hRisetime1090 = hRisetime1090_ext;
	hRisetime3090 = hRisetime3090_ext;
}
void THypGeMWD::Connect2DEnergyRisetimeHistograms(TH2D** hEnergyRise1090Corr_ext, TH2D** hEnergyRise3090Corr_ext)
{
	hEnergyRise1090Corr = hEnergyRise1090Corr_ext;
	hEnergyRise3090Corr = hEnergyRise3090Corr_ext;
}
void THypGeMWD::Connect2DEnergyTimeSinceLastPulseHistograms(TH2D** hEnergyTimeSinceLastPulse_ext, TH2D** hEnergyTimeSinceLastPulseCorr_ext ,TH2D** hEnergyTimeSinceLastPulse_WithCuts_ext, TH2D** hEnergyTimeSinceLastPulseCorr_WithCuts_ext, Int_t NumberOfCuts)
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
void THypGeMWD::ConnectTestHistograms(TH1D* hDeri1_ext, TH1D* hDeri2_ext, TH1D* hDeri3_ext, TH1D* hDeri4_ext)
{
	hTraceDeri1 = hDeri1_ext;
	hTraceDeri2 = hDeri2_ext;
	hTraceDeri3 = hDeri3_ext;
	hTraceDeri4 = hDeri4_ext;
}
Int_t THypGeMWD::AnaStep_FillEnergySpectrumWithPileUpCut(Double_t CutValue)
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		for (Int_t i =1; i <Int_t(SignalTime[ChanNumber].size()); i++)
		{
			if (SignalTime[ChanNumber][i]- SignalTime[ChanNumber][i-1] > CutValue )
				hEnergySpectrumWithCut[ChanNumber]->Fill(energy[ChanNumber][i]);
			//cout << i << "\t" <<energy[ChanNumber][i]<< "\t" << SignalTime[ChanNumber][i]- SignalTime[ChanNumber][i-1] << endl;
		}
	}
	return 0;
}

void THypGeMWD::CalculateDerivatives(TH1D* hInput,Int_t ChanNumber)
{
	Derivative1array[ChanNumber][0]=0;
	Derivative2array[ChanNumber][0]=0;
	Derivative3array[ChanNumber][0]=0;
	Derivative4array[ChanNumber][0]=0;
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Derivative1array[ChanNumber][i] = hInput->GetBinContent(i+1) - hInput->GetBinContent(i);
	}
	for(Int_t i=1;i<=TraceLength-1;i++)
	{
		Derivative2array[ChanNumber][i] = Derivative1array[ChanNumber][i+1] - Derivative1array[ChanNumber][i];
	}
	Derivative2array[ChanNumber][TraceLength]=0;
	for(Int_t i=1;i<=TraceLength-1;i++)
	{
		Derivative3array[ChanNumber][i] = Derivative2array[ChanNumber][i+1] - Derivative2array[ChanNumber][i];
	}
	Derivative3array[ChanNumber][TraceLength]=0;
	for(Int_t i=1;i<=TraceLength-1;i++)
	{
		Derivative4array[ChanNumber][i] = Derivative3array[ChanNumber][i+1] - Derivative3array[ChanNumber][i];
	}
	Derivative4array[ChanNumber][TraceLength]=0;
}
