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
#include "TTree.h"
//#include "TBranch.h"

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

	if (UseTimer)
	{ 
		timer.Reset();
		timer.Start();
	}

	if ( SmoothingMethod)
	{
		AnaStep_Smoothing();
		PossiblePrintTimer("Smoothing step used ");
	}
	else
	{
		for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
		{
			for(Int_t i=1;i<=TraceLength;i++)				//loop over every point of the trace
			{
				hSmoothedTrace[ChanNumber]->SetBinContent(i,hTrace[ChanNumber]->GetBinContent (i));
			}
		}
	}
	AnaStep_BaselineCorrection();
	PossiblePrintTimer("Baseline correction step used");
	
	if (AnaStep_DoMovingWindowDeconvolution() == -1)
		return -1;
	PossiblePrintTimer("Deconvolution and MWD step used");
	
	if (EnableMA)
	{
		AnaStep_DoMovingAverageFilter();
		PossiblePrintTimer("MA step used");
	}
	
	AnaStep_FillEnergyspectrum();				//needs MWD and (if EnableMA) MA too
	PossiblePrintTimer("Energyspectrum step used");
	
	AnaStep_ExtractRisetime();		//needs Energyspectrum
	PossiblePrintTimer("Risetime step used ");

	AnaStep_DoDirectFilter();
	PossiblePrintTimer("Direct filter step used ");

	//cout << "finished direct filter"<< endl;
	AnaStep_FillHistograms();
	PossiblePrintTimer("Filling of Histograms step used ");
	if (UseTimer)
			timer.Stop();
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
			hTrace_Direct[ChanNumber]->SetBinContent(n,Sarray[ChanNumber][n]*2/1.242/1000000);				// factor is used to have the same output as the MA filter, factor is gained by trial and error, steinen, 21.1.15
		}
	}
	return 0;
}	
	//PileupCompensation and Energyspectrum
		
Int_t THypGeMWD::AnaStep_FillEnergyspectrum()
{
	if (EnableMA)
		EvaluateMA();
	else
	{
			if (useMWD)
			{
				EvaluateMWD();
			}
			else
			{
				cout << "No energy spectrum! O'rly?" << endl;
				return -1;
			}
	}
	return 0;
}
	
	
//Risetime
Int_t THypGeMWD::AnaStep_ExtractRisetime()
{
	//cout << "new Trace" << endl;
		
	Double_t RiseX10, RiseX30,RiseX90;
	Double_t RiseX1090, RiseX3090, RiseT1090, RiseT3090;
	Double_t threshold_Risetime = 50;		//threshold for signals to be identified as useful for Risetime calculation

	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{

		mTraceposRt10[ChanNumber].clear();
		mTraceposRt30[ChanNumber].clear();
		mTraceposRt90[ChanNumber].clear();
		for (std::map<Int_t,Double_t>::iterator it=mTraceposEnergy[ChanNumber].begin(); it!=mTraceposEnergy[ChanNumber].end(); ++it)
		{
			//it->first is bin number of the beginning of the rising edge

			Int_t MinimumPosition = FindLocalMinimumBin(hTrace[ChanNumber], it->first-10,it->first+10);
			Int_t MaximumPosition = FindLocalMaximumBin(hTrace[ChanNumber],it->first,it->first+60);
			//cout << "lb "<< it->first<<endl;
			if(hTrace[ChanNumber]->GetBinContent(MaximumPosition) > threshold_Risetime+hTrace[ChanNumber]->GetBinContent(MinimumPosition))
			{
				Double_t Amp10Threshold = (    hTrace[ChanNumber]->GetBinContent(MaximumPosition)+ 9 * hTrace[ChanNumber]->GetBinContent(MinimumPosition))/10;		//Calculates 10% threshold of pulse
				Double_t Amp30Threshold = (3 * hTrace[ChanNumber]->GetBinContent(MaximumPosition)+ 7 * hTrace[ChanNumber]->GetBinContent(MinimumPosition))/10;		//Calculates 30% threshold of pulse
				Double_t Amp90Threshold = (9 * hTrace[ChanNumber]->GetBinContent(MaximumPosition)+     hTrace[ChanNumber]->GetBinContent(MinimumPosition))/10;		//Calculates 90% threshold of pulse
				//cout << "Max " << hTrace[ChanNumber]->GetBinContent(MaximumPosition) << ", Min " << hTrace[ChanNumber]->GetBinContent(MinimumPosition) << ", Amp10Threshold " << Amp10Threshold << endl;
				//cout << "10t "<< Amp10Threshold << endl;
				RiseX10 = Double_t(FindFirstBinAbove(hTrace[ChanNumber],Amp10Threshold,MinimumPosition,MaximumPosition));
				RiseX30 = Double_t(FindFirstBinAbove(hTrace[ChanNumber],Amp30Threshold,MinimumPosition,MaximumPosition));
				RiseX90 = Double_t(FindFirstBinAbove(hTrace[ChanNumber],Amp90Threshold,MinimumPosition,MaximumPosition));

				//cout << "R10" << RiseX10 << endl;
				//cout << "R30" << RiseX30 << endl;
				//cout << "R90" << RiseX90 << endl;

				RiseX1090 = RiseX90 - RiseX10;
				RiseX3090 = RiseX90 - RiseX30;
				
				RiseT1090 = RiseX1090 * 10;		//Conversion into nanoseconds (FADC 100 MS/s)
				RiseT3090 = RiseX3090 * 10;
				//cout << RiseT1090 << endl;

				mTraceposRt10[ChanNumber][it->first] = RiseX10*10;
				mTraceposRt30[ChanNumber][it->first] = RiseX30*10;
				mTraceposRt90[ChanNumber][it->first] = RiseX90*10;
			}
		}
	}
	return 0;
}

Int_t		THypGeMWD::AnaStep_FillHistograms()
{
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		Int_t PulseNumber = 0;

		for (std::map<Int_t,Double_t>::iterator it=mTraceposEnergy[ChanNumber].begin(); it!=mTraceposEnergy[ChanNumber].end(); ++it)
		{
			std::map<Int_t,Double_t>::iterator it2 = mTraceposRt10[ChanNumber].find(it->first);		// check if a risetime was calculated for a found pulse with calculated energy, only pulses where both is found should be considered
			if (it2->first)
			{
				PulseStartingTime = it->first;
				PulseEnergy = it->second;
				//1D histos
				//fill energy spectrum after MA
				hEnergySpectrumMA[ChanNumber]->Fill(PulseEnergy);
				//fill energy spectrum after MA -- rise time correction
				hEnergySpectrumMACorr[ChanNumber]->Fill(EnergyRtCorrection(PulseEnergy,mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime]));
				//fill rt1090 distribution
				hRisetime1090[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime]);
				//fill rt3090 distribution
				hRisetime3090[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt30[ChanNumber][PulseStartingTime]);
				//fill rt1090 distribution - Co60 1332 keV line only
				if (PulseEnergy > Co60PeakRangeMin && PulseEnergy < Co60PeakRangeMax )
				{
					hRisetime1090Co1332Only[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime]);
				}

				//2D histos
				//fill rt1090 - energy correlation
				hEnergyRt1090[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime],PulseEnergy);
				if (PulseEnergy > Co60PeakRangeMin && PulseEnergy < Co60PeakRangeMax )
					hEnergyRt1090Co1332Only[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime],PulseEnergy);
				//fill rt3090 - energy correlation
				hEnergyRt3090[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt30[ChanNumber][PulseStartingTime],PulseEnergy);
				//fill rt1090 - energy correlation -- risetime correction
				if (isSR)
				{
					hEnergyRt1090CorrectionRt[ChanNumber]->Fill(mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime],EnergyRtCorrection(PulseEnergy,mTraceposRt90[ChanNumber][PulseStartingTime] - mTraceposRt10[ChanNumber][PulseStartingTime]));
					//cout <<EnergyRtCorrection(energyMA[ChanNumber][n],risetime1090[ChanNumber][n])<<endl;
				}

				// all the pile up time stuff, not really useful atm but adapted to work if needed
				if (PulseNumber >0)
				{
					PileUpTime= mSignalTime[ChanNumber][PulseStartingTime]- mSignalTime[ChanNumber][PulseStartingTimeBefore];
					hEnergyTimeSinceLastPulse[ChanNumber]->Fill(PileUpTime,PulseEnergy);
					hEnergyTimeSinceLastPulseCorr[ChanNumber]->Fill(PileUpTime,EnergyPileUpTimeCorrection(PulseEnergy,PileUpTime));
					if (PileUpTime > PileUpTimeThreshold )
						hEnergySpectrumWithCut[ChanNumber]->Fill(PulseEnergy);
					if (ChanNumber ==1)		// Cut version of analysis is at the moment only done for the first channel
					{
						for(Int_t j = 0; j < NumberOfPileUpTimeHistograms; j++)
						{
							if (PulseEnergyBefore > j*100 && PulseEnergyBefore < (j+1)*100-1 )
								hEnergyTimeSinceLastPulse_WithCuts[j]->Fill(PileUpTime,PulseEnergy);
						}

					}
				}

				PulseStartingTimeBefore = PulseStartingTime;
				PulseEnergyBefore = PulseEnergy;
				PulseNumber++;
			} // end of if (it2->first)
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

void THypGeMWD::IsSecondRun(Bool_t isSR_ext)
{
	isSR = isSR_ext;
}
void THypGeMWD::SetSecondRunParametersFileName(TString SecondRunParametersFileName_ext)
{
	SecondRunParametersFileName = SecondRunParametersFileName_ext;
}

void THypGeMWD::Init()
{
	UseTimer = 0;
	if (isSR)
	{
		TFile *InParameterFile = new TFile(SecondRunParametersFileName);


		EnergyRtCorrFuncPol = (TF1*) InParameterFile->Get("fProfileFitFuncPol");
		EnergyRtCorrFuncConst= (TF1*) InParameterFile->Get("fProfileFitFuncConstant");
		TTree *Tree = (TTree*) InParameterFile->Get("tNtupleD");
		Double_t Co60PeakRangeMinBuf, Co60PeakRangeMaxBuf;
		Tree->SetBranchAddress("RangeMin",&Co60PeakRangeMinBuf);
		Tree->SetBranchAddress("RangeMax",&Co60PeakRangeMaxBuf);
		Tree->GetEntry(0);
		Co60PeakRangeMin = Co60PeakRangeMinBuf;		// values copied from buffers -> file can be closed
		Co60PeakRangeMax = Co60PeakRangeMaxBuf;
			Tree->Delete();
			Tree = 0;
		InParameterFile->Close();
	}
	else
	{
		Co60PeakRangeMin = -1;
		Co60PeakRangeMax = -1;
	}
		// PileUpTimeCorr to be reworked! (maybe not neccessary) -- steinen, 13.1.15
	EnergyPileUpTimeCorrFunc = new TF1("EnergyPileUpTimeCorrFunc","[0]*(1-[1]/pow(x,[2]))",0,2000);
		EnergyPileUpTimeCorrFunc->SetParameters(1,3.74495e-02,1.95113);
	cout << "Init of THypGeMWD successful!"<<endl;
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
		mSignalTime[ChanNumber].clear();
		mTraceposEnergy[ChanNumber].clear();
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

			/*	for(i=rightborder_low+1;i<=rightborder_low+10;i++)				// calculate average value of base after pulse
				{
					if(hMWD[ChanNumber]->GetBinContent(i) < threshold_MWD)
					{
						Sumnoise_right = Sumnoise_right + hMWD[ChanNumber]->GetBinContent(i);
						k++;
					}
				}
			*/
				Noise_av = (Sumnoise_left + Sumnoise_right)/k;					// average of complete noise (left + right)/NoAllBaselinePoints

				Energy_av = Energy_av - Noise_av; //+ M*mCorrection/50;			// added correction for inclined amplitude signal --> most likely not working



				if (Energy_av > threshold_MWD)
				{
					mTraceposEnergy[ChanNumber][leftborder_low] = Energy_av;
					mSignalTime[ChanNumber][leftborder_low]=(leftborder_low+leftborder_high)/2/100;	// used for Energy - time before pulse Correlation, time value in µs
				}
			}
		}
	}
}

void THypGeMWD::EvaluateMA()
{

	Double_t EnergyThreshold = 100;		// threshold to fill the energy histogram (get rid of low energy gammas)

	EvalMAThreshold = 1;
	Int_t posMax1 = 0;
	Int_t posMax2 = 0;
	Int_t SignalCenter = 0;
	Int_t SumMax = 0;
	Int_t SumOffsetLeft, SumOffsetRight, OffsetWidth;
	Int_t a = M-L;				// width of top of trapezoid, feet width is M+L
	Int_t b = M;
	Double_t Energy;
	for (Int_t ChanNumber = 0; ChanNumber < NumberOfChannels; ChanNumber++)
	{
		mTraceposEnergy[ChanNumber].clear();
		if (a<0) break;
		CalculateDerivatives(hMWDMA[ChanNumber],ChanNumber);
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hTraceDeri1 ->SetBinContent(i, Derivative1array[ChanNumber][i]*100);
			hTraceDeri2 ->SetBinContent(i, Derivative2array[ChanNumber][i]*1000);
			hTraceDeri3 ->SetBinContent(i, Derivative3array[ChanNumber][i]*10000);
			hTraceDeri4 ->SetBinContent(i, Derivative4array[ChanNumber][i]*100000);
		}
		for(Int_t i=1;i<=TraceLength;i++)
		{
			if (Derivative1array[ChanNumber][i] > EvalMAThreshold)
			{
				SumMax = 0;
				SumOffsetLeft = 0;
				SumOffsetRight = 0;
				OffsetWidth = 0;
				posMax1 = FindLocalMaximumBin(Derivative2array[ChanNumber], i,i+ a/2.);
				if (posMax1 == TraceLength)
					continue;
				i = i + b + a;
				posMax2 = FindLocalMaximumBin(Derivative2array[ChanNumber], i,i+ a/2.);
				i = posMax2 + a/2.;
				SignalCenter = (posMax1 + posMax2)/2;
				//cout << "PM " << posMax1 << "\t" << posMax2  << endl;
				//cout << "SC " << SignalCenter  << endl;
				for (Int_t nBin = SignalCenter-a/4 ; nBin <= SignalCenter + a/4 ; nBin++)
				{
					SumMax +=  hMWDMA[ChanNumber]->GetBinContent(nBin);
				}
				for (Int_t nBin = posMax1-a-20 ; nBin <= posMax1-a ; nBin++)
				{
					if (nBin >0)
					{
						SumOffsetLeft +=  hMWDMA[ChanNumber]->GetBinContent(nBin);
						OffsetWidth++;
					}
				}
				/*for (Int_t nBin = posMax2+a ; nBin <= posMax2+a+20 ; nBin++)
				{
					if (nBin <= TraceLength)
					{
						SumOffsetRight +=  hMWDMA[ChanNumber]->GetBinContent(nBin);
						OffsetWidth++;
					}
				}
				*/

				Energy = Double_t(SumMax)/(a/2+1) - Double_t(SumOffsetRight + SumOffsetLeft)/OffsetWidth;
				if (Energy > EnergyThreshold)
				{
					int leftborder_low = hMWDMA[ChanNumber]->GetBin(SignalCenter-1.2*(M+L)/2);			// start of MA signal
					mTraceposEnergy[ChanNumber][leftborder_low]=Energy;
					mSignalTime[ChanNumber][leftborder_low]= SignalCenter;
				}
			}
		}// end of loop over trace

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

Double_t THypGeMWD::EnergyRtCorrection( Double_t EnergyUncorr, Double_t Rt)
{
	//cout << 1/EnergyRtCorrFuncPol->Eval(Rt) * EnergyRtCorrFuncConst->Eval(Rt) << endl;
	Double_t EnergyCorr = EnergyUncorr / EnergyRtCorrFuncPol->Eval(Rt) * EnergyRtCorrFuncConst->Eval(Rt);
	return EnergyCorr;
}


Double_t	THypGeMWD::EnergyPileUpTimeCorrection(Double_t EnergyUncorr,Double_t PileUpTime)
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
void THypGeMWD::Connect1DEnergySpectraHistograms(TH1D **hEnergySpectrumMA_ext, TH1D **hEnergySpectrumMACorr_ext ,TH1D **hEnergySpectrumWithCut_ext)
{

	hEnergySpectrumMA = hEnergySpectrumMA_ext;
	hEnergySpectrumMACorr = hEnergySpectrumMACorr_ext;
	hEnergySpectrumWithCut = hEnergySpectrumWithCut_ext;
}
void THypGeMWD::Connect1DRisetimeHistograms(TH1D** hRisetime1090_ext, TH1D** hRisetime3090_ext,TH1D** hRisetime1090Co1332Only_ext)
{
	hRisetime1090 = hRisetime1090_ext;
	hRisetime3090 = hRisetime3090_ext;
	hRisetime1090Co1332Only =hRisetime1090Co1332Only_ext;
}
void THypGeMWD::Connect2DEnergyRisetimeHistograms(TH2D** hEnergyRise1090Corr_ext,TH2D** fhEnergyRt1090CorrectionRt_ext, TH2D **hEnergyRt1090Co1332Only_ext, TH2D** hEnergyRise3090Corr_ext, TH2D** hEnergyRise1090CorrBallistic_ext)
{
	hEnergyRt1090 = hEnergyRise1090Corr_ext;
	hEnergyRt1090CorrectionRt= fhEnergyRt1090CorrectionRt_ext;
	hEnergyRt3090 = hEnergyRise3090Corr_ext;
	hEnergyRt1090Co1332Only = hEnergyRt1090Co1332Only_ext;
	//hEnergyRt1090CorrectionRt = hEnergyRise1090CorrBallistic_ext;
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

void THypGeMWD::PossiblePrintTimer(const char *text )
{
	if (UseTimer)
	{
		timer.Stop();
		cout << text << timer.RealTime() << " seconds" << endl;
		timer.Reset();
		timer.Start();
	}
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
