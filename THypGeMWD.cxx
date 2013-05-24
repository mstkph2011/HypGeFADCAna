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
#include <vector>

using namespace std;

//default constructor DON'T use!
THypGeMWD::THypGeMWD()
{
}

THypGeMWD::THypGeMWD(Int_t TraceLength_ext)             //constructor
{
	TraceLength=TraceLength_ext;
	cout << "Analysis created"<<endl;
}

THypGeMWD::~THypGeMWD()            //destructor
{
	
}

Double_t THypGeMWD::FullAnalysis (TH1D* hTrace_ext, TH1D* hSmoothedTrace, TH1D* hTrace_bc, TH1D* hAmplitude,TH1D* hMWD, TH1D* hEnergy, TH1D* hRisetime1090, TH1D* hRisetime3090, TH1D* hMWDMA, TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr)
{
	Smoothing(hTrace_ext, hSmoothedTrace,100);
	Baseline(hTrace_ext, hTrace_bc);
	MWD(hTrace_ext, hMWD, hAmplitude);
	//MA(hTrace_ext, hMWDMA, hMWD);
	Energyspectrum(hTrace_ext, hEnergy, hMWD);				//needs MWD
	Risetime(hTrace_ext, hRisetime1090, hRisetime3090);		//needs Energyspectrum
	ERC(hEnergyRise1090Corr, hEnergyRise3090Corr);			//needs Energyspectrum and Risetime
		
	return 0;
}

	
	//Smoothing
	
Int_t THypGeMWD::Smoothing(TH1D* hTrace_ext, TH1D* hSmoothedTrace, Int_t NoOfSmothing)
{
	SetTrace(hTrace_ext);
	
	for(Int_t k=0;k<NoOfSmothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)	
		{
			if (i == 1)	
				hTrace_internal->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/2);
			else
				if (i == TraceLength)
					hTrace_internal->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/2 + Double_t(hTrace_internal->GetBinContent(i))/2);
				else
					hTrace_internal->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/4 + Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/4);
		}
	}

	*hSmoothedTrace = *hTrace_internal;
	hSmoothedTrace->SetName("SmoothedTrace");
	return 0;
}
	
	//Baseline
	
Int_t THypGeMWD::Baseline(TH1D* hTrace_ext, TH1D* hTrace_bc)
{
	if (hTrace_internal == 0)
		SetTrace(hTrace_ext);
		
	if (hTrace_internal->GetMaximum() == 0)
		return -1;
	Double_t sumoffset = 0;
	Int_t BaseLineCorrectionBins = 100;		//building average of first 100 bins
	
		for(Int_t i=1;i<=BaseLineCorrectionBins;i++)
		{
			sumoffset = sumoffset + hTrace_internal->GetBinContent(i);
		}

		Double_t offset_av=sumoffset/BaseLineCorrectionBins;
		
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hTrace_internal->SetBinContent(i,hTrace_internal->GetBinContent(i) - offset_av);
		}
		
	*hTrace_bc = *hTrace_internal;
	hTrace_bc->SetName("Trace_bc_1");
	return 0;
}
	
	//Moving-Window-Deconvolution
	
Int_t THypGeMWD::MWD(TH1D* hTrace_ext, TH1D* hMWD, TH1D* hAmplitude)
{
	if (hTrace_internal == 0)
	SetTrace(hTrace_ext);
	
	Double_t A[TraceLength+1], MWD[TraceLength+1];
		
	Double_t tau = 5294;		//tau is time-constant of Preamlifier in 1/100 µs
	Double_t M = 200;		//Window broadness
		
	for(Int_t i=1;i<=TraceLength;i++)
	{
		A[0] = 0;
		A[i] = hTrace_internal->GetBinContent(i) - hTrace_internal->GetBinContent(i-1) * (1.-(1./tau)) + A[i-1];
		
		if (i > Int_t(M))
			MWD[i] = A[i] - A[i-Int_t(M)];
		else
			MWD[i] = A[i];
			
		hAmplitude->SetBinContent(i,A[i]);
		
		hMWD->SetBinContent(i,MWD[i]);
	}
	
	*hTrace_internal = *hMWD;
	return 0;
}
	
	//Moving-Average
	
Int_t THypGeMWD::MA(TH1D* hTrace_ext, TH1D* hMWDMA, TH1D* hMWD)
{
	if (hTrace_internal == 0)
	SetTrace(hTrace_ext);
	
	Double_t MWD[TraceLength+1], MWDMA[TraceLength+1];
	
	Double_t L;	
		
	L = 100;		//Plateau broadness
	
	for(Int_t i=1;i<=TraceLength;i++)
		{
			MWDMA[i] = MWD[i];
		}
		
	for(Int_t n=1;n<=TraceLength;n++)
	{
		if (n > Int_t(L))
			MWDMA[n] = MWDMA[n-1] + 1./L * (Double_t(hTrace_internal->GetBinContent(n)) - Double_t(hMWD->GetBinContent(n-Int_t(L))));
		else
			MWDMA[n] = MWDMA[n-1] + 1./L * Double_t(hMWD->GetBinContent(n));
		
		hMWDMA->SetBinContent(n,MWDMA[n]);
	}
	*hTrace_internal = *hMWDMA;
	return 0;
}


	//PileupCompensation and Energyspectrum
		
Int_t THypGeMWD::Energyspectrum(TH1D* hTrace_ext, TH1D* hEnergy, TH1D* hMWD)
{
	if (hTrace_internal == 0)
	SetTrace(hTrace_ext);

	/*
	//test

	Double_t n;
	n = hTrace_internal->GetMaximum();
	hEnergy->Fill(n);
	
	//test end
	*/
	
	Double_t threshold_MWD = 50;		//threshold up where signals were identified as useful for MWD
	Double_t grad_MWD = 1;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high, rightborder_high;
	Double_t rightborder_low = -9;
	Int_t abort1 = 0;
	leftborder.clear();
	energy.clear();

	do
	{
		Double_t Sumenergy = 0;
		Double_t Sumnoise_left = 0;
		Double_t Sumnoise_right = 0;
		Double_t Energy_av = 0;
		Double_t Noise_av = 0;
	
		Int_t i = rightborder_low+10;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) < grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) < grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) < grad_MWD);
		
		leftborder_low = hTrace_internal->GetBin(i);
						
		leftborder.push_back(leftborder_low+10);		//used for Risetime (+10 results from the discrepance between smoothed and unsmoothed MWD-shape)

		i = leftborder_low+20;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) >= grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) >= grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) >= grad_MWD);
		
		leftborder_high = hTrace_internal->GetBin(i);
				
		i = leftborder_high+(M/4);
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) > -grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) > -grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) > -grad_MWD);
	
		rightborder_high = hTrace_internal->GetBin(i);
		
		i = rightborder_high+20;
		do
		{
			if(i+10 >= TraceLength)
			{
				abort1 = 1;
				break;
			}
			i++;
		}	
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) <= -grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) <= -grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) <= -grad_MWD);
	
		rightborder_low = hTrace_internal->GetBin(i);
		
		if(hTrace_internal->GetBinContent(leftborder_low) < threshold_MWD && hTrace_internal->GetBinContent(rightborder_low) < threshold_MWD && hTrace_internal->GetBinContent(leftborder_high) > threshold_MWD && hTrace_internal->GetBinContent(rightborder_high) > threshold_MWD)
		{
			
			for(i=leftborder_high+10;i<=rightborder_high-10;i++)
			{
				Sumenergy = Sumenergy + hTrace_internal->GetBinContent(i);
			}
			
			Energy_av = Sumenergy/(rightborder_high-leftborder_high-20);
			
			Int_t k = 0;
			
			for(i=leftborder_low-10;i<leftborder_low;i++)
			{	
				if(hTrace_internal->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_left = Sumnoise_left + hTrace_internal->GetBinContent(i);
					k++;
				}
			}
			
			for(i=rightborder_low+1;i<=rightborder_low+10;i++)
			{	
				if(hTrace_internal->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_right = Sumnoise_right + hTrace_internal->GetBinContent(i);
					k++;
				}
			}
			
			Noise_av = (Sumnoise_left + Sumnoise_right)/k;
		
			Energy_av = Energy_av - Noise_av;
			if (Energy_av >10)
				hEnergy->Fill(Energy_av);
			
			energy.push_back(Energy_av);		//used for Energy-Risetime-Correlation 
		}
	}
	while(abort1==0);
	return 0;
}
	
	//Risetime
	
Int_t THypGeMWD::Risetime(TH1D* hTrace_ext, TH1D* hRisetime1090, TH1D* hRisetime3090)
{
	if (hTrace_ext == 0)
	SetTrace(hTrace_ext);
		
	risetime1090.clear();
	risetime3090.clear();

	Double_t threshold_Risetime = 50;		//threshold up where signals were identified as useful for Risetime calculation

	for(Int_t k=0;k<Int_t(leftborder.size());k++)
	{		
		Double_t rightborder = FindLocalMaximumBin(hTrace_ext,leftborder[k],leftborder[k]+100);
		
		if(hTrace_ext->GetBinContent(rightborder) > threshold_Risetime+hTrace_ext->GetBinContent(leftborder[k]))
		{
			Double_t l = (hTrace_ext->GetBinContent(rightborder)-hTrace_ext->GetBinContent(leftborder[k]))/10;		//Calculates 10% of height of the rising signal
						
			Double_t RiseX1090 = Double_t(FindFirstBinAbove(hTrace_ext,hTrace_ext->GetBinContent(rightborder)-l,leftborder[k],rightborder)) - Double_t(FindFirstBinAbove(hTrace_ext,hTrace_ext->GetBinContent(rightborder)-9*l,leftborder[k],rightborder));
			Double_t RiseX3090 = Double_t(FindFirstBinAbove(hTrace_ext,hTrace_ext->GetBinContent(rightborder)-l,leftborder[k],rightborder)) - Double_t(FindFirstBinAbove(hTrace_ext,hTrace_ext->GetBinContent(rightborder)-7*l,leftborder[k],rightborder));
			
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
	
Int_t THypGeMWD::ERC(TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr)
{
		
	if(Int_t(energy.size()) == Int_t(risetime1090.size()))
	{
		for(Int_t n=0;n < Int_t(energy.size());n++)
		{
			hEnergyRise1090Corr->Fill(energy[n],risetime1090[n]);
		}
	}
	
	if(Int_t(energy.size()) == Int_t(risetime3090.size()))
	{
		for(Int_t n=0;n < Int_t(energy.size());n++)
		{
			hEnergyRise3090Corr->Fill(energy[n],risetime3090[n]);
		}
	}
	return 0; 
}
	
//delete hTrace_internal;


	

void THypGeMWD::DoMWD()
{
	cout << "Doing Analysis" << endl;
}

TH1D* THypGeMWD::GetTrace()
{
	return hTrace_internal;
}

void THypGeMWD::SetTrace(TH1D* hTrace_ext)
{
	hTrace_internal = (TH1D*) hTrace_ext->Clone();
	hTrace_internal->SetName("hTrace_internal");
}
void THypGeMWD::DoAna()
{

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
