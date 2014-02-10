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
	hTrace_internal =new TH1D();
	M = 200;					// window width for MWD
	L = 100;					// top width for MA
	NoOfSmoothing = 100;
	Width = 3;
	Sigma =10;					// sigma for gaussian smoothing
	Sigma2 = 10;
	tau = 50000;
	EnableMA = 0;			// Switch for second moving average filter
	EnableSmoothing = 1;	// Switch smoothing on or off
	EnableBaselineCorrection = 1; 	//Switch baseline correction on or off
	
	useMWD= 1;
	CalculateGausCoeff();
	//CalculateSecondGausCoeff();
}

THypGeMWD::~THypGeMWD()            //destructor
{
	
}

Double_t THypGeMWD::FullAnalysis (TH1D* hTrace_ext, TH1D* hSmoothedTrace, TH1D* hTrace_bc, TH1D* hAmplitude,TH1D* hMWD, TH1D* hEnergy, TH1D* hRisetime1090, TH1D* hRisetime3090, TH1D* hMWDMA, TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr)
{
	//cout << EnableSmoothing << endl;
	if ( EnableSmoothing)
		Smoothing(hTrace_ext, hSmoothedTrace);
	Baseline(hTrace_ext, hTrace_bc);
	if (MWD(hTrace_ext, hMWD, hAmplitude) == -1)
		return -1;
	if (EnableMA)
		MA(hTrace_ext, hMWDMA, hMWD);
	
	if ( useMWD)
		Energyspectrum(hTrace_ext, hEnergy, hMWD);				//needs MWD
	else
		Energyspectrum(hTrace_ext, hEnergy, hAmplitude);		//needs MWD
	Risetime(hTrace_ext, hRisetime1090, hRisetime3090);		//needs Energyspectrum
	ERC(hEnergyRise1090Corr, hEnergyRise3090Corr);			//needs Energyspectrum and Risetime
		
	return 0;
}

	
	//Smoothing
	
Int_t THypGeMWD::Smoothing(TH1D* hTrace_ext, TH1D* hSmoothedTrace)
{
	
	SetTrace(hTrace_ext);
	
	/*


	Double_t lower_frequency = 1;
	Double_t upper_frequency = 10;
	Double_t attenuation = 1;
    //Double_t M_PI = 3.1415;

    
	for(Int_t k=1; k<=TraceLength; k++)
	{
		std::complex<double> sum(0.0,0.0);
		for(Int_t j=1; j<=TraceLength; j++)
		{
			Int_t factor1 = -2*j*k;
			std::complex<Double_t> exponent1(0.0, M_PI/TraceLength*(Double_t)factor1);
			sum += hTrace_internal->GetBinContent(j) * std::exp(exponent1);
		}
		hTrace_fourier->SetBinContent(k,sum);
	}

	for(Int_t k=lower_frequency; k<=upper_frequency; k++)
	{ 
		hTrace_fourier->SetBinContent(k,hTrace_fourier->GetBinContent(k)/attenuation);
	}

	for(Int_t k=1; k<=TraceLength; k++)
	{
		std::complex<double> sum2(0.0,0.0);
		for(int j=1; j<=TraceLength; j++)
		{
			Int_t factor2 = 2*j*k;
			std::complex<Double_t> exponent2(0.0, M_PI/TraceLength*(Double_t)factor2);
			sum2 += hTrace_fourier->GetBinContent(j) * std::exp(exponent2);
		}
		hTrace_internal->SetBinContent(k,sum2);
	}

	*/
	
	/*
	cout << "test" <<endl;
	
	w = new TVectorD(TraceLength);
	v = new TVectorD(TraceLength);
	for(Int_t i=0;i<=TraceLength-1;i++)
	{
		(*v)(i) = hTrace_internal->GetBinContent(i+1);
	}
	cout << "test2" <<endl;
	
	//TMatrixD *S = TFile::Open("local/data/Go4_Analysen/Kai/S_Matrices/SmoothingMatrices/2048/MatrixTr2048_Sm100.root");
	//gFile->GetObject("local/data/Go4_Analysen/Kai/S_Matrices/SmoothingMatrices/2048/MatrixTr2048_Sm100.root",S);
	//S = S_Matrices.Get;
	
	
	(*w) = (*S) * (*v);
	for(Int_t i=0;i<=TraceLength-1;i++)
	{
		hTrace_internal->SetBinContent(i+1,(*w)(i));
	}
	*/
/*
	// weighted average
	for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)	
		{
			if (i == 1)	
				hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/2);
			else
				if (i == TraceLength)
					hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/2 + Double_t(hTrace_internal->GetBinContent(i))/2);
				else
					hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/4 + Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/4);
		}
	}
	
*/



/*
	//Rectangular
	Double_t width = 3;
	
	for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Double_t Content = 0;
			if(i <= width)
			{
				for(Int_t j=1-i;j<=width;j++)
				{
					Content += hTrace_internal->GetBinContent(i+j);
				}
				Content = Content/(width+i);
			}
			else
			{
				if(i >= TraceLength-width)
				{
					for(Int_t j=-width;j<=TraceLength-i;j++)
					{
						Content += hTrace_internal->GetBinContent(i+j);
					}
					Content = Content/(width+TraceLength-i+1);
				}	
				else
				{
					for(Int_t j=-width;j<=width;j++)
					{
						Content += hTrace_internal->GetBinContent(i+j);
					}
					Content = Content/(2*width+1);
				}
			}
			hSmoothedTrace->SetBinContent(i,Content);
		}
	}
*/	

	// gauss
	//Double_t Sigma = 7;		// changed to external parameter
	
		
		/*for(Int_t i=1;i<=TraceLength;i++)
		{	
			Double_t u = 0;
			for(Int_t k=0;k<GausBreakUp;k++)
			{
				if (i < GausBreakUp)	
				{
					u += g[k] * (2*hTrace_internal->GetBinContent(i+k));
					if(k==0)
					{
						u = u - g[0] * hTrace_internal->GetBinContent(i);
					}
				}
				else
				{
					if (i > TraceLength-GausBreakUp-1)
					{
						u += g[k] * (2*hTrace_internal->GetBinContent(i-k));
						if(k==0)
						{
							u = u - g[0] * hTrace_internal->GetBinContent(i);
						}
					}
					else
					{					
						u += g[k] * (hTrace_internal->GetBinContent(i-k) + hTrace_internal->GetBinContent(i+k));
						if(k==0)
						{
							u = u - g[0] * hTrace_internal->GetBinContent(i);
						}
					}
				}
			}
			hSmoothedTrace->SetBinContent(i,u);
		}*/
	if( EnableSmoothing == 1)
		DoMeanFilter(hSmoothedTrace);
	if( EnableSmoothing == 2)
		DoWeightedAverageFilter(hSmoothedTrace);
	if( EnableSmoothing == 3)
		DoGaussianFilter(hSmoothedTrace);
	if( EnableSmoothing == 4)	
		DoBilateralFilter(hSmoothedTrace);
	//// bilateral
	//Double_t value;
	//Double_t norm;
	//for (Int_t i=1;i<=TraceLength;i++)
	//{
		//value = 0;
		//norm = 0;
		//for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)
		//{
			////cout << j << endl;
			//if (j <1 )
				//continue;
			//if (j > TraceLength)
				//continue;
			//value += g[abs(j-i)]* Gaus(hTrace_internal->GetBinContent (i)-hTrace_internal->GetBinContent (j)) *hTrace_internal->GetBinContent (j);
			//norm += g[abs(j-i)]* Gaus(hTrace_internal->GetBinContent (i)-hTrace_internal->GetBinContent (j));
			////cout << g[j] << endl;
		//}
		//value = value/norm;
		
		////cout << value << "\t" << norm << endl;
		//hSmoothedTrace->SetBinContent(i,value);
	//}
	
	
//	*hSmoothedTrace = *hTrace_internal;			//use for weighted average and rectangular only
	*hTrace_internal = *hSmoothedTrace;
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
	Int_t BaseLineCorrectionBins = 20;//100;		//building average of first 100 bins
	
		for(Int_t i=1;i<=BaseLineCorrectionBins;i++)
		{
			sumoffset = sumoffset + hTrace_internal->GetBinContent(i);
		}

		 offset_av=hTrace_internal->GetMinimum();//sumoffset/BaseLineCorrectionBins;			// changed in Jülich, trace must be positive to avoid a strong negative drift of the deconvoluted signal!, remaining offset is corrected after deconvolution
		
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
		
	//tau =5000;//5383; //5294; //400; //1158;		//tau is time-constant of Preamlifier in 1/100 µs			// changed to external parameter
	//Double_t M = 200;		//Window broadness				// changed to external parameter
		
	A[0] = hTrace_internal->GetBinContent(1);			// first value of baseline corrected trace
	
	for(Int_t i=1;i<=TraceLength;i++)
	{
		//A[0] = 0;
		
		A[i] = hTrace_internal->GetBinContent(i) - hTrace_internal->GetBinContent(i-1) * (1.-(1./tau)) + A[i-1];
		
		
	}
	//bc correction of amp signal
	A[0]=0;
	Double_t SumAmpOffset=0;
	Int_t LengthSumOffsetAverage=20;
	for(Int_t i = 1; i <= LengthSumOffsetAverage; i++)
		SumAmpOffset+= A[i];
	SumAmpOffset= SumAmpOffset/LengthSumOffsetAverage;
		
	// try to correct slope in baseline - not working (yet)
	Double_t mCorrA[301];			// mCorrA[0] = sum and later on the mean
	mCorrA[0] = 0;
	for (Int_t i=1; i <= 300; i++)
	{
		mCorrA[i] = (A[1]-A[i+1])/i;				// get gradient of first 300 slope triangles
		mCorrA[0] += mCorrA[i];
	}
	mCorrection = mCorrA[0]/300;					// mean of slope
	
	for(Int_t i=1;i<=TraceLength;i++)
	{
		A[i]= A[i]-SumAmpOffset; //+ mCorrection*i;						// bc correction of amplitude signal

		if (i > Int_t(M))
			MWD[i] = A[i] - A[i-Int_t(M)];
		else
			MWD[i] = A[i]-A[1];				// added to cover edge efects 
			
		
		//hAmplitude->SetBinContent(i,A[i]);
		
		hMWD->SetBinContent(i,MWD[i]);
	}
	
			Double_t GradMWD1[TraceLength];
		Double_t GradMWD2[TraceLength];
		for(Int_t i=1;i<=TraceLength;i++)
		{
			GradMWD1[i] = hMWD->GetBinContent(i+5) - hMWD->GetBinContent(i);
		}
		for(Int_t i=1;i<=TraceLength;i++)
		{
			GradMWD2[i] = GradMWD1[i+5] - GradMWD1[i];
		}
		for(Int_t i=1;i<=TraceLength;i++)
		{
			hAmplitude->SetBinContent(i,GradMWD2[i]);
		}

	//if ((hAmplitude->GetBinContent(1)-hAmplitude->GetBinContent(301))/300 > 0.01)
		//return -1;
	
	return 0;
}
	
	//Moving-Average
	
Int_t THypGeMWD::MA(TH1D* hTrace_ext, TH1D* hMWDMA, TH1D* hMWD)
{
	if (hTrace_internal == 0)
		SetTrace(hTrace_ext);
	
	Double_t  MWDMA[TraceLength+1];
	
	MWDMA[0] = hMWD->GetBinContent(1);
	for(Int_t n=1;n<=TraceLength;n++)
	{
		if (n > Int_t(L))
			MWDMA[n] = MWDMA[n-1] + 1./L * (hMWD->GetBinContent(n) - hMWD->GetBinContent(n-Int_t(L)));
		else
			MWDMA[n] = MWDMA[n-1] + 1./L * (hMWD->GetBinContent(n) - hMWD->GetBinContent(1));
		
		hMWDMA->SetBinContent(n,MWDMA[n]);
	}

	*hTrace_internal = *hMWDMA;
	return 0;
}


	//PileupCompensation and Energyspectrum
		
Int_t THypGeMWD::Energyspectrum(TH1D* hTrace_ext, TH1D* hEnergy, TH1D* hMWD)
{
	
	
	SetTrace(hMWD);
	if (useMWD)
	{
		//if ( !EnableMA)
			//SetTrace(hMWD);
		EvaluateMWD( hEnergy);
	}
	else
	{
		//SetTrace(hMWD);
		EvaluateAmplitude( hEnergy);
	}
	
	return 0;
}
	
	//Risetime
	
Int_t THypGeMWD::Risetime(TH1D* hTrace_ext, TH1D* hRisetime1090, TH1D* hRisetime3090)
{
	//cout << "new Trace" << endl;
	if (hTrace_ext == 0)
		SetTrace(hTrace_ext);
	/*
	risetime1090.clear();
	risetime3090.clear();
	
	if(leftborder.size() == rightborder.size())
	{	
		Double_t Grad1[TraceLength];
		Double_t Grad2[TraceLength];
		for(Int_t i=1;i<=TraceLength;i++)
		{
			Grad1[i] = hTrace_ext->GetBinContent(i+5) - hTrace_ext->GetBinContent(i);
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
			
			Double_t leftborderrise = FindLocalMinimumBin(hTrace_ext,max(leftborderrisegrad-5,1),min(leftborderrisegrad+5,TraceLength));
			Double_t rightborderrise = FindLocalMaximumBin(hTrace_ext,max(rightborderrisegrad-5,1),min(rightborderrisegrad+5, TraceLength));
			
			//Double_t rightborderrise = FindFirstBinAboveInterpolated(Grad2,0,FindLocalMinimumBin(Grad2,max(rightborder[k]-40,1),min(rightborder[k]+20, TraceLength)),min(FindLocalMinimumBin(Grad2,max(rightborder[k]-40,1),min(rightborder[k]+20, TraceLength))+20,TraceLength));
			//cout << "lfbr" << hTrace_ext->GetBinContent(leftborderrise) << " rfbr "  << hTrace_ext->GetBinContent(rightborderrise) <<endl;
			Double_t mean_low = 0;
			for(Int_t i=leftborderrise-9;i<=leftborderrise;i++)
			{
				mean_low += hTrace_ext->GetBinContent(i);
			}
			mean_low = mean_low/10.;
			//cout << "mv1 " << mean_low << endl;
			Double_t mean_up = 0;
			for(Int_t i=rightborderrise;i<=rightborderrise+9;i++)
			{
				mean_up += hTrace_ext->GetBinContent(i);
			}
			mean_up = mean_up/10.;
			//cout << "mv2 " << mean_up << endl;
			if(mean_up > mean_low)
			{
				Double_t l = (mean_up-mean_low)/10;		//Calculates 10% of height of the rising signal
				//cout << "llllll " << mean_low+9*l << endl;
				//Double_t RiseX1090 = FindFirstBinAboveInterpolated(hTrace_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace_ext,mean_low+1*l,leftborderrise,rightborderrise);
				//Double_t RiseX1090 = Double_t(FindFirstBinAbove(hRisetime1090,leftborderrise+mean_low+9*l,leftborderrise,rightborderrise)) - Double_t(FindFirstBinAbove(hRisetime1090,leftborderrise+mean_low+l,leftborderrise,rightborderrise));
				//Double_t RiseX3090 = FindFirstBinAbove(hTrace_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAbove(hTrace_ext,mean_low+3*l,leftborderrise,rightborderrise);
				//cout << "bla" <<RiseX3090 << endl;	
				Double_t RiseX1090 = FindFirstBinAboveInterpolated(hTrace_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace_ext,mean_low+l,leftborderrise,rightborderrise);
				Double_t RiseX3090 = FindFirstBinAboveInterpolated(hTrace_ext,mean_low+9*l,leftborderrise,rightborderrise) - FindFirstBinAboveInterpolated(hTrace_ext,mean_low+3*l,leftborderrise,rightborderrise);
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
		while(hTrace_ext->GetBinContent(i+1) - hTrace_ext->GetBinContent(i) < gradRisetime && hTrace_ext->GetBinContent(i+2) - hTrace_ext->GetBinContent(i+1) < gradRisetime && hTrace_ext->GetBinContent(i+3) - hTrace_ext->GetBinContent(i+2) < gradRisetime);
		Double_t leftborder = hTrace_ext->GetBinContent(i);
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
		while(hTrace_ext->GetBinContent(i+1) - hTrace_ext->GetBinContent(i) >= 0 && hTrace_ext->GetBinContent(i+2) - hTrace_ext->GetBinContent(i+1) >= 0 && hTrace_ext->GetBinContent(i+3) - hTrace_ext->GetBinContent(i+2) >= 0);
		Double_t rightborder = hTrace_ext->GetBinContent(i);
		
		if(hTrace_ext->GetBinContent(rightborder) - hTrace_ext->GetBinContent(leftborder) > tresholdRisetime && rightborder - leftborder < 100)
		{
			Double_t mean_low = 0;
			for(Int_t i=leftborder-29;i<=leftborder;i++)
			{
				mean_low += hTrace_ext->GetBinContent(i);
			}
			mean_low = mean_low/30;
		
			Double_t mean_up = 0;
			for(Int_t i=rightborder;i<=rightborder+29;i++)
			{
				mean_up += hTrace_ext->GetBinContent(i);
			}
			mean_up = mean_up/30;
		
			cout << leftborder << endl;
			cout << rightborder << endl;
		
			Double_t l = (mean_up-mean_low)/10;		//Calculates 10% of height of the rising signal
		
			Double_t RiseX1090 = Double_t(FindFirstBinAbove(hTrace_ext,leftborder+mean_low+9*l,leftborder,rightborder)) - Double_t(FindFirstBinAbove(hTrace_ext,leftborder+mean_low+l,leftborder,rightborder));
			Double_t RiseX3090 = Double_t(FindFirstBinAbove(hTrace_ext,leftborder+mean_low+l,leftborder,rightborder)) - Double_t(FindFirstBinAbove(hTrace_ext,leftborder+mean_low+3*l,leftborder,rightborder));
						
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
		while(hTrace_ext->GetBinContent(d+1) - hTrace_ext->GetBinContent(d) >= 0 && abort==0);
		if(abort==0)
		{
		Double_t rightborder = hTrace_ext->GetBin(d);
		

		}
	}
	*/
	if (hTrace_ext == 0)
	SetTrace(hTrace_ext);
		
	risetime1090.clear();
	risetime3090.clear();

	Double_t threshold_Risetime = 50;		//threshold up where signals were identified as useful for Risetime calculation

	for(Int_t k=0;k<Int_t(leftborder.size());k++)
	{		
		Double_t rightborder = FindLocalMaximumBin(hTrace_ext,leftborder[k],leftborder[k]+40);
						
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

	if (hTrace_internal )
		delete hTrace_internal;
	hTrace_internal = (TH1D*) hTrace_ext->Clone("hTrace_internal");
	
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

Int_t THypGeMWD::FindLocalMinimumBin(TH1D* fHisto2,Int_t low2, Int_t high2)
{
	//Find local minimum bin between low and high
   Double_t Min = 0;
   Int_t MinBin = low2;
   
   for (Int_t bin=high2;bin<=low2;bin-1) 
   {
      if (fHisto2->GetBinContent(bin) < Min)
      {
				Min = fHisto2->GetBinContent(bin);
				MinBin = bin;
			}
   }
   return MinBin;
}

void THypGeMWD::GetS(TMatrixD* S_ext)
{
	S = S_ext;
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
	EnableSmoothing = EnaSmo;	// Switch smoothing on or off
	EnableBaselineCorrection = EnaBC; 	//Switch baseline correction on or off

}

void THypGeMWD::EvaluateAmplitude(TH1D* hEnergy)
{
	
	//for MWD without MA
	Double_t threshold_MWD = 50;		//threshold up where signals were identified as useful for MWD
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high, rightborder_high;
	Double_t rightborder_low = -9;
	Int_t abort1 = 0;
	leftborder.clear();
	energy.clear();
	Double_t energysum = 0;

	
	do
	{
		Double_t Sumenergy = 0;
		Double_t Sumnoise_left = 0;
		Double_t Sumnoise_right = 0;
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
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) < grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) < grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) < grad_MWD);
		
		leftborder_low = hTrace_internal->GetBin(i);
						
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
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) >= grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) >= grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) >= grad_MWD);
		
		leftborder_high = hTrace_internal->GetBin(i);
				
				
		if(hTrace_internal->GetBinContent(leftborder_low) < threshold_MWD && hTrace_internal->GetBinContent(leftborder_high) > threshold_MWD)
		{
			
			for(i=leftborder_high+5;i<=leftborder_high+15;i++)
			{
				Sumenergy = Sumenergy + hTrace_internal->GetBinContent(i);
			}
			
			Energy_av = Sumenergy/11;
			
			Int_t k = 0;
			
			for(i=leftborder_low-10;i<leftborder_low;i++)
			{	
				if(hTrace_internal->GetBinContent(i) < threshold_MWD)
				{
					Sumnoise_left = Sumnoise_left + hTrace_internal->GetBinContent(i);
					k++;
				}
			}
									
			Noise_av = Sumnoise_left/k;
		
			Energy_av = Energy_av - Noise_av;
			if (Energy_av > threshold_MWD)
			{
				hEnergy->Fill(Energy_av);
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
	if(hTrace_internal->GetBinContent(i) > threshold_MWD)
	{
		Double_t leftborder_low = FindLocalMinimumBin(Grad2_MWD,Grad2_MWD->GetBin(i)-30,hTrace_internal->GetBin(i)+30);
		Double_t leftborder_high = FindLocalMaximumBin(Grad2_MWD,Grad2_MWD->GetBin(i)-30,hTrace_internal->GetBin(i)+30);

		Double_t leftborder_low_base = FindLocalMinimumBin(Grad2_MWD,leftborder_low-30,leftborder_low);
		Double_t leftborder_high_base = FindLocalMaximumBin(Grad2_MWD,leftborder_high,leftborder_high+30);

		Double_t rightborder_high = FindLocalMaximumBin(Grad2_MWD,leftborder_high,leftborder_high+M);
		Double_t rightborder_low = FindLocalMinimumBin(Grad2_MWD,rightborder_high,rightborder_high+30);
		
		Double_t rightborder_high_base = FindLocalMaximumBin(Grad2_MWD,rightborder_high-30,rightborder_high);
		Double_t rightborder_low_base = FindLocalMinimumBin(Grad2_MWD,rightborder_low,rightborder_low+30);
*/		
		
		
		

void THypGeMWD::EvaluateMWD(TH1D* hEnergy)
{
	//for MWD
 
	Double_t threshold_MWD = 150;		//threshold up where signals were identified as useful for MWD
	Double_t grad_MWD = 2;			//gradient threshold to identify the borders of MWD-Signal
	
	Double_t leftborder_low, leftborder_high, rightborder_high;
	/*
	Double_t Grad3[TraceLength];
	Double_t Grad4[TraceLength];
	for(Int_t i=1;i<=TraceLength;i++)
	{
		Grad3[i] = hTrace_internal->GetBinContent(i+5) - hTrace_internal->GetBinContent(i);
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
					Sumenergy += hTrace_internal->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_min[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-20;i<=border_max[k]-10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1]+10;i<=border_max[k+1]+20;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergy->Fill(Energy);
				}
			}
			if(border_max[k] < border_min[k] && border_max[k+1] < border_min[k+1] && border_max[k+2] > border_min[k+2] && border_max[k+3] > border_min[k+3])
			{
				
				Double_t Sumnoise = 0;
				Int_t p = 0;
				for(i=border_max[k]-20;i<=border_max[k]-10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+3]+10;i<=border_max[k+3]+20;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Double_t Sumenergy = 0;
				for(Int_t i=border_min[k]+10;i<=border_max[k+1]-10;i++)
				{
					Sumenergy += hTrace_internal->GetBinContent(i);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergy->Fill(Energy);
				}
				Sumenergy = 0;
				for(Int_t i=border_max[k+2]+10;i<=border_min[k+3]-10;i++)
				{
					Sumenergy += hTrace_internal->GetBinContent(i);
				}
				Energy = Sumenergy/(border_min[k+3]-border_max[k+2]);
				Energy = Energy - Noise;
				if(Energy>20 && Energy<3500)
				{
					hEnergy->Fill(Energy);
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
					Sumenergy += hTrace_internal->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_min[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-10;i<=border_max[k];i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1];i<=border_max[k+1]+10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergy->Fill(Energy);
			}
			if(border_max[k] > border_min[k] && border_max[k+1] < border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_max[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTrace_internal->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_max[k]);
				Int_t p = 0;
				for(i=border_min[k]-10;i<=border_min[k];i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_min[k+1];i<=border_min[k+1]+10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergy->Fill(Energy);
			}
			if(border_max[k] < border_min[k] && border_max[k+1] < border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_min[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTrace_internal->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_max[k]-10;i<=border_max[k];i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_min[k+1];i<=border_min[k+1]+10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergy->Fill(Energy);
			}
			if(border_max[k] > border_min[k] && border_max[k+1] > border_min[k+1])
			{
				Double_t Sumenergy = 0;
				//Double_t Sumnoise = 0;
				for(Int_t l=border_min[k];l<=border_max[k+1];l++)
				{
					Sumenergy += hTrace_internal->GetBinContent(l);
				}
				Double_t Energy = Sumenergy/(border_max[k+1]-border_min[k]);
				Int_t p = 0;
				for(i=border_min[k]-10;i<=border_min[k];i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				for(i=border_max[k+1];i<=border_max[k+1]+10;i++)
				{
					if(hTrace_internal->GetBinContent(i)<treshold)
					{
						Sumnoise += hTrace_internal->GetBinContent(i);
						p++;
					}
				}
				Double_t Noise = Sumnoise/p;
				Energy = Energy - Noise;
				hEnergy->Fill(Energy);
			}
		}
	}
	*/
	/*
	for(Int_t i=0;i<=TraceLength;i++)
	{
		hEnergy->SetBinContent(i,Grad4[i]);
	}
	*/
	
	Double_t rightborder_low = -9;
	Int_t abort1 = 0;
	leftborder.clear();
	rightborder.clear();
	energy.clear();
	Double_t energysum = 0;

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
						
		//leftborder.push_back(leftborder_low);		//used for Risetime (add +10 to leftborder_low a high smoothing rates -> results from the discrepance between smoothed and unsmoothed MWD-shape)
		//cout << leftborder_low << endl;
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
		while(hTrace_internal->GetBinContent(i+1) - hTrace_internal->GetBinContent(i) >= grad_MWD && hTrace_internal->GetBinContent(i+2) - hTrace_internal->GetBinContent(i+1) >= grad_MWD && hTrace_internal->GetBinContent(i+3) - hTrace_internal->GetBinContent(i+2) >= grad_MWD);
		
		leftborder_high = hTrace_internal->GetBin(i);
				
		i = leftborder_high+5;
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
		
		i = rightborder_high+5;
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
			
			Energy_av = Energy_av - Noise_av; //+ M*mCorrection/50;			// added correction for inclined amplitude signal --> most likely not working
			
			
			
			if (Energy_av > threshold_MWD)
			{
				hEnergy->Fill(Energy_av);
				energy.push_back(Energy_av);		//used for Energy-Risetime-Correlation 
				leftborder.push_back(leftborder_low);		//used for Risetime
				rightborder.push_back(leftborder_high);		//used for Risetime 
			}
		}
	}
	while(abort1==0);
}

void THypGeMWD::SetUseMWD(Bool_t useMWD_ext)
{
	useMWD = useMWD_ext;
}

void THypGeMWD::CalculateGausCoeff()
{
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
	//cout << "GC calcd" << endl;
}

Double_t THypGeMWD::EnergyRtCorrection(Double_t Rt, Double_t EnergyUncorr )
{
	Double_t EnergyCorr = EnergyUncorr / CorrFunc(Rt);
	
	return EnergyCorr;
}

Double_t THypGeMWD::Gaus(Double_t x)
{
	Double_t Norm =sqrt(2*TMath::Pi());
	return (TMath::Exp(-pow((x/(sqrt(2)*Sigma2)),2)))/Double_t(Sigma2)/Norm;
}

Double_t THypGeMWD::CorrFunc(Double_t x)
{
	Double_t par[5];
	par[0] =1;
	par[1] = -9.1237e-05;
	par[2] = 7.95504e-07;
	par[3] =  -3.0778e-09;
	par[4] = 4.39764e-12;
	
	
	Double_t CorrValue = par[0] + par[1] * x + par[2] * x *x + par[3] * x * x * x +  par[4] * x * x * x * x ;
	return CorrValue;
}
void THypGeMWD::DoMeanFilter(TH1D* hSmoothedTrace)
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
					Content += hTrace_internal->GetBinContent(i+j);
				}
				Content = Content/(Width+i);
			}
			else
			{
				if(i >= TraceLength-Width)
				{
					for(Int_t j=-Width;j<=TraceLength-i;j++)
					{
						Content += hTrace_internal->GetBinContent(i+j);
					}
					Content = Content/(Width+TraceLength-i+1);
				}	
				else
				{
					for(Int_t j=-Width;j<=Width;j++)
					{
						Content += hTrace_internal->GetBinContent(i+j);
					}
					Content = Content/(2*Width+1);
				}
			}
			hSmoothedTrace->SetBinContent(i,Content);
		}
	}
}
void THypGeMWD::DoWeightedAverageFilter(TH1D* hSmoothedTrace)
{
	// weighted average
	for(Int_t k=0;k<NoOfSmoothing;k++)		//Smoothing k times
	{
		for(Int_t i=1;i<=TraceLength;i++)	
		{
			if (i == 1)	
				hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/2);
			else
				if (i == TraceLength)
					hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/2 + Double_t(hTrace_internal->GetBinContent(i))/2);
				else
					hSmoothedTrace->SetBinContent(i,Double_t(hTrace_internal->GetBinContent(i-1))/4 + Double_t(hTrace_internal->GetBinContent(i))/2 + Double_t(hTrace_internal->GetBinContent(i+1))/4);
		}
	}
}
void THypGeMWD::DoGaussianFilter(TH1D* hSmoothedTrace)
{
	// gauss
	for(Int_t i=1;i<=TraceLength;i++)
	{	
		Double_t u = 0;
		for(Int_t k=0;k<GausBreakUp;k++)
		{
			if (i < GausBreakUp)	
			{
				u += g[k] * (2*hTrace_internal->GetBinContent(i+k));
				if(k==0)
				{
					u = u - g[0] * hTrace_internal->GetBinContent(i);
				}
			}
			else
			{
				if (i > TraceLength-GausBreakUp-1)
				{
					u += g[k] * (2*hTrace_internal->GetBinContent(i-k));
					if(k==0)
					{
						u = u - g[0] * hTrace_internal->GetBinContent(i);
					}
				}
				else
				{					
					u += g[k] * (hTrace_internal->GetBinContent(i-k) + hTrace_internal->GetBinContent(i+k));
					if(k==0)
					{
						u = u - g[0] * hTrace_internal->GetBinContent(i);
					}
				}
			}
		}
		hSmoothedTrace->SetBinContent(i,u);
	}
}
void THypGeMWD::DoBilateralFilter(TH1D* hSmoothedTrace)
{
// bilateral
	Double_t value;
	Double_t norm;
	for (Int_t i=1;i<=TraceLength;i++)
	{
		value = 0;
		norm = 0;
		for (Int_t j = i-GausBreakUp; j <= i+GausBreakUp; j++)
		{
			//cout << j << endl;
			if (j <1 )
				continue;
			if (j > TraceLength)
				continue;
			value += g[abs(j-i)]* Gaus(hTrace_internal->GetBinContent (i)-hTrace_internal->GetBinContent (j)) *hTrace_internal->GetBinContent (j);
			norm += g[abs(j-i)]* Gaus(hTrace_internal->GetBinContent (i)-hTrace_internal->GetBinContent (j));
			//cout << g[j] << endl;
		}
		value = value/norm;
		
		//cout << value << "\t" << norm << endl;
		hSmoothedTrace->SetBinContent(i,value);
	}
}
