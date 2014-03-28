//      THypGeMWD.h
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


#ifndef THYPGEMWD_H
#define THYPGEMWD_H
#define htrace
#define h_CalculateRisetime
#define h_CalculateRisetime90
#define h_CalculateRisetime30
#define h_CalculateRisetime9030
#define fxT90
#define fxT30
#define fxT90T30

#include "TH2D.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"

class THypGeMWD
{
	public:
		THypGeMWD();			//default constructor DON'T use!
		THypGeMWD(Int_t TraceLength_ext);
		virtual ~THypGeMWD();
	
		//Double_t CalculateRisetime(TH1D* hTrace_ext,Double_t* times);
		
		Double_t 	FullAnalysis (TH1D* hTrace_ext, TH1D* hSmoothedTrace, TH1D* hTrace_bc, TH1D* hAmplitude,TH1D* hMWD, TH1D* hEnergy, TH1D* hRisetime1090, TH1D* hRisetime3090, TH1D* hMWDMA, TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr, TH2D* hEnergyTimeSinceLastPulse);
		
		void 				ConnectHistograms(TH1D *hEnergySpectrumWithCut_ext);
		
		
		TH1D*			GetTrace();
		
		void 			SetTrace(TH1D* hTrace_ext);
		
		//void			ConnectHistograms() // NYI 25.03.14  restructuring of analysis is needed for this. internal histograms necessary 
		
		//Double_t	EnergyRiseCorr(TH2D* fHisto);
		Int_t		Smoothing(TH1D* hTrace_ext, TH1D* hSmoothedTrace);
		Int_t		Baseline(TH1D* hTrace_ext, TH1D* hTrace_bc);
		Int_t		MWD(TH1D* hTrace_ext, TH1D* hMWD, TH1D* hAmplitude);
		Int_t		MA(TH1D* hTrace_ext, TH1D* hMWDMA, TH1D* hMWD);
		Int_t		Energyspectrum(TH1D* hTrace_ext, TH1D* hEnergy,TH1D* hMWD);
		Int_t		Risetime(TH1D* hTrace_ext, TH1D* hRisetime1090, TH1D* hRisetime3090);
		Int_t		ERC(TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr);
		Int_t 	EnergyTimeSinceLastPulseCorrelation(TH2D* hEnergyTimeSinceLastPulse);
		Int_t		FillEnergySpectrumWithPileUpCut(Double_t CutValue);
		
		private:
		Int_t 		FindFirstBinAbove(TH1D* fHisto,Double_t threshold,Int_t low, Int_t high);
		Double_t 	FindFirstBinAboveInterpolated(TH1D* fHisto, Double_t threshold,Int_t low, Int_t high);
		Double_t 	FindFirstBinAboveInterpolated(Double_t *Array, Double_t threshold,Int_t low, Int_t high);
		Double_t 	FindLocalMaximum(TH1D* fHisto,Int_t low, Int_t high);
		Double_t 	FindLocalMaximum(Double_t *Array,Int_t low, Int_t high);
		Double_t 	FindLocalMinimum(Double_t *Array,Int_t low, Int_t high);
		Double_t 	FindLocalMinimum(TH1D* fHisto2,Int_t low2, Int_t high2);
		Int_t 		FindLocalMaximumBin(TH1D* fHisto,Int_t low, Int_t high);
		Int_t 		FindLocalMaximumBin(Double_t *Array,Int_t low, Int_t high);
		Int_t 		FindLocalMinimumBin(Double_t *Array,Int_t low, Int_t high);
		Int_t 		FindLocalMinimumBin(TH1D* fHisto2,Int_t low2, Int_t high2);
		
		public:
		void 		SetParameters( Int_t M_ext, Int_t L_ext,Int_t NoS_ext, Int_t Width_ext, Int_t Sigma_ext, Int_t SigmaBil_ext, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC);
		
		void 		EvaluateAmplitude(TH1D* hEnergy);
		void		EvaluateMWD(TH1D* hEnergy);
		
		void		SetUseMWD(Bool_t useMWD_ext);
		
		void		CalculateGausCoeff();
		Double_t 	Gaus(Double_t x);
		
		Double_t 	CorrFunc(Double_t x);
		Double_t	EnergyRtCorrection(Double_t Rt, Double_t EnergyUncorr );
		
		void 		DoMeanFilter(TH1D* hSmoothedTrace);
		void 		DoWeightedAverageFilter(TH1D* hSmoothedTrace);
		void 		DoGaussianFilter(TH1D* hSmoothedTrace);
		void	 	DoBilateralFilter(TH1D* hSmoothedTrace);
		
		private : void 		DoFourierTransformation();
		private : void 		DoBandStopFilter();
		private : void 		DoFourierBackTransformation();
		
		
		
	private:
		Int_t 		TraceLength;
		Int_t 		PeakCounter;
		TH1D			*hTrace_internal;

		TH1D			*hEnergySpectrumWithCut;

		Double_t offset_av;
		
		
		std::vector<Int_t> border_max;
		std::vector<Int_t> border_min;
		std::vector<Double_t> m;
		std::vector<Double_t> x;
		std::vector<Double_t> y;
		std::vector<Int_t> leftborder;
		std::vector<Int_t> rightborder;
		std::vector<Double_t> energy;
		std::vector<Double_t> risetime1090;
		std::vector<Double_t> risetime3090;
		
		std::vector<Int_t> SignalTime;
	
			//external parameters
		Double_t 	M;					// window width for MWD
		Double_t 	L;					// top width for MA
		Int_t		NoOfSmoothing;		// Number of smoothings of mean and WA filter
		Int_t		Width;					// Width of mean filter
		Double_t 	Sigma;					// sigma for gaussian smoothing
		Double_t 	Sigma2;					// sigma for second gaussian of bilateral filter
		Double_t 	tau;
		Int_t 		EnableMA;			// Switch for second moving average filter
		Int_t 		SmoothingMethod;	// Switch smoothing on or off
		Int_t 		EnableBaselineCorrection; 	//Switch baseline correction on or off
		Int_t 		PileUpTimeCut;
		
		Bool_t 		useMWD;				// Switch between MWD and Amplitude evaluation for energy spetrum
		
		Double_t g[10000];
		
		Double_t *Aarray;
		Double_t *MWDarray;
		Double_t *GradMWD1array;
		Double_t *GradMWD2array;

		Int_t 		GausBreakUp;
		Double_t 	*GausNorm;
		
		TStopwatch timer;
		
		
		Double_t mCorrection;
		//TH1D* Baseline;
		/* add your private declarations */
		
	ClassDef(THypGeMWD,1)			// doesn't compile with this line
};

#endif /* THYPGEMWD_H */ 
