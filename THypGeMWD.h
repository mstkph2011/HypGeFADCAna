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
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"

class THypGeMWD
{
	public:
		THypGeMWD();			//default constructor DON'T use!
		THypGeMWD(Int_t TraceLength_ext,Int_t NumberOfChannels_ext);
		virtual ~THypGeMWD();
	
		//Double_t CalculateRisetime(TH1D* hTrace_ext,Double_t* times);
		
		Double_t FullAnalysis ();
		void				ConnectTraceHistograms(TH1D** hTrace_ext, TH1D** hSmoothedTrace_ext, TH1D** hTrace_bc_ext, TH1D** hAmplitude_ext,TH1D** hMWD_ext,TH1D** hMWDMA_ext,TH1D** hTrace_Direct_ext);
		void 				Connect1DEnergySpectraHistograms(TH1D **hEnergySpectrum_ext,TH1D **hEnergySpectrumWithCut_ext);
		void 				Connect1DRisetimeHistograms(TH1D** hRisetime1090_ext, TH1D** hRisetime3090_ext);
		void 				Connect2DEnergyRisetimeHistograms(TH2D** hEnergyRise1090Corr_ext, TH2D** hEnergyRise3090Corr_ext);
		void 				Connect2DEnergyTimeSinceLastPulseHistograms(TH2D** hEnergyTimeSinceLastPulse_ext,TH2D** hEnergyTimeSinceLastPulseCorr_ext, TH2D** hEnergyTimeSinceLastPulse_WithCuts_ext, TH2D** hEnergyTimeSinceLastPulseCorr_WithCuts_ext, Int_t NumberOfCuts);
		void				ConnectTestHistograms(TH1D* hDeri1_ext, TH1D* hDeri2_ext, TH1D* hDeri3_ext, TH1D* hDeri4_ext);
		private :
		
		//single steps of the analysis
		Int_t		AnaStep_Smoothing();
		Int_t		AnaStep_BaselineCorrection();
		Int_t		AnaStep_DoMovingWindowDeconvolution();
		Int_t		AnaStep_DoMovingAverageFilter();
		Int_t		AnaStep_DoDirectFilter();
		Int_t		AnaStep_FillEnergyspectrum();
		Int_t		AnaStep_ExtractRisetime();
		Int_t		AnaStep_DoEnergyRisetimeCorrelation();
		Int_t 	AnaStep_EnergyTimeSinceLastPulseCorrelation();
		Int_t		AnaStep_FillEnergySpectrumWithPileUpCut(Double_t CutValue);
		
		//functions used in the Energyspectrum step
		void 		EvaluateAmplitude();
		void		EvaluateMWD();
		void		EvaluateMA();
		
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
		
		
		
		void		SetUseMWD(Bool_t useMWD_ext);
		private :
		
		void		CalculateGausCoeff();
		Double_t 	Gaus(Double_t x);
		
		TF1				*EnergyRtCorrFunc;
		Double_t	EnergyRtCorrection(Double_t Rt, Double_t EnergyUncorr );
		
		TF1				*EnergyPileUpTimeCorrFunc;
		Double_t	EnergyPileUpTimeCorrection(Double_t PileUpTime, Double_t EnergyUncorr);						// energy correction for high rates (close signals)
		
		void 		DoMeanFilter();
		void 		DoWeightedAverageFilter();
		void 		DoGaussianFilter();
		void	 	DoBilateralFilter();
		
		private : void 		DoFourierTransformation();
		private : void 		DoBandStopFilter();
		private : void 		DoFourierBackTransformation();
		
		void CalculateDerivatives(TH1D* hInput,Int_t ChanNumber);
		
	private:
		Int_t 		TraceLength;
		Int_t			NumberOfChannels;
		Int_t 		PeakCounter;
		TH1D			*hTraceBuffer;
		
		//internal pointers to histograms
		
		TH1D			**hTrace; 
		TH1D			**hSmoothedTrace;
		TH1D			**hTrace_bc;
		TH1D			**hAmplitude;
		TH1D			**hMWD;
		TH1D			**hMWDMA;
		TH1D			**hTrace_Direct;
		
		TH1D			**hEnergySpectrum;
		TH1D			**hEnergySpectrumWithCut;
		
		TH1D			**hRisetime1090;
		TH1D			**hRisetime3090;
		
		TH2D			**hEnergyRise1090Corr;
		TH2D			**hEnergyRise3090Corr;
		
		TH2D			**hEnergyTimeSinceLastPulse;
		TH2D			**hEnergyTimeSinceLastPulse_WithCuts;
		TH2D			**hEnergyTimeSinceLastPulseCorr;
		TH2D			**hEnergyTimeSinceLastPulseCorr_WithCuts;

		TH1D			*hTraceDeri1;
		TH1D			*hTraceDeri2;
		TH1D			*hTraceDeri3;
		TH1D			*hTraceDeri4;

		Double_t offset_av;
		
		
		std::vector<Int_t> border_max;
		std::vector<Int_t> border_min;
		std::vector<Double_t> m;
		std::vector<Double_t> x;
		std::vector<Double_t> y;
		std::vector<Int_t> leftborder[8];
		std::vector<Int_t> rightborder[8];
		std::vector<Double_t> energy[8];
		std::vector<Double_t> risetime1090[8];
		std::vector<Double_t> risetime3090[8];
		
		std::vector<Int_t> SignalTime[8];
	
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
		Int_t 		PileUpTimeThreshold;
		Bool_t 		useMWD;				// Switch between MWD and Amplitude evaluation for energy spetrum
		
		Double_t g[10000];			// array for coefficients of gaussian filter
		
		Double_t **Aarray;
		Double_t **MWDarray;
		Double_t **Derivative1array;
		Double_t **Derivative2array;
		Double_t **Derivative3array;
		Double_t **Derivative4array;
		Double_t **MWDMAarray;

		Double_t **Sarray;		// array to store the result of the direct filter
		Double_t **Parray;		// array to store an intermediate result of the direct filter
		
		Int_t 		GausBreakUp;
		Double_t 	*GausNorm;
		
		Int_t 		NumberOfPileUpTimeHistograms;

		TStopwatch timer;
		
		
		Double_t mCorrection;
		//TH1D* Baseline;
		/* add your private declarations */
		
	ClassDef(THypGeMWD,1)			// doesn't compile with this line
};

#endif /* THYPGEMWD_H */ 
