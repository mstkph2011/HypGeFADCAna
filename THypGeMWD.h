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

#include <map>
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"
#include <fstream>


class THypGeMWD
{
	public:
		THypGeMWD();			//default constructor DON'T use!
		THypGeMWD(Int_t TraceLength_ext,Int_t NumberOfChannels_ext);
		virtual ~THypGeMWD();
	
		//Double_t CalculateRisetime(TH1D* hTrace_ext,Double_t* times);
		
		Double_t 					FullAnalysis ();
		void							ConnectTraceHistograms(TH1D** hTrace_ext,TH1D** hTraceDeri_ext, TH1D** hSmoothedTrace_ext, TH1D** hTrace_bc_ext, TH1D** hAmplitude_ext,TH1D** hMWD_ext,TH1D** hMWDMA_ext,TH1D** hTrace_Direct_ext);
		void 							Connect1DEnergySpectraHistograms(TH1D **hEnergySpectrumMA_ext, TH1D **hEnergySpectrumMACorr_ext ,TH1D **hEnergySpectrumWithCut_ext);
		void 							Connect1DRisetimeHistograms(TH1D** hRisetime1090_ext, TH1D** hRisetime3090_ext,TH1D** hRisetime1090Co1332Only_ext);
		void 							Connect2DEnergyRisetimeHistograms(TH2D** hEnergyRise1090Corr_ext,TH2D** fhEnergyRt1090CorrectionRt_ext, TH2D **hEnergyRt1090Co1332Only_ext, TH2D** hEnergyRise3090Corr_ext, TH2D** hEnergyRise1090CorrBallistic_ext);
		void							ConnectRratioHistograms(TH1D** fhRratioCo1332Only_ext, TH2D** fhEnergyRratio_ext, TH2D** fhEnergyRratioCorrectionR_ext, TH2D** fhEnergyRratioCo1332Only_ext);
		void 							Connect2DEnergyTimeSinceLastPulseHistograms(TH2D** hEnergyTimeSinceLastPulse_ext,TH2D** hEnergyTimeSinceLastPulseCorr_ext, TH2D** hEnergyTimeSinceLastPulse_WithCuts_ext, TH2D** hEnergyTimeSinceLastPulseCorr_WithCuts_ext, Int_t NumberOfCuts);
		void							ConnectTestHistograms(TH1D* hDeri1_ext, TH1D* hDeri2_ext, TH1D* hDeri3_ext, TH1D* hDeri4_ext);
		void							ConnectOutputTraces(TH1D **hTraceOutput_ext);
		void							ConnectPulseFilteredHistograms(TH1D **hEnergyFirst_ext, TH1D **hEnergyOther_ext, TH1D **hRt1090First_ext, TH1D **hRt1090Other_ext);
		void							ConnectRtCorrelationHistograms(TH2D **hRt1030Rt1090Co1332Only_ext, TH2D **hRt1030Rt80100Co1332Only_ext);
		void							ConnectTraceDeriMaximumHistograms(TH1D **hTraceDeriMaximum_ext);

	private :
		
		//single steps of the analysis
		Int_t							AnaStep_Smoothing();
		Int_t							AnaStep_Derivative();
		Int_t							AnaStep_CheckTrace();
		Int_t							AnaStep_BaselineCorrection();
		Int_t							AnaStep_DoMovingWindowDeconvolution();
		Int_t							AnaStep_DoMovingAverageFilter();
		Int_t							AnaStep_DoDirectFilter();
		Int_t							AnaStep_FillEnergyspectrum();
		Int_t							AnaStep_ExtractRisetime();
		Int_t							AnaStep_ExtractRisetimeDeconv();
		Int_t							AnaStep_FillHistograms();
		
		//functions used in the Energyspectrum step
		void							EvaluateMWD();
		void							EvaluateMA();
		
		private:
		Int_t 						FindFirstBinAbove(TH1D* fHisto,Double_t threshold,Int_t low, Int_t high, Int_t NumberOfConsequentPoints = 0);
		Double_t 					FindFirstBinAboveInterpolated(TH1D* fHisto, Double_t threshold,Int_t low, Int_t high, Int_t NumberOfConsequentPoints = 0);
		Double_t 					FindFirstBinAboveInterpolated(Double_t *Array, Double_t threshold,Int_t low, Int_t high, Int_t NumberOfConsequentPoints = 0);
		Double_t 					FindLocalMaximum(TH1D* fHisto,Int_t low, Int_t high);
		Double_t 					FindLocalMaximum(Double_t *Array,Int_t low, Int_t high);
		Double_t 					FindLocalMinimum(Double_t *Array,Int_t low, Int_t high);
		Double_t 					FindLocalMinimum(TH1D* fHisto2,Int_t low2, Int_t high2);
		Int_t 						FindLocalMaximumBin(TH1D* fHisto,Int_t low, Int_t high);
		Int_t 						FindLocalMaximumBin(Double_t *Array,Int_t low, Int_t high);
		Int_t 						FindLocalMinimumBin(Double_t *Array,Int_t low, Int_t high);
		Int_t 						FindLocalMinimumBin(TH1D* fHisto2,Int_t low2, Int_t high2);
		
		public:
		void 							SetParameters( Int_t M_ext, Int_t L_ext,Int_t NoS_ext, Int_t Width_ext, Int_t Sigma_ext, Int_t SigmaBil_ext, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC);
		void							IsSecondRun(Bool_t isSR_ext = 1);
		void							SetSecondRunParametersFileName(TString SecondRunParametersFileName_ext);
		
		void							Init();
		void							Finish();
		
		void							SetUseMWD(Bool_t useMWD_ext);
		private :
		
		void							CalculateGausCoeff();
		Double_t 					Gaus(Double_t x);
		
		TF1								*EnergyRtCorrFuncPol;
		TF1								*EnergyRtCorrFuncConst;
		Double_t					EnergyRtCorrection(Double_t EnergyUncorr,Double_t Rt );
		TF1								*EnergyRratioCorrFuncPol;
		TF1								*EnergyRratioCorrFuncConst;
		Double_t					EnergyRratioCorrection(Double_t EnergyUncorr,Double_t Rratio );
		TF1								*EnergyPileUpTimeCorrFunc;
		Double_t					EnergyPileUpTimeCorrection(Double_t EnergyUncorr,Double_t PileUpTime);						// energy correction for high rates (close signals)
		
		void 							DoMeanFilter();
		void 							DoWeightedAverageFilter();
		void 							DoGaussianFilter();
		void	 						DoBilateralFilter();
		
		void 							DoFourierTransformation();
		void 							DoBandStopFilter();
		void 							DoFourierBackTransformation();
		
		void 							PossiblePrintTimer(const char *text );
		void 							CalculateDerivatives(TH1D* hInput,Int_t ChanNumber);
		void							WriteTestTraces(Int_t IntervalBetweenOutputs,Int_t MaxTraces);
		
		Int_t							FindDerivativeMaximum(Int_t ChanNumber, Int_t StartPosition);

	private:
		Int_t 						TraceLength;
		Int_t							NumberOfChannels;
		Int_t 						PulseCounter;
		TH1D							*hTraceBuffer;
		
		//internal pointers to histograms
		
		TH1D							**hTrace;
		TH1D							**hTraceDeri;
		TH1D							**hSmoothedTrace;
		TH1D							**hTrace_bc;
		TH1D							**hAmplitude;
		TH1D							**hMWD;
		TH1D							**hMWDMA;
		TH1D							**hTrace_Direct;
		

		TH1D							**hEnergySpectrumMA;
		TH1D							**hEnergySpectrumMACorr;
		TH1D							**hEnergySpectrumWithCut;
		TH1D							**fhEnergySpectrumFirstPulse;
		TH1D							**fhEnergySpectrumOtherPulses;
		
		TH1D							**hRisetime1090;
		TH1D							**hRisetime1090Co1332Only;
		TH1D							**hRisetime3090;
		TH1D							**fhRisetime1090FirstPulse;
		TH1D							**fhRisetime1090OtherPulses;
		
		TH2D							**hEnergyRt1090;
		TH2D							**hEnergyRt1090CorrectionRt;
		TH2D							**hEnergyRt1090Co1332Only;
		TH2D							**hEnergyRt3090;
		
		TH2D							**hRt1030Rt1090Co1332Only;
		TH2D 							**hRt1030Rt80100Co1332Only;

		TH1D							**fhRratioCo1332Only;			// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
		TH2D							**fhEnergyRratio;								// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
		TH2D							**fhEnergyRratioCorrectionR;		// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
		TH2D							**fhEnergyRratioCo1332Only;	// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)


		TH2D							**hEnergyTimeSinceLastPulse;
		TH2D							**hEnergyTimeSinceLastPulse_WithCuts;
		TH2D							**hEnergyTimeSinceLastPulseCorr;
		TH2D							**hEnergyTimeSinceLastPulseCorr_WithCuts;

		TH1D							*hTraceDeri1;
		TH1D							*hTraceDeri2;
		TH1D							*hTraceDeri3;
		TH1D							*hTraceDeri4;

		TH1D							**hTraceOutput;

		TH1D 							**hTraceDeriMaximum;
		Int_t							OutputTraceNumber;
		ofstream					TxtOutputFile;

		Double_t					offset_av;
		
		
		Double_t 					EvalMAThreshold;

		//std::vector<Int_t> border_max;
		//std::vector<Int_t> border_min;
		//std::vector<Double_t> m;
		//std::vector<Double_t> x;
		//std::vector<Double_t> y;

		std::map<Int_t,Double_t> 	mTraceposEnergy[8];
		std::map<Int_t,Double_t> 	mTraceposRt10[8];
		std::map<Int_t,Double_t> 	mTraceposRt30[8];
		std::map<Int_t,Double_t> 	mTraceposRt80[8];
		std::map<Int_t,Double_t> 	mTraceposRt90[8];
		std::map<Int_t,Double_t> 	mTraceposRratio[8];
		std::map<Int_t,Double_t> 	mSignalTime[8];
		std::map<Int_t,Int_t> 	mTraceposMinPositionBin[8];
		std::map<Int_t,Int_t> 	mTraceposMaxPositionBin[8];
		std::map<Int_t,Int_t> 	mPulseNumber[8];



			//external parameters
		Double_t 	M;					// window width for MWD
		Double_t 	L;					// top width for MA
		Int_t			NoOfSmoothing;		// Number of smoothings of mean and WA filter
		Int_t			Width;					// Width of mean filter
		Double_t 	Sigma;					// sigma for gaussian smoothing
		Double_t 	Sigma2;					// sigma for second gaussian of bilateral filter
		Double_t 	tau;
		Int_t 		EnableMA;			// Switch for second moving average filter
		Int_t 		SmoothingMethod;	// Switch smoothing on or off
		Int_t 		EnableBaselineCorrection; 	//Switch baseline correction on or off
		Int_t 		PileUpTimeThreshold;
		Bool_t 		useMWD;				// Switch between MWD and Amplitude evaluation for energy spetrum
		Bool_t 		isSR;					// Switch between first and second analysis round (second uses analysis results from first to make corrections)
		TString		SecondRunParametersFileName;	// Name of Filename for parameters from gained in first analysis round


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
		
		Int_t PulseStartingTime, PulseStartingTimeBefore;
		Double_t PulseEnergy,PulseEnergyBefore;
		Double_t PileUpTime;

		Int_t 		NumberOfPileUpTimeHistograms;

		Double_t Co60PeakRangeMin;
		Double_t Co60PeakRangeMax;


		TStopwatch timer;
		Bool_t UseTimer;
		
		Int_t EventNumber;

		Double_t mCorrection;
		
	ClassDef(THypGeMWD,2)			// doesn't compile with this line
};

#endif /* THYPGEMWD_H */ 
