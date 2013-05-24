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

class THypGeMWD
{
	public:
		THypGeMWD();			//default constructor DON'T use!
		THypGeMWD(Int_t TraceLength_ext);
		virtual ~THypGeMWD();
		
	
		void DoMWD();
		//Double_t CalculateRisetime(TH1D* hTrace_ext,Double_t* times);
		
		Double_t 	FullAnalysis (TH1D* hTrace_ext, TH1D* hSmoothedTrace, TH1D* hTrace_bc, TH1D* hAmplitude,TH1D* hMWD, TH1D* hEnergy, TH1D* hRisetime1090, TH1D* hRisetime3090, TH1D* hMWDMA, TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr);
		//TH1D*			MWD(TH1D* hTrace_ext);
		TH1D*			GetTrace();
		
		void 			SetTrace(TH1D* hTrace_ext);
		void			DoAna();
		
		//Double_t	EnergyRiseCorr(TH2D* fHisto);
		Int_t		Smoothing(TH1D* hTrace_ext, TH1D* hSmoothedTrace, Int_t NoOfSmothing);
		Int_t		Baseline(TH1D* hTrace_ext, TH1D* hTrace_bc);
		Int_t		MWD(TH1D* hTrace_ext, TH1D* hMWD, TH1D* hAmplitude);
		Int_t		MA(TH1D* hTrace_ext, TH1D* hMWDMA, TH1D* hMWD);
		Int_t		Energyspectrum(TH1D* hTrace_ext, TH1D* hEnergy,TH1D* hMWD);
		Int_t		Risetime(TH1D* hTrace_ext, TH1D* hRisetime1090, TH1D* hRisetime3090);
		Int_t		ERC(TH2D* hEnergyRise1090Corr, TH2D* hEnergyRise3090Corr);
		Int_t 		FindFirstBinAbove(TH1D* fHisto,Double_t threshold,Int_t low, Int_t high);
		Double_t 	FindLocalMaximum(TH1D* fHisto,Int_t low, Int_t high);
		Int_t 		FindLocalMaximumBin(TH1D* fHisto,Int_t low, Int_t high);
		
	
	private:
		Int_t 		TraceLength;
		Int_t 		PeakCounter;
		TH1D		*hTrace_internal;
		TH1D		*hTrace_bc;
		Double_t 	M;
		Double_t 	L;
		std::vector<Double_t> leftborder;
		std::vector<Double_t> energy;
		std::vector<Double_t> risetime1090;
		std::vector<Double_t> risetime3090;
		//TH1D* Baseline;
		/* add your private declarations */
		
	ClassDef(THypGeMWD,1)			// doesn't compile with this line
};

#endif /* THYPGEMWD_H */ 
