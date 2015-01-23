// $Id: THypGeParameter.h 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef THypGeParameter_H
#define THypGeParameter_H

#include "TGo4Parameter.h"
#include "TString.h"
#include "defines.h"

class THypGeParameter : public TGo4Parameter {
   public:
      THypGeParameter(const char* name = 0);
      virtual ~THypGeParameter() {}
	Int_t 							MWDm;			// M of MWD
	Int_t 							MAl;				// L of MA
	Int_t 							NoOfSmoothing;		// Number of smoothings of mean and WA filter
	Int_t 							Width;				// Width of mean filter
	Int_t 							sigmaGaus;			// sigma of gaussian shaper
	Int_t 							sigmaBil;			// second sigma of bilateral shaper
	Double_t 						tau;					//tau of MWD, time constant of preamp
	Int_t 							EnableMA;			// Switch for second moving average filter
	Int_t 							SmoothingMethod;	// Choose Smoothing Filter: 0 = Off, 1 = Mean, 2 = WA, 3 = Gaus, 4 = Bil
	Int_t 							EnableBaselineCorrection; 	//Switch baseline correction on or off
	Int_t								SecondAnalysisRound;	// parameter file from first analysis run with parameters for corrections from first analysis run exists and should be read
	TString							ParameterFileName;		// name and path of parameters file


	Bool_t 							ParametersChanged;

		
	void 								SetParameters( Int_t M, Int_t L, Int_t NOS_ext, Int_t Width_ext, Int_t Sigma, Int_t SigmaBil, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC);
	void								SetSecondAnaRoundParameters(Int_t SecondAnalysisRound_ext, TString ParameterFileName_ext);

	void 								PrintParameters();
	Int_t 							GetMWDm();
	Int_t 							GetMAl();
	Int_t 							GetNoOfSmoothing();
	Int_t 							GetWidth();
	Int_t 							GetSigmaGaus();
	Int_t 							GetSigmaBil();
	Double_t 						GetTau();
	Int_t 							GetEnableMA();
	Int_t 							GetSmoothingMethod();
	Int_t 							GetEnableBaselineCorrection();
	Int_t 							GetSecondAnalysisRound();
	TString							GetParameterFileName();

	Bool_t 							GetParametersChanged();
	void 								SetParametersChanged(Bool_t ParamaterChangedValue = 1);

	
   ClassDef(THypGeParameter,5)
};

#endif
