// $Id: THypGeParameter.h 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef THypGeParameter_H
#define THypGeParameter_H

#include "TGo4Parameter.h"

#include "defines.h"

class THypGeParameter : public TGo4Parameter {
   public:
      THypGeParameter(const char* name = 0);
      virtual ~THypGeParameter() {}
#ifdef EXA_CODE
      Float_t frP1; // Offset for calibration
      Float_t frP2; // Factor for Calibration
      Bool_t fbHisto; // Enable Histogramming
#endif

	Int_t MWDm;			// M of MWD
	Int_t MAl;				// L of MA
	Int_t NoOfSmoothing;		// Number of smoothings of mean and WA filter
	Int_t Width;				// Width of mean filter
	Int_t sigmaGaus;			// sigma of gaussian shaper
	Int_t sigmaBil;			// second sigma of bilateral shaper
	Double_t tau;					//tau of MWD, time constant of preamp
	Int_t EnableMA;			// Switch for second moving average filter
	Int_t EnableSmoothing;	// Choose Smoothing Filter: 0 = Off, 1 = Mean, 2 = WA, 3 = Gaus, 4 = Bil
	Int_t EnableBaselineCorrection; 	//Switch baseline correction on or off
		
		
	void SetParameters( Int_t M, Int_t L, Int_t NOS_ext, Int_t Width_ext, Int_t Sigma, Int_t SigmaBil, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC);
	void PrintParameters();
	Int_t GetMWDm();
	Int_t GetMAl();
	Int_t GetNoOfSmoothing();		
	Int_t GetWidth();				
	Int_t GetSigmaGaus();
	Int_t GetSigmaBil();
	Double_t GetTau();
	Int_t GetEnableMA();
	Int_t GetEnableSmoothing();
	Int_t GetEnableBaselineCorrection();
	

	
   ClassDef(THypGeParameter,3)
};

#endif
