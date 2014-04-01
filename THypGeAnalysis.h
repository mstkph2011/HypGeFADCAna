// $Id: THypGeAnalysis.h 524 2009-11-11 09:53:34Z adamczew $
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

#ifndef THypGeANALYSIS_H
#define THypGeANALYSIS_H

#include "TGo4Analysis.h"
#include "defines.h"						// defines (~globals) are in this file, add this line to every file

#include "THypGeMWD.h"

class TH1D;
class TH2D;
class TGo4MbsEvent;
class TGo4WinCond;
class THypGeUnpackEvent;
class THypGeAnlEvent;
class THypGeParameter;

class THypGeAnalysis : public TGo4Analysis  {
   public:
      THypGeAnalysis();
      THypGeAnalysis(int argc, char** argv);
      virtual ~THypGeAnalysis() ;
      virtual Int_t UserPreLoop();
      virtual Int_t UserEventFunc();
      virtual Int_t UserPostLoop();
   private:
      TGo4MbsEvent       *fMbsEvent;
      THypGeUnpackEvent    *fRawEvent;
      THypGeAnlEvent       *fCalEvent;
#ifdef EXA_CODE
      
      TGo4WinCond        	*fWinCon1;
#endif

			THypGeParameter      *fPar;

      Int_t               fEvents;
      Int_t               fLastEvent;
      
      //put elements here, that have to be accesed in more than one step or in UserPostLoop()
      TH1D							 	*fhTrace[FADC_CHAN];



			int 								MWDm ;			// M of MWD
			int 								MAl ;				// L of MA
			int 								NoS;				// Number of smoothings of mean and WA filter
			int 								Width;				// Width of mean filter
			int 								sigmaGaus ;			// sigma of gaussian shaper
			int 								sigmaBil ;			// second sigma of bil shaper
			Double_t 						tau;					//tau of MWD, time constant of preamp
			int 								EnableMA ;			// Switch for second moving average filter
			int 								SmoothingMethod ;	// Switch smoothing on or off
			int 								EnableBaselineCorrection ; 	//Switch baseline correction on or off
			
   ClassDef(THypGeAnalysis,1)
};

#endif //THypGeANALYSIS_H



