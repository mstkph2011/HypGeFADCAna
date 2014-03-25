// $Id: THypGeAnlProc.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeANLPROCESSOR_H
#define THypGeANLPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "THypGeUnpackEvent.h"

#include "defines.h"						// defines (~globals) are in this file, add this line to every file

#include "THypGeMWD.h"

#include "TMatrixD.h"
#include "TVectorD.h"

class THypGeParameter;
class THypGeAnlEvent;

class THypGeAnlProc : public TGo4EventProcessor {
   public:
      THypGeAnlProc();
      THypGeAnlProc(const char * name);
      virtual ~THypGeAnlProc();

      virtual Bool_t BuildEvent(TGo4EventElement* dest);

			THypGeParameter  	*fHypPar;
			Int_t 						EventCounter;
			//Add items (histograms, conditions, ...) used only  in this step here
			TH1D							*fhTrace[FADC_CHAN];
			TH1D							*fhTrace_Smoothed[FADC_CHAN];
			TH1D							*fhTrace_BaseCorr[FADC_CHAN];
			TH1D							*fhTrace_deconv[FADC_CHAN];
			TH1D							*fhTrace_MWD[FADC_CHAN];
			TH1D							*fhTrace_MA[FADC_CHAN];
			
			TH1D							*fhEnergySpectrum;
			TH1D							*fhRisetime1090;
			TH1D							*fhRisetime3090;
			TH2D							*fhEnergyRise1090Corr;
			TH2D							*fhEnergyRise3090Corr;
					
			TH2D							*fhEnergyTimeSinceLastPulse;
			
			TH1D							*fhAmplBaselinegradient;
			

			//real analysis object
			THypGeMWD					*fMWDAna;
			

   ClassDef(THypGeAnlProc, 1)
};
#endif //THypGeANLPROCESSOR_H
