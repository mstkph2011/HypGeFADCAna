// $Id: THypGeNaIAnalysisProc.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeNaIAnalysisPROCESSOR_H
#define THypGeNaIAnalysisPROCESSOR_H

#include "defines.h"						// defines (~globals) are in this file, add this line to every file
#include "TGo4EventProcessor.h"

#include "TLatex.h"

//class THypGeParameter;
//class THypGeNaIAnalysisEvent;

#include "THypGeNaIAnalysisEvent.h"

class THypGeNaIAnalysisProc : public TGo4EventProcessor {
   public:
      THypGeNaIAnalysisProc();
      THypGeNaIAnalysisProc(const char* name);
      virtual ~THypGeNaIAnalysisProc();

      virtual Bool_t BuildEvent(TGo4EventElement* dest);

   protected:
     


      
      TH1D					*fhNaISpectrum;
      TH1D					*fhNaITrace;


   ClassDef(THypGeNaIAnalysisProc,1)
};

#endif //THypGeNaIAnalysisPROCESSOR_H
