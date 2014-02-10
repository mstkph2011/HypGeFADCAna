// $Id: THypGeUnpackProc.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeUNPACKPROCESSOR_H
#define THypGeUNPACKPROCESSOR_H

#include "defines.h"						// defines (~globals) are in this file, add this line to every file
#include "TGo4EventProcessor.h"

#include "TLatex.h"

class THypGeParameter;
class THypGeUnpackEvent;

#include "THypGeUnpackEvent.h"

class THypGeUnpackProc : public TGo4EventProcessor {
   public:
      THypGeUnpackProc();
      THypGeUnpackProc(const char* name);
      virtual ~THypGeUnpackProc();

      virtual Bool_t BuildEvent(TGo4EventElement* dest);

   protected:
#ifdef EXA_CODE
      TH1           *fCr1Ch[HypGe_NUM_CHAN];
      TH1           *fCr2Ch[HypGe_NUM_CHAN];
      TH2           *fCr1Ch1x2;
      TH1           *fHis1;
      TH1           *fHis1gate;
      TH1           *fHis2;
      TH1           *fHis2gate;
      TGo4WinCond   *fconHis1;
      TGo4WinCond   *fconHis2;
      TGo4WinCond   *fWinCon1;
      TGo4WinCond   *fWinCon2;
      TGo4PolyCond  *fPolyCon1;
      TGo4CondArray *fConArr1;
      TGo4CondArray *fConArr2;
      THypGeParameter *fParam1;
      TGo4Picture   *fPicture1;
      TGo4Picture   *fcondSet;
      TLatex* fLaText;
#endif
     
      THypGeParameter *fParam1; 
      Int_t 				UnpackCounter;
      
      TH1D					*fhTrace[FADC_CHAN];


   ClassDef(THypGeUnpackProc,1)
};

#endif //THypGeUNPACKPROCESSOR_H
