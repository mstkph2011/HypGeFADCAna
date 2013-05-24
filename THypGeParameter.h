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

		Int_t M;
		Int_t L;
		
   ClassDef(THypGeParameter,1)
};

#endif
