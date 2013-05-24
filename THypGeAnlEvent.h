// $Id: THypGeAnlEvent.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeANLEVENT_H
#define THypGeANLEVENT_H

#include "TGo4EventElement.h"
#include "THypGeUnpackEvent.h"

class THypGeAnlEvent : public TGo4EventElement {
   public:
      THypGeAnlEvent() : TGo4EventElement() {}
      THypGeAnlEvent(const char* name) : TGo4EventElement(name) {}
      virtual ~THypGeAnlEvent() {}

      virtual void  Clear(Option_t *t="");

      Float_t frData[HypGe_NUM_CHAN];

   ClassDef(THypGeAnlEvent,1)
};
#endif //THypGeANLEVENT_H



