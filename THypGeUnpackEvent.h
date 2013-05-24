// $Id: THypGeUnpackEvent.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeEVENT_H
#define THypGeEVENT_H

#define HypGe_NUM_CHAN 8

#include "TGo4EventElement.h"
#include "THypGeUnpackProc.h"
#include "TGo4FileSource.h"

class THypGeUnpackProc;

class THypGeUnpackEvent : public TGo4EventElement {
   public:
      THypGeUnpackEvent() : TGo4EventElement() {}
      THypGeUnpackEvent(const char* name) : TGo4EventElement(name) {}
      virtual ~THypGeUnpackEvent() {}

			//Int_t Init();				// doesn't work yet 28.03.13
			//Int_t Fill();				// doesn't work yet 28.03.13
			
      /**
       * Method called by the event owner (analysis step) to clear the
       * event element.
       */
      virtual void Clear(Option_t *t="");

      Int_t fiCrate1[HypGe_NUM_CHAN];
      Int_t fiCrate2[HypGe_NUM_CHAN];
      Int_t fiCrate3[HypGe_NUM_CHAN];
      Int_t fiCrate4[HypGe_NUM_CHAN];


	private:
		THypGeUnpackProc* 		fHypGeUP;  		// Don't put this to file
    TGo4FileSource* 			fHypGeFS; 		// Don't put this to file
    
  public:
   ClassDef(THypGeUnpackEvent,1)
};
#endif //THypGeEVENT_H



