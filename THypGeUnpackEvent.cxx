// $Id: THypGeUnpackEvent.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeUnpackEvent.h"
#include <iostream>

using namespace std;


//Int_t THypGeUnpackEvent::Init()				// doesn't work yet 28.03.13
//{
  //Int_t rev=0;
  ////cout << "+++ Init event" << endl;
  //Clear();
  //// is it used by Unpack step as output?
  //if(CheckEventSource("THypGeUnpackProc")){
    //fHypGeUP = (THypGeUnpackProc*)GetEventSource();
    //cout << "**** HypGe_UnpackEvent init for Unpack step"<< endl;
  //}
  //// or is it used from Analysis step as input
  //else if(CheckEventSource("TGo4FileSource")){
    //fHypGeFS = (TGo4FileSource*)GetEventSource();
    //cout << "**** HypGe_UnpackEvent init for Analysis step"<< endl;
  //}
  //else          rev=1;
  //return rev;
//}
////-----------------------------------------------------------
//Int_t THypGeUnpackEvent::Fill()				// doesn't work yet 28.03.13
//{
   //Int_t rev=0;
   //Clear();
   //if(fHypGeUP)fHypGeUP->BuildEvent(this);  // user event processing method
   //if(fHypGeFS)fHypGeFS->BuildEvent(this); // method from framework to restore event from file
   //return rev;
//}

void THypGeUnpackEvent::Clear(Option_t *t)
{
   for (int i=0;i<HypGe_NUM_CHAN;i++) {
      fiCrate1[i] = 0;
      fiCrate2[i] = 0;
      fiCrate3[i] = 0;
      fiCrate4[i] = 0;
   }
}
