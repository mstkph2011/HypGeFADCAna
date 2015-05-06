// $Id: THypGeNaIAnalysisEvent.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeNaIAnalysisEvent.h"
#include <iostream>

using namespace std;


void THypGeNaIAnalysisEvent::Clear(Option_t *t)
{
   for (int i=0;i<HypGe_NUM_CHAN;i++) {
      fiCrate1[i] = 0;
      fiCrate2[i] = 0;
      fiCrate3[i] = 0;
      fiCrate4[i] = 0;
   }
}
