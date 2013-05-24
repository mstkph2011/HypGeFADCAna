// $Id: THypGeParameter.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeParameter.h"

THypGeParameter::THypGeParameter(const char* name) :
   TGo4Parameter(name)
{
#ifdef EXA_CODE
   frP1=10;
   frP2=20;
   fbHisto=kTRUE;
#endif

	M = 300;
	L = 100;
}

