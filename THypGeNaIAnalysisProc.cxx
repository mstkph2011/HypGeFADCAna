// $Id: THypGeNaIAnalysisProc.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeNaIAnalysisProc.h"
#include "TGo4UserException.h"

#include "defines.h"						// defines (~globals) are in this file, add this line to every file

#include "Riostream.h"
#include <time.h>
#include "TROOT.h"

#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TMath.h"

#include "s_filhe_swap.h"
#include "s_bufhe_swap.h"

#include "TGo4MbsEvent.h"
#include "TGo4WinCond.h"
#include "TGo4PolyCond.h"
#include "TGo4CondArray.h"
#include "TGo4Picture.h"

#include "THypGeParameter.h"
#include "THypGeNaIAnalysisEvent.h"

//-----------------------------------------------------------
THypGeNaIAnalysisProc::THypGeNaIAnalysisProc() :
   TGo4EventProcessor()
{
}
//-----------------------------------------------------------
THypGeNaIAnalysisProc::THypGeNaIAnalysisProc(const char* name) :
   TGo4EventProcessor(name)
{
	
  fhNaISpectrum = new TH1D("fhNaISpectrum","fhNaISpectrum",5000,0,5000);
  AddHistogram(fhNaISpectrum,"NaI");

  fhNaITrace = (TH1D*) GetHistogram("Traces/Trace_02");
	if (fhNaITrace)
		cout << "fhNaITrace not found"<< endl;
}
//-----------------------------------------------------------
THypGeNaIAnalysisProc::~THypGeNaIAnalysisProc()
{

}
//-----------------------------------------------------------
Bool_t THypGeNaIAnalysisProc::BuildEvent(TGo4EventElement* dest)
{
	
	//peak is inverted!!!
	double Min = fhNaITrace->GetMinimum();
	int MinBin = fhNaITrace->GetMinimumBin();
	double Max = fhNaITrace->GetBinContent(MinBin-30);
	
	fhNaISpectrum->Fill(TMath::Abs(Max-Min));
	return 1;
}
