// $Id: THypGeParameter.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeParameter.h"
#include <iostream>
using namespace std;
THypGeParameter::THypGeParameter(const char* name) :
   TGo4Parameter(name)
{
	ParametersChanged = 0;


}

void THypGeParameter::SetParameters( Int_t M, Int_t L, Int_t NOS_ext, Int_t Width_ext, Int_t Sigma, Int_t SigmaBil, Double_t tau_ext, Int_t EnaMA, Int_t EnaSmo, Int_t EnaBC)
{
	if (MWDm != M || MAl != L || NoOfSmoothing != NOS_ext || Width != Width_ext || sigmaGaus != Sigma || sigmaBil != SigmaBil || tau != tau_ext || EnableMA != EnaMA || SmoothingMethod != EnaSmo || EnableBaselineCorrection != EnaBC )
	{
		ParametersChanged = 1;
		MWDm = M;			// M of MWD
		MAl = L;				// L of MA
		NoOfSmoothing = NOS_ext;
		Width = Width_ext;
		sigmaGaus = Sigma;			// sigma of gaussian shaper
		sigmaBil = SigmaBil;			// second sigma of bil shaper
		tau = tau_ext;
		EnableMA = EnaMA;			// Switch for second moving average filter
		SmoothingMethod = EnaSmo;	// Switch smoothing on or off
		EnableBaselineCorrection = EnaBC; 	//Switch baseline correction on or off
		//std::cout << "MWDm is \t" << MWDm << "\t" << M << std::endl;
	}
	else
		ParametersChanged = 0;
}
void THypGeParameter::SetSecondAnaRoundParameters(Int_t SecondAnalysisRound_ext, TString ParameterFileName_ext)
{
	SecondAnalysisRound = SecondAnalysisRound_ext;
	ParameterFileName = ParameterFileName_ext;
}

void THypGeParameter::PrintParameters()
{
	std::cout << "MWDm is \t" << MWDm << endl;
	std::cout << "MAl is \t" << MAl << endl;
	std::cout << "NoS is \t" << NoOfSmoothing << endl;
	std::cout << "Width is \t" << Width << endl;
	std::cout << "sG is \t" << sigmaGaus << endl;
	std::cout << "sB is \t" << sigmaBil << endl;
	std::cout << "tau is \t" << tau << endl;
	std::cout << "EnaMA is \t" << EnableMA << endl;
	std::cout << "EnaSmo is \t" << SmoothingMethod << endl;
	std::cout << "EnaBC is \t" << EnableBaselineCorrection << endl;
}

Int_t THypGeParameter::GetMWDm()
{
	return MWDm;
}
Int_t THypGeParameter::GetMAl()
{
	return MAl;
}
Int_t THypGeParameter::GetNoOfSmoothing()
{
	return NoOfSmoothing;
}
Int_t THypGeParameter::GetWidth()
{
	return Width;
}

Int_t THypGeParameter::GetSigmaGaus()
{
	return sigmaGaus;
}

Int_t THypGeParameter::GetSigmaBil()
{
	return sigmaBil;
}
Double_t THypGeParameter::GetTau()
{
	return tau;
}
Int_t THypGeParameter::GetEnableMA()
{
	return EnableMA;
}
Int_t THypGeParameter::GetSmoothingMethod()
{
	return SmoothingMethod;
}
Int_t THypGeParameter::GetEnableBaselineCorrection()
{
	return EnableBaselineCorrection;
}

Int_t THypGeParameter::GetSecondAnalysisRound()
{
	return SecondAnalysisRound;
}
TString THypGeParameter::GetParameterFileName()
{
	return ParameterFileName;
}
Bool_t THypGeParameter::GetParametersChanged()
{
	return ParametersChanged;
}
void THypGeParameter::SetParametersChanged(Bool_t ParamaterChangedValue)
{
	ParametersChanged = ParamaterChangedValue;
}

