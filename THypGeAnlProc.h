// $Id: THypGeAnlProc.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef THypGeANLPROCESSOR_H
#define THypGeANLPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "THypGeUnpackEvent.h"

#include "defines.h"						// defines (~globals) are in this file, add this line to every file

#include "THypGeMWD.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TVectorD.h"


#define FADC_CHAN 1
// value from defines.h overwritten !!!!! in cxx too (for safety)


class THypGeParameter;
class THypGeAnlEvent;
class TTree;
//typedef struct EVENTN {double energy, t1090,t3090,t10max,tmax90;};

class THypGeAnlProc : public TGo4EventProcessor {
   public:
      THypGeAnlProc();
      THypGeAnlProc(const char * name);
      //void ConnectTree(TTree *tree_ext);
      virtual ~THypGeAnlProc();
		void UserPostLoop();
      virtual Bool_t BuildEvent(TGo4EventElement* dest);

			THypGeParameter  	*fHypPar;
			Int_t 						EventCounter;
			//Add items (histograms, conditions, ...) used only  in this step here
			TH1D							*fhTrace[FADC_CHAN];
			TH1D							*fhTrace_Deri[FADC_CHAN];
			TH1D							*fhTrace_Smoothed[FADC_CHAN];
			TH1D							*fhTrace_BaseCorr[FADC_CHAN];
			TH1D							*fhTrace_deconv[FADC_CHAN];
			TH1D							*fhTrace_MWD[FADC_CHAN];
			TH1D							*fhTrace_MA[FADC_CHAN];
			TH1D							*fhTrace_MA2[FADC_CHAN];
			TH1D							*fhTrace_Direct[FADC_CHAN];
			

			TH1D							*fhEnergySpectrum[FADC_CHAN];
			TH1D							*fhEnergySpectrum_MA2[FADC_CHAN];
			TH1D							*fhEnergySpectrumCorr[FADC_CHAN];
			TH1D							*fhEnergySpectrumFirstPulse[FADC_CHAN];
			TH1D							*fhEnergySpectrumOtherPulses[FADC_CHAN];
			TH1D							*fhRisetime1090[FADC_CHAN];
			TH1D							*fhRisetime1090Co1332FirstPulse[FADC_CHAN];
			TH1D							*fhRisetime1090Co1332OtherPulses[FADC_CHAN];
			TH1D							*fhRisetime3090[FADC_CHAN];
			TH1D							*fhRisetime1090Co1332Only[FADC_CHAN];
			TH2D							*fhEnergyRt1090[FADC_CHAN];
			TH2D							*fhEnergyRt1090CorrectionRt[FADC_CHAN];
			TH2D							*fhEnergyRt1090Co1332Only[FADC_CHAN];

			TH2D							*fhRt1030Rt1090Co1332Only[FADC_CHAN];
			TH2D							*fhRt1030Rt80100Co1332Only[FADC_CHAN];

			TH1D							*fhRratioCo1332Only[FADC_CHAN];			// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
			TH2D							*fhEnergyRratio[FADC_CHAN];								// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
			TH2D							*fhEnergyRratioCorrectionR[FADC_CHAN];		// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)
			TH2D							*fhEnergyRratioCo1332Only[FADC_CHAN];	// Rratio = (t_90 - t_30)/(t_90+t_30-2*t_10)

			TH2D							*fhEnergyRt1090_ballistic[FADC_CHAN];
			TH2D							*fhEnergyRt3090[FADC_CHAN];
					
			TH2D							*fhEnergyTimeSinceLastPulse[FADC_CHAN];
			TH2D							*fhEnergyTimeSinceLastPulse_withCut[20];
					
			TH2D							*fhEnergyTimeSinceLastPulseCorr[FADC_CHAN];
			TH2D							*fhEnergyTimeSinceLastPulseCorr_withCut[20];
			
			TH1D							*fhTraceDeriMaximumPosition[FADC_CHAN];
			//TH1D							*fh[FADC_CHAN];
			
			TH1D							*fhEnergySpectrum_withCut[FADC_CHAN];
			TH1D							*fhPreAmpTauFit[FADC_CHAN];
			
			TH1D 							*fhTrace_deri1;
			TH1D 							*fhTrace_deri2;
			TH1D 							*fhTrace_deri3;
			TH1D 							*fhTrace_deri4;

			TH1D							*hTraceOutput[1000];

			TH1D							*fhDeriMaxT90[FADC_CHAN];
			TH1D							*fhDeriMaxT90Rel[FADC_CHAN];
			TH2D							*fhEnergy_DeriMaxT90[FADC_CHAN];
			TH2D							*fhEnergy_DeriMaxT90Rel[FADC_CHAN];
			TH2D							*fDeriMaxT90Rel_T1090[FADC_CHAN];
			TH2D							*fDeriMaxT90_T1090[FADC_CHAN];
			TH2D							*fDeriMaxT90_DeriMaxT90Rel[FADC_CHAN];

			TH2D							*fhEnergyCorr_DeriMaxT90Rel[FADC_CHAN];
			
			TH1D							*fhT10DeriMax[FADC_CHAN];
			TH1D							*fhT10DeriMaxRel[FADC_CHAN];
			TH2D							*fhEnergy_T10DeriMax[FADC_CHAN];
			TH2D							*fhEnergy_T10DeriMaxRel[FADC_CHAN];
			TH2D							*fhEnergyCorr_T10DeriMaxRel[FADC_CHAN];
			TH1D							*fhMin[FADC_CHAN];
			
			TH1D							*fhEnergySpectrumCorrCorr[FADC_CHAN];
			TH2D							*fhEnergyRt1090CorrCorr[FADC_CHAN];

			//TTree							*ftDataTree;
			//real analysis object
			THypGeMWD					*fMWDAna;
			
			double 						fEnergy;
			double 						ft1090;
			double 						ft3090;
			double 						ft10max;
			double 						ftmax90;
			public:
			//struct EVENTN 
			//{
				//double energy;
				//double t1090;
				//double t3090; 
				//double t10max;
				//double tmax90;
				//};
			//EVENTN eventn;

   ClassDef(THypGeAnlProc, 1)
};
#endif //THypGeANLPROCESSOR_H
