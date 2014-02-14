// $Id: THypGeUnpackProc.cxx 755 2011-05-20 08:04:11Z linev $
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

#include "THypGeUnpackProc.h"
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

#include "s_filhe_swap.h"
#include "s_bufhe_swap.h"

#include "TGo4MbsEvent.h"
#include "TGo4WinCond.h"
#include "TGo4PolyCond.h"
#include "TGo4CondArray.h"
#include "TGo4Picture.h"

#include "THypGeParameter.h"
#include "THypGeUnpackEvent.h"

//-----------------------------------------------------------
THypGeUnpackProc::THypGeUnpackProc() :
   TGo4EventProcessor()
{
}
//-----------------------------------------------------------
THypGeUnpackProc::THypGeUnpackProc(const char* name) :
   TGo4EventProcessor(name)
{
  cout << "**** THypGeUnpackProc: Create" << endl;
	UnpackCounter = 0;

	fParam1   = (THypGeParameter *)   GetParameter("HypGeParameter");
	
	char chis[100];
	for (Int_t i =0;i < FADC_CHAN; i++)
	{
		snprintf(chis,63,"V1724/Trace%02d",i+1);
		fhTrace[i] = (TH1D*) GetHistogram(chis);
	}
	if (!fhTrace[0])
		cout << "fhTrace[0] not found"<< endl;
}
//-----------------------------------------------------------
THypGeUnpackProc::~THypGeUnpackProc()
{

}
//-----------------------------------------------------------
Bool_t THypGeUnpackProc::BuildEvent(TGo4EventElement* dest)
{
   Bool_t isValid=kFALSE; // validity of output event

   TGo4MbsEvent* inp_evt = (TGo4MbsEvent* ) GetInputEvent(); // from this
   THypGeUnpackEvent* out_evt = (THypGeUnpackEvent*) dest;

   if (inp_evt==0) {
      cout << "HypGeUnpackProc: no input event !"<< endl;
      out_evt->SetValid(isValid); // to store or not to store
      // default calling Fill method will set validity of out_evt to return value!
      return isValid;
   }
   isValid=kTRUE;
/*
	////////////////////////////////////////////////////
	// Some examples how to skip event processing or stop analysis by exception
	// for convenience, we provide  GO4_ macros to invoke appropriate exception throws
	// NOTE: You need to #include "TGo4UserException.h" for this feature
    //	 static UInt_t count=0;
	//	 if((count++ % 100000)==0 && count>1) // user may put a real analysis condition here
	//		 {
	//			 // this macro will skip event and subsequent analysis steps and send specified message to gui log window:
	//			 // GO4_SKIP_EVENT_MESSAGE("Skipped Event %d",count-1)
	//
	//			 // this macro will skip event and subsequent analysis steps without message:
	//			 GO4_SKIP_EVENT
	//
	//			// this macro will stop analysis and send specified message to gui log window:
	//			 //GO4_STOP_ANALYSIS_MESSAGE("Stopped after Event %d",count-1)
	//
	//			// this macro will stop analysis without message
	//			//GO4_STOP_ANALYSIS
	//		 }
	////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////
   ////// use this if you want access to the mbs file header data:
   //      s_filhe* head=inp_evt->GetMbsSourceHeader();
   //      if(head)
   //         {
   //            cout <<"found filhe structure:" << endl;
   //            cout <<"\tdatalen: "<<head->filhe_dlen << endl;
   //            cout <<"\tfilename_l: "<<head->filhe_file_l << endl;
   //            cout <<"\tfilename: "<<head->filhe_file << endl;
   //            cout <<"\ttype: "<<head->filhe_type << endl;
   //            cout <<"\tsubtype: "<<head->filhe_subtype << endl;
   //            cout <<"\t#commentlines: "<<head->filhe_lines << endl;
   //         }
   //      else
   //         {
   //            cout <<"zero file header" << endl;
   //         }
   //////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////
   ////// use this if you want access to the mbs buffer header data:
   //      s_bufhe* head=inp_evt->GetMbsBufferHeader();
   //      if(head) {
   //            cout <<"\nfound bufhe structure:" << endl;
   //            cout <<"\tbuffernumber: "<<head->l_buf << endl;
   //            cout <<"\tdatalen: "<<head->l_dlen << endl;
   //            cout <<"\ttime lo: "<<head->l_time[0] << endl; // seconds since epoch 1970
   //            cout <<"\ttime hi: "<<head->l_time[1] << endl; // microseconds since time lo
   //            char sbuf[1000]; f_ut_utime(head->l_time[0], head->l_time[1], sbuf);
   //            cout <<"\ttimestring: " << sbuf <<endl;
   //            cout <<"\ttype: "<<head->i_type << endl;
   //            cout <<"\tsubtype: "<<head->i_subtype << endl;
   //         }
   //      else
   //         {
   //            cout <<"zero buffer header" << endl;
   //         }
   //////////////////////////////////////////////////////////////////
   */ 
	Int_t lwords;
	Int_t *pdata;

	// HS, NEW variables

	Int_t *pl_data;
	short *pl_data16;

	Int_t j =0,
			 header,     // header contains event size
			 ev_size,
			 ev_cnt,
			 ch_cnt,
			 ch_mask,
			 mask,
			 time_tag,
			 samples,    // no of samples per trace
			 chan = 0;       // channel number

	 inp_evt->ResetIterator();
	 TGo4MbsSubEvent* psubevt(0);
	 while ((psubevt = inp_evt->NextSubEvent()) != 0) // subevent loop
	{

			pdata=psubevt->GetDataField(); //get data from subevent
      lwords= psubevt->GetIntLen();	//get length of subevent (number of 32 bit long words)

      if (psubevt->GetProcid()==25 && psubevt->GetIntLen()>2)     // V1724 subevent; 4 words header: words event size,1 word board ID (24 bit of board info, 1 word event counter (8 bit Event counter),1 word trigger time tag
      {
	//	sleep(1);
	//        printf("inside V1724 SubEvent!\n");

			pl_data=pdata;

			// read four header longwords
	// read size of event
			header=*pl_data++;
//          ev_size=header&0x0FFFFFFF;
			ev_size=header&0x00FFFFFF;         // NEW 12-06-08 bit 24..27 used for decimation factor of sampling rate
	// extract channel mask from 2nd word (see manual 3.3.5.1) 
			ch_mask=*pl_data++;
			ch_mask=ch_mask&0xFF;

			if (!ch_mask)     // empty event
				return kTRUE;
	// read event count
			ev_cnt=*pl_data++;
	// read trigger time tag
			time_tag=*pl_data++;

			// count number of active channels

			mask=1;
			ch_cnt=0;
			for(Int_t i=0;i<8;i++)			//FADC_CHAN max 8
			{
				if (mask&ch_mask)	// bitwise comparison, check 2^i'th bit, if set increase channel count by 1
					ch_cnt++;

				mask=2*mask;		// shift to next bit 2^n*2=2^(n+1)
			}

//printf("Number of channels: %d\n",ch_cnt);

			samples=2*((ev_size-4)/ch_cnt);	// get number of samples per event; 2 samples per 32 bit word, header must be substracted and divided by the number of channels

//	  printf("samples/chan: %d\n\n",samples);

			// fill spectra

			pl_data16 = (INTS2 *) pl_data;	// typecast to INTS2 (16 bit short, defined in <go4dir>/include/typedefs.h)
									// first 14 bits (0:13) used for data
			
			mask=1;
			for(;;)  // find next channel from channel mask // could also be a while, breaks if j == ch_cnt
			{
				if(mask&ch_mask)
				{                
//      printf("chan: %d\n",chan);
					for (Int_t i=0;i<samples;i++)   // Channel 
					{
						fhTrace[chan]->SetBinContent(i+1,(*pl_data16++)&0x3FFF);// fill first 14 bits to histogram
						//fTrace[chan]->GetBinContent(i+1)
						//(*pl_data16-1)&0x3FFF 
					}

						// length of measured trace shorter than size of spectrum
					if (samples<TRACE_LENGTH)     // NEW 11-06-08, fill rest of spectra with last measured sample
					{
						for (Int_t i=0;i<(TRACE_LENGTH-samples);i++)
						{
							fhTrace[chan]->SetBinContent(i+1+samples,(*(pl_data16-1))&0x3FFF);  
						}
					}
					j++;
					if (j==ch_cnt)		// ch_cnt = <activate channels> -1
						break;
				}

				chan++;							// just for histograms; chan = j
				mask=2*mask; 				// for bitwise comparison
			} // end for(;;)
	  //          printf("\n");
			//UnpackCounter++;
			//cout << "UnpackCounter " << UnpackCounter << endl;
		}  //  end V1724 subevent
	}  // end while

#ifdef EXA_CODE
   TString latext;
   latext.Form("#scale[3.0]{#color[2]{Event number:%d}}",inp_evt->GetCount());
   //latext.Form("Event number:%d",inp_evt->GetCount());
   fLaText->SetText(0.5,0.5,latext.Data());
#endif
   out_evt->SetValid(isValid); // to store or not to store
   // default calling Fill method will set validity of out_evt to return value!


   return isValid; // this will overwrite out_evt->SetValid
   // except one would have a Fill method implemented in the output event class!
   // In this case the return value can be handled differently there.
}
