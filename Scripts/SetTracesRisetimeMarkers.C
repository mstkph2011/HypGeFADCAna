#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMarker.h"
#include "TString.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

void makeMarker(TMarker *Min, Double_t MinimumBin, TH1D *hTrace, Int_t kColor)
{
	Min = new TMarker(hTrace->GetBinCenter(MinimumBin), hTrace->GetBinContent(MinimumBin),3);
	//cout << hTrace->GetBinCenter(MinimumBin) << "\t" << hTrace->GetBinContent(MinimumBin) << endl;
	Min->SetMarkerSize(10);
	Min->SetMarkerColor(kColor);
	Min->Draw();
}

int SetTracesRisetimeMarkers(TString TraceRootFileName = "../MyAnalysisASF.root", TString TraceInfoFileName = "../TraceOutput.txt")
{
	char buf[256];
	Double_t PeakRangeMin,PeakRangeMax, ADCValue, Rt10Pos, Rt30Pos, Rt90Pos;
	Int_t NumberOfTraces,NumberOfEventsInTrace, MinimumBin, MaximumBin ;
	ifstream TraceInfoFile;
	TraceInfoFile.open(TraceInfoFileName.Data());

	if (!TraceInfoFile.is_open())
	{
		cout << "File not found!" << endl;
		return -1;
	}
	for (int i = 0; i <8; i++)
	{
		TraceInfoFile.getline(buf,256);// just header
		//cout << buf << endl;
	}
	TraceInfoFile.getline(buf,256);
	//cout << buf << endl;
	sscanf (buf,"%*s %lf\t%lf",&PeakRangeMin,&PeakRangeMax);
	//cout << PeakRangeMin << "\t" <<PeakRangeMax<< endl;
	TraceInfoFile.getline(buf,256);
	sscanf (buf,"%*s %d",&NumberOfTraces);
	const int NoT = NumberOfTraces;
	const int NoTen = NumberOfTraces*10;


	TFile *RootInput = new TFile(TraceRootFileName);
	TH1D *hTrace[NoT];
	for (int i = 0; i <NumberOfTraces; i++)
	{
		sprintf(buf,"/Histograms/OutputTraces/OutPutTrace_%d;1",i);
		RootInput->GetObject(buf,hTrace[i]);
	}

	TH1D *hCoCompare[NoT];

	TFile *RootOutput = new TFile("Output.root","RECREATE");
	TMarker *Min[NoTen], *Max[NoTen],*Rt10[NoTen],*Rt30[NoTen],*Rt90[NoTen];
	int counter = 0;
	int bämcounter =0;
	TCanvas *cPeakCompare = cTrace = new TCanvas("cPeakCompare","cPeakCompare",800,600);
	TCanvas *cTrace = cTrace = new TCanvas("cTrace","cTrace",800,600);


	for (int i = 0; i <NumberOfTraces; i++)
	{
		TraceInfoFile.getline(buf,256);	// Trace begin indicator
		TraceInfoFile.getline(buf,256);
		NumberOfEventsInTrace = atoi(buf);


		cTrace->Clear();
		hTrace[i]->Draw();

		for(int j = 0; j < NumberOfEventsInTrace;j++)
		{
			TraceInfoFile.getline(buf,256);
			sscanf (buf,"%lf\t%d\t%lf\t%lf\t%lf\t%d",&ADCValue, &MinimumBin,&Rt10Pos, &Rt30Pos, &Rt90Pos, &MaximumBin );
			//cout << ADCValue << "\t" << MinimumBin <<"\t" << Rt10Pos <<"\t" <<Rt30Pos <<"\t" <<Rt90Pos<<"\t" <<MaximumBin << endl;
			makeMarker(Min[counter], MinimumBin, hTrace[i], 2);
			makeMarker(Rt10[counter], Rt10Pos, hTrace[i], 3);
			makeMarker(Rt30[counter], Rt30Pos, hTrace[i], 3);
			makeMarker(Rt90[counter], Rt90Pos, hTrace[i], 3);
			makeMarker(Max[counter], MaximumBin, hTrace[i], 2);
			counter++;
			if (ADCValue > PeakRangeMin&& ADCValue < PeakRangeMax)
			{
				cout << "bäm" << endl;
				cPeakCompare->cd();
				sprintf(buf, "PeakCompare_%d",bämcounter);
				hCoCompare[bämcounter] = new TH1D(buf,buf,100,0,100);
				for (int k = 1; k<100;k++)
				{
					hCoCompare[bämcounter]->SetBinContent(k, hTrace[i]->GetBinContent(MinimumBin-10+k)-hTrace[i]->GetBinContent(MinimumBin)+100);
				}
				hCoCompare[bämcounter]->Draw("same");
				cTrace->cd();
				bämcounter++;
			}
		}
		sprintf(buf,"OutputTrace_%d;1",i);
		cTrace->Write(buf);

	}
	cPeakCompare->Write("c1");
	RootInput->Close();
	RootOutput->Close();

	return 0;
}
