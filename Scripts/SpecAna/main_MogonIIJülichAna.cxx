#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TNtupleD.h"
#include "TRint.h"							// starts the interactive mode of ROOT (needed for canvas etc.)
#include "TGraphErrors.h"

#include "THypGeSpectrumAnalyser.h"
#include "THypGePeakFitFunction.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

using namespace std;




int main(int argc, char* argv[] )
{
	TString InputFile,OutPutFileRoot,OutPutFileTxt;
	//TString SpektrumModel= "co60";
	TString SpektrumModel= "jülich2";
	
	
	bool UseFreeSkewedFitting=1;
	int SecondRun=0;
	if (argc==1)
	{
		cout << "Please give a InputFile" << endl;
		return -1;
	}
	if (argc==2)
	{
		InputFile= argv[1];
		cout << InputFile.Data() << endl;
		cout << "Please give a OutputFile!"<<endl;
		return 0;
	}
	if (argc==3)
	{
		 InputFile= argv[1];
		cout << InputFile.Data() << endl;
		if(InputFile.Contains("SR0"))
			SecondRun=0;
		if(InputFile.Contains("SR1"))
			SecondRun=1;
		if(InputFile.Contains("SR2"))
			SecondRun=2;
		OutPutFileRoot= argv[2];
		cout << OutPutFileRoot.Data() << endl;
		OutPutFileTxt=OutPutFileRoot;
		OutPutFileTxt.ReplaceAll(".root",".txt");
		cout << OutPutFileTxt.Data() << endl;
	}
	if (argc==4)
	{
		SpektrumModel=argv[3];
	}
	
	TH1D* hEnergyMA;
	TH2D* hEnergyRt1090;
	
	TFile *Input= new TFile(InputFile.Data());
		if (!SecondRun)
		{
			Input->GetObject("/Histograms/Energyspectrum/Energy_01;1",hEnergyMA);									// get histrogram
			Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090_01;1",hEnergyRt1090);									// get histrogram
		}
		else
		{
			if(SecondRun==1)
			{
				Input->GetObject("/Histograms/Energyspectrum/EnergyCorr_01;1",hEnergyMA);									// get corrected histrogram
				Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01;1",hEnergyRt1090);									// get corrected histrogram
			}
			else
			{
				if(SecondRun==2)
				{
					Input->GetObject("/Histograms/Energyspectrum/EnergyCorrCorr_01;1",hEnergyMA);									// get corrected histrogram
					Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01;1",hEnergyRt1090);									// get corrected histrogram
				}
			}
		}
		hEnergyMA->SetDirectory(0);
		hEnergyRt1090->SetDirectory(0);
	Input->Close();
	//TRint *App = new TRint("ROOT",0,0,0,0,kTRUE); //&argc,argv);
	THypGeSpectrumAnalyser *Ana;

	Ana= new THypGeSpectrumAnalyser(hEnergyMA,SpektrumModel.Data(), 35 );

		Ana->SetSearchRange(500,4000);
		//Ana->SetOutputPath(Path);
		Ana->SetTxtFileOutputName(OutPutFileTxt);
		Ana->SetRootFileOutputName(OutPutFileRoot);
		if(UseFreeSkewedFitting==1)
		{
			Ana->SetNewFunctionFitting();
		}
		else{
			Ana->SetGaussianFitting();
		}
		Ana->AnalyseSpectrum();
		cout << "FWHM 1332 kEv:\t\t"<< Ana->GetFWHMCo() << endl;

	return 0;
}
