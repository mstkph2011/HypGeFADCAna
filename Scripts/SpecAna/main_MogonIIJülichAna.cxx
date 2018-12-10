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
#include "TSystem.h"

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
	TString PathInRootFileCorr1="";
	TString PathInRootFileCorr2="";
	
	
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
	if (argc>=3)
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
	if (argc>=4)
	{
		SpektrumModel=argv[3];
	}
	if (argc>=5)
	{
		PathInRootFileCorr1=argv[4];
		cout << PathInRootFileCorr1.Data()<<endl;
	}
	if (argc==6)
	{
		PathInRootFileCorr2=argv[5];
		cout << PathInRootFileCorr2.Data()<<endl;
	}
	TH1D* hEnergyMA;
	//TH2D* hEnergyRt1090;
	
	TString EnergyHistoString="hEnergy;1";
	EnergyHistoString.Form("hEnergy%s%s;1",PathInRootFileCorr1.Data(),PathInRootFileCorr2.Data());
	//if(!PathInRootFileCorr1.IsNull() && PathInRootFileCorr2.IsNull())
	//{
		////EnergyHistoString.Form("%s;1/hEnergy%s;1",PathInRootFileCorr1.Data(),PathInRootFileCorr1.Data());
		//EnergyHistoString.Form("hEnergy%s;1",PathInRootFileCorr1.Data());
	//}
	//if(!PathInRootFileCorr1.IsNull()&& !PathInRootFileCorr2.IsNull())
	//{
		////EnergyHistoString.Form("%s;1/%s;1/hEnergy%s%s;1",PathInRootFileCorr1.Data(),PathInRootFileCorr2.Data(),PathInRootFileCorr1.Data(),PathInRootFileCorr2.Data());
		//EnergyHistoString.Form("hEnergy%s%s;1",PathInRootFileCorr1.Data(),PathInRootFileCorr2.Data());
	//}
	cout << EnergyHistoString.Data() << endl;
	//return 0;
	//changes here for the proper input
	TFile *Input= new TFile(InputFile.Data());
		if(!PathInRootFileCorr1.IsNull())
		{
			TDirectory *Subdir = Input->GetDirectory(PathInRootFileCorr1.Data());
			gDirectory->cd(PathInRootFileCorr1.Data());
			if(!PathInRootFileCorr2.IsNull())
			{
				TDirectory *SubSubdir = Subdir->GetDirectory(PathInRootFileCorr2.Data());
				gDirectory->cd(PathInRootFileCorr1.Data());
			}
			
		}
		cout << "directory?"<< endl;
		cout <<gDirectory->GetPath()<<endl;
		gDirectory->ls();
		cout << "directory?"<< endl;
		//Input->GetObject(EnergyHistoString.Data(),hEnergyMA);									// get histrogram
		hEnergyMA= (TH1D*) gDirectory->Get(EnergyHistoString.Data());
		//hEnergyMA= (TH1D*) gDirectory->Get("hEnergyCorr1_T10DeriMaxEnergyNorm_2;1");
		//Input->GetObject(EnergyHistoString.Data(),hEnergyMA);									// get histrogram
	
		//if (!SecondRun)
		//{
			//Input->GetObject("/Histograms/Energyspectrum/Energy_01;1",hEnergyMA);									// get histrogram
			//Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090_01;1",hEnergyRt1090);									// get histrogram
		//}
		//else
		//{
			//if(SecondRun==1)
			//{
				//Input->GetObject("/Histograms/Energyspectrum/EnergyCorr_01;1",hEnergyMA);									// get corrected histrogram
				//Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01;1",hEnergyRt1090);									// get corrected histrogram
			//}
			//else
			//{
				//if(SecondRun==2)
				//{
					//Input->GetObject("/Histograms/CorrCorr/EnergyCorrCorr_01;1",hEnergyMA);									// get corrected histrogram
					//Input->GetObject("/Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01;1",hEnergyRt1090);									// get corrected histrogram
				//}
			//}
		//}
		cout << "histo" << hEnergyMA << endl;
		//till here
		hEnergyMA->SetDirectory(0);
		//hEnergyRt1090->SetDirectory(0);
	Input->Close();
	//TRint *App = new TRint("ROOT",0,0,0,0,kTRUE); //&argc,argv);
	THypGeSpectrumAnalyser *Ana;
	cout << "histo" << hEnergyMA->GetEntries() << endl;
	Ana= new THypGeSpectrumAnalyser(hEnergyMA,SpektrumModel.Data(), 35 );

		Ana->SetSearchRange(500,4000);
		//Ana->SetOutputPath(Path);
		Ana->SetTxtFileOutputName(OutPutFileTxt);
		Ana->SetRootFileOutputName(OutPutFileRoot,PathInRootFileCorr1,PathInRootFileCorr2);
		if(UseFreeSkewedFitting==1)
		{
			Ana->SetNewFunctionFitting();
		}
		else{
			Ana->SetGaussianFitting();
		}
		Ana->AnalyseSpectrum();
		cout << "FWHM 1332 kEv:\t\t"<< Ana->GetFWHMCo() << endl;
	//proper output would be nice
	
	TSystem *sys=new TSystem();
	TString TxtCollectionName= sys->DirName(OutPutFileRoot.Data());
	double RadiationPerStep;			//in Panda days
	
		if(OutPutFileRoot.Contains("June2014"))
		{
			TxtCollectionName+="/CollectionFileJune2014";
			RadiationPerStep=93.8/14.;
		}
		if(OutPutFileRoot.Contains("July2014"))
		{
			TxtCollectionName+="/CollectionFileJuly2014";
			RadiationPerStep=59/10.;
		}
		TxtCollectionName+=PathInRootFileCorr1;
		TxtCollectionName+=PathInRootFileCorr2;
		TxtCollectionName+=".txt";
		
	cout <<TxtCollectionName.Data()<< endl;
	
	TString OutputRootfileBasename= sys->BaseName(OutPutFileRoot.Data());
	OutputRootfileBasename.ReplaceAll("Dataset","Dataset ");
	OutputRootfileBasename.ReplaceAll("_"," _");
	
	int iDataset=-1;
		cout <<sscanf(OutputRootfileBasename.Data(),"%*s %i %*s.root",&iDataset) << endl;
						//sTreeCOSYJune2014Dataset0_200,100,0,5339_SR0.root
	
	cout << iDataset << endl;
	cout << iDataset* RadiationPerStep<< endl;
	
	
	ofstream TxtOutputFile(TxtCollectionName.Data(),std::ofstream::out|std::ofstream::app);
	
	cout	<< iDataset << "\t" << iDataset* RadiationPerStep << "\t" << Ana->GetFWHMCo() <<"\t"<< Ana->GetFWTMCo() << "\t"<<	Ana->GetFWHM511() <<"\t"<< Ana->GetFWTM511()<<endl;
	TxtOutputFile	<< iDataset << "\t" << iDataset* RadiationPerStep << "\t" << Ana->GetFWHMCo() <<"\t"<< Ana->GetFWTMCo() << "\t"<<	Ana->GetFWHM511() <<"\t"<< Ana->GetFWTM511()<<endl;
	TxtOutputFile.close();
	return 0;
}
