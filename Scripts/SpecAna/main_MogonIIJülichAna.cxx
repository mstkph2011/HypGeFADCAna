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
	TString UseT10DMCut="";
	
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
	if (argc>=6)
	{
		PathInRootFileCorr2=argv[5];
		cout << PathInRootFileCorr2.Data()<<endl;
	}
	if (argc>=7)
	{
		UseT10DMCut=argv[6];
		cout << UseT10DMCut.Data()<<endl;
	}
	TH1D* hEnergyMA;
	TH2D* hEnergyT10DM;
	
	TString EnergyHistoString="hEnergy;1";
	EnergyHistoString.Form("hEnergy%s%s;1",PathInRootFileCorr2.Data(),PathInRootFileCorr1.Data());
	TString EnergyT10DMHistoString="hT10DeriMaxEnergy;1";
	if(!UseT10DMCut.IsNull())
	{
		EnergyT10DMHistoString.Form("hT10DeriMaxEnergy%s%s;1",PathInRootFileCorr2.Data(),PathInRootFileCorr1.Data());
	}
	cout << EnergyHistoString.Data() << endl;
	cout << EnergyT10DMHistoString.Data() << endl;
	//return 0;
	//change here for the proper input
	TFile *Input= new TFile(InputFile.Data());
		if(!PathInRootFileCorr1.IsNull())
		{
			TDirectory *Subdir = Input->GetDirectory(PathInRootFileCorr1.Data());
			gDirectory->cd(PathInRootFileCorr1.Data());
			if(!PathInRootFileCorr2.IsNull())
			{
				cout << "test" << endl;
				TDirectory *SubSubdir = gDirectory->GetDirectory(PathInRootFileCorr2.Data());
				cout << "SSd" << SubSubdir<< endl;
				gDirectory->cd(PathInRootFileCorr2.Data());
			}
			
		}
		//cout << "directory?"<< endl;
		//cout <<gDirectory->GetPath()<<endl;

		//gDirectory->ls();
		//cout << "directory?"<< endl;

		hEnergyMA= (TH1D*) gDirectory->Get(EnergyHistoString.Data());
		if(UseT10DMCut.IsNull())
		{
			hEnergyMA= (TH1D*) gDirectory->Get(EnergyHistoString.Data());
		}
		else
		{
			hEnergyT10DM=(TH2D*) gDirectory->Get(EnergyT10DMHistoString.Data());
			hEnergyMA= (TH1D*) hEnergyT10DM->ProjectionY("_py",35,-1);
		}
		cout << "histo" << hEnergyMA << endl;
		cout << "histo2d" << hEnergyT10DM << endl;
	
		//till here
		hEnergyMA->SetDirectory(0);
		hEnergyT10DM->SetDirectory(0);
		//hEnergyRt1090->SetDirectory(0);
	Input->Close();
	//TRint *App = new TRint("ROOT",0,0,0,0,kTRUE); //&argc,argv);
	
	cout <<hEnergyMA->GetEntries()<<endl;
		
	
	THypGeSpectrumAnalyser *Ana;
	cout << "histo" << hEnergyMA->GetEntries() << endl;
	Ana= new THypGeSpectrumAnalyser(hEnergyMA,SpektrumModel.Data(), 35 );

		Ana->SetSearchRange(500,4000);
		//Ana->SetOutputPath(Path);
		Ana->SetTxtFileOutputName(OutPutFileTxt);
		if(!UseT10DMCut.IsNull())
		{
			PathInRootFileCorr2+=UseT10DMCut;
		}
		Ana->SetRootFileOutputName(OutPutFileRoot,PathInRootFileCorr1,PathInRootFileCorr2);
		if(UseFreeSkewedFitting==1)
		{
			Ana->SetNewFunctionFitting();
		}
		else{
			Ana->SetGaussianFitting();
		}
		//Ana->SetSecondGausianFitting();
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
						
	
	cout << iDataset << endl;
	cout << iDataset* RadiationPerStep<< endl;
	
	
	ofstream TxtOutputFile(TxtCollectionName.Data(),std::ofstream::out|std::ofstream::app);
	
	cout	<< iDataset << "\t" << iDataset* RadiationPerStep << "\t" << Ana->GetFWHMCo() <<"\t"<< Ana->GetFWTMCo() << "\t"<<	Ana->GetFWHM511() <<"\t"<< Ana->GetFWTM511()<<endl;
	TxtOutputFile	<< iDataset << "\t" << iDataset* RadiationPerStep << "\t" << Ana->GetFWHMCo() <<"\t"<< Ana->GetFWTMCo() << "\t"<<	Ana->GetFWHM511() <<"\t"<< Ana->GetFWTM511()<<endl;
	TxtOutputFile.close();
	return 0;
}
