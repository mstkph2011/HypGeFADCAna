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

#include "THypGeSpectrumAnalyser.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

int main(int argc, char* argv[] )
{
	Bool_t WriteEnabled;
	if (argc >1)
	{
		WriteEnabled = std::atoi(argv[1]);
	}
	else
	{
		WriteEnabled = 1;		// Default is first run
	}
	//return 0;
	TRint *App = new TRint("ROOT",0,0,0,0,kTRUE); //&argc,argv);
	int ArrayNumber=100;
	TFile *InFile[ArrayNumber];
	TH1D* hRiseTimePeak1332[ArrayNumber];
	TProfile* hXProfileRiseTimePeak1332[ArrayNumber];
	TF1* fProfileFitFuncPol[ArrayNumber];
	TF1* fProfileFitFuncConstant[ArrayNumber];
	TFile *OutputParametersFile[ArrayNumber];
  int i = 0;
  int UseFreeSkewedFitting=1;			//Variable ob FreeSkewedFitting benutzt wird 1=JA, 0=Nein @Torben Rathmann
  int UseMovingAv=1;					//Variable/Einstellung ob MA genutz wird 1=JA , 0=Nein  @Torben Rathmann
  float tauArray[ArrayNumber];				//Arrays zum eintragen der tau's bzw FWHM(tau) @Torben Rathmann

  float FWHMArray511[ArrayNumber];				//Array für FWHM 511 annihilation
  float FWHMArray511ERR[ArrayNumber];
  float FWTMArray511[ArrayNumber];			//Array für FWTM
  float FWTMArray511ERR[ArrayNumber];
  float RatioArray511[ArrayNumber];
  float RatioArray511ERR[ArrayNumber];

  float FWHMArrayCo1[ArrayNumber];				//Array für FWHM 60Co peak1
  float FWHMArrayCo1ERR[ArrayNumber];
  float FWTMArrayCo1[ArrayNumber];
  float FWTMArrayCo1ERR[ArrayNumber];
  float RatioArrayCo1[ArrayNumber];
  float RatioArrayCo1ERR[ArrayNumber];

  float FWHMArrayCo[ArrayNumber];				//Array für FWHM 60Co Peak 2
  float FWHMArrayCoERR[ArrayNumber];			//Array für FWHM
  float FWTMArrayCo[ArrayNumber];
  float FWTMArrayCoERR[ArrayNumber];
  float RatioArrayCo[ArrayNumber];
  float RatioArrayCoERR[ArrayNumber];

  float FWHMArrayAl[ArrayNumber];				//Array für FWHM 27Al
  float FWHMArrayAlERR[ArrayNumber];
  float FWTMArrayAl[ArrayNumber];
  float FWTMArrayAlERR[ArrayNumber];
  float RatioArrayAl[ArrayNumber];
  float RatioArrayAlERR[ArrayNumber];

  float MWDArray[ArrayNumber];
  float MALArray[ArrayNumber];
  float iArray[ArrayNumber];
  float zeroArray[ArrayNumber];
  float geraten[ArrayNumber];
  
  
  TString COSYTESTANADIR= getenv("COSYTESTANADIR");
  TString InputListPath = COSYTESTANADIR;
  InputListPath +="/june2014/txtfiles/";
  TString InputFile = InputListPath + "Filestofit.txt";
  TString CorruptedFileList = InputListPath +"CorruptedFiles.txt";

	ifstream InputList(InputFile.Data());

	ofstream CorruptedList (CorruptedFileList.Data(),std::fstream::trunc);
	
	Int_t MRange [13] = {60,80,100,120,140,160,180,200,220,240,260,280,300};
	Int_t SigmaBilRange [10] = {300,700,900,1100,1300,1500,1700,1900,2100};
	Int_t SigmaGausRange [6] = {1,3,5,7,9,11};
	
	//create maps here!!!		130 for bil filter + 13 for gaus filter	, first 10 for first M value and running SigmaBil
	std::map<Int_t,Double_t> ResultMapBil[130];
	std::map<Int_t,Double_t> ResultMapGaus[13];
	
	if(!InputList)
	{
		cout << "ERROR - file does NOT exist" << endl;
	}
	else
	{
	  while (InputList.good())																//loop over all lines
		{
			char buffer[200] = "";
			InputList.getline(buffer,200);												//read name of a ROOT file
			cout << "File:\t\t\t\t"<<buffer << endl;
			TString InFileName = buffer; 
			cout << "BufferLength: \t" << InFileName.Length() << endl;
			if (InFileName.Length() == 0)
				continue;
			//extract parameter values from file name
			int M, FilterType,SigmaGaus,SigmaBil, tau, MAL, Date, runNo;
			TString ComparisonString = COSYTESTANADIR + "/june2014/CombinedData/COSY_Ana_%i_run%i___%i,%i,%i,%i,%i.root";
			cout << "File:\t\t\t\t"<< ComparisonString << endl;
			sscanf(InFileName.Data(),ComparisonString.Data(),&Date,&runNo,&M,&MAL,&FilterType,&SigmaGaus,&tau);
			InFile[i] = new TFile(InFileName);								//open ROOT file
			if (InFile[i]->GetSize() < 600)
			{
				CorruptedList << InFileName.Data() << endl;
				cout << "corrupted file" << endl;
				continue;
			}

			TH1D* hEnergyMA;
			TH2D* hEnergyRt1090;

			InFile[i]->GetObject("/Histograms/Energyspectrum/Energy_01;1",hEnergyMA);									// get histrogram
			InFile[i]->GetObject("/Histograms/EnergyRise1090Corr/EnergyRise1090Corr_01;1",hEnergyRt1090);									// get histrogram

			if (UseMovingAv)
				if (!hEnergyMA)
					continue;
	
			TString Path, InfileName,RootFilename, TxtFilename;
			Int_t LengthOfPath;
			
			Path = COSYTESTANADIR;
				Path += "/june2014/CombinedData/Fit/";
			char buf[100];
				sprintf(buf,"Fitted_COSY_Ana_%i_%i___%i,%i,%i,%i,%i,%i,MA",Date,runNo,M,MAL,FilterType,SigmaGaus,SigmaBil,tau);
			RootFilename = buf;
				if(UseFreeSkewedFitting==1){	//Anfügen "FSFit" vom Dateinamen wenn FreeSkewed genutzt wird @Torben Rathmann
					RootFilename+= ",";
					RootFilename+= "FSFit";
				}
				RootFilename += ".root";
				cout << RootFilename.Data() << endl;
			TxtFilename = buf;
				if(UseFreeSkewedFitting==1){	//Anfügen "FSFit" vom Dateinamen wenn FreeSkewed genutzt wird @Torben Rathmann
					TxtFilename+= ",";
					TxtFilename+= "FSFit";
				}
				TxtFilename+= ".txt";

			THypGeSpectrumAnalyser *Ana;
			Ana= new THypGeSpectrumAnalyser(hEnergyMA,"jülich2", 35 );
			Ana->SetSearchRange(500,4000);
			Ana->SetOutputPath(Path);
			Ana->SetTxtFileOutputName(TxtFilename);
			Ana->SetRootFileOutputName(RootFilename);
			if(UseFreeSkewedFitting==1)
			{
				Ana->SetNewFunctionFitting();
			}
			else{
				Ana->SetGaussianFitting();
			}
			Ana->AnalyseSpectrum();
			cout << "FWHM 1332 kEv:\t\t"<< Ana->GetFWHMCo() << endl;

			Double_t* PeakRange1332 = Ana->GetPeakRangeChannels(2);
			TNtupleD *tNtupleD = new TNtupleD("tNtupleD","Tree with vectors","RangeMin:RangeMax");
				tNtupleD->Fill(PeakRange1332);
									//cout << "PeakLow "<< PeakRange1332[0] << "\t\t" << "PeakHigh" << PeakRange1332[1] << endl;
			sprintf(buf,"hRiseTimePeak1332_run%i",i);
			// extract risetime histogram of only 1332 keV Peak range from 2d histogram
			hRiseTimePeak1332[i]= hEnergyRt1090->ProjectionX(buf,PeakRange1332[0],PeakRange1332[1],"");

			// do XProfile of Correlation histogram to get correction functions
			sprintf(buf,"hXProfileRiseTimePeak1332_run%i",i);
			hXProfileRiseTimePeak1332[i] = hEnergyRt1090->ProfileX(buf,PeakRange1332[0],PeakRange1332[1],"");
			sprintf(buf,"fProfileFitFuncPol_run%i",i);
			fProfileFitFuncPol[i] = new TF1(buf, "pol9", 0,400);
			hXProfileRiseTimePeak1332[i]->Fit(fProfileFitFuncPol[i],"R");
			sprintf(buf,"fProfileFitFuncConstant_run%i",i);
			fProfileFitFuncConstant[i]= new TF1(buf, "pol0", 0,400);
			fProfileFitFuncConstant[i]->SetLineColor(kBlue);
			hXProfileRiseTimePeak1332[i]->Fit(fProfileFitFuncConstant[i],"R+");



				if (WriteEnabled)
				{
				// write extracted parameters to file
				sprintf(buf,"ParametersFirstAnaStepCOSY_Ana_%i_%i___%i,%i,%i,%i,%i,%i,MA.root",Date,runNo,M,MAL,FilterType,SigmaGaus,SigmaBil,tau);
				TString OutputParametersFileName = COSYTESTANADIR + "/june2014/DatabaseFirstAnalysisStep/" +buf;
				OutputParametersFile[i] = new TFile(OutputParametersFileName,"RECREATE");
					hEnergyRt1090->Write();
					hRiseTimePeak1332[i]->Write("hRiseTimePeak1332");
					hXProfileRiseTimePeak1332[i]->Write("hXProfileRiseTimePeak1332");
					fProfileFitFuncPol[i]->Write("fProfileFitFuncPol");
					fProfileFitFuncConstant[i]->Write("fProfileFitFuncConstant");
					tNtupleD->Write();
				OutputParametersFile[i]->Close();
				}

			tauArray[i]=tau;								//Daten werden in Array geschrieben
			MWDArray[i]=M;
			MALArray[i]=MAL;
			iArray[i]=i*4E12;
			zeroArray[i]=0;
			geraten[i]=0.03;

			//peak 511annihi
			FWHMArray511[i]=Ana->GetFWHM511();
			//FWHMArray511ERR[i]=Ana->	;
			FWTMArray511[i]=Ana->GetFWTM511();
			//FWTMArray511ERR[i]=Ana->	;
			RatioArray511[i]=Ana->GetFWTM511()/Ana->GetFWHM511();
			//RatioArray511ERR[i]= ;


			//peak1 60Co
			FWHMArrayCo1[i]=Ana->GetFWHMCo1();
			//FWHMArrayCo1ERR[i]=Ana->	;
			FWTMArrayCo1[i]=Ana->GetFWTMCo1();
			//FWTMArrayCo1ERR[i]=Ana-> ;
			RatioArrayCo1[i]=Ana->GetFWTMCo1()/Ana->GetFWHMCo1();
			//RatioArrayCo1ERR[i]= ;

			//peak2 60Co
			FWHMArrayCo[i]=Ana->GetFWHMCo();
			//FWHMArrayCoERR[i]=Ana-> ;
			FWTMArrayCo[i]=Ana->GetFWTMCo();
			//FWTMArrayCoERR[i]=Ana-> ;
			RatioArrayCo[i]=Ana->GetFWTMCo()/Ana->GetFWHMCo();
			//RatioArrayCoErr[i]=;

			//peak 27Al
			FWHMArrayAl[i]=Ana->GetFWHMAl();
			//FWHMArrayAlERR[i]=Ana->	;
			FWTMArrayAl[i]=Ana->GetFWTMAl();
			//FWTMArrayAlERR[i]=Ana->	;
			RatioArrayAl[i]=Ana->GetFWTMAl()/Ana->GetFWHMAl();
			//RatioArrayAlERR[i]=;

			//write to maps

			if (SigmaBil != 3)	// bil
			{
				ResultMapBil[(M /20 -3 )*10 + ((SigmaBil -100)/200 )-1][SigmaGaus] = Ana->GetFWHMCo();
			}
			if (SigmaBil == 3)	// gaus
			{
				ResultMapGaus[M /20 -3][SigmaGaus] = Ana->GetFWHMCo();
			}
			delete Ana;
			delete hEnergyMA;
			InFile[i]->Close();
			delete InFile[i];
			i++;
		}
	}

	TGraphErrors *FWHM511= new TGraphErrors(i,&iArray[0],&FWHMArray511[0],&zeroArray[0],&geraten[0]);//&FWHMArray511ERR[0]);				//Graphen von FWHM über Protonen
	TGraphErrors *FWHMCo1= new TGraphErrors(i,&iArray[0],&FWHMArrayCo1[0],&zeroArray[0],&geraten[0]);//&FWHMArrayCo1ERR[0]);
	TGraphErrors *FWHMCo = new TGraphErrors(i,&iArray[0],&FWHMArrayCo[0] ,&zeroArray[0],&geraten[0]);//&FWHMArrayCoERR[0] );
	TGraphErrors *FWHMAl = new TGraphErrors(i,&iArray[0],&FWHMArrayAl[0] ,&zeroArray[0],&geraten[0] );//&FWHMArrayAlERR[0] );

	TGraphErrors *FWTMHM511= new TGraphErrors(i,&iArray[0],&RatioArray511[0],&zeroArray[0],&geraten[0]);//&RatioArray511ERR[0]);	//Graphen von FWTM/FWHM über Protonen
	TGraphErrors *FWTMHMCo1= new TGraphErrors(i,&iArray[0],&RatioArrayCo1[0],&zeroArray[0],&geraten[0]);//&RatioArrayCo1ERR[0]);
	TGraphErrors *FWTMHMCo = new TGraphErrors(i,&iArray[0],&RatioArrayCo[0] ,&zeroArray[0],&geraten[0]);//&RatioArrayCoERR[0]);
	TGraphErrors *FWTMHMAl = new TGraphErrors(i,&iArray[0],&RatioArrayAl[0] ,&zeroArray[0],&geraten[0]);//&RatioArrayAlERR[0]);

	/*TGraph *verschlechtFWHM = new TGraph(i,&iArray[0],&FWHMArrayCo[0]);
	TGraph *verschlechtFWTM = new TGraph(i,&iArray[0],&FWTMArrayCo[0]);
	TGraph *verschlecht = new TGraph(i,&iArray[0],&RatioArrayCo[0]);
	*/
	TString VerschName= "FWTMFWHM";
	if(UseMovingAv==1){
			VerschName+= "MA";
		}
	VerschName+=".root";
	TFile *verhaeltout=new TFile(VerschName,"RECREATE");
	/*verschlecht->Write("FWTMFWHM");
	verschlechtFWHM->Write("FWHM");
	verschlechtFWTM->Write("FWTM");
	*/
	FWHM511->Write("FWHM511");
	FWHMCo1->Write("FWHMCo1");
	FWHMCo ->Write("FWHMCo2");
	FWHMAl ->Write("FWHMAl");
	cout << i << endl;
	for (int ii=0; ii < i; ii++)
		hRiseTimePeak1332[ii]->Write();

	FWTMHM511->Write("FWTM /FWHM 511");
	FWTMHMCo1->Write("FWTM /FWHM Co1");
	FWTMHMCo ->Write("FWTM /FWHM Co2");
	FWTMHMAl ->Write("FWTM /FWHM Al");

	verhaeltout->Close();

	//for (int j = 0; j < i; j++)
	////	cout << tauArray[j][0] << "\t" << tauArray[j][1] << endl;
	//}
	//@Torben Rathmann
	/* TGraph *tauAuf = new TGraph(i,&tauArray[0],&FWHMArray[0]);	//Graph zur Bestimmung FWHM_min(tau)
	TString tauAufloesName ="tauAufloesung";
	if(UseMovingAv==1){
		tauAufloesName+= "MA";
	}
	tauAufloesName+=".root";
	TFile *rootoutput = new TFile(tauAufloesName,"RECREATE"); //zum direkten Vergleich ob tau bei MA besser oder schlechter wird
	tauAuf->Write("tauAufloesung");
	rootoutput->Close();
	*/

	/*
	TGraph2D *MLAuf = new TGraph2D(i, &MWDArray[0], &MALArray[0], &FWHMArray[0]);			//2dim graph für M,L,FWHM
	TString MLAufloesName="MLAufloesung_tau6210";
	MLAuf->SetNameTitle("FWHM (M,L) [keV]","FWHM (M,L) [keV]");
	if(UseMovingAv==1){
		MLAufloesName+= "MA";
	}
	MLAufloesName+=".root";
	TFile *outi = new TFile(MLAufloesName,"RECREATE");
	MLAuf->Write("MLAufloesung_tau6210");
	outi->Close();
	
	
	
	TGraph2D *MLTenthAuf = new TGraph2D(i, &MWDArray[0], &MALArray[0], &TFWHMArray[0]);		//2dim graph für M,L,FWTM
	TString MLTenthAufloesName="MLTenthAufloesung_tau6210";
	if(UseMovingAv==1){
		MLTenthAufloesName+= "MA";
	}
	MLTenthAufloesName+=".root";
	TFile *TenthOut = new TFile(MLTenthAufloesName,"RECREATE");
	MLTenthAuf->Write("MLTenthAufloesung_tau6210");
	TenthOut->Close();
	*/
	

	/*for(int j=0;j<i;j++){					//Ausgabe diverser Werte zum identifizieren guter M,L Werte
		if(FWHMArray[j]<4.3 && TFWHMArray[j]<7.9){
			cout<<"tau:"<<tauArray[j]<<"; M:"<<MWDArray[j]<<"; L:"<<MALArray[j]<<"; aufloesung:"<<FWHMArray[j]<<"; TenthAufloesung:"<<TFWHMArray[j]<<endl;
		}	
	}*/
	
	//

	
	TString OutputFile = COSYTESTANADIR;
				OutputFile += "/COSY/CombinedData/Fit/Output.txt";
	ofstream Output (OutputFile.Data(),std::fstream::trunc);
	
	/*Output << "SigmaGausRange" << endl;
	for (int i = 0; i<6;i++)
	{
		cout << SigmaGausRange[i] << "\t";
		Output << SigmaGausRange[i] << "\t";
	}
	cout << endl;
	Output << endl;
	cout << "bil" <<endl;
	Output << "bil" <<endl;
	for (int i = 0; i<130; i++)
	{
		cout << "M\t SigmaBil: \t\t\t\t" <<  (int(i/10)+3)*20 <<"\t" << ((i % 10)+1)*200+100 << endl;
		Output << "M\t SigmaBil: \t\t\t\t" << (int(i/10)+3)*20 <<"\t" << ((i % 10)+1)*200+100 << endl;
		
		for (std::map<Int_t,Double_t>::iterator it=ResultMapBil[i].begin(); it!=ResultMapBil[i].end(); ++it)
    {
			cout << it->second << "\t";
			Output << it->second << "\t";
		}
		cout << endl;
		Output << endl;
	}
	cout << "gaus" <<endl;
	Output << "gaus" <<endl;
	
	for (int i = 0; i<13; i++)
	{
		cout << "M: \t\t\t\t" <<  (i+3)*20  << endl;
		Output << "M: \t\t\t\t" << (i+3)*20  << endl;
		
		for (std::map<Int_t,Double_t>::iterator it=ResultMapGaus[i].begin(); it!=ResultMapGaus[i].end(); ++it)
    {
			cout << it->second << "\t";
			Output << it->second << "\t";
		}
		cout << endl;
		Output << endl;
	}*/
	Output.close();
	InputList.close();
	CorruptedList.close();
	App->Run();
	return 0;
}
