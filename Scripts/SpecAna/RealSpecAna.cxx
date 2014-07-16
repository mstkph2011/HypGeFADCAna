#include "TH1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TRint.h"							// starts the interactive mode of ROOT (needed for canvas etc.)

#include "THypGeSpectrumAnalyser.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

int main(int argc, char* argv[] )
{
	
	//TRint *App = new TRint("ROOT",0,0); //&argc,argv);

	TFile *InFile[1000];
  int i = 0;
  int UseMovingAv=1;					//Variable/Einstellung ob MA genutz wird 1=JA , 0=Nein  @Torben Rathmann
  float tauArray[1000];				//Arrays zum eintragen der tau's bzw FWHM(tau) @Torben Rathmann
  float aufArray[1000];
  
  TString COSYTESTANADIR= getenv("COSYTESTANADIR");
  TString InputListPath = COSYTESTANADIR;
  //InputListPath +="/COSY/txtfiles/";
  InputListPath +="/june2014/txtfiles/";
  TString InputFile = InputListPath + "Filestofit.txt";
  TString CorruptedFileList = InputListPath +"CorruptedFiles.txt";
	//ifstream InputList("files.txt");
	cout << InputFile.Data() << endl;
	ifstream InputList(InputFile.Data());
	//ifstream InputList("files2.txt");
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
			int M, FilterType,SigmaGaus,SigmaBil, tau;
			TString ComparisonString = COSYTESTANADIR + "/june2014/CombinedData/COSY_Ana%i,%i,%i,%i.root";
			cout << "File:\t\t\t\t"<< ComparisonString << endl;
			sscanf(InFileName.Data(),ComparisonString.Data(),&M,&FilterType,&SigmaGaus,&tau);
			InFile[i] = new TFile(InFileName);								//open ROOT file
			//InFile[i] = new TFile("COSY_Ana200,4,11,300.root");								//open ROOT file
			cout << InFile[i]->GetSize() << endl;
			if (InFile[i]->GetSize() < 600)
			{
				CorruptedList << InFileName.Data() << endl;
				cout << "corrupted file" << endl;
				continue;
			}
			TH1D* hEnergy;
			TH1D* hEnergyMA;
			InFile[i]->GetObject("/Histograms/Energyspectrum/Energy_01;1",hEnergy);									// get histrogram
			InFile[i]->GetObject("/Histograms/Energyspectrum/EnergyMA_01;1",hEnergyMA);									// get histrogram
			cout <<hEnergy->GetMaximum() << endl;
	
			TString Path, InfileName,RootFilename, TxtFilename;
			Int_t LengthOfPath;
			
			Path = COSYTESTANADIR;
				Path += "/june2014/CombinedData/Fit/";
			RootFilename = "Fitted_COSY_Ana";
				RootFilename+=M;
				RootFilename+= ",";
				RootFilename+= FilterType;
				RootFilename+= ","; 
				RootFilename += SigmaGaus; 
				RootFilename+= ",";
				RootFilename+= SigmaBil; 
				RootFilename+= ",";
				RootFilename+= tau;
				if(UseMovingAv==1){				//Anfügen "MA" vom Dateinamen wenn Moving Average genutzt wird @Torben Rathmann
					RootFilename+= ",";
					RootFilename+= "MA";
				}
				RootFilename += ".root";
				cout << RootFilename.Data() << endl;
			TxtFilename = "Fitted_COSY_Ana";
				TxtFilename += ",";
				TxtFilename += FilterType;
				TxtFilename += ","; 
				TxtFilename += SigmaGaus; 
				TxtFilename+= ",";
				TxtFilename+= SigmaBil;
				TxtFilename+= ",";
				TxtFilename+= tau;
				if(UseMovingAv==1){				//Anfügen "MA" vom Dateinamen wenn Moving Average genutzt wird @Torben Rathmann
					TxtFilename+= ",";
					TxtFilename+= "MA";
				}
				TxtFilename+= ".txt";

				THypGeSpectrumAnalyser *Ana;
			if(UseMovingAv==1){ 						//Überschreibt *Ana um Moving Average zu benutzen @Torben Rathmann
				Ana= new THypGeSpectrumAnalyser(hEnergyMA,"co60", 35 );
			}
			else
				Ana = new THypGeSpectrumAnalyser(hEnergy,"co60", 35 );
			
			Ana->SetSearchRange(1000,2000);
			Ana->SetOutputPath(Path);
			Ana->SetTxtFileOutputName(TxtFilename);
			Ana->SetRootFileOutputName(RootFilename);
			Ana->SetGaussianFitting();
			//Ana->SetFreeSkewedFitting();
			Ana->SetSecondGausianFitting();
			Ana->AnalyseSpectrum();
			cout << Ana->GetFWHMCo() << endl;


			//@Torben Rathmann
			/*
			if(i+1==sizeof(tauArray)/sizeof(int)/2){						//Test ob Array voll ist  !!!Achtung bei /2 <- 2 ist die anzahl Spalten(2)
				float tmpTauEffi[sizeof(tauArray)/sizeof(int)/2+50][2];	//tmp array erstellen, Arraygröße wird um 50 erhöht
				for(int k=0; k<int(sizeof(tauArray)/sizeof(int)/2.);k++){		//zwischenspeichern der daten
					tmpTauEffi[k][0]=tauArray[k][0];
					tmpTauEffi[k][1]=tauArray[k][1];
				}
				tauArray=tmpTauEffi;						//altes Array wird benutzt
			}*/
			//cout << tau << endl << endl << endl;
			tauArray[i]=tau;								//Daten werden in Array geschrieben
			aufArray[i]=Ana->GetFWHMCo();



			//write to maps
			if (SigmaBil != 3)	// bil
			{
				ResultMapBil[(M /20 -3 )*10 + ((SigmaBil -100)/200 )-1][SigmaGaus] = Ana->GetFWHMCo();
				
			}
			if (SigmaBil == 3)	// gaus
			{
				ResultMapGaus[M /20 -3][SigmaGaus] = Ana->GetFWHMCo();
			}
			cout << "test1" << endl;
			delete Ana;
			delete hEnergy;
			InFile[i]->Close();
			delete InFile[i];
			i++;
			cout << "test2" << endl;
		}
	}
	for (int j = 0; j < i; j++)
	{
		cout << tauArray[j][0] << "\t" << tauArray[j][1] << endl;
	}
	//@Torben Rathmann
	TGraph *tauAuf = new TGraph(i,&tauArray[0],&aufArray[0]);
	TFile *rootoutput = new TFile("test.root","RECREATE");
	tauAuf->Write("tauAufloesung");
	rootoutput->Close();
	delete tauAuf;
	//

	
	TString OutputFile = COSYTESTANADIR;
				OutputFile += "/COSY/CombinedData/Fit/Output.txt";
	ofstream Output (OutputFile.Data(),std::fstream::trunc);
	
	Output << "SigmaGausRange" << endl;
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
	}
	Output.close();
	InputList.close();
	CorruptedList.close();
	//App->Run();
	return 0;
}
