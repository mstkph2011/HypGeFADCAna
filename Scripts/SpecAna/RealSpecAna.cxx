#include "TH1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
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
	//TString BeamTimeMonth = "july2014";
	TString BeamTimeMonth = "COSYnewMogon";
	if (argc == 1)
		BeamTimeMonth=argv[0];
		
	cout <<BeamTimeMonth.Data() << endl;
	return 0;
	TFile *InFile[1000];
  int i = 0;
  int UseFreeSkewedFitting=1;			//Variable ob FreeSkewedFitting benutzt wird 1=JA, 0=Nein @Torben Rathmann
  int UseMovingAv=1;					//Variable/Einstellung ob MA genutz wird 1=JA , 0=Nein  @Torben Rathmann
  float tauArray[1000];				//Arrays zum eintragen der tau's bzw FWHM(tau) @Torben Rathmann

  float aufArray511[1000];				//Array für FWHM 511 annihilation
  float aufArray511ERR[1000];
  float TenthArray511[1000];			//Array für FWTM
  float TenthArray511ERR[1000];
  float VerschlArray511[1000];
  float VerschlArray511ERR[1000];

  float aufArrayCo1[1000];				//Array für FWHM 60Co peak1
  float aufArrayCo1ERR[1000];
  float TenthArrayCo1[1000];
  float TenthArrayCo1ERR[1000];
  float VerschlArrayCo1[1000];
  float VerschlArrayCo1ERR[1000];

  float aufArrayCo[1000];				//Array für FWHM 60Co Peak 2
  float aufArrayCoERR[1000];			//Array für FWHM
  float TenthArrayCo[1000];
  float TenthArrayCoERR[1000];
  float VerschlArrayCo[1000];
  float VerschlArrayCoERR[1000];

  float aufArrayAl[1000];				//Array für FWHM 27Al
  float aufArrayAlERR[1000];
  float TenthArrayAl[1000];
  float TenthArrayAlERR[1000];
  float VerschlArrayAl[1000];
  float VerschlArrayAlERR[1000];

  float MWDArray[1000];
  float MALArray[1000];
  float iArray[1000];
  float zeroArray[1000];
  float geraten[1000];
  
  
  TString COSYTESTANADIR= getenv("COSYTESTANADIR");
  TString InputListPath = COSYTESTANADIR;
  //InputListPath +="/COSY/txtfiles/";
  InputListPath +="/";
  InputListPath +=BeamTimeMonth;
  InputListPath +="/txtfiles/";
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
		//for(int iii= 0; iii <10 ; iii++)
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
			TString ComparisonString = COSYTESTANADIR +"/" + BeamTimeMonth "/CombinedData/COSY_Ana_%i_run%i___%i,%i,%i,%i,%i.root"; //XXXXXXXXXXXXXXXXXXXXXXXXXHier fehlt was
			cout << "File:\t\t\t\t"<< ComparisonString << endl;
			sscanf(InFileName.Data(),ComparisonString.Data(),&Date,&runNo,&M,&MAL,&FilterType,&SigmaGaus,&tau);
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

			if (!UseMovingAv)
				if (!hEnergy)
					continue;
			if (UseMovingAv)
				if (!hEnergyMA)
					continue;
			cout <<hEnergy->GetMaximum() << endl;
	
			TString Path, InfileName,RootFilename, TxtFilename;
			Int_t LengthOfPath;
			
			Path = COSYTESTANADIR;
				Path += "/";
				Path += BeamTimeMonth;
				Path +="/CombinedData/Fit/";
			RootFilename = "Fitted_COSY_Ana_";
				RootFilename+=Date;
				RootFilename+= "_";
				RootFilename+=runNo;
				RootFilename+= "___";
				RootFilename+=M;
				RootFilename+= ",";
				RootFilename+=MAL;
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
				if(UseFreeSkewedFitting==1){	//Anfügen "FSFit" vom Dateinamen wenn FreeSkewed genutzt wird @Torben Rathmann
					RootFilename+= ",";
					RootFilename+= "FSFit";
				}
				RootFilename += ".root";
				cout << RootFilename.Data() << endl;
			TxtFilename = "Fitted_COSY_Ana";
				TxtFilename+=M;
				TxtFilename+= ",";
				TxtFilename+=MAL;
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
				if(UseFreeSkewedFitting==1){	//Anfügen "FSFit" vom Dateinamen wenn FreeSkewed genutzt wird @Torben Rathmann
					TxtFilename+= ",";
					TxtFilename+= "FSFit";
				}
				TxtFilename+= ".txt";

				THypGeSpectrumAnalyser *Ana;
			if(UseMovingAv==1){ 						//Überschreibt *Ana um Moving Average zu benutzen @Torben Rathmann
				Ana= new THypGeSpectrumAnalyser(hEnergyMA,"jülich2", 35 );
			}
			else
				Ana = new THypGeSpectrumAnalyser(hEnergy,"jülich2", 35 );
			
			Ana->SetSearchRange(500,4000);
			Ana->SetOutputPath(Path);
			Ana->SetTxtFileOutputName(TxtFilename);
			Ana->SetRootFileOutputName(RootFilename);
			if(UseFreeSkewedFitting==1){
				//Ana->SetFreeSkewedFitting();
				Ana->SetNewFunctionFitting();
			}
			else{
				Ana->SetGaussianFitting();
			}
			//Ana->SetFreeSkewedFitting();
			//Ana->SetSecondGausianFitting();
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
			MWDArray[i]=M;
			MALArray[i]=MAL;
			iArray[i]=i*4E12;
			zeroArray[i]=0;
			geraten[i]=0.03;

			//peak 511annihi
			aufArray511[i]=Ana->GetFWHM511();
			//aufArray511ERR[i]=Ana->	;
			TenthArray511[i]=Ana->GetFWTM511();
			//TenthArray511ERR[i]=Ana->	;
			VerschlArray511[i]=Ana->GetFWTM511()/Ana->GetFWHM511();
			//VerschlArray511ERR[i]= ;


			//peak1 60Co
			aufArrayCo1[i]=Ana->GetFWHMCo1();
			//aufArrayCo1ERR[i]=Ana->	;
			TenthArrayCo1[i]=Ana->GetFWTMCo1();
			//TenthArrayCo1ERR[i]=Ana-> ;
			VerschlArrayCo1[i]=Ana->GetFWTMCo1()/Ana->GetFWHMCo1();
			//VerschlArrayCo1ERR[i]= ;

			//peak2 60Co
			aufArrayCo[i]=Ana->GetFWHMCo();
			//aufArrayCoERR[i]=Ana-> ;
			TenthArrayCo[i]=Ana->GetFWTMCo();
			//TenthArrayCoERR[i]=Ana-> ;
			VerschlArrayCo[i]=Ana->GetFWTMCo()/Ana->GetFWHMCo();
			//VerschlArrayCoErr[i]=;

			//peak 27Al
			aufArrayAl[i]=Ana->GetFWHMAl();
			//aufArrayAlERR[i]=Ana->	;
			TenthArrayAl[i]=Ana->GetFWTMAl();
			//TenthArrayAlERR[i]=Ana->	;
			VerschlArrayAl[i]=Ana->GetFWTMAl()/Ana->GetFWHMAl();
			//VerschlArrayAlERR[i]=;


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

	TGraphErrors *FWHM511= new TGraphErrors(i,&iArray[0],&aufArray511[0],&zeroArray[0],&geraten[0]);//&aufArray511ERR[0]);				//Graphen von FWHM über Protonen
	TGraphErrors *FWHMCo1= new TGraphErrors(i,&iArray[0],&aufArrayCo1[0],&zeroArray[0],&geraten[0]);//&aufArrayCo1ERR[0]);
	TGraphErrors *FWHMCo = new TGraphErrors(i,&iArray[0],&aufArrayCo[0] ,&zeroArray[0],&geraten[0]);//&aufArrayCoERR[0] );
	TGraphErrors *FWHMAl = new TGraphErrors(i,&iArray[0],&aufArrayAl[0] ,&zeroArray[0],&geraten[0] );//&aufArrayAlERR[0] );

	TGraphErrors *FWTMHM511= new TGraphErrors(i,&iArray[0],&VerschlArray511[0],&zeroArray[0],&geraten[0]);//&VerschlArray511ERR[0]);	//Graphen von FWTM/FWHM über Protonen
	TGraphErrors *FWTMHMCo1= new TGraphErrors(i,&iArray[0],&VerschlArrayCo1[0],&zeroArray[0],&geraten[0]);//&VerschlArrayCo1ERR[0]);
	TGraphErrors *FWTMHMCo = new TGraphErrors(i,&iArray[0],&VerschlArrayCo[0] ,&zeroArray[0],&geraten[0]);//&VerschlArrayCoERR[0]);
	TGraphErrors *FWTMHMAl = new TGraphErrors(i,&iArray[0],&VerschlArrayAl[0] ,&zeroArray[0],&geraten[0]);//&VerschlArrayAlERR[0]);

	/*TGraph *verschlechtFWHM = new TGraph(i,&iArray[0],&aufArrayCo[0]);
	TGraph *verschlechtFWTM = new TGraph(i,&iArray[0],&TenthArrayCo[0]);
	TGraph *verschlecht = new TGraph(i,&iArray[0],&VerschlArrayCo[0]);
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

	FWTMHM511->Write("FWTM /FWHM 511");
	FWTMHMCo1->Write("FWTM /FWHM Co1");
	FWTMHMCo ->Write("FWTM /FWHM Co2");
	FWTMHMAl ->Write("FWTM /FWHM Al");

	verhaeltout->Close();

	//for (int j = 0; j < i; j++)
	////	cout << tauArray[j][0] << "\t" << tauArray[j][1] << endl;
	//}
	//@Torben Rathmann
	/* TGraph *tauAuf = new TGraph(i,&tauArray[0],&aufArray[0]);	//Graph zur Bestimmung FWHM_min(tau)
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
	TGraph2D *MLAuf = new TGraph2D(i, &MWDArray[0], &MALArray[0], &aufArray[0]);			//2dim graph für M,L,FWHM
	TString MLAufloesName="MLAufloesung_tau6210";
	MLAuf->SetNameTitle("FWHM (M,L) [keV]","FWHM (M,L) [keV]");
	if(UseMovingAv==1){
		MLAufloesName+= "MA";
	}
	MLAufloesName+=".root";
	TFile *outi = new TFile(MLAufloesName,"RECREATE");
	MLAuf->Write("MLAufloesung_tau6210");
	outi->Close();
	
	
	
	TGraph2D *MLTenthAuf = new TGraph2D(i, &MWDArray[0], &MALArray[0], &TaufArray[0]);		//2dim graph für M,L,FWTM
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
		if(aufArray[j]<4.3 && TaufArray[j]<7.9){
			cout<<"tau:"<<tauArray[j]<<"; M:"<<MWDArray[j]<<"; L:"<<MALArray[j]<<"; aufloesung:"<<aufArray[j]<<"; TenthAufloesung:"<<TaufArray[j]<<endl;
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
	//App->Run();
	return 0;
}
