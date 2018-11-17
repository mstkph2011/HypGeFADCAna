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
	Bool_t WriteEnabled;
	cout << "argc"<< argc << endl;
	if (argc >1)
	{
		WriteEnabled = std::atoi(argv[1]);
	}
	else
	{
		WriteEnabled = 1;		// Default is first run
	}
	//TRint *App = new TRint("ROOT",0,0); //&argc,argv);
		TString BeamTimeMonth = "july2014";
		//if (argc == 1)
			//BeamTimeMonth=argv[0];

		cout <<BeamTimeMonth.Data() << endl;
	bool SR = 0;// second run? corrected histograms will be used if true
	//return 0;
	TRint *App = new TRint("ROOT",0,0,0,0,kTRUE); //&argc,argv);
	int ArrayNumber=100;
	int PeakNumber =4;
	int PolOrder = 5;		//Order of Polynom for energy-rt-corr fitting
	TFile *InFile[ArrayNumber];
	TH1D* hRiseTimePeak[ArrayNumber][PeakNumber];
	TProfile* hXProfileRiseTimePeak[ArrayNumber][PeakNumber];
	TF1* fProfileFitFuncPol[ArrayNumber][PeakNumber];
	TF1* fProfileFitFuncConstant[ArrayNumber][PeakNumber];
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
  InputListPath +="/";
	InputListPath +=BeamTimeMonth;
	InputListPath +="/txtfiles/";
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
		int incr = 0;
	  while (InputList.good())																//loop over all lines
		{
	  	incr++;
	  	cout << "blaaaaaa" << endl;
			char buffer[200] = "";
			char buf2[10] = "";
			InputList.getline(buffer,200);												//read name of a ROOT file
			cout << "File:\t\t\t\t"<<buffer << endl;
			TString InFileName = buffer; 
			cout << "BufferLength: \t" << InFileName.Length() << endl;
			if (InFileName.Length() == 0)
				continue;
			//extract parameter values from file name
			int M, FilterType,SigmaGaus,SigmaBil, tau, MAL,  runNo;
			int StartFile, EndFile;
			TString Date;


			//TString ComparisonString = COSYTESTANADIR +"/" + BeamTimeMonth + "/COSY_Ana_%*i_run%*i_%*i_%*i___%*i,%*i,%*i,%*i,%*i,%*i/Ana_COSY_Ana_%s_run%i_%i_%i___%i,%i,%i,%i,%i,%i.root";
			TString ComparisonString = COSYTESTANADIR +"/" + BeamTimeMonth + "/CombinedData/Ana_COSY_Ana_%s_run%i_%i_%i___%i,%i,%i,%i,%i,%i.root";

			cout << InFileName.Data() << endl;
			cout << "File:\t\t\t\t"<< ComparisonString << endl;

			//sscanf(InFileName.Data(),ComparisonString.Data(),&Date,&runNo,&StartFile, &EndFile,&M,&MAL,&FilterType,&SigmaGaus,&SigmaBil,&tau);
			//cout <<Date<< "\t"<<runNo<< "\t"<<StartFile<< "\t" <<EndFile<< "\t"<<M<< "\t"<<MAL<< "\t"<<FilterType<< "\t"<<SigmaGaus<< "\t"<<tau << endl;
			if (incr == 1)
				Date="0108";
			if (incr == 2)
				Date="2407";

			runNo = 1;
			StartFile = 1;
			EndFile = 1;
			M = 200;
			MAL = 100;
			FilterType= 0;
			SigmaGaus=3;
			SigmaBil=3;
			tau=6210;

			InFile[i] = new TFile(InFileName);								//open ROOT file
			if (InFile[i]->GetSize() < 600)
			{
				CorruptedList << InFileName.Data() << endl;
				cout << "corrupted file" << endl;
				continue;
			}

			TH1D* hEnergyMA;
			TH2D* hEnergyRt1090;
			if (!SR)
			{
				InFile[i]->GetObject("/Histograms/Energyspectrum/Energy_01;1",hEnergyMA);									// get histrogram
				InFile[i]->GetObject("/Histograms/EnergyRt1090/EnergyRt1090_01;1",hEnergyRt1090);									// get histrogram
			}
			else
			{
				InFile[i]->GetObject("/Histograms/Energyspectrum/EnergyCorr_01;1",hEnergyMA);									// get corrected histrogram
				InFile[i]->GetObject("/Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01;1",hEnergyRt1090);									// get corrected histrogram

			}

			if (UseMovingAv)
				if (!hEnergyMA)
					continue;
	
			TString Path, InfileName,RootFilename, TxtFilename;
			Int_t LengthOfPath;
			
			Path = COSYTESTANADIR;
			Path+="/";
			Path+=BeamTimeMonth;
				Path += "/CombinedData/Fit/";
			char buf[100];
			if (!SR)
				sprintf(buf,"Fitted_COSY_Ana_%s_run%i_%i_%i___%i,%i,%i,%i,%i.root",Date.Data(),runNo,StartFile, EndFile,M,MAL,FilterType,SigmaGaus,SigmaBil,tau);
			else
				sprintf(buf,"Fitted_COSY_Ana_%s_run%i_%i_%i___%i,%i,%i,%i,%iSR.root",Date.Data(),runNo,StartFile, EndFile,M,MAL,FilterType,SigmaGaus,SigmaBil,tau);
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
			if (incr == 1)
				Ana= new THypGeSpectrumAnalyser(hEnergyMA,"jülich2", 35 );
			if (incr == 2)
				Ana= new THypGeSpectrumAnalyser(hEnergyMA,"co60", 35 );
			//Ana= new THypGeSpectrumAnalyser(hEnergyMA,"jülich2", 35 );
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

			double PeakRange[PeakNumber][3];
			double* PeakRangeBuffer;
			double PeakPositionX[PeakNumber][2];
			double* PeakPositionXBuffer;
			TNtupleD *tNtupleD = new TNtupleD("tNtupleD","Tree with vectors","Energy:RangeMin:RangeMax");
			if (incr ==2)
				PeakNumber=2;
			for (int j = 0; j < PeakNumber; j++)
			{
				PeakRangeBuffer= Ana->GetPeakRangeChannels(j);
				PeakRange[j][0] = PeakRangeBuffer[0];
				PeakRange[j][1] = PeakRangeBuffer[1];
				PeakRange[j][2] = PeakRangeBuffer[2];
				tNtupleD->Fill(PeakRange[j]);
				PeakPositionXBuffer = Ana->GetPeakMaximumX(j);
				PeakPositionX[j][0] = PeakPositionXBuffer[0];
				PeakPositionX[j][1] = PeakPositionXBuffer[1];
									//cout << "PeakLow "<< PeakRange[0] << "\t\t" << "PeakHigh" << PeakRange[1] << endl;
				sprintf(buf,"hRiseTimePeak%i_run%i",j,i);
				// extract risetime histogram of only the 4 Peaks from 2d histogram
				hRiseTimePeak[i][j]= hEnergyRt1090->ProjectionX(buf,PeakRange[j][1],PeakRange[j][2],"");

				// do XProfile of Correlation histogram to get correction functions
				sprintf(buf,"hXProfileRiseTimePeak%i_run%i",j,i);
				hXProfileRiseTimePeak[i][j] = hEnergyRt1090->ProfileX(buf,PeakRange[j][1],PeakRange[j][2],"");
				sprintf(buf,"fProfileFitFuncPolPeak%i_run%i",j,i);
				sprintf(buf2,"pol%i",PolOrder);
				fProfileFitFuncPol[i][j] = new TF1(buf, buf2, 80,220);
				hXProfileRiseTimePeak[i][j]->Fit(fProfileFitFuncPol[i][j],"R");
				sprintf(buf,"fProfileFitFuncConstantPeak%i_run%i",j,i);
				fProfileFitFuncConstant[i][j]= new TF1(buf, "pol0", 80,220);
				fProfileFitFuncConstant[i][j]->SetLineColor(kBlue);
				hXProfileRiseTimePeak[i][j]->Fit(fProfileFitFuncConstant[i][j],"R+");

				//normalization of pol fit by division of each parameter by the value of the constant fit
//				for(Int_t iPara =0; iPara <= PolOrder; iPara++)
//				{
//					cout << fProfileFitFuncPol[i][j]->GetParameter(iPara) << "\t";
//					fProfileFitFuncPol[i][j]->SetParameter(iPara,fProfileFitFuncPol[i][j]->GetParameter(iPara)/fProfileFitFuncConstant[i][j]->GetParameter(0));
//					cout << fProfileFitFuncPol[i][j]->GetParameter(iPara) << endl;
//				}


			}	// end of loop over Peaks (j)

			cout << "starting the slice fitting"<< endl;
			Int_t XBinStart =14;
			Int_t XBinEnd = 28;
			Int_t nXBins = XBinEnd - XBinStart;
			Double_t *X= new Double_t[nXBins];
			Double_t *YSimple= new Double_t[nXBins];
			Double_t *Y= new Double_t[nXBins];
			Double_t *Yerror= new Double_t[nXBins];
			Int_t YRangeMin = PeakRange[2][1];
			Int_t YRangeMax = PeakRange[2][2];

			Double_t NormalizationFactor = PeakPositionX[2][2];
			TH1D *hBinScan[nXBins];
			TF1 *fFit[nXBins];
			//cout << YRangeMin <<"\t"<< YRangeMax<<endl;

			// loop over Xbins to create slices and fit them
			for (int iXBin = XBinStart; iXBin < XBinEnd; iXBin ++)
			{
				sprintf(buf,"nBinScan%i",iXBin);

				hBinScan[iXBin] = new TH1D(buf,"Scan", YRangeMax-YRangeMin+1, YRangeMin, YRangeMax);
				for (int iYBin = YRangeMin; iYBin <= YRangeMax; iYBin++)
				{
					hBinScan[iXBin]->SetBinContent(hBinScan[iXBin]->FindBin(iYBin),hEnergyRt1090->GetBinContent(iXBin, iYBin));
					//cout << hEnergyRt1090->GetBinContent(iXBin, iYBin)<< endl;
				}

				sprintf(buf,"fFit%i",iXBin);
				fFit[iXBin] = new TF1(buf,"gaus(0)+pol2(3)", hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()-2), hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()+15));
//				fFit[iXBin]->SetParameter(0, 10);
//				fFit[iXBin]->SetParameter(1, hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()));
//				fFit[iXBin]->SetParLimits(1,fFit[iXBin]->GetParameter(1)-3,fFit[iXBin]->GetParameter(1)+3);
//
				THypGePeakFitFunction *func = new THypGePeakFitFunction();

	//			fFit[iXBin] = new TF1(buf,func,&THypGePeakFitFunction::NewFunction, hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()-5), hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()+5),8,"THypGePeakFitFunction","NewFunction");
				fFit[iXBin]->SetParameter(0, 10);
				fFit[iXBin]->SetParameter(1, hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()));
				fFit[iXBin]->SetParLimits(1,fFit[iXBin]->GetParameter(1)-1,fFit[iXBin]->GetParameter(1)+5);

				cout <<hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin()) << endl;
				fFit[iXBin]->SetParameter(2, 3);
				//hBinScan->Draw();
				hBinScan[iXBin]->Fit(fFit[iXBin],"0");
				X[iXBin-XBinStart]= iXBin*10;
				YSimple[iXBin-XBinStart]= hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin());
				Y[iXBin-XBinStart]= fFit[iXBin]->GetParameter(1);///2 +hBinScan[iXBin]->GetBinCenter(hBinScan[iXBin]->GetMaximumBin())/2;
				Yerror[iXBin-XBinStart]=fFit[iXBin]->GetParError(1);
				cout <<X[iXBin-XBinStart]<<" " << Y[iXBin-XBinStart] << " " << Yerror[iXBin-XBinStart]<< endl;
				//hBinScan[iXBin]->Delete();
				//fFit[iXBin]->Delete();
			}
			for (int iXBin = XBinStart; iXBin < XBinEnd; iXBin ++)
				cout <<X[iXBin-XBinStart]<<" " << Y[iXBin-XBinStart] << " " << Yerror[iXBin-XBinStart]<< endl;

			TGraphErrors *gBinScanSimple = new TGraphErrors(nXBins,X,YSimple,0,0);
			TGraphErrors *gBinScan = new TGraphErrors(nXBins,X,Y,0,Yerror);
			//TGraphErrors *gBinScan = new TGraphErrors(nXBins,X,Y,0,0);
			TCanvas *cCan = new TCanvas("cCan","cCan",800,600);
			cCan->cd();

			//gBinScan->Draw("same");
			TF1 *fFitBinScanConst = new TF1("fFitBinScanConst","pol0", XBinStart*10, XBinEnd*10);
			TF1 *fFitBinScan = new TF1("fFitBinScan","pol3", XBinStart*10, XBinEnd*10);
			//gBinScan->Fit(fFitBinScanConst,"r");

			gBinScan->Fit(fFitBinScan,"r");
			TF1 *fFitBinScanNorm = new TF1("fFitBinScanNorm","pol4", XBinStart*10, XBinEnd*10);
			for(Int_t iPara =0; iPara <= 4; iPara++)
			{
				//fFitBinScanNorm->SetParameter(iPara,fFitBinScan->GetParameter(iPara)/fFitBinScanConst->GetParameter(0));
				fFitBinScanNorm->SetParameter(iPara,fFitBinScan->GetParameter(iPara)/NormalizationFactor);
			}
			gBinScan->Draw();
			gBinScanSimple->SetLineColor(kRed);
			//gBinScanSimple->Draw("same");
			fFitBinScan->Draw("same");


				if (WriteEnabled)
				{
				// write extracted parameters to file
				sprintf(buf,"ParametersFirstAnaStepCOSY_Ana_%s_run%i_%i_%i___%i,%i,%i,%i,%i,%i,MA.root",Date.Data(),runNo,StartFile, EndFile,M,MAL,FilterType,SigmaGaus,SigmaBil,tau);

				cout << "Writing Parameter File" << endl;
				TString OutputParametersFileName = COSYTESTANADIR + "/"+BeamTimeMonth+"/DatabaseFirstAnalysisStep/" +buf;

				OutputParametersFile[i] = new TFile(OutputParametersFileName,"RECREATE");
				for (int iXBin = XBinStart; iXBin < XBinEnd; iXBin ++)
				{
					hBinScan[iXBin]->Write();
					fFit[iXBin]->Write();
				}
				gBinScan->Write();
					hEnergyRt1090->Write();
					fFitBinScan->Write();
					fFitBinScanConst->Write();
					fFitBinScanNorm->Write();

					for (int j = 0; j < PeakNumber; j ++)
					{
						sprintf(buf,"hRiseTimePeak%i",j);
						hRiseTimePeak[i][j]->Write(buf);
						sprintf(buf,"hXProfileRiseTimePeak%i",j);
						hXProfileRiseTimePeak[i][j]->Write(buf);
						sprintf(buf,"fProfileFitFuncPolPeak%i",j);
						fProfileFitFuncPol[i][j]->Write(buf);
						sprintf(buf,"fProfileFitFuncConstantPeak%i",j);
						fProfileFitFuncConstant[i][j]->Write(buf);

					}
					tNtupleD->Write();
				OutputParametersFile[i]->Close();
				cout << "Parameter File finished!" << endl;
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
			cout << SigmaBil << endl;
			if (SigmaBil != 3)	// bil
			{
				cout << "resultmap" << i<< endl;
				ResultMapBil[(M /20 -3 )*10 + ((SigmaBil -100)/200 )-1][SigmaGaus] = Ana->GetFWHMCo();
			}
			if (SigmaBil == 3)	// gaus
			{
				ResultMapGaus[M /20 -3][SigmaGaus] = Ana->GetFWHMCo();
			}
			cout << "blaaaarrrr" << endl;
			if (Ana)
			delete Ana;

			delete hEnergyMA;
			cout << "blaaaa" << endl;
			InFile[i]->Close();
			if (InFile[i])
			{
				cout << "InFile" << i <<endl;
				InFile[i]->Delete();
				delete InFile[i];
			}
			i++;
			cout << "blaaaa" << endl;
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
//	for (int ii=0; ii < i; ii++)
//	{
//		hRiseTimePeak[ii][2]->Write();
//	}
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

	
	//TString OutputFile = COSYTESTANADIR;
	//			OutputFile += "/COSY/CombinedData/Fit/Output.txt";
	//ofstream Output (OutputFile.Data(),std::fstream::trunc);
	
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
	//Output.close();
	InputList.close();
	CorruptedList.close();
	App->Run();
	return 0;
}
