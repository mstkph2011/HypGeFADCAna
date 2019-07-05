TF1* GeneralCorrectionFunction(
								double ChannelRangeMin, 
								double ChannelRangeMax,
								double ChannelPeakPos ,
								TH2D *h2DInputForCorrection ,
								int LineIndex,
								TString InputType="T10DeriMaxEnergy",
								TString CorrNumber="1",
								double XRangeMin=0, 
								double XRangeMax=300,
								TString FitFuncCorr="pol2",
								double FitCorrRangeMin=10, 
								double FitCorrRangeMax= 270,
								double TresholdForCorrection=10,
								TString FitFuncSlicesString="gaus(0)+[3]+gaus(4)"
								
								
)
{
	h2DInputForCorrection->GetYaxis()->SetRangeUser(ChannelRangeMin,ChannelRangeMax);
	char buf[60];
	//cout << "histo" << h2DInputForCorrection <<endl;
	sprintf(buf, "hMaxCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	TH1D *hMaxPosManually=new TH1D(buf,"",h2DInputForCorrection->GetNbinsX(),h2DInputForCorrection->GetXaxis()->GetXmin(),h2DInputForCorrection->GetXaxis()->GetXmax());
	sprintf(buf, "hMaxFitCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	TH1D *hMaxPosManuallyFit=new TH1D(buf,"",h2DInputForCorrection->GetNbinsX(),h2DInputForCorrection->GetXaxis()->GetXmin(),h2DInputForCorrection->GetXaxis()->GetXmax());
	sprintf(buf, "hGausSigmaCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	TH1D *hGausSigmaManually=new TH1D(buf,"",h2DInputForCorrection->GetNbinsX(),h2DInputForCorrection->GetXaxis()->GetXmin(),h2DInputForCorrection->GetXaxis()->GetXmax());
	
	//cout << "histo again " << h2DInputForCorrection <<endl;
	for(int binX = h2DInputForCorrection->GetXaxis()->FindBin(XRangeMin);binX <= h2DInputForCorrection->GetXaxis()->FindBin(XRangeMax);binX++)
	//for(int binX = h2DInputForCorrection->GetXaxis()->FindBin(200);binX <= h2DInputForCorrection->GetXaxis()->FindBin(200);binX++)
	{
		
		TH1D *hProfileY =h2DInputForCorrection->ProjectionY("_py",binX,binX);	
		double MaxValue=hProfileY->GetBinCenter(hProfileY->GetMaximumBin());
		
		//hMaxPosManually->SetBinContent(binX, MaxValue);
		//h2DInputForCorrection
		//cout <<hProfileY->GetEntries()<<endl;
		//TF1* FitFuncSlices = new TF1("FitFuncSlices","gaus(0)+[3]",MaxValue-20,MaxValue+20);
		//cout << TMath::Max(MaxValue-20,double(ChannelRangeMin)) << "\t" << TMath::Min(MaxValue+20,double(ChannelRangeMax)) << "\t"<<endl;
		TF1* FitFuncGausSlices = new TF1("FitFuncGausSlices","gaus(0)",TMath::Max(MaxValue-20,double(ChannelRangeMin)),TMath::Min(MaxValue+20,double(ChannelRangeMax)));
		FitFuncGausSlices->SetParameters(hProfileY->GetBinContent(hProfileY->GetMaximumBin()),MaxValue,4);
		
		hProfileY->Fit(FitFuncGausSlices,"RNIQ");
		TF1* FitFuncSlices = new TF1("FitFuncSlices",FitFuncSlicesString.Data(),TMath::Max(MaxValue-3*FitFuncGausSlices->GetParameter(2),double(ChannelRangeMin)),TMath::Min(MaxValue+3*FitFuncGausSlices->GetParameter(2),double(ChannelRangeMax)));
		FitFuncSlices->SetParameters(FitFuncGausSlices->GetParameter(0),FitFuncGausSlices->GetParameter(1),FitFuncGausSlices->GetParameter(2),10,10,FitFuncGausSlices->GetParameter(1)-5,5);
		
		FitFuncSlices->SetParLimits(0,FitFuncGausSlices->GetParameter(0)*0.8,FitFuncGausSlices->GetParameter(0)*1.5);
		FitFuncSlices->SetParLimits(1,TMath::Max(FitFuncGausSlices->GetParameter(1)-10,double(ChannelRangeMin)),TMath::Min(FitFuncGausSlices->GetParameter(1)+10,double(ChannelRangeMax)));
		FitFuncSlices->SetParLimits(2,0,FitFuncGausSlices->GetParameter(2)*2);
		FitFuncSlices->SetParLimits(3,0,500);
		FitFuncSlices->SetParLimits(4,0,FitFuncGausSlices->GetParameter(0)*0.3);
		
		FitFuncSlices->SetParLimits(5,TMath::Max(FitFuncGausSlices->GetParameter(1)-10,double(ChannelRangeMin)),TMath::Min(MaxValue-1,double(ChannelRangeMax)));
		FitFuncSlices->SetParLimits(6,0,10);
		hProfileY->Fit(FitFuncSlices,"RINQ");
		//hProfileY->DrawCopy();
		
		//cout <<MaxValue<<"  " << FitFuncSlices->GetParameter(1) << "   " << FitFuncSlices->GetParError(1) <<endl;
		//cout <<MaxValue<<"  " << FitFuncSlices->GetParameter(1) << "   " << FitFuncSlices->GetMaximumX() <<endl;
		hMaxPosManually->SetBinContent(binX, (FitFuncSlices->GetParameter(1))/ChannelPeakPos);
		hMaxPosManually->SetBinError(binX, FitFuncSlices->GetParError(1)/ChannelPeakPos);
		hGausSigmaManually->SetBinContent(binX, FitFuncSlices->GetParameter(2));
		hGausSigmaManually->SetBinError(binX, FitFuncSlices->GetParError(2));
		if(FitFuncSlices->GetParameter(2)<TresholdForCorrection && FitFuncSlices->GetParError(2)<5)
		{
			hMaxPosManuallyFit->SetBinContent(binX, (FitFuncSlices->GetParameter(1))/ChannelPeakPos);
			hMaxPosManuallyFit->SetBinError(binX, FitFuncSlices->GetParError(1)/ChannelPeakPos);
		}
		//hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->SetBinError(binX, FitFuncSlices->GetParameter(2)/ChannelPeakPos);
		hProfileY->Delete();
		//cin.ignore();
	}
	//write histos to file
		
	//sprintf(buf, "hMaxCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	hMaxPosManually->Write(0,TObject::kOverwrite);
	//sprintf(buf, "hGausSigmaCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	hGausSigmaManually->Write(0,TObject::kOverwrite);
	hMaxPosManuallyFit->Write(0,TObject::kOverwrite);
	sprintf(buf, "funcCorr%s_%sNorm_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	//fit corr function and write it to file
	TF1 *funcCorrNorm = new TF1(buf,FitFuncCorr.Data(),FitCorrRangeMin,FitCorrRangeMax);
	funcCorrNorm->SetParameters(1,0,-0);
	
	funcCorrNorm->SetParLimits(0,0.8,1);
	//funcCorrNorm->SetParLimits(2,-1E5,0);
	//if(LineIndex==2)
	//	hMaxPosManuallyFit->Fit(funcCorrNorm,"RBI");
	//else
		hMaxPosManuallyFit->Fit(funcCorrNorm,"RBQI");
	sprintf(buf, "funcCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	TF1 *funcCorr = new TF1(buf,FitFuncCorr.Data(),FitCorrRangeMin,FitCorrRangeMax);
	for(int i= 0; i<funcCorr->GetNpar();i++)
	{
		funcCorr->SetParameter(i,funcCorrNorm->GetParameter(i)*ChannelPeakPos);
	}
	//sprintf(buf, "funcCorr%s_%sNorm_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	funcCorrNorm->Write(0,TObject::kOverwrite);
	//sprintf(buf, "funcCorr%s_%s_%d",CorrNumber.Data(),InputType.Data(),LineIndex);
	funcCorr->Write(0,TObject::kOverwrite);
	h2DInputForCorrection->GetYaxis()->UnZoom();
	return funcCorr;
}



void DoTreeBasedT10DeriMaxCorrection(TString SpectrumFileInput = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/June2014/HistosTreeCOSYJune2014Dataset11_200,100,0,5339_SR0.root", 
							  TString FitFileInput = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJune2014Dataset11_200,100,0,5339_SR0.root",
							  int PeakNumber=2)
{
	TH2D *hT10DeriMaxEnergy;
	TH2D *hT1090Energy;
	TH2D *hTDeriMax90Energy;

	TH2D	*hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2;
	TH2D	*hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2;

	TH2D	*hTDeriMax90EnergyCorr1_T1090EnergyNorm_2;	
	TH2D	*hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2;	
 

	TH2D	*hT1090EnergyCorr1_TDeriMax90EnergyNorm_2;
	TH2D	*hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2;
 
char buf[200];
	//directories
	//			Corr1_T10DeriMaxEnergyNorm_2;1
	//				hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2;1
	//				hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2;1
  	//			Corr1_T1090EnergyNorm_2;1	
	//				hTDeriMax90EnergyCorr1_T1090EnergyNorm_2;1	
	//				hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2;1	
  	//			Corr1_TDeriMax90EnergyNorm_2;1	
  	//				hT1090EnergyCorr1_TDeriMax90EnergyNorm_2;1
  	//				hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2;1
  	

	
	TFile *SpectrumInput = new TFile(SpectrumFileInput.Data());
		hT10DeriMaxEnergy=(TH2D*) SpectrumInput->Get("hT10DeriMaxEnergy");
		hT10DeriMaxEnergy->SetDirectory(0);
		hT1090Energy=(TH2D*) SpectrumInput->Get("hT1090Energy");
		hT1090Energy->SetDirectory(0);
		hTDeriMax90Energy=(TH2D*) SpectrumInput->Get("hTDeriMax90Energy");
		hTDeriMax90Energy->SetDirectory(0);
		
		sprintf(buf,"Corr1_T10DeriMaxEnergyNorm_%i",PeakNumber);
		gDirectory->cd(buf);
			sprintf(buf,"hT1090EnergyCorr1_T10DeriMaxEnergyNorm_%i;1",PeakNumber);
			hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if (hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2)
				hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2->SetDirectory(0);
			//cout<<hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2<<endl;
			//return 0;
			sprintf(buf,"hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_%i;1",PeakNumber);
			hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if(hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2)
				hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2->SetDirectory(0);
		
		sprintf(buf,"../Corr1_T1090EnergyNorm_%i;1",PeakNumber);
		gDirectory->cd(buf);
			sprintf(buf,"hTDeriMax90EnergyCorr1_T1090EnergyNorm_%i;1",PeakNumber);
			hTDeriMax90EnergyCorr1_T1090EnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if(hTDeriMax90EnergyCorr1_T1090EnergyNorm_2)
				hTDeriMax90EnergyCorr1_T1090EnergyNorm_2->SetDirectory(0);
			sprintf(buf,"hT10DeriMaxEnergyCorr1_T1090EnergyNorm_%i;1",PeakNumber);
			hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if(hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2)
				hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2->SetDirectory(0);
				
		sprintf(buf,"../Corr1_TDeriMax90EnergyNorm_%i;1",PeakNumber);
		gDirectory->cd(buf);
			sprintf(buf,"hT1090EnergyCorr1_TDeriMax90EnergyNorm_%i;1",PeakNumber);
			hT1090EnergyCorr1_TDeriMax90EnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if(hT1090EnergyCorr1_TDeriMax90EnergyNorm_2)
				hT1090EnergyCorr1_TDeriMax90EnergyNorm_2->SetDirectory(0);
			sprintf(buf,"hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_%i;1",PeakNumber);
			hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2=(TH2D*) gDirectory->Get(buf);
			if(hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2)
				hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2->SetDirectory(0);
			
	SpectrumInput->Close();
	
	TFile *FitInput = new TFile(FitFileInput.Data(),"UPDATE");
	DataNTuple = (TNtuple*)FitInput->Get("DataNTuple");
	DataNTuple->Scan();
	
	Int_t entries = (Int_t)DataNTuple->GetEntries();
	cout<<"Number of Entries: "<<entries<< endl;
	const int entriesArrayValue =entries;
	TF1 *FitFunc[entriesArrayValue];
	float Energy=0;
	float ChannelPeakPos=0;
	float ChannelRangeMin=0;
	float ChannelRangeMax=0;
	DataNTuple->SetBranchAddress("Energy",&Energy);
	DataNTuple->SetBranchAddress("ChannelPeakPos",&ChannelPeakPos);
	DataNTuple->SetBranchAddress("ChannelRangeMin",&ChannelRangeMin);
	DataNTuple->SetBranchAddress("ChannelRangeMax",&ChannelRangeMax);
	
	TF1 * funcCorr1T10DeriMax90Norm ;
	double ChannelPeakPos1;
	for (Int_t ki=0;ki<entries;ki++)
	{
		DataNTuple->GetEntry(ki);
		cout << ki << endl;
		//if (int(Energy) == 1332)
		//if (int(Energy) == 510)
		//if (ki == entries-1)
		{
			//cout << ChannelRangeMin << " " << ChannelRangeMax << endl;
			if(hT10DeriMaxEnergy)
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT10DeriMaxEnergy,ki);
			if(hT1090Energy)
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT1090Energy,ki,"T1090Energy","1",50,300,"pol2",80,250,10);
			if(hTDeriMax90Energy)
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hTDeriMax90Energy,ki,"TDeriMax90Energy","1",-20,300,"pol2",0,250,15);
			cout << "first order corrections done"<<endl;
			//second round of corrections
			if (hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2)
			{
				sprintf(buf,"T1090EnergyNorm_%iCorr1_T10DeriMaxEnergy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2,ki,buf,"2",50,300,"pol2",80,250,10);
			}
			if(hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2)
			{
				sprintf(buf,"TDeriMax90EnergyNorm_%iCorr1_T10DeriMaxEnergy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2,ki,buf,"2",-20,300,"pol2",0,250,15);
			}
			if(hTDeriMax90EnergyCorr1_T1090EnergyNorm_2)
			{
				cout << "test" << hTDeriMax90EnergyCorr1_T1090EnergyNorm_2<<endl;
				sprintf(buf,"TDeriMax90EnergyNorm_%iCorr1_T1090Energy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hTDeriMax90EnergyCorr1_T1090EnergyNorm_2,ki,buf,"2",-20,300,"pol2",0,250,15);
			}
			if(hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2)
			{
				sprintf(buf,"T10DeriMaxEnergyNorm_%iCorr1_T1090Energy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2,ki,buf,"2");
			}	
			if(hT1090EnergyCorr1_TDeriMax90EnergyNorm_2)
			{
				sprintf(buf,"T1090EnergyNorm_%iCorr1_TDeriMax90Energy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT1090EnergyCorr1_TDeriMax90EnergyNorm_2,ki,buf,"2",50,300,"pol2",80,250,10);
			}
			if(hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2)
			{
				sprintf(buf,"T10DeriMaxEnergyNorm_%iCorr1_TDeriMax90Energy",PeakNumber);
				GeneralCorrectionFunction(ChannelRangeMin,ChannelRangeMax,ChannelPeakPos,hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2,ki,buf,"2");
			}
		}
	}
	
	
	TCanvas* can=new TCanvas();
	//hT10DeriMaxEnergy->Draw("colz");
	hT10DeriMaxEnergy->Write(0,TObject::kOverwrite);
	hT1090Energy->Draw("colz");
	hT1090Energy->Write(0,TObject::kOverwrite);
	hTDeriMax90Energy->Write(0,TObject::kOverwrite);

	hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
	hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);

	hTDeriMax90EnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
	hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);


	hT1090EnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
	hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
	
	
	
	//funcCorr1T10DeriMax90Norm->SetParameters(funcCorr1T10DeriMax90Norm->GetParameter(0)*ChannelPeakPos1,funcCorr1T10DeriMax90Norm->GetParameter(1)*ChannelPeakPos1,funcCorr1T10DeriMax90Norm->GetParameter(2)*ChannelPeakPos1);
	//funcCorr1T10DeriMax90Norm->Draw("same");
	//gPad->SetLogz();
	//hT10DeriMaxEnergy->GetXaxis()->SetRangeUser(-30,250);
	FitInput->Close();
	
}

