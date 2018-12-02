double PolPiecewise(double *x, double *par)
{
	double value;
	if (x[0] < par[0])
		value=par[1];
	else
		value = par[1]+par[2]*(x[0]-par[0])+par[3]*(x[0]-par[0])*(x[0]-par[0])+par[4]*(x[0]-par[0])*(x[0]-par[0])*(x[0]-par[0])+par[5]*(x[0]-par[0])*(x[0]-par[0])*(x[0]-par[0])*(x[0]-par[0]);
	return value;
}


void DoTDeriMax1090Correction(TString SpectrumFileInput = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/June2014/COSYJune2014Dataset11_200,100,0,5339_SR1.root", TString FitFileInput = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJune2014Dataset11_200,100,0,5339_SR1.root")
{
	TH2D *hSpectrumTDeriMax1090_EnergyChannel;
	TH2D *hSpectrumTDeriMax1090Rel_EnergyChannel;
	
	TH2D *hSpectrumT1090_EnergyChannelCorr1;
	TNtuple *DataNTuple;
	TFile *SpectrumInput = new TFile(SpectrumFileInput.Data());
	
	hSpectrumTDeriMax1090_EnergyChannel = (TH2D*) SpectrumInput->Get("Histograms/Energy_DeriMaxT90/Energy_DeriMaxT90_01");
	hSpectrumTDeriMax1090_EnergyChannel->SetDirectory(0);
	hSpectrumTDeriMax1090Rel_EnergyChannel = (TH2D*) SpectrumInput->Get("Histograms/Energy_DeriMaxT90Rel/Energy_DeriMaxT90Rel_01");
	hSpectrumTDeriMax1090Rel_EnergyChannel->SetDirectory(0);
	hSpectrumT1090_EnergyChannelCorr1 = (TH2D*) SpectrumInput->Get("Histograms/EnergyRt1090/EnergyRt1090CorrectionRt_01");
	hSpectrumT1090_EnergyChannelCorr1->SetDirectory(0);
	SpectrumInput->Close();
	//hSpectrumTDeriMax1090_EnergyChannel->Draw("colz");
	
	
	TFile *FitInput = new TFile(FitFileInput.Data(),"Update");
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
	
	TCanvas* can=new TCanvas();
	gPad->SetLogz();
	hSpectrumTDeriMax1090_EnergyChannel->GetXaxis()->SetRangeUser(-30,250);
	hSpectrumTDeriMax1090Rel_EnergyChannel->GetXaxis()->SetRangeUser(-0.1,0.95);
	
	for (Int_t ki=0;ki<entries;ki++)
	{
		DataNTuple->GetEntry(ki);
		//if (int(Energy) == 1332)
		//if (int(Energy) == 510)
		//if (ki == entries-1)
		{
			cout << ChannelRangeMin << " " << ChannelRangeMax << endl;
			
			
			//first correction via TDeriMaxT90Rel
			hSpectrumTDeriMax1090_EnergyChannel->GetYaxis()->SetRangeUser(ChannelRangeMin,ChannelRangeMax);
			hSpectrumTDeriMax1090Rel_EnergyChannel->GetYaxis()->SetRangeUser(ChannelRangeMin,ChannelRangeMax);
			//TF1* FitFuncSlices = new TF1("FitFuncSlices","gaus(0)+pol0(3)",ChannelRangeMin,ChannelRangeMax);
			//FitFuncSlices->SetParameters(1000,ChannelPeakPos-10,4,10);
			//FitFuncSlices->SetParLimits(1,ChannelRangeMin,ChannelRangeMax);
			//FitFuncSlices->SetParLimits(2,0,10);
			//FitFuncSlices->SetParLimits(3,0,100);
	
			//gDirectory->ls();
			
			TH1D *hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually=new TH1D("hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually","",hSpectrumTDeriMax1090Rel_EnergyChannel->GetNbinsX(),-0.3,1.3);
			//cout <<hSpectrumTDeriMax1090_EnergyChannel_MaxPos->GetEntries()<< endl;
			for(int binX = hSpectrumTDeriMax1090Rel_EnergyChannel->GetXaxis()->FindBin(-0.1);binX <= hSpectrumTDeriMax1090Rel_EnergyChannel->GetXaxis()->FindBin(0.90);binX++)
			{
				cout << "binx " << binX << endl;
				TH1D *hProfileY =hSpectrumTDeriMax1090Rel_EnergyChannel->ProjectionY("_py",binX,binX);	
				double MaxValue=hProfileY->GetBinCenter(hProfileY->GetMaximumBin());
				TF1* FitFuncSlices = new TF1("FitFuncSlices","gaus(0)+[3]",MaxValue-20,MaxValue+20);
				TF1* FitFuncGausSlices = new TF1("FitFuncGausSlices","gaus(0)",MaxValue-20,MaxValue+20);
				FitFuncGausSlices->SetParameters(hProfileY->GetBinContent(hProfileY->GetMaximumBin()),MaxValue,4);
				
				hProfileY->Fit(FitFuncGausSlices,"RNQ");
				
				FitFuncSlices->SetParameters(FitFuncGausSlices->GetParameter(0),FitFuncGausSlices->GetParameter(1),FitFuncGausSlices->GetParameter(2),10);
				
				FitFuncSlices->SetParLimits(0,0,10000);
				FitFuncSlices->SetParLimits(1,MaxValue-10,MaxValue+10);
				FitFuncSlices->SetParLimits(2,0,10);
				FitFuncSlices->SetParLimits(3,0,100);
				hProfileY->Fit(FitFuncSlices,"RNQ");
				
				cout <<MaxValue<<"  " << FitFuncSlices->GetParameter(1) << "   " << FitFuncSlices->GetParError(1) <<endl;
				hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->SetBinContent(binX, FitFuncSlices->GetParameter(1));
				hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->SetBinError(binX, FitFuncSlices->GetParError(1));
				
				
			}
			hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->GetYaxis()->SetRangeUser(ChannelPeakPos-100,ChannelPeakPos+50);
			hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->GetXaxis()->SetRangeUser(-0.05,0.9);
			TF1 * funcCorrConst = new TF1("funcCorrConst","pol0",-0.03,0.03);
			funcCorrConst->SetLineColor(kRed);
			TF1 * funcCorrPol = new TF1("funcCorrPol",PolPiecewise,-0.05,0.9,6);
			
			funcCorrPol->SetLineColor(kBlue);
			funcCorrPol->SetParameter(0,0.04);
			funcCorrPol->SetParameter(1,ChannelPeakPos);
			hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->Fit(funcCorrConst,"R");
			hSpectrumTDeriMax1090Rel_EnergyChannel_MaxPosManually->Fit(funcCorrPol,"R+");
			TCanvas* can2=new TCanvas();
			gPad->SetLogz();
			
			hSpectrumTDeriMax1090Rel_EnergyChannel->Draw("colz");
			funcCorrPol->DrawCopy("same");
			for(int i = 7;i>0;i--)
				funcCorrPol->SetParameter(i,funcCorrPol->GetParameter(i)/funcCorrPol->GetParameter(1));
			funcCorrPol->Print();
			char buf[20];
			sprintf(buf, "funcCorrNorm_%d",ki);
			funcCorrPol->Write(buf,TObject::kOverwrite);
			
			funcCorrConst->Delete();
			funcCorrPol->Delete();
			
			
			//second correction via T1090 (after first)
			
			hSpectrumT1090_EnergyChannelCorr1->GetYaxis()->SetRangeUser(ChannelRangeMin,ChannelRangeMax);
			
			//TF1* FitFuncSlices = new TF1("FitFuncSlices","gaus(0)+pol0(3)",ChannelRangeMin,ChannelRangeMax);
			//FitFuncSlices->SetParameters(1000,ChannelPeakPos-10,4,10);
			//FitFuncSlices->SetParLimits(1,ChannelRangeMin,ChannelRangeMax);
			//FitFuncSlices->SetParLimits(2,0,10);
			//FitFuncSlices->SetParLimits(3,0,100);
	
			//gDirectory->ls();
			
			TH1D *hSpectrumT1090_EnergyChannelCorr1_MaxPosManually=new TH1D("hSpectrumT1090_EnergyChannelCorr1_MaxPosManually","",hSpectrumT1090_EnergyChannelCorr1->GetNbinsX(),0,1000);
			//cout <<hSpectrumTDeriMax1090_EnergyChannel_MaxPos->GetEntries()<< endl;
			for(int binX = hSpectrumT1090_EnergyChannelCorr1->GetXaxis()->FindBin(30);binX <= hSpectrumT1090_EnergyChannelCorr1->GetXaxis()->FindBin(500);binX++)
			{
				cout << "binx " << binX << endl;
				TH1D *hProfileY =hSpectrumT1090_EnergyChannelCorr1->ProjectionY("_py",binX,binX);
				//hProfileY->Draw();
				//return 0;
				double MaxValue=hProfileY->GetBinCenter(hProfileY->GetMaximumBin());
				FitFuncSlices = new TF1("FitFuncSlices","gaus(0)+[3]",MaxValue-20,MaxValue+20);
				FitFuncGausSlices = new TF1("FitFuncGausSlices","gaus(0)",MaxValue-20,MaxValue+20);
				FitFuncGausSlices->SetParameters(hProfileY->GetBinContent(hProfileY->GetMaximumBin()),MaxValue,4);
				
				hProfileY->Fit(FitFuncGausSlices,"RNQ");
				
				FitFuncSlices->SetParameters(FitFuncGausSlices->GetParameter(0),FitFuncGausSlices->GetParameter(1),FitFuncGausSlices->GetParameter(2),10);
				
				FitFuncSlices->SetParLimits(0,0,10000);
				FitFuncSlices->SetParLimits(1,MaxValue-10,MaxValue+10);
				FitFuncSlices->SetParLimits(2,0,10);
				FitFuncSlices->SetParLimits(3,0,100);
				hProfileY->Fit(FitFuncSlices,"RNQ");
				
				cout <<MaxValue<<"  " << FitFuncSlices->GetParameter(1) << "   " << FitFuncSlices->GetParError(1) <<endl;
				hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->SetBinContent(binX, FitFuncSlices->GetParameter(1));
				hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->SetBinError(binX, FitFuncSlices->GetParError(1));
				
			}
			hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->GetYaxis()->SetRangeUser(ChannelPeakPos-100,ChannelPeakPos+50);
			hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->GetXaxis()->SetRangeUser(30,500);
			new TCanvas();
			hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->Draw();
			
			TF1 * funcCorrConst = new TF1("funcCorrConst","pol0",50,350);
			funcCorrConst->SetLineColor(kRed);
			TF1 * funcCorrPol = new TF1("funcCorrPol","pol6",50,300);
			
			funcCorrPol->SetLineColor(kBlue);
			funcCorrPol->SetParameter(0,0.04);
			funcCorrPol->SetParameter(1,ChannelPeakPos);
			hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->Fit(funcCorrConst,"R");
			hSpectrumT1090_EnergyChannelCorr1_MaxPosManually->Fit(funcCorrPol,"R+");
			
			TCanvas* can2=new TCanvas();
			gPad->SetLogz();
			
			hSpectrumT1090_EnergyChannelCorr1->Draw("colz");
			funcCorrPol->DrawCopy("same");
			for(int i = 7;i>=0;i--)
				funcCorrPol->SetParameter(i,funcCorrPol->GetParameter(i)/funcCorrConst->GetParameter(0));
			funcCorrPol->Print();
			char buf[20];
			sprintf(buf, "funcCorr2Norm_%d",ki);
			funcCorrPol->Write(buf,TObject::kOverwrite);
			
			funcCorrConst->Delete();
			funcCorrPol->Delete();
			//return 0;
			
		}
	}
	//can->ls();
	gDirectory->Delete("ncCorrPol;1");
	
	FitInput->ls();
	//TCanvas* can3=new TCanvas();
	//hSpectrumTDeriMax1090_EnergyChannel->ProfileX()->Draw();
	//hSpectrumTDeriMax1090Rel_EnergyChannel->ProfileX()->Draw();
	//hSpectrumTDeriMax1090Rel_EnergyChannel->Fit("pol3");
	FitInput->Close();
	
	
}
