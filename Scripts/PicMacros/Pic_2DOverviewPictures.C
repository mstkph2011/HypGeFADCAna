void CleanHisto(TH1D *histo,int NFunc=4)
{
	TPolyMarker * pm =(TPolyMarker*)histo->GetListOfFunctions()->FindObject("TPolyMarker");

	if (pm) 
	{
		histo->GetListOfFunctions()->Remove(pm);
		pm->Delete();
	}
	char buf[20];
	for(int i= 1;i<=NFunc;i++)
	{
		sprintf(buf,"FitFunc_%i",i);
		TF1 * func =(TF1*)histo->GetListOfFunctions()->FindObject(buf);
		if (func) 
		{
			histo->GetListOfFunctions()->Remove(func);
			func->Delete();
		}

	}
}

void Pic_2DOverviewPictures(
	TString InputName = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJune2014Dataset14_200,100,0,5339_SR0.root")
	//TString InputName = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJuly2014Dataset9_200,100,0,5339_SR0.root")
{
	andi::setCustomStyle(0,0);
	 gROOT->ForceStyle();
	 
	bool SaveImages=0;
	
	
	TH1D *hEnergy, *hEnergyCorr2TDM90_2Corr1T10DM_2, *hEnergyCorr2TDM90_2Corr1T10DM_2_cut;
	
	
	TFile *Input = new TFile(InputName.Data());
	hEnergy=(TH1D*) Input->Get("Energy_01;1");
		hEnergy->SetDirectory(0);
	gDirectory->cd("Corr1_T10DeriMaxEnergyNorm_2;1/Corr2_TDeriMax90EnergyNorm_2;1");
	//cout<<gDirectory->GetPath()<<endl;
	gDirectory->ls();
	hEnergyCorr2TDM90_2Corr1T10DM_2=(TH1D*) gDirectory->Get("hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2;1");
		hEnergyCorr2TDM90_2Corr1T10DM_2->SetDirectory(0);
	gDirectory->cd("../Corr2_TDeriMax90EnergyNorm_2cut;1");
	gDirectory->ls();
	hEnergyCorr2TDM90_2Corr1T10DM_2_cut=(TH1D*) gDirectory->Get("hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2_py;1");
		hEnergyCorr2TDM90_2Corr1T10DM_2_cut->SetDirectory(0);
	Input->Close();
	hEnergyCorr2TDM90_2Corr1T10DM_2->Draw();
	
	hEnergy->SetLineColor(kRed);
	hEnergyCorr2TDM90_2Corr1T10DM_2_cut->SetLineColor(kGreen+2);
	hEnergyCorr2TDM90_2Corr1T10DM_2_cut->Scale(0.5);
	 CleanHisto(hEnergy);
	 CleanHisto(hEnergyCorr2TDM90_2Corr1T10DM_2);
	 CleanHisto(hEnergyCorr2TDM90_2Corr1T10DM_2_cut);
    gPad->Update();
    hEnergyCorr2TDM90_2Corr1T10DM_2->GetXaxis()->SetRangeUser(1800,1900);
    
    TCanvas *Can=new TCanvas("Can","Can",1200,900);
    hEnergyCorr2TDM90_2Corr1T10DM_2->SetStats(0);
	hEnergyCorr2TDM90_2Corr1T10DM_2->Draw();
	hEnergyCorr2TDM90_2Corr1T10DM_2_cut->Draw("same");
	hEnergy->Draw("same");
	
	TLegend *leg = new TLegend(0.15,0.9,0.55,0.7);
	leg->SetLineColor(0);
	leg->SetFillStyle(2);
	leg->SetHeader("93 days #bar{P}anda");
	leg->AddEntry(hEnergy,"uncorrected","l");
	leg->AddEntry(hEnergyCorr2TDM90_2Corr1T10DM_2,"C1 T_{10DM} C2 T_{DM90}","l");
	leg->AddEntry(hEnergyCorr2TDM90_2Corr1T10DM_2_cut,"C1 T_{10DM} C2 T_{DM90}, T_{10DM} > 50 ns","l");
	leg->Draw();
	
	if (SaveImages) andi::saveCanvas_allFileNames(Can, "~/pictures/JuelichResolution/Mar_COSY_PeakShapeJuneEnd" );
	
	
}
