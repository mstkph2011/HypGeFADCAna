void HistoMakeup(TH2D *histo,Color_t color, Style_t MarkerStyle,double MarkerSize=1.5,bool useDottedLines=0)
{
		histo->SetMarkerStyle(MarkerStyle);
		histo->SetMarkerColor(color);
		histo->SetMarkerSize(MarkerSize);
		histo->SetLineColor(color);
		histo->SetLineWidth(2);
		//histo->SetFillColor(color);
		if(useDottedLines)
			histo->SetLineStyle(1);
		histo->GetXaxis()->SetTitle("t_{10Imax} [#mus]");
		histo->GetYaxis()->SetTitle("Pulse height [ADC channels]");
		histo->GetXaxis()->CenterTitle(0);
		histo->GetYaxis()->CenterTitle(0);
		histo->GetYaxis()->SetTitleOffset(1.4);
		histo->GetXaxis()->SetRangeUser(-10,350);
		histo->GetYaxis()->SetRangeUser(1550,1900);
		//histo->GetYaxis()->SetTitleOffset(1.35);
		histo->GetXaxis()->SetNdivisions(507);
		histo->GetYaxis()->SetNdivisions(509);
		histo->SetStats(0);
		
		

}

void Pic_2DHistoExamples(
	TString InputName = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/June2014/HistosTreeCOSYJune2014Dataset-1_200,100,0,5339_SR0.root",
	TString InputName2 = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/June2014/HistosTreeCOSYJune2014Dataset14_200,100,0,5339_SR0.root")
{
	andi::setCustomStyle(0,0);
	 gROOT->ForceStyle();
	 
	bool SaveImages=1;
	TString PicRootDir = gSystem->Getenv("PICTUREDIR");
	int CanvasScaler=1;
	if (SaveImages)
	{
		gROOT->SetBatch(kTRUE);
		CanvasScaler=1;
	}
	TH2D* hT10ImaxEnergy0;
	TH2D* hT10ImaxEnergy14;
	
	TFile *Input1=new TFile(InputName.Data());
	hT10ImaxEnergy0=(TH2D*) Input1->Get("hT10DeriMaxEnergy");
	hT10ImaxEnergy0->SetDirectory(0);
	HistoMakeup(hT10ImaxEnergy0,1,1);
	Input1->Close();
	TFile *Input2=new TFile(InputName2.Data());
	hT10ImaxEnergy14=(TH2D*) Input2->Get("hT10DeriMaxEnergy");
	hT10ImaxEnergy14->SetDirectory(0);
	HistoMakeup(hT10ImaxEnergy14,1,1);
	Input2->Close();
	
	TCanvas *Can=new TCanvas("Can","",1400*CanvasScaler,600*CanvasScaler);
	//Can->UseCurrentStyle();
	//gROOT->ForceStyle();
	Can->Divide(2,1);
	Can->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.13);
	hT10ImaxEnergy0->Draw("colz");
	Can->cd(2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.13);
	hT10ImaxEnergy14->Draw("colz");
	
	if(SaveImages)
		andi::saveCanvas_allFileNames(Can,TString::Format("%s/Mar_COSY_Corr2DHistoExamples",PicRootDir.Data()));
	return 0;
}
