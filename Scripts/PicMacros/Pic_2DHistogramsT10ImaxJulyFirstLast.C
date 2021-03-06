void HistoMakeup(TH1 *histo,Color_t color, Style_t MarkerStyle,double MarkerSize=1.5,bool useDottedLines=0)
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
		//histo->GetXaxis()->SetLabelOffset(0.01);
		//histo->GetYaxis()->SetLabelOffset(0.01);
		
		//histo->GetXaxis()->SetTitleSize(0.07);
		//histo->GetYaxis()->SetTitleSize(0.07);
		
		//histo->GetXaxis()->SetLabelSize(0.06);
		//histo->GetYaxis()->SetLabelSize(0.06);
		//histo->GetZaxis()->SetLabelSize(0.06);
		
		//histo->GetXaxis()->SetTitleOffset(1);
		histo->GetYaxis()->SetTitleOffset(1.4);
		histo->GetXaxis()->SetRangeUser(-10,350);
		histo->GetYaxis()->SetRangeUser(2460,2530);
		//histo->GetYaxis()->SetTitleOffset(1.35);
		histo->GetXaxis()->SetNdivisions(507);
		histo->GetYaxis()->SetNdivisions(506);
		histo->GetZaxis()->SetNdivisions(506);
		histo->SetStats(0);
		
		

}

void Pic_2DHistogramsT10ImaxJulyFirstLast(
	TString InputName = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJune2014Dataset14_200,100,0,5339_SR0.root",
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
		CanvasScaler=2;
	}
	
	
	TString InputNameArray[15];
	TString InputName2Array[15];
	for (int i= 0;i<=9;i+=9) 
	{
		InputNameArray[i]="/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJuly2014Dataset";
		InputNameArray[i]+=i;
		InputNameArray[i]+="_200,100,0,5339_SR0.root";
		InputName2Array[i]="/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/July2014/HistosTreeCOSYJuly2014Dataset";
		InputName2Array[i]+=i;
		InputName2Array[i]+="_200,100,0,5339_SR0.root";
	}
	
	double pd=5.6/14.;
	//TH1D* hMaxCorr1,*hSigmaCorr1;
	TH2F* hT10ImaxEnergy14[15];
	//TH2D* hT10ImaxEnergy14Corr;
	TF1* funcCorr[15];
	TString Inset1[15];
	for (int i= 0;i<=9;i+=9) 
	{
		TFile *Input1=new TFile(InputNameArray[i].Data());
		if(i==0)
			funcCorr[i]=(TF1*) Input1->Get("funcCorr1_T10DeriMaxEnergy_2;1");
		else
			funcCorr[i]=(TF1*) Input1->Get("funcCorr1_T10DeriMaxEnergy_3;1");
		funcCorr[i]->SetLineWidth(3);
		Input1->Close();
		Input1->Delete();
		
		TFile *Input2=new TFile(InputName2Array[i].Data());
		hT10ImaxEnergy14[i]=(TH2F*) Input2->Get("hT10DeriMaxEnergy;1");
		hT10ImaxEnergy14[i]->SetDirectory(0);
		HistoMakeup(hT10ImaxEnergy14[i],1,1);
		Input2->Close();
		Input2->Delete();
		
		Inset1[i].Form( "%.1lf x 10^{9} neutrons/cm^{2}",i*pd);
		cout << Inset1[i].Data()<< endl;
	}
	
	double cw =1400*CanvasScaler;
	double ch =600*CanvasScaler;
	TCanvas *Can=new TCanvas("Can","Can",cw,ch);
	
	Can->SetCanvasSize(cw,ch);
	Can->Divide(2);
	
	TLatex *Text1=new TLatex();
	Text1->SetTextSize(0.07);
	Text1->SetTextColor(0);
	
	for (int i= 0;i<=9;i+=9) 
	{
		if(i==0)
			Can->cd(1);
		else
			Can->cd(2);
		gPad->SetLeftMargin(0.13);
		if(i%3 ==0)
			gPad->SetLeftMargin(0.2);	
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.13);
		hT10ImaxEnergy14[i]->Draw("colz");
		funcCorr[i]->Draw("same");
		if(i==1||i==2 ||i==3)
		{
			Text1->DrawLatex(180,1870,"Detector issues");
		}
		Text1->DrawLatex(50,1810,Inset1[i].Data());
		cout << Text1->GetTextFont()<<endl;
	}
	
	if(SaveImages)
		andi::saveCanvas_allFileNames(Can,TString::Format("%s/Mar_COSY_2DHistogramsT10ImaxJulyFirstLast",PicRootDir.Data()));
		return 0;
		
	TFile *Input1=new TFile(InputName.Data());
	hMaxCorr1=(TH1D*) Input1->Get("hMaxCorr1_T10DeriMaxEnergy_2;1");
	hMaxCorr1->SetDirectory(0);
	hSigmaCorr1=(TH1D*) Input1->Get("hGausSigmaCorr1_T10DeriMaxEnergy_2;1");
	hSigmaCorr1->SetDirectory(0);
	//funcCorr1_T10DeriMaxEnergy2;1
	funcCorr=(TF1*) Input1->Get("funcCorr1_T10DeriMaxEnergy_2;1");
	
	//funcCorr->SetDirectory(0);
	HistoMakeup(hMaxCorr1,1,8,1);
	Input1->Close();
	for (int i=0; i<=hMaxCorr1->GetNbinsX();i++)
	{
		hMaxCorr1->SetBinContent(i,hMaxCorr1->GetBinContent(i)*1859.5579);
		hMaxCorr1->SetBinError(i,hSigmaCorr1->GetBinContent(i));
	}
	TFile *Input2=new TFile(InputName2.Data());
	//Input2->ls();
	gDirectory->cd("Corr1_T10DeriMaxEnergyNorm_2;1");
	//gDirectory->ls();
	hT10ImaxEnergy14Corr=(TH2D*) gDirectory->Get("hT10DeriMaxEnergyCorr1_T10DeriMaxEnergyNorm_2;1");
	hT10ImaxEnergy14Corr->SetDirectory(0);
	HistoMakeup(hT10ImaxEnergy14Corr,1,1);
	hT10ImaxEnergy14=(TH2F*) Input2->Get("hT10DeriMaxEnergy;1");
	hT10ImaxEnergy14->SetDirectory(0);
	HistoMakeup(hT10ImaxEnergy14,1,1);
	Input2->Close();
	
	TCanvas *Can=new TCanvas("Can","",1400*CanvasScaler,1200*CanvasScaler);
	//Can->UseCurrentStyle();
	//gROOT->ForceStyle();
	Can->Divide(2,2);
	Can->cd(1);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.13);
	hT10ImaxEnergy14->Draw("colz");
	Can->cd(2);
	
	gPad->SetLeftMargin(0.12);
	//gPad->SetRightMargin(0.13);
	hMaxCorr1->GetYaxis()->SetTitle("Slice fit result [ADC channels]");
	hMaxCorr1->Draw("E");
	funcCorr->Draw("same");
	Can->cd(3);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.13);
	hT10ImaxEnergy14Corr->GetYaxis()->SetTitle("Corrected pulse height [ADC channels]");
	hT10ImaxEnergy14Corr->Draw("colz");
	
	Can->cd(4);
	gPad->SetLeftMargin(0.12);
	//gPad->SetRightMargin(0.13);
	TH1D *hPeakShape =(TH1D*) hT10ImaxEnergy14 ->ProjectionY();
	//hPeakShape->SetLineColor(kRed);
	TH1D *hPeakShapeCorr =(TH1D*) hT10ImaxEnergy14Corr ->ProjectionY();
	hPeakShapeCorr->SetLineColor(kRed);
	hPeakShapeCorr->SetStats(0);
	hPeakShapeCorr->GetXaxis()->SetTitle("Pulse height [ADC channels]");
	hPeakShapeCorr->GetXaxis()->SetTitleOffset(1);
	hPeakShapeCorr->GetYaxis()->SetTitle("Counts [a.u.]");
	hPeakShapeCorr->GetYaxis()->SetTitleOffset(1.4);
	hPeakShapeCorr->Draw();
	hPeakShape->Draw("same");
	hPeakShapeCorr->Draw("same");
	
	TLegend *leg=new TLegend(0.2,0.4,0.4,0.57);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.05);
	leg->SetHeader("Peak shape:");
	leg->AddEntry(hPeakShape,"uncorrected","l");
	leg->AddEntry(hPeakShapeCorr,"corrected","l");
	
	leg->Draw();
	
	if(SaveImages)
		andi::saveCanvas_allFileNames(Can,TString::Format("%s/Mar_COSY_CorrectionExamplesJuneDataset14",PicRootDir.Data()));
	return 0;
}
