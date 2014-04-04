
//Double_t fitf(Double_t *x, Double_t *par)
//{
  ////Double_t fitval = par[0]*x[0] + par[1];
  
  ////Double_t fitval = par[0]+par[1]*(*x)+par[2]*(*x)*(*x);
  //Double_t fitval = par[0]+par[1]*(*x)+par[2]*(*x)*(*x)+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
   //return fitval;
//}

void CalculatePileUpTimeCorr(TString Filename = "/data/work/kpha1/steinen/COSYBeamtestAna/COSYnew/CombinedData/COSY_Ana200,3,9.root")
{
   
		TFile* File = new TFile(Filename.Data());
		// set edges for cuts
		Int_t LowerRtCut = 1;																
		Int_t UpperRtCut = 120;
		Int_t LowerChannelCut = 1810;
		Int_t UpperChannelCut = 1870;
	
	File->ls();
	
	// get histogram with correlation between pile up time and ADC value from file
	TH2D* hEPTC;
	File->GetObject("Histograms/EnergyTimeSinceLastPulse/EnergyTimeSinceLastPulse",hEPTC);
	TCanvas *c123 = new TCanvas("c123");
		hEPTC->Draw("colz");
	TH2D *hEPTCcut = new TH2D("hEPTCcut", "EPTCcut",hEPTC->GetNbinsX(),hEPTC->GetXaxis()->GetXmin() ,hEPTC->GetXaxis()->GetXmax() , hEPTC->GetNbinsY(),hEPTC->GetYaxis()->GetXmin() ,hEPTC->GetYaxis()->GetXmax() ); // make histogram for cut spectrum with same bins (number, min max; x and y dir)
	 
	 //define cut 
   const Int_t n = 5;
   Float_t x[n] = {LowerRtCut,UpperRtCut,UpperRtCut,LowerRtCut,LowerRtCut};
   Float_t y[n] = {LowerChannelCut,LowerChannelCut,UpperChannelCut,UpperChannelCut,LowerChannelCut};
   TCutG *cut = new TCutG("cut",n,x,y);
   
   //loop over all bins and fill cut histo if inside cut
   for (Int_t iX=1;iX<hEPTC->GetNbinsX();iX++) 
   {
		for (Int_t iY=1;iY<hEPTC->GetNbinsY();iY++) 
		{
			if (cut->IsInside(hEPTC->GetXaxis()->GetBinCenter(iX),hEPTC->GetYaxis()->GetBinCenter(iY))) 			// GetBinCenter to be able to use bin widths != 1
			{
				hEPTCcut->Fill(hEPTC->GetXaxis()->GetBinCenter(iX),hEPTC->GetYaxis()->GetBinCenter(iY), hEPTC->GetBinContent(iX,iY));
			}
		}
   }
   
	TCanvas *c1 = new TCanvas("c1","show profile",600,900);
	c1->Divide(3,2);
	c1->cd(1);
		hEPTC->Draw("colz");
		cut->Draw("same");
	c1->cd(2);
		hEPTCcut->Draw("colz");
   //use a TProfile to convert the 2-d to 1-d problem
	TProfile *prof = hEPTC->ProfileX();
	TProfile *profcut = hEPTCcut->ProfileX();
   //prof->Fit("pol1");
 
		// fit of both histograms, tuning of fit (function, range, ... necessary)
		Int_t NoOfPar = 5;
		TF1 *funcProf = new TF1("funcProf","pol4",hEPTC->GetXaxis()->GetXmin(),hEPTC->GetXaxis()->GetXmax());
		TF1 *funcProfcut = new TF1("funcProfcut","[0]*(1-[1]/pow(x,[2]))",LowerRtCut,UpperRtCut);
		//TF1 *funcProfcut = new TF1("funcProfcut","[0]*(1-TMath::Exp(-[1]*pow(x,[2])))",LowerRtCut,UpperRtCut);
		funcProfcut->SetParameters((LowerChannelCut+UpperChannelCut)/2, 1);
		//funcProfcut->SetParameter(3,0.1);
		//funcProfcut->SetParLimits(3,0,1);
		
	
	cout << "test" << endl;
	c1->cd(4);
		prof->Fit("funcProf","rB");
	c1->cd(5);
		profcut->Fit("funcProfcut","rB");
	cout << "test" << endl;
	c1->cd(3);
			TH2D *hEPTCcutCorr = new TH2D("hEPTCcutCorr", "EPTCcutCorr",hEPTC->GetNbinsX(),hEPTC->GetXaxis()->GetXmin() ,hEPTC->GetXaxis()->GetXmax() , hEPTC->GetNbinsY(),hEPTC->GetYaxis()->GetXmin() ,hEPTC->GetYaxis()->GetXmax() ); // make histogram for cut spectrum with same bins (number, min max; x and y dir)

		//Double_t norm = funcProfcut->GetParameter(0);
		////Double_t norm = funcProf->GetParameter(0);
		
		//for (Int_t i = 0; i < 2; i++)
		//{
			//funcProfcut->SetParameter(i,  funcProfcut->GetParameter(i)/norm);
			////funcProf->SetParameter(i,  funcProf->GetParameter(i)/norm);
			//cout << funcProfcut->GetParameter(i) << endl;
		//}
		funcProfcut->SetParameter(0,1);
		for (Int_t iX=1;iX<hEPTC->GetNbinsX();iX++) 
   {
			for (Int_t iY=1;iY<hEPTC->GetNbinsY();iY++) 
			{
				//hEPTCcutCorr->Fill(hEPTC->GetXaxis()->GetBinCenter(iX),hEPTC->GetYaxis()->GetBinCenter(iY)/funcProfcut->Eval(hEPTC->GetXaxis()->GetBinCenter(iX)), hEPTC->GetBinContent(iX,iY));
				hEPTCcutCorr->Fill(hEPTCcut->GetXaxis()->GetBinCenter(iX) , hEPTCcut->GetYaxis()->GetBinCenter(iY)/funcProfcut->Eval(hEPTCcut->GetXaxis()->GetBinCenter(iX)), hEPTCcut->GetBinContent(iX,iY));
				//hEPTCcutCorr->Fill(hEPTC->GetXaxis()->GetBinCenter(iX),hEPTC->GetYaxis()->GetBinCenter(iY)/funcProf->Eval(hEPTC->GetXaxis()->GetBinCenter(iX)), hEPTC->GetBinContent(iX,iY));
			}
		}
		hEPTCcutCorr->Draw("colz");
		
	c1->cd(6);
	//TProfile *profcutECorr =hEPTCcutCorr->ProjectionY();
	//TProfile *profcutE =hEPTCcut->ProjectionY();
	
	//profcutECorr->Fit("gausn");
	//profcutE->Draw("SAME");
	//profcutE->SetLineColor(kBlue);
	
	TProfile *profcutCorr = hEPTCcutCorr->ProfileX();
	profcutCorr->Draw();
}
   
