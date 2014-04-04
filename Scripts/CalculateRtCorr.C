
//Double_t fitf(Double_t *x, Double_t *par)
//{
  ////Double_t fitval = par[0]*x[0] + par[1];
  
  ////Double_t fitval = par[0]+par[1]*(*x)+par[2]*(*x)*(*x);
  //Double_t fitval = par[0]+par[1]*(*x)+par[2]*(*x)*(*x)+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
   //return fitval;
//}

void CalculateRtCorr()
{
   
	Bool_t firstFile = 1;
	if (firstFile)
	{
		TFile* File = new TFile("ERCorrTest.root");
		// set edges for cuts
		Int_t LowerRtCut = 1;																
		Int_t UpperRtCut = 300;
		Int_t LowerChannelCut = 1900;
		Int_t UpperChannelCut = 1925;
	}
	else
	{
		TFile* File = new TFile("ERCorrTest2.root");
		// set edges for cuts
		Int_t LowerRtCut = 1;
		Int_t UpperRtCut = 300;
		Int_t LowerChannelCut = 1835;
		Int_t UpperChannelCut = 1860;
   }
		//File->ls();
	TCanvas* Can = (TCanvas*) File->Get("Panel3;1");					// get Canvas from File
	//Can->ls();
	TH2D* hERC = (TH2D*) Can->FindObject("EnergyRise1090Corr");		// get histogram from Canvas
	TH2D *hERCcut = new TH2D("hERCcut", "ERCcut",hERC->GetNbinsX(),hERC->GetXaxis()->GetXmin() ,hERC->GetXaxis()->GetXmax() , hERC->GetNbinsY(),hERC->GetYaxis()->GetXmin() ,hERC->GetYaxis()->GetXmax() ); // make histogram for cut spectrum with same bins (number, min max; x and y dir)
	 
	 //define cut 
   const Int_t n = 5;
   Float_t x[n] = {LowerRtCut,UpperRtCut,UpperRtCut,LowerRtCut,LowerRtCut};
   Float_t y[n] = {LowerChannelCut,LowerChannelCut,UpperChannelCut,UpperChannelCut,LowerChannelCut};
   TCutG *cut = new TCutG("cut",n,x,y);
   
   //loop over all bins and fill cut histo if inside cut
   for (Int_t iX=1;iX<hERC->GetNbinsX();iX++) 
   {
		for (Int_t iY=1;iY<hERC->GetNbinsY();iY++) 
		{
			if (cut->IsInside(hERC->GetXaxis()->GetBinCenter(iX),hERC->GetYaxis()->GetBinCenter(iY))) 			// GetBinCenter to be able to use bin widths != 1
			{
				hERCcut->Fill(hERC->GetXaxis()->GetBinCenter(iX),hERC->GetYaxis()->GetBinCenter(iY), hERC->GetBinContent(iX,iY));
			}
		}
   }
   
	TCanvas *c1 = new TCanvas("c1","show profile",600,900);
	c1->Divide(3,2);
	c1->cd(1);
		hERC->Draw("colz");
	cut->Draw("");
	c1->cd(2);
		hERCcut->Draw("colz");
   //use a TProfile to convert the 2-d to 1-d problem
	TProfile *prof = hERC->ProfileX();
	TProfile *profcut = hERCcut->ProfileX();
   //prof->Fit("pol1");
 
		// fit of both histograms, tuning of fit (function, range, ... necessary)
		Int_t NoOfPar = 5;
		TF1 *funcProf = new TF1("funcProf","pol4",hERC->GetXaxis()->GetXmin(),hERC->GetXaxis()->GetXmax(),NoOfPar);
		TF1 *funcProfcut = new TF1("funcProfcut","pol4",LowerRtCut,UpperRtCut,NoOfPar);
	c1->cd(4);
		prof->Fit("funcProf","rB");
	c1->cd(5);
		profcut->Fit("funcProfcut","rB");
		
	c1->cd(3);
			TH2D *hERCcutCorr = new TH2D("hERCcutCorr", "ERCcutCorr",hERC->GetNbinsX(),hERC->GetXaxis()->GetXmin() ,hERC->GetXaxis()->GetXmax() , hERC->GetNbinsY(),hERC->GetYaxis()->GetXmin() ,hERC->GetYaxis()->GetXmax() ); // make histogram for cut spectrum with same bins (number, min max; x and y dir)

		Double_t norm = funcProfcut->GetParameter(0);
		//Double_t norm = funcProf->GetParameter(0);
		for (Int_t i = 0; i < NoOfPar; i++)
		{
			funcProfcut->SetParameter(i,  funcProfcut->GetParameter(i)/norm);
			//funcProf->SetParameter(i,  funcProf->GetParameter(i)/norm);
			cout << funcProfcut->GetParameter(i) << endl;
		}
			
		for (Int_t iX=1;iX<hERC->GetNbinsX();iX++) 
   {
			for (Int_t iY=1;iY<hERC->GetNbinsY();iY++) 
			{
				hERCcutCorr->Fill(hERC->GetXaxis()->GetBinCenter(iX),hERC->GetYaxis()->GetBinCenter(iY)/funcProfcut->Eval(hERC->GetXaxis()->GetBinCenter(iX)), hERC->GetBinContent(iX,iY));
				//hERCcutCorr->Fill(hERC->GetXaxis()->GetBinCenter(iX),hERC->GetYaxis()->GetBinCenter(iY)/funcProf->Eval(hERC->GetXaxis()->GetBinCenter(iX)), hERC->GetBinContent(iX,iY));
			}
		}
		hERCcutCorr->Draw("colz");
		
	c1->cd(6);
	TProfile *profcutECorr =hERCcutCorr->ProjectionY();
	TProfile *profcutE =hERCcut->ProjectionY();
	
	profcutECorr->Fit("gausn");
	profcutE->Draw("SAME");
	profcutE->SetLineColor(kBlue);
	
}
   
