void HistoMakeup(TGraph *histo[6],Color_t color, Style_t MarkerStyle,double MarkerSize=2.5,bool useDottedLines=1)
{
	for(int i =0;i<6;i++)
	{
		histo[i]->SetMarkerStyle(MarkerStyle);
		histo[i]->SetMarkerColor(color);
		histo[i]->SetMarkerSize(MarkerSize);
		histo[i]->SetLineColor(color);
		histo[i]->SetFillColor(color);
		if(useDottedLines)
			histo[i]->SetLineStyle(2);
	}
}
void HistoMakeupSingle(TGraph *histo,Color_t color, Style_t MarkerStyle,double MarkerSize=2.5,bool useDottedLines=1)
{
	
		histo->SetMarkerStyle(MarkerStyle);
		histo->SetMarkerColor(color);
		histo->SetMarkerSize(MarkerSize);
		histo->SetLineColor(color);
		histo->SetFillColor(color);
		if(useDottedLines)
			histo->SetLineStyle(2);
	
}

void FillGraphs(TGraph *Graphs[6],TString Filename,const int nPoints=15,bool inNumberOfNeutrons=1)
{
	//TGraph *Graphs[4];
	//**Graphs=new TGraph*[4];
	double PandaDays[nPoints];
	double FWHMCo[nPoints],FWTMCo[nPoints],FWTMHMCo[nPoints],FWHM511[nPoints],FWTM511[nPoints],FWTMHM511[nPoints];
	double pd,fhCo,ftCo,fh511, ft511;
	int dataset;
	char buf[60];
	
	ifstream file(Filename.Data());
	if(!file.is_open())
	{
		cout << "Error opening File: " << Filename.Data() << endl;
		return -1;
	}
	for(int i =0; i<nPoints;i++)
	{
		file.getline(buf,60);
		sscanf(buf,"%d\t%lf\t%lf\t%lf\t%lf\t%lf",&dataset,&pd,&fhCo,&ftCo,&fh511,&ft511);
		//cout << buf << endl;
		//cout << "\t" << pd << " " << fhCo << " " << ftCo << " " << fh511 << " " << ft511 << endl;
				
		if(inNumberOfNeutrons)
			pd=pd*5.6/93.8;
		//cout << pd << endl;
		PandaDays[i]=pd;
		FWHMCo[i]=fhCo;
		FWTMCo[i]=ftCo;
		FWTMHMCo[i]=ftCo/fhCo;
		FWHM511[i]=fh511;
		FWTM511[i]=ft511;
		FWTMHM511[i]=ftCo/fhCo;
	}
	
	Graphs[0]=new TGraph(nPoints,PandaDays, FWHMCo);
	Graphs[0]->SetTitle("FWHMCo;Panda days; FWHM [keV]");
	Graphs[1]=new TGraph(nPoints,PandaDays, FWTMCo);
	Graphs[1]->SetTitle("FWTMCo;Panda days; FWTM [keV]");
	Graphs[2]=new TGraph(nPoints,PandaDays, FWTMHMCo);
	Graphs[2]->SetTitle("FWTM/FWHMCo;Panda days; FWTM/FWHM ");
	Graphs[3]=new TGraph(nPoints,PandaDays, FWHM511);
	Graphs[3]->SetTitle("FWHM511;Panda days; FWHM [kev]");
	Graphs[4]=new TGraph(nPoints,PandaDays, FWTM511);
	Graphs[4]->SetTitle("FWTM511;Panda days; FWTM [kev]");
	Graphs[5]=new TGraph(nPoints,PandaDays, FWTMHM511);
	Graphs[5]->SetTitle("FWTM/FWHM511;Panda days; FWTM/FWHM ");
	
	if(inNumberOfNeutrons)
	{
		Graphs[0]->SetTitle("FWHMCo;Number of neutrons #[]{#frac{10^{9}}{cm^{2}}}; FWHM [keV]");
		Graphs[1]->SetTitle("FWTMCo;Number of neutrons  #[]{#frac{10^{9}}{cm^{2}}}; FWTM [keV]");
		Graphs[2]->SetTitle("FWTM/FWHMCo;Number of neutrons  #[]{#frac{10^{9}}{cm^{2}}}; FWTM/FWHM");
		Graphs[3]->SetTitle("FWHM511;Number of neutrons  #[]{#frac{10^{9}}{cm^{2}}}; FWHM [kev]");
		Graphs[4]->SetTitle("FWTM511;Number of neutrons  #[]{#frac{10^{9}}{cm^{2}}}; FWTM [kev]");
		Graphs[5]->SetTitle("FWTM/FWHM511;Number of neutrons  #[]{#frac{10^{9}}{cm^{2}}}; FWTM/FWHM");
	}
	
	for(int i=0;i<6;i++)
		Graphs[i]->Sort();
	cout << "asdasd" <<endl;
	//return Graphs;
}



void Pic_Mar_COSY_IrradiationSummaryPicure()
{
	andi::setCustomStyle(0,0);
	 gROOT->ForceStyle();
	 
	bool SaveImages=0;
	TString PicRootDir = gSystem->Getenv("PICTUREDIR");
	int CanvasScaler=1;
	if (SaveImages)
	{
		gROOT->SetBatch(kTRUE);
		CanvasScaler=1;
	}
	
	const int nFilesJune=12;
	const int nFilesJuly=10;
	//June2014
	TGraph *gUncorrectedJune[6];
		FillGraphs(gUncorrectedJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014NoIssued.txt",nFilesJune);
		HistoMakeup( gUncorrectedJune,1, 24);
	//first corrections
	//~ TGraph *gCorr1T1090June[6];
		//~ FillGraphs(gCorr1T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2.txt",nFilesJune);
		//~ HistoMakeup( gCorr1T1090June,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJune[6];
		FillGraphs(gCorr1T10DeriMaxJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2NoIssued.txt",nFilesJune);
		HistoMakeup( gCorr1T10DeriMaxJune,1, 20);
	TGraph *gCorr1T10DeriMaxJuneSmall[6];
		FillGraphs(gCorr1T10DeriMaxJuneSmall, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2NoIssued.txt",nFilesJune);
		HistoMakeup( gCorr1T10DeriMaxJuneSmall,0, 20);
		gCorr1T10DeriMaxJuneSmall[0]->SetMarkerSize(gCorr1T10DeriMaxJuneSmall[0]->GetMarkerSize()*0.5);
	//TGraph *gCorr1TDeriMax90June[6];
		//FillGraphs(gCorr1TDeriMax90June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2.txt",nFilesJune);
		//HistoMakeup( gCorr1TDeriMax90June,2, 23);
	
	
	
	//June2014 Lab Results (HEX 158)
	
	double JuneLabEnergies[2] = {2.25,3.52};
	double JuneLabEnergiesCorr[2] = {2.19,2.67};
	double JuneLabX[2] = {0,5.6};
	
	double JuneLabEnergiesCorrShort[1] = {2.67};
	double JuneLabXShort[1] = {5.6}; 
	
	TGraph *gJuneLabEnergies=new TGraph(2,JuneLabX,JuneLabEnergies);
	TGraph *gJuneLabEnergiesCorr=new TGraph(2,JuneLabX,JuneLabEnergiesCorr);
	TGraph *gJuneLabEnergiesCorrSmall=new TGraph(2,JuneLabX,JuneLabEnergiesCorr);
	//TGraph *gJuneLabEnergiesCorr=new TGraph(1,JuneLabXShort,JuneLabEnergiesCorrShort);
	//TGraph *gJuneLabEnergiesCorrSmall=new TGraph(1,JuneLabXShort,JuneLabEnergiesCorrShort);
		
		HistoMakeupSingle( gJuneLabEnergies,4, 25);
		HistoMakeupSingle( gJuneLabEnergiesCorr,4, 21);
		HistoMakeupSingle( gJuneLabEnergiesCorrSmall,0, 21);
		gJuneLabEnergiesCorrSmall->SetMarkerSize(gJuneLabEnergiesCorrSmall->GetMarkerSize()*0.5);
		
	double JuneLabAnnealingEnergies[2]={2.6,2.42};
	//double JuneLabAnnealingX[2]={1.5,2.5};
	double JuneLabAnnealingX[2]={1.5,2.5};
	TGraph *gJuneLabEnergiesAnnealing=new TGraph(2,JuneLabAnnealingX,JuneLabAnnealingEnergies);
		HistoMakeupSingle( gJuneLabEnergiesAnnealing,4, 25);
	
	//July2014 Lab Results (HEX 146)
	
	double JulyLabEnergies[2] = {1.96,3.44};
	double JulyLabX[2] = {0,3.5};
	
	TGraph *gJulyLabEnergies=new TGraph(2,JulyLabX,JulyLabEnergies);
		HistoMakeupSingle( gJulyLabEnergies,2, 26);
	//gJulyLabEnergies->SetMarkerStyle(32);
	double JulyLabAnnealingEnergies[1]={2.09};
	//double JulyLabAnnealingX[1]={0.5};
	double JulyLabAnnealingX[1]={0.5};
	TGraph *gJulyLabEnergiesAnnealing=new TGraph(1,JulyLabAnnealingX,JulyLabAnnealingEnergies);
		HistoMakeupSingle( gJulyLabEnergiesAnnealing,2, 26);
	
	//return ;
	//TCanvas *Can=new TCanvas("Can","Can",1600*CanvasScaler,500*CanvasScaler);
	double cw =1100*CanvasScaler;
	double ch =800*CanvasScaler;
	TCanvas *Can=new TCanvas("Can","Can",cw,ch);
	//Can->SetCanvasSize(cw,ch);
	Can->SetWindowSize(cw + (cw - Can->GetWw()), ch + (ch - Can->GetWh()));
	
	//Can->SetWindowSize(1066, 1000);
	
	//Can->Divide(3);
	double SplitHorizontal=0.75;
	double SplitVertical=0.12;
	
	TPad *Pads[3];
	Pads[0]=new TPad("TopLeftPad","TopLeftPad",		0,					SplitVertical,	SplitHorizontal,	1);
	Pads[1]=new TPad("TopRightPad","TopRightPad",	SplitHorizontal,	0,	1,					1);
	Pads[2]=new TPad("BottomPad","BottomPad",		0.0,				0,				SplitHorizontal,	SplitVertical);
	//Pads[2]->SetFillColor(kYellow);
	TBox *ProblemArea = new TBox();
	ProblemArea->SetFillColorAlpha(kBlack,.20);
	ProblemArea->SetLineColorAlpha(0,1);
	
	gUncorrectedJune[0]->GetYaxis()->SetRangeUser(0,10);
	gUncorrectedJune[1]->GetYaxis()->SetRangeUser(0,35);
	gUncorrectedJune[2]->GetYaxis()->SetRangeUser(0,4);
	
	for(int i =0;i<3;i++)
	{
		Can->cd();
		Pads[i]->Draw();
	}
	
	double bottommargin=0.12;
	double XMin=-0.2;
	double XMax=5.8;
	
	TLine *TwokeVLine=new TLine();
	TwokeVLine->SetLineColor(12);
	TwokeVLine->SetLineStyle(9);
	double *Xarray =gUncorrectedJune[0]->GetX();
	//for(int i= 0;i<15;i++)
	//	cout << i << " " << Xarray[i]<< endl;
	Pads[0]->cd();
	gPad->SetRightMargin(0);
	gPad->SetBottomMargin(bottommargin);
	
	//gUncorrectedJune[0]->Print();
	gUncorrectedJune[0]->Draw("ap");
	gUncorrectedJune[0]->GetXaxis()->SetLimits(XMin,XMax);
	gUncorrectedJune[0]->GetYaxis()->SetNdivisions(509);
	gUncorrectedJune[0]->GetXaxis()->SetTitleOffset(1.035);
	gUncorrectedJune[0]->GetXaxis()->SetRangeUser(XMin,XMax);
	gCorr1T10DeriMaxJune[0]->Draw("p");
	//gCorr1T10DeriMaxJuneSmall[0]->Draw("p");
	gJuneLabEnergies->Draw("p");
	gJuneLabEnergiesCorr->Draw("p");
	//gJuneLabEnergiesCorrSmall->Draw("p");
	
	gJulyLabEnergies->Draw("p");
	
	//TwokeVLine->DrawLine(0,2,10,2);
	
	
	Pads[1]->cd();
		gPad->SetLeftMargin(0);
		//gPad->SetBottomMargin(bottommargin);
		gPad->SetBottomMargin(SplitVertical+bottommargin*(1-SplitVertical));
		gPad->SetTopMargin(gPad->GetTopMargin()*(1-SplitVertical));
		gJuneLabEnergiesAnnealing->GetYaxis()->SetRangeUser(0,10);
		gJuneLabEnergiesAnnealing->GetYaxis()->SetTickLength(gJuneLabEnergiesAnnealing->GetYaxis()->GetTickLength()*SplitHorizontal/(1.-SplitHorizontal) );
		
		gJuneLabEnergiesAnnealing->GetXaxis()->Set(3,0,3);
		//gJuneLabEnergiesAnnealing->GetXaxis()->SetLimits(0.1,2.9);
		gJuneLabEnergiesAnnealing->GetXaxis()->SetRangeUser(0,3);
		gJuneLabEnergiesAnnealing->Draw("ap");
		gJuneLabEnergiesAnnealing->GetYaxis()->SetNdivisions(509);
		gJulyLabEnergiesAnnealing->Draw("p");
		
		//TwokeVLine->DrawLine(0,2,10,2);
		
		
		
		gJuneLabEnergiesAnnealing->GetXaxis()->SetLabelSize(0);
		gJuneLabEnergiesAnnealing->GetYaxis()->SetLabelSize(0);
		gJuneLabEnergiesAnnealing->GetXaxis()->SetLabelOffset(0.004);
		gJuneLabEnergiesAnnealing->GetXaxis()->SetTickLength(0.03*(1-SplitVertical));
		//gJuneLabEnergiesAnnealing->GetXaxis()->CenterLabels();
		//gJuneLabEnergiesAnnealing->GetXaxis()->SetBinLabel(1,"#splitline{Prompt annealing}{(4 months)}");
		gJuneLabEnergiesAnnealing->GetXaxis()->SetBinLabel(1,"#splitline{  Prompt}{annealing}");
		//gJuneLabEnergiesAnnealing->GetXaxis()->SetBinLabel(2,"#splitline{Delayed annealing}{(4 years)}");
		gJuneLabEnergiesAnnealing->GetXaxis()->SetBinLabel(2,"#splitline{ Delayed}{annealing}");
		gJuneLabEnergiesAnnealing->GetXaxis()->SetBinLabel(3,"#splitline{ Second}{annealing}");
		
		TLatex* BinLabels=new TLatex();
		BinLabels->SetTextSize(0.1);
		BinLabels->SetTextAngle(-70);
		BinLabels->SetTextAlign(22);
		BinLabels->DrawLatexNDC(0.24,0.14,"#splitline{Prompt}{annealing}");
		BinLabels->DrawLatexNDC(0.54,0.14,"#splitline{Delayed}{annealing}");
		BinLabels->DrawLatexNDC(0.84,0.14,"#splitline{Second}{annealing}");
		
		//for(int i = 1;i<4;i++)
		//	gJuneLabEnergiesAnnealing->GetXaxis()->ChangeLabel();//i,90,0.05);
		
		TGraph *gDummy1=new TGraph();
		gDummy1->SetMarkerStyle(24);
		gDummy1->SetMarkerColor(kBlack);
		gDummy1->SetMarkerSize(2.5);
		TGraph *gDummy2=new TGraph();
		gDummy2->SetMarkerStyle(25);
		gDummy2->SetMarkerColor(kRed);
		gDummy2->SetMarkerSize(2.5);
		
		TLegend *leg = new TLegend(0.08,0.66,0.58,0.85);
		leg->SetTextSize(0.07);
		leg->AddEntry(gUncorrectedJune[0],"HEX 158 uncorrected","p");
		//leg->AddEntry(gDummy1,"HEX 158 corrected","p");
		leg->AddEntry(gCorr1T10DeriMaxJune[0],"HEX 158 corrected","p");
		leg->AddEntry(gJuneLabEnergies,"HEX 158 lab uncorr.","p");
		leg->AddEntry(gJuneLabEnergiesCorr,"HEX 158 lab corr.","p");
		leg->AddEntry(gJulyLabEnergies,"HEX 146 lab uncorr.","p");
		
		leg->Draw();
	
	
	Pads[2]->cd();
		gPad->SetRightMargin(0);
		
		TGaxis *PandaAxis = new TGaxis(0.13,0.7,
									96.1/(93*XMax/5.6.),0.7,
									0,96.1,//93*XMax/5.6,
									510,"+S");
		PandaAxis->SetName("PandaAxis");
		PandaAxis->SetTitle("Days of PANDA beam time");
		PandaAxis->SetLabelSize(gUncorrectedJune[0]->GetXaxis()->GetLabelSize()*(1-SplitVertical)/SplitVertical);
		PandaAxis->SetTitleSize(gUncorrectedJune[0]->GetXaxis()->GetTitleSize()*(1-SplitVertical)/SplitVertical);
	   
		PandaAxis->SetTickLength(gUncorrectedJune[0]->GetXaxis()->GetTickLength()*(1-SplitVertical)/SplitVertical);
	   
		PandaAxis->Draw();
	
	
	if (1) 
		andi::saveCanvas_allFileNames(Can, "~/pictures/Mar_COSY_IrradiationSummaryPicure" );
	return 0;
	
}
	
	
	
	
	
	//for(int i =0;i<3;i++)
	//{
		//Can->cd();
		//Pads[i]->Draw();
		
		//Pads[i]->cd();
		////Can->cd(i+1);
		//gPad->SetBottomMargin(0.15);
		////gCorr1TDeriMax90June[i]->Draw("apl");
		//gUncorrectedJune[i]->Draw("apl");
		//gUncorrectedJune[i]->GetXaxis()->SetTitleOffset(1.35);
		//gPad->Update();
		//gPad->Modified();
		//ProblemArea->DrawBox(0.3,gPad->GetUymin(),1.3,gPad->GetUymax());
		//gUncorrectedJune[i]->Draw("pl");
							
		//gCorr1T10DeriMaxJune[i]->Draw("pl");
		//gCorr1T1090Corr2T10DMJune[i]->Draw("pl");
		//gCorr1T1090Corr2TDM90June[i]->Draw("pl");
		//gCorr1T10DMCorr2T1090June[i]->Draw("pl");
		//gCorr1T10DMCorr2T90DMJune[i]->Draw("pl");
		//gCorr1TDM90Corr2T1090June[i]->Draw("pl");
		//gCorr1TDM90Corr2T10DMJune[i]->Draw("pl");
		//ProblemText->DrawLatexNDC(0.205,0.89,"#splitline{Detector}{issues}");
	//}
	
	//Pads[0]->cd();
	////Can->cd(1);
		////gUncorrectedJune[0]->GetYaxis()->SetRangeUser(0,10);
		
		////ProblemText->DrawLatexNDC(0.25,0.67,"#splitline{Detector}{problems}");
		//gPad->Update();
		//gPad->Modified();
		
		//TLegend *leg=new TLegend(0.28,0.2,0.9,0.5);
			//leg->SetFillColor(0);
			//leg->SetFillStyle(0);
			//leg->SetTextSize(0.05);
			//leg->SetHeader("Corr 1, Corr 2");
			//leg->SetNColumns(2);
			//leg->AddEntry(gUncorrectedJune[0],"uncorr.","p");
			//leg->AddEntry(gCorr1T10DeriMaxJune[0],"t_{10Imax}","p");
				//leg->AddEntry(gCorr1T1090Corr2T10DMJune[0],"t_{1090} t_{10Imax}","p");
				//leg->AddEntry(gCorr1T1090Corr2TDM90June[0],"t_{1090} t_{Imax90}","p");
				//leg->AddEntry(gCorr1T10DMCorr2T1090June[0],"t_{10Imax} t_{1090}","p");
				//leg->AddEntry(gCorr1T10DMCorr2T90DMJune[0],"t_{10Imax} t_{Imax90}","p");
				//leg->AddEntry(gCorr1TDM90Corr2T1090June[0],"t_{Imax90} t_{1090}","p");
				//leg->AddEntry(gCorr1TDM90Corr2T10DMJune[0],"t_{Imax90} t_{10Imax}","p");
			////leg->AddEntry(gCorr1T1090June[0],"t_{1090}","p");
			////leg->AddEntry(gCorr1TDeriMax90June[0],"t_{Imax90}","p");
		//leg->Draw();
		
	//Pads[1]->cd();
	////Can->cd(2);
		////gUncorrectedJune[1]->GetYaxis()->SetRangeUser(0,35);
		////ProblemText->DrawLatexNDC(0.25,0.67,"#splitline{Detector}{problems}");
		//gPad->Update();
		//gPad->Modified();
		
	//Pads[2]->cd();
	////Can->cd(3);
	//gPad->Update();
		//gPad->Modified();
		////gUncorrectedJune[2]->GetYaxis()->SetRangeUser(0,4);
		////gUncorrectedJune[2]
		//TLine *GausLine=new TLine();
		//GausLine->SetLineColor(14);
		//GausLine->SetLineStyle(9);
		//gPad->Update();
		//gPad->Modified();
		//GausLine->DrawLine(0,1.8226,gPad->GetUxmax(),1.8226);
		//TText *GausText=new TText();
		//GausText->SetTextSize(0.05);
		//GausText->SetTextColor(14);
		//GausText->DrawText(2.3,1.6,"FWTM/FWHM gaussian");
		//gPad->Update();
		//gPad->Modified();
		//Can->SaveAs("Mar_APP_SecondOrderCorrectionJune.png");
		////if (1) andi::saveCanvas_allFileNames(Can, "~/pictures/Mar_COSY_IrradiationSummaryPicure" );
		//if (SaveImages) andi::saveCanvas_allFileNames(Can, "~/pictures/Mar_COSY_IrradiationSummaryPicure" );
	//return 0;
//}






	////TCanvas *Can=new TCanvas("Can","Can",1600,500*CanvasS);
	////Can->Divide(3,2);
	//////June
	////for(int i =0;i<6;i++)
	////{
		////Can->cd(i+1);
	
		////gUncorrectedJune[i]->Draw("apl");
		////gCorr1T1090June[i]->Draw("pl");
		////gCorr1T10DeriMaxJune[i]->Draw("pl");
		////gCorr1TDeriMax90June[i]->Draw("pl");
		////gCorr1T1090Corr2T10DMJune[i]->Draw("pl");
		////gCorr1T1090Corr2TDM90June[i]->Draw("pl");
		////gCorr1T10DMCorr2T1090June[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJune[i]->Draw("pl");
		////gCorr1TDM90Corr2T1090June[i]->Draw("pl");
		////gCorr1TDM90Corr2T10DMJune[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJuneCut[i]->Draw("pl");
	////}
	
	////TLegend *leg=new TLegend(0.3,0.1,0.9,0.3);
	////leg->SetFillColor(0);
	////leg->SetHeader("June2014");
	////leg->SetNColumns(2);
	////leg->AddEntry(gUncorrectedJune[0],"uncorrected","p");
	////leg->AddEntry(gCorr1T1090June[0],"Corr1 T1090","p");
	////leg->AddEntry(gCorr1T10DeriMaxJune[0],"Corr1 T10DM","p");
	////leg->AddEntry(gCorr1TDeriMax90June[0],"Corr1 TDM90","p");
	////leg->AddEntry(gCorr1T1090Corr2T10DMJune[0],"Corr1 T1090 Corr2 T10DM","p");
	////leg->AddEntry(gCorr1T1090Corr2TDM90June[0],"Corr1 T1090 Corr2 TDM90","p");
	////leg->AddEntry(gCorr1T10DMCorr2T1090June[0],"Corr1 T10DM Corr2 T1090","p");
	////leg->AddEntry(gCorr1T10DMCorr2T90DMJune[0],"Corr1 T10DM Corr2 T90DM","p");
	////leg->AddEntry(gCorr1TDM90Corr2T1090June[0],"Corr1 TDM90 Corr2 T1090","p");
	////leg->AddEntry(gCorr1TDM90Corr2T10DMJune[0],"Corr1 TDM90 Corr2 T10DM","p");
	////leg->AddEntry(gCorr1T10DMCorr2T90DMJuneCut[0],"Corr1 T10DM Corr2 T90DM (T10DM>50ns)","p");
	////Can->cd(1);
	////leg->Draw();
	////if (SaveImages) andi::saveCanvas_allFileNames(Can, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_2June2014" );
	////TCanvas *Can_0=new TCanvas("Can_0","Can_0",1600,1000);
	////Can_0->Divide(3,2);
	////for(int i =0;i<6;i++)
	////{
		////Can_0->cd(i+1);
	
		////gUncorrectedJune[i]->Draw("apl");
		////gCorr1T1090June_0[i]->Draw("pl");
		////gCorr1T10DeriMaxJune_0[i]->Draw("pl");
		////gCorr1TDeriMax90June_0[i]->Draw("pl");
		////gCorr1T1090Corr2T10DMJune_0[i]->Draw("pl");
		////gCorr1T1090Corr2TDM90June_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T1090June_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJune_0[i]->Draw("pl");
		////gCorr1TDM90Corr2T1090June_0[i]->Draw("pl");
		////gCorr1TDM90Corr2T10DMJune_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJune_0Cut[i]->Draw("pl");
	////}
	////Can_0->cd(1);
	////leg->Draw();
	////if (SaveImages) andi::saveCanvas_allFileNames(Can_0, "~/pictures/Mar_COSY_ThreeFirstOrderCorrectionJune" );
	
	//////July
	////TCanvas *CanJuly=new TCanvas("CanJuly","CanJuly",1600,1000);
	////CanJuly->Divide(3,2);
	////for(int i =0;i<6;i++)
	////{
		////CanJuly->cd(i+1);
	
		////gUncorrectedJuly[i]->Draw("apl");
		////gCorr1T1090July[i]->Draw("pl");
		////gCorr1T10DeriMaxJuly[i]->Draw("pl");
		////gCorr1TDeriMax90July[i]->Draw("pl");
		////gCorr1T1090Corr2T10DMJuly[i]->Draw("pl");
		////gCorr1T1090Corr2TDM90July[i]->Draw("pl");
		////gCorr1T10DMCorr2T1090July[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJuly[i]->Draw("pl");
		////gCorr1TDM90Corr2T1090July[i]->Draw("pl");
		////gCorr1TDM90Corr2T10DMJuly[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJulyCut[i]->Draw("pl");
	////}
	
	////TLegend *leg=new TLegend(0.3,0.1,0.9,0.3);
	////leg->SetFillColor(0);
	////leg->SetHeader("July2014");
	////leg->SetNColumns(2);
	////leg->AddEntry(gUncorrectedJuly[0],"uncorrected","p");
	////leg->AddEntry(gCorr1T1090July[0],"Corr1 T1090","p");
	////leg->AddEntry(gCorr1T10DeriMaxJuly[0],"Corr1 T10DM","p");
	////leg->AddEntry(gCorr1TDeriMax90July[0],"Corr1 TDM90","p");
	////leg->AddEntry(gCorr1T1090Corr2T10DMJuly[0],"Corr1 T1090 Corr2 T10DM","p");
	////leg->AddEntry(gCorr1T1090Corr2TDM90July[0],"Corr1 T1090 Corr2 TDM90","p");
	////leg->AddEntry(gCorr1T10DMCorr2T1090July[0],"Corr1 T10DM Corr2 T1090","p");
	////leg->AddEntry(gCorr1T10DMCorr2T90DMJuly[0],"Corr1 T10DM Corr2 T90DM","p");
	////leg->AddEntry(gCorr1TDM90Corr2T1090July[0],"Corr1 TDM90 Corr2 T1090","p");
	////leg->AddEntry(gCorr1TDM90Corr2T10DMJuly[0],"Corr1 TDM90 Corr2 T10DM","p");
	////leg->AddEntry(gCorr1T10DMCorr2T90DMJulyCut[0],"Corr1 T10DM Corr2 T90DM (T10DM>50ns)","p");
	////CanJuly->cd(1);
	////leg->Draw();
	////if (SaveImages) andi::saveCanvas_allFileNames(CanJuly, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_2July2014" );
	////TCanvas *CanJuly_0=new TCanvas("CanJuly_0","CanJuly_0",1600,1000);
	////CanJuly_0->Divide(3,2);
	////for(int i =0;i<6;i++)
	////{
		////CanJuly_0->cd(i+1);
	
		////gUncorrectedJuly[i]->Draw("apl");
		////gCorr1T1090July_0[i]->Draw("pl");
		////gCorr1T10DeriMaxJuly_0[i]->Draw("pl");
		////gCorr1TDeriMax90July_0[i]->Draw("pl");
		////gCorr1T1090Corr2T10DMJuly_0[i]->Draw("pl");
		////gCorr1T1090Corr2TDM90July_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T1090July_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJuly_0[i]->Draw("pl");
		////gCorr1TDM90Corr2T1090July_0[i]->Draw("pl");
		////gCorr1TDM90Corr2T10DMJuly_0[i]->Draw("pl");
		////gCorr1T10DMCorr2T90DMJuly_0Cut[i]->Draw("pl");
	////}
	////CanJuly_0->cd(1);
	////leg->Draw();
	////if (SaveImages) andi::saveCanvas_allFileNames(CanJuly_0, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_0July2014" );
	
	
////}
