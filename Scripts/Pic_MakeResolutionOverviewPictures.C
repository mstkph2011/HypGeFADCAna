void HistoMakeup(TGraph *histo[6],Color_t color, Style_t MarkerStyle,double MarkerSize=1.5,bool useDottedLines=1)
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


void FillGraphs(TGraph *Graphs[6],TString Filename,const int nPoints=15)
{
	//TGraph *Graphs[4];
	//**Graphs=new TGraph*[4];
	double PandaDays[nPoints];
	double FWHMCo[nPoints],FWTMCo[nPoints],FWTMHMCo[nPoints],FWHM511[nPoints],FWTM511[nPoints],FWTMHM511[nPoints];
	double pd,fhCo,ftCo,fh511, ft511;
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
		sscanf(buf,"%*lf\t%lf\t%lf\t%lf\t%lf\t%lf",&pd,&fhCo,&ftCo,&fh511,&ft511);
		cout << buf << endl;
		cout << "\t" << pd << " " << fhCo << " " << ftCo << " " << fh511 << " " << ft511 << endl;
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
	for(int i=0;i<6;i++)
		Graphs[i]->Sort();
	cout << "asdasd" <<endl;
	//return Graphs;
}



void Pic_MakeResolutionOverviewPictures()
{
	andi::setCustomStyle(0,0);
	 gROOT->ForceStyle();
	
	//June2014
	TGraph *gUncorrectedJune[6];
	gUncorrectedJune=FillGraphs(gUncorrectedJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014.txt",15);
	HistoMakeup( gUncorrectedJune,1, 20);
	//first corrections
	TGraph *gCorr1T1090June[6];
	gUncorrectedJune=FillGraphs(gCorr1T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T1090June,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJune[6];
	gUncorrectedJune=FillGraphs(gCorr1T10DeriMaxJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T10DeriMaxJune,kGreen+2, 22);
	TGraph *gCorr1TDeriMax90June[6];
	gUncorrectedJune=FillGraphs(gCorr1TDeriMax90June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1TDeriMax90June,2, 23);
	//second corrections
	TGraph *gCorr1T1090Corr2T10DMJune[6];
	gUncorrectedJune=FillGraphs(gCorr1T1090Corr2T10DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T1090Corr2T10DMJune,9, 24);
	TGraph *gCorr1T1090Corr2TDM90June[6];
	gUncorrectedJune=FillGraphs(gCorr1T1090Corr2TDM90June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T1090Corr2TDM90June,38, 25);
	
	TGraph *gCorr1T10DMCorr2T1090June[6];
	gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T10DMCorr2T1090June,8, 26);
	TGraph *gCorr1T10DMCorr2T90DMJune[6];
	gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1T10DMCorr2T90DMJune,30, 27);
	
	TGraph *gCorr1TDM90Corr2T1090June[6];
	gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2.txt",15);
	HistoMakeup( gCorr1TDM90Corr2T1090June,46, 28);
	TGraph *gCorr1TDM90Corr2T10DMJune[6];
	gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T10DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",15);
	HistoMakeup( gCorr1TDM90Corr2T10DMJune,6, 29);
	
	
	//T10DM T90DM, cut on T10DM >50 ns
	TGraph *gCorr1T10DMCorr2T90DMJuneCut[6];
	gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJuneCut, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2cut.txt",15);
	HistoMakeup(gCorr1T10DMCorr2T90DMJuneCut,3, 47);
	
	
	TCanvas *Can=new TCanvas("Can","Can",1600,1000);
	Can->Divide(3,2);
	
	for(int i =0;i<6;i++)
	{
		Can->cd(i+1);
	
		gUncorrectedJune[i]->Draw("apl");
		gCorr1T1090June[i]->Draw("pl");
		gCorr1T10DeriMaxJune[i]->Draw("pl");
		gCorr1TDeriMax90June[i]->Draw("pl");
		gCorr1T1090Corr2T10DMJune[i]->Draw("pl");
		gCorr1T1090Corr2TDM90June[i]->Draw("pl");
		gCorr1T10DMCorr2T1090June[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJune[i]->Draw("pl");
		gCorr1TDM90Corr2T1090June[i]->Draw("pl");
		gCorr1TDM90Corr2T10DMJune[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJuneCut[i]->Draw("pl");
	}
	TLegend *leg=new TLegend(0.3,0.1,0.9,0.3);
	leg->SetFillColor(0);
	leg->SetHeader("June2014");
	leg->SetNColumns(2);
	leg->AddEntry(gUncorrectedJune[0],"uncorrected","p");
	leg->AddEntry(gCorr1T1090June[0],"Corr1 T1090","p");
	leg->AddEntry(gCorr1T10DeriMaxJune[0],"Corr1 T10DM","p");
	leg->AddEntry(gCorr1TDeriMax90June[0],"Corr1 TDM90","p");
	leg->AddEntry(gCorr1T1090Corr2T10DMJune[0],"Corr1 T1090 Corr2 T10DM","p");
	leg->AddEntry(gCorr1T1090Corr2TDM90June[0],"Corr1 T1090 Corr2 TDM90","p");
	leg->AddEntry(gCorr1T10DMCorr2T1090June[0],"Corr1 T10DM Corr2 T1090","p");
	leg->AddEntry(gCorr1T10DMCorr2T90DMJune[0],"Corr1 T10DM Corr2 T90DM","p");
	leg->AddEntry(gCorr1TDM90Corr2T1090June[0],"Corr1 TDM90 Corr2 T1090","p");
	leg->AddEntry(gCorr1TDM90Corr2T10DMJune[0],"Corr1 TDM90 Corr2 T10DM","p");
	leg->AddEntry(gCorr1T10DMCorr2T90DMJuneCut[0],"Corr1 T10DM Corr2 T90DM (T10DM>50ns)","p");
	Can->cd(1);
	leg->Draw();
}
