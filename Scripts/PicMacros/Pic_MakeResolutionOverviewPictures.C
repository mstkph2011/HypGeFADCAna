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
	 
	bool SaveImages=0;
	
	const int nFilesJune=15;
	const int nFilesJuly=10;
	//June2014
	TGraph *gUncorrectedJune[6];
		gUncorrectedJune=FillGraphs(gUncorrectedJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014.txt",nFilesJune);
		HistoMakeup( gUncorrectedJune,1, 20);
	//first corrections
	TGraph *gCorr1T1090June[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T1090June,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJune[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DeriMaxJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T10DeriMaxJune,kGreen+2, 22);
	TGraph *gCorr1TDeriMax90June[6];
		gUncorrectedJune=FillGraphs(gCorr1TDeriMax90June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1TDeriMax90June,2, 23);
	//second corrections
	TGraph *gCorr1T1090Corr2T10DMJune[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090Corr2T10DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T1090Corr2T10DMJune,9, 24);
	TGraph *gCorr1T1090Corr2TDM90June[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090Corr2TDM90June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T1090Corr2TDM90June,38, 25);
	
	TGraph *gCorr1T10DMCorr2T1090June[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T10DMCorr2T1090June,8, 26);
	TGraph *gCorr1T10DMCorr2T90DMJune[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1T10DMCorr2T90DMJune,30, 27);
	
	TGraph *gCorr1TDM90Corr2T1090June[6];
		gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T1090June, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1TDM90Corr2T1090June,46, 28);
	TGraph *gCorr1TDM90Corr2T10DMJune[6];
		gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T10DMJune, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",nFilesJune);
		HistoMakeup( gCorr1TDM90Corr2T10DMJune,6, 29);
	//T10DM T90DM, cut on T10DM >50 ns
	TGraph *gCorr1T10DMCorr2T90DMJuneCut[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJuneCut, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2cut.txt",nFilesJune);
		HistoMakeup(gCorr1T10DMCorr2T90DMJuneCut,3, 47);
	
	//corrections on 511 keV line
	//first corrections
	TGraph *gCorr1T1090June_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090June_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T1090June_0,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJune_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DeriMaxJune_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T10DeriMaxJune_0,kGreen+2, 22);
	TGraph *gCorr1TDeriMax90June_0[6];
		gUncorrectedJune=FillGraphs(gCorr1TDeriMax90June_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1TDeriMax90June_0,2, 23);
	//second corrections
	TGraph *gCorr1T1090Corr2T10DMJune_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090Corr2T10DMJune_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_0Corr2_T10DeriMaxEnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T1090Corr2T10DMJune_0,9, 24);
	TGraph *gCorr1T1090Corr2TDM90June_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T1090Corr2TDM90June_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T1090EnergyNorm_0Corr2_TDeriMax90EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T1090Corr2TDM90June_0,38, 25);
	
	TGraph *gCorr1T10DMCorr2T1090June_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T1090June_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_0Corr2_T1090EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T10DMCorr2T1090June_0,8, 26);
	TGraph *gCorr1T10DMCorr2T90DMJune_0[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJune_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_0Corr2_TDeriMax90EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1T10DMCorr2T90DMJune_0,30, 27);
	
	TGraph *gCorr1TDM90Corr2T1090June_0[6];
		gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T1090June_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_0Corr2_T1090EnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1TDM90Corr2T1090June_0,46, 28);
	TGraph *gCorr1TDM90Corr2T10DMJune_0[6];
		gUncorrectedJune=FillGraphs(gCorr1TDM90Corr2T10DMJune_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_TDeriMax90EnergyNorm_0Corr2_T10DeriMaxEnergyNorm_0.txt",nFilesJune);
		HistoMakeup( gCorr1TDM90Corr2T10DMJune_0,6, 29);
	//T10DM T90DM, cut on T10DM >50 ns
	TGraph *gCorr1T10DMCorr2T90DMJune_0Cut[6];
		gUncorrectedJune=FillGraphs(gCorr1T10DMCorr2T90DMJune_0Cut, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJune2014Corr1_T10DeriMaxEnergyNorm_0Corr2_TDeriMax90EnergyNorm_0cut.txt",nFilesJune);
		HistoMakeup(gCorr1T10DMCorr2T90DMJune_0Cut,3, 47);
	
	//July2014
	TGraph *gUncorrectedJuly[6];
		gUncorrectedJuly=FillGraphs(gUncorrectedJuly, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014.txt",nFilesJuly);
		HistoMakeup( gUncorrectedJuly,1, 20);
	//first corrections
	TGraph *gCorr1T1090July[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090July, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090July,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJuly[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DeriMaxJuly, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DeriMaxJuly,kGreen+2, 22);
	TGraph *gCorr1TDeriMax90July[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDeriMax90July, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1TDeriMax90July,2, 23);
	//second corrections
	TGraph *gCorr1T1090Corr2T10DMJuly[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090Corr2T10DMJuly, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090Corr2T10DMJuly,9, 24);
	TGraph *gCorr1T1090Corr2TDM90July[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090Corr2TDM90July, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090Corr2TDM90July,38, 25);
	
	TGraph *gCorr1T10DMCorr2T1090July[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T1090July, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DMCorr2T1090July,8, 26);
	TGraph *gCorr1T10DMCorr2T90DMJuly[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T90DMJuly, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DMCorr2T90DMJuly,30, 27);
	
	TGraph *gCorr1TDM90Corr2T1090July[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDM90Corr2T1090July, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1TDM90Corr2T1090July,46, 28);
	TGraph *gCorr1TDM90Corr2T10DMJuly[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDM90Corr2T10DMJuly, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2.txt",nFilesJuly);
		HistoMakeup( gCorr1TDM90Corr2T10DMJuly,6, 29);
	//T10DM T90DM, cut on T10DM >50 ns
	TGraph *gCorr1T10DMCorr2T90DMJulyCut[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T90DMJulyCut, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2cut.txt",nFilesJuly);
		HistoMakeup(gCorr1T10DMCorr2T90DMJulyCut,3, 47);
	
	//corrections on 511 keV line
	//first corrections
	TGraph *gCorr1T1090July_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090July_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090July_0,kBlue, 21);
	TGraph *gCorr1T10DeriMaxJuly_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DeriMaxJuly_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DeriMaxJuly_0,kGreen+2, 22);
	TGraph *gCorr1TDeriMax90July_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDeriMax90July_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1TDeriMax90July_0,2, 23);
	//second corrections
	TGraph *gCorr1T1090Corr2T10DMJuly_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090Corr2T10DMJuly_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_0Corr2_T10DeriMaxEnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090Corr2T10DMJuly_0,9, 24);
	TGraph *gCorr1T1090Corr2TDM90July_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T1090Corr2TDM90July_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T1090EnergyNorm_0Corr2_TDeriMax90EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T1090Corr2TDM90July_0,38, 25);
	
	TGraph *gCorr1T10DMCorr2T1090July_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T1090July_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_0Corr2_T1090EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DMCorr2T1090July_0,8, 26);
	TGraph *gCorr1T10DMCorr2T90DMJuly_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T90DMJuly_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_0Corr2_TDeriMax90EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1T10DMCorr2T90DMJuly_0,30, 27);
	
	TGraph *gCorr1TDM90Corr2T1090July_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDM90Corr2T1090July_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_0Corr2_T1090EnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1TDM90Corr2T1090July_0,46, 28);
	TGraph *gCorr1TDM90Corr2T10DMJuly_0[6];
		gUncorrectedJuly=FillGraphs(gCorr1TDM90Corr2T10DMJuly_0, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_TDeriMax90EnergyNorm_0Corr2_T10DeriMaxEnergyNorm_0.txt",nFilesJuly);
		HistoMakeup( gCorr1TDM90Corr2T10DMJuly_0,6, 29);
	//T10DM T90DM, cut on T10DM >50 ns
	TGraph *gCorr1T10DMCorr2T90DMJuly_0Cut[6];
		gUncorrectedJuly=FillGraphs(gCorr1T10DMCorr2T90DMJuly_0Cut, "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/CollectionFileJuly2014Corr1_T10DeriMaxEnergyNorm_0Corr2_TDeriMax90EnergyNorm_0cut.txt",nFilesJuly);
		HistoMakeup(gCorr1T10DMCorr2T90DMJuly_0Cut,3, 47);
	
	
	
	TCanvas *Can=new TCanvas("Can","Can",1600,1000);
	Can->Divide(3,2);
	//June
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
	if (SaveImages) andi::saveCanvas_allFileNames(Can, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_2June2014" );
	TCanvas *Can_0=new TCanvas("Can_0","Can_0",1600,1000);
	Can_0->Divide(3,2);
	for(int i =0;i<6;i++)
	{
		Can_0->cd(i+1);
	
		gUncorrectedJune[i]->Draw("apl");
		gCorr1T1090June_0[i]->Draw("pl");
		gCorr1T10DeriMaxJune_0[i]->Draw("pl");
		gCorr1TDeriMax90June_0[i]->Draw("pl");
		gCorr1T1090Corr2T10DMJune_0[i]->Draw("pl");
		gCorr1T1090Corr2TDM90June_0[i]->Draw("pl");
		gCorr1T10DMCorr2T1090June_0[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJune_0[i]->Draw("pl");
		gCorr1TDM90Corr2T1090June_0[i]->Draw("pl");
		gCorr1TDM90Corr2T10DMJune_0[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJune_0Cut[i]->Draw("pl");
	}
	Can_0->cd(1);
	leg->Draw();
	if (SaveImages) andi::saveCanvas_allFileNames(Can_0, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_0June2014" );
	
	//July
	TCanvas *CanJuly=new TCanvas("CanJuly","CanJuly",1600,1000);
	CanJuly->Divide(3,2);
	for(int i =0;i<6;i++)
	{
		CanJuly->cd(i+1);
	
		gUncorrectedJuly[i]->Draw("apl");
		gCorr1T1090July[i]->Draw("pl");
		gCorr1T10DeriMaxJuly[i]->Draw("pl");
		gCorr1TDeriMax90July[i]->Draw("pl");
		gCorr1T1090Corr2T10DMJuly[i]->Draw("pl");
		gCorr1T1090Corr2TDM90July[i]->Draw("pl");
		gCorr1T10DMCorr2T1090July[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJuly[i]->Draw("pl");
		gCorr1TDM90Corr2T1090July[i]->Draw("pl");
		gCorr1TDM90Corr2T10DMJuly[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJulyCut[i]->Draw("pl");
	}
	
	TLegend *leg=new TLegend(0.3,0.1,0.9,0.3);
	leg->SetFillColor(0);
	leg->SetHeader("July2014");
	leg->SetNColumns(2);
	leg->AddEntry(gUncorrectedJuly[0],"uncorrected","p");
	leg->AddEntry(gCorr1T1090July[0],"Corr1 T1090","p");
	leg->AddEntry(gCorr1T10DeriMaxJuly[0],"Corr1 T10DM","p");
	leg->AddEntry(gCorr1TDeriMax90July[0],"Corr1 TDM90","p");
	leg->AddEntry(gCorr1T1090Corr2T10DMJuly[0],"Corr1 T1090 Corr2 T10DM","p");
	leg->AddEntry(gCorr1T1090Corr2TDM90July[0],"Corr1 T1090 Corr2 TDM90","p");
	leg->AddEntry(gCorr1T10DMCorr2T1090July[0],"Corr1 T10DM Corr2 T1090","p");
	leg->AddEntry(gCorr1T10DMCorr2T90DMJuly[0],"Corr1 T10DM Corr2 T90DM","p");
	leg->AddEntry(gCorr1TDM90Corr2T1090July[0],"Corr1 TDM90 Corr2 T1090","p");
	leg->AddEntry(gCorr1TDM90Corr2T10DMJuly[0],"Corr1 TDM90 Corr2 T10DM","p");
	leg->AddEntry(gCorr1T10DMCorr2T90DMJulyCut[0],"Corr1 T10DM Corr2 T90DM (T10DM>50ns)","p");
	CanJuly->cd(1);
	leg->Draw();
	if (SaveImages) andi::saveCanvas_allFileNames(CanJuly, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_2July2014" );
	TCanvas *CanJuly_0=new TCanvas("CanJuly_0","CanJuly_0",1600,1000);
	CanJuly_0->Divide(3,2);
	for(int i =0;i<6;i++)
	{
		CanJuly_0->cd(i+1);
	
		gUncorrectedJuly[i]->Draw("apl");
		gCorr1T1090July_0[i]->Draw("pl");
		gCorr1T10DeriMaxJuly_0[i]->Draw("pl");
		gCorr1TDeriMax90July_0[i]->Draw("pl");
		gCorr1T1090Corr2T10DMJuly_0[i]->Draw("pl");
		gCorr1T1090Corr2TDM90July_0[i]->Draw("pl");
		gCorr1T10DMCorr2T1090July_0[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJuly_0[i]->Draw("pl");
		gCorr1TDM90Corr2T1090July_0[i]->Draw("pl");
		gCorr1TDM90Corr2T10DMJuly_0[i]->Draw("pl");
		gCorr1T10DMCorr2T90DMJuly_0Cut[i]->Draw("pl");
	}
	CanJuly_0->cd(1);
	leg->Draw();
	if (SaveImages) andi::saveCanvas_allFileNames(CanJuly_0, "~/pictures/JuelichResolution/Mar_COSY_AllCorrections_0July2014" );
	
	
}
