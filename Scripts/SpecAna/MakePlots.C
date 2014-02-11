void GraphCosmetics(TGraph *graph, int color, int MarkerStyle)
{
	graph->SetMarkerStyle(MarkerStyle);
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
}

void FirstGraphSettings(TGraph *graph, char* Title, char* XTitle,Bool_t XCentered, Float_t XOffset, Float_t XSize, char* YTitle, Bool_t YCentered, Float_t YOffset, Float_t YSize, Double_t YMin, Double_t YMax)
{
	graph->SetTitle(Title);
	graph->GetHistogram()->SetXTitle(XTitle);
	graph->GetHistogram()->GetXaxis()->CenterTitle(XCentered);
	graph->GetHistogram()->GetXaxis()->SetTitleOffset(XOffset);
	graph->GetHistogram()->GetXaxis()->SetTitleSize(XSize);
	
	graph->GetHistogram()->SetYTitle(YTitle);
	graph->GetHistogram()->GetYaxis()->CenterTitle(YCentered);
	graph->GetHistogram()->GetYaxis()->SetTitleOffset(YOffset);
	graph->GetHistogram()->GetYaxis()->SetTitleSize(YSize);
	graph->GetHistogram()->SetMinimum(YMin);
	graph->GetHistogram()->SetMaximum(YMax);
}
void MakePlots()
{
	ifstream InputFile("Output.txt");
	Double_t FWHMValuesBil[130][6];
	Double_t FWHMValuesGaus[130][6];
	Int_t MRange [13] = {60,80,100,120,140,160,180,200,220,240,260,280,300};
	Int_t SigmaBilRange [10] = {300,700,900,1100,1300,1500,1700,1900,2100};
	Double_t SigmaGausRange [6] = {1,3,5,7,9,11};
	
	TGraph * GraphBil[130];
	TGraph * GraphGaus[13];
	
	
	char buffer[200] = "";
	InputFile.getline(buffer,200);	//
	InputFile.getline(buffer,200);	// line with SigmaGaus values
	InputFile.getline(buffer,200);	// line "bil"
	for (int i = 0; i<130; i++)
	{
		InputFile.getline(buffer,200);
		for (int j = 0;j<5;j++)
		{
			FWHMValuesBil[i][j] = 0;			//init value array
		}
		InputFile.getline(buffer,200);
		TString buf = buffer;
		switch (buf.CountChar('\t'))		//find number of values and extract
		{
			case 1: sscanf(buffer,"%lf\t",&FWHMValuesBil[i][0]);
						break;
			case 2: sscanf(buffer,"%lf\t%lf\t",&FWHMValuesBil[i][0],&FWHMValuesBil[i][1]);
						break;
			case 3: sscanf(buffer,"%lf\t%lf\t%lf\t",&FWHMValuesBil[i][0],&FWHMValuesBil[i][1],&FWHMValuesBil[i][2]);
						break;
			case 4: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t",&FWHMValuesBil[i][0],&FWHMValuesBil[i][1],&FWHMValuesBil[i][2],&FWHMValuesBil[i][3]);
						break;
			case 5: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t",&FWHMValuesBil[i][0],&FWHMValuesBil[i][1],&FWHMValuesBil[i][2],&FWHMValuesBil[i][3],&FWHMValuesBil[i][4]);
						break;
			case 6: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",&FWHMValuesBil[i][0],&FWHMValuesBil[i][1],&FWHMValuesBil[i][2],&FWHMValuesBil[i][3],&FWHMValuesBil[i][4],&FWHMValuesBil[i][5]);
		}
		for (int j = 0; j <buf.CountChar('\t'); j++)
			cout << FWHMValuesBil[i][j] << ", ";
		cout << endl;
		GraphBil[i]= new TGraph(6, SigmaGausRange,FWHMValuesBil[i]);
	}
	InputFile.getline(buffer,200);	// line "gaus"
	for (int i = 0; i<13; i++)
	{
		InputFile.getline(buffer,200);
		for (int j = 0;j<5;j++)
		{
			FWHMValuesGaus[i][j]=0;
		}
		InputFile.getline(buffer,200);
		TString buf = buffer;
		switch (buf.CountChar('\t'))		//find number of values and extract
		{
			case 1: sscanf(buffer,"%lf\t",&FWHMValuesGaus[i][0]);
						break;
			case 2: sscanf(buffer,"%lf\t%lf\t",&FWHMValuesGaus[i][0],&FWHMValuesGaus[i][1]);
						break;
			case 3: sscanf(buffer,"%lf\t%lf\t%lf\t",&FWHMValuesGaus[i][0],&FWHMValuesGaus[i][1],&FWHMValuesGaus[i][2]);
						break;
			case 4: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t",&FWHMValuesGaus[i][0],&FWHMValuesGaus[i][1],&FWHMValuesGaus[i][2],&FWHMValuesGaus[i][3]);
						break;
			case 5: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t",&FWHMValuesGaus[i][0],&FWHMValuesGaus[i][1],&FWHMValuesGaus[i][2],&FWHMValuesGaus[i][3],&FWHMValuesGaus[i][4]);
						break;
			case 6: sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",&FWHMValuesGaus[i][0],&FWHMValuesGaus[i][1],&FWHMValuesGaus[i][2],&FWHMValuesGaus[i][3],&FWHMValuesGaus[i][4],&FWHMValuesGaus[i][5]);
		}
		for (int j = 0; j <buf.CountChar('\t'); j++)
			cout << FWHMValuesGaus[i][j] << ", ";
		cout << endl;
		GraphGaus[i]= new TGraph(6, SigmaGausRange,FWHMValuesGaus[i]);
	}
	
	TCanvas *Canvas[13];
	TLegend *Legend[13];
	char CanvasName[70];
	char PictureName[70];
	
	Int_t YRangeMax [13] = {25,25,25,25,25,25,25,25,25,25,25,25,25};
	for (int i = 0; i < 13; i ++)
	{
		sprintf(CanvasName,"Results with Bilateral Filter, M = %d",MRange[i]); 
		Canvas[i] = new TCanvas(CanvasName,CanvasName,800,600);
		
		GraphCosmetics(GraphBil[10*i],1,2);
		
		FirstGraphSettings(GraphBil[10*i],CanvasName,"#sigma_{gaus}",true, 1,0.04,"FWHM [keV]",true, 1,0.04,0,YRangeMax[i]);
		
		
		GraphCosmetics(GraphBil[10*i+1],2,3);
		GraphCosmetics(GraphBil[10*i+2],3,4);
		GraphCosmetics(GraphBil[10*i+3],4,5);
		GraphCosmetics(GraphBil[10*i+4],5,8);
		GraphCosmetics(GraphBil[10*i+5],6,21);
		GraphCosmetics(GraphBil[10*i+6],7,22);
		GraphCosmetics(GraphBil[10*i+7],8,23);
		GraphCosmetics(GraphBil[10*i+8],9,28);
		GraphCosmetics(GraphBil[10*i+9],14,29);
		
		GraphBil[10*i]->Draw("APL");
		GraphBil[10*i+1]->Draw("PL");
		GraphBil[10*i+2]->Draw("PL");
		GraphBil[10*i+3]->Draw("PL");
		GraphBil[10*i+4]->Draw("PL");
		GraphBil[10*i+5]->Draw("PL");
		GraphBil[10*i+6]->Draw("PL");
		GraphBil[10*i+7]->Draw("PL");
		GraphBil[10*i+8]->Draw("PL");
		GraphBil[10*i+9]->Draw("PL");
		
		Legend[i]= new TLegend(0.75,0.6,0.9,0.9);
		Legend[i]->SetFillColor(0);
		Legend[i]->AddEntry(GraphBil[10*i],"#sigma_{Bil} = 300","lp");
		Legend[i]->AddEntry(GraphBil[10*i+1],"#sigma_{Bil} = 500","lp");
		Legend[i]->AddEntry(GraphBil[10*i+2],"#sigma_{Bil} = 700","lp");
		Legend[i]->AddEntry(GraphBil[10*i+3],"#sigma_{Bil} = 900","lp");
		Legend[i]->AddEntry(GraphBil[10*i+4],"#sigma_{Bil} = 1100","lp");
		Legend[i]->AddEntry(GraphBil[10*i+5],"#sigma_{Bil} = 1300","lp");
		Legend[i]->AddEntry(GraphBil[10*i+6],"#sigma_{Bil} = 1500","lp");
		Legend[i]->AddEntry(GraphBil[10*i+7],"#sigma_{Bil} = 1700","lp");
		Legend[i]->AddEntry(GraphBil[10*i+8],"#sigma_{Bil} = 1900","lp");
		Legend[i]->AddEntry(GraphBil[10*i+9],"#sigma_{Bil} = 2100","lp");
		Legend[i]->Draw();
		
		sprintf(PictureName,"pics/BilFilterWithM%d.png",MRange[i]);
		Canvas[i]->Print(PictureName,"png");
	}
	
	TCanvas *CanvasGaus;
	TLegend *LegendGaus;
		CanvasGaus = new TCanvas("Gaussian Filter","Gaussian Filter",800,600);
		
		GraphCosmetics(GraphGaus[0],1,2);
		GraphCosmetics(GraphGaus[1],2,3);
		GraphCosmetics(GraphGaus[2],3,4);
		GraphCosmetics(GraphGaus[3],4,5);
		GraphCosmetics(GraphGaus[4],5,8);
		GraphCosmetics(GraphGaus[5],6,21);
		GraphCosmetics(GraphGaus[6],7,22);
		GraphCosmetics(GraphGaus[7],8,23);
		GraphCosmetics(GraphGaus[8],9,28);
		GraphCosmetics(GraphGaus[9],14,29);
		GraphCosmetics(GraphGaus[10],28,27);
		GraphCosmetics(GraphGaus[11],30,30);
		GraphCosmetics(GraphGaus[12],44,26);
		
		
		FirstGraphSettings(GraphGaus[0],"Results with gaussian filter","#sigma_{gaus}",true, 1,0.04,"FWHM [keV]",true, 1,0.04,0,25);
		
		
		GraphGaus[0]->Draw("APL");
		GraphGaus[1]->Draw("PL");
		GraphGaus[2]->Draw("PL");
		GraphGaus[3]->Draw("PL");
		GraphGaus[4]->Draw("PL");
		GraphGaus[5]->Draw("PL");
		GraphGaus[6]->Draw("PL");
		GraphGaus[7]->Draw("PL");
		GraphGaus[8]->Draw("PL");
		GraphGaus[9]->Draw("PL");
		GraphGaus[10]->Draw("PL");
		GraphGaus[11]->Draw("PL");
		GraphGaus[12]->Draw("PL");
		
		LegendGaus= new TLegend(0.8,0.6,0.9,0.9);
		LegendGaus->SetFillColor(0);
		LegendGaus->AddEntry(GraphGaus[0],"M = 60","lp");
		LegendGaus->AddEntry(GraphGaus[1],"M = 80","lp");
		LegendGaus->AddEntry(GraphGaus[2],"M = 100","lp");
		LegendGaus->AddEntry(GraphGaus[3],"M = 120","lp");
		LegendGaus->AddEntry(GraphGaus[4],"M = 140","lp");
		LegendGaus->AddEntry(GraphGaus[5],"M = 160","lp");
		LegendGaus->AddEntry(GraphGaus[6],"M = 180","lp");
		LegendGaus->AddEntry(GraphGaus[7],"M = 200","lp");
		LegendGaus->AddEntry(GraphGaus[8],"M = 220","lp");
		LegendGaus->AddEntry(GraphGaus[9],"M = 240","lp");
		LegendGaus->AddEntry(GraphGaus[10],"M = 260","lp");
		LegendGaus->AddEntry(GraphGaus[11],"M = 280","lp");
		LegendGaus->AddEntry(GraphGaus[12],"M = 300","lp");
		
		LegendGaus->Draw();
		
		CanvasGaus->Print("pics/GausFilter.png","png");
		//CanvasGaus->Print("pics/GausFilter.pdf","pdf");					//not working

}
