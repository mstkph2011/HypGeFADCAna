void HypGeAnaCreateHistogramsFromTree(TString TreeInputFilename ="/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/June2014/TreeCOSYJune2014Dataset11_200,100,0,5339_SR0.root", 
									  TString OutputFilename = "Treetest.root", 
									  TString ParameterFilename = "/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSYJune2014Dataset11_200,100,0,5339_SR0.root")
{

TF1*funcCorr1_T10DeriMaxEnergyNorm_2;
	TF1*funcCorr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2;
	TF1*funcCorr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2;
TF1*funcCorr1_T1090EnergyNorm_2;
	TF1*funcCorr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2;
	TF1*funcCorr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2;
TF1*funcCorr1_TDeriMax90EnergyNorm_2;
	TF1*funcCorr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2;
	TF1*funcCorr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2;
//get peak range(s) and correction functions from file
if(ParameterFilename.Contains(".root"))
{
	TFile *ParameterFile= new TFile(ParameterFilename.Data());
	funcCorr1_T10DeriMaxEnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr1_T10DeriMaxEnergyNorm_2");
		funcCorr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2");
		funcCorr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2");
	funcCorr1_T1090EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr1_T1090EnergyNorm_2");
		funcCorr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2");
		funcCorr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2");

	funcCorr1_TDeriMax90EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr1_TDeriMax90EnergyNorm_2");
		funcCorr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2");
		funcCorr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2=(TF1*) ParameterFile->Get("funcCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2");
	

	ParameterFile->Close();
}

cout << funcCorr1_T10DeriMaxEnergyNorm_2->GetParameter(0) << endl;
cout << funcCorr1_T1090EnergyNorm_2->GetParameter(0) << endl;
cout << funcCorr1_TDeriMax90EnergyNorm_2->GetParameter(0) << endl;

//init histograms


TH1F * hEnergy ;
TH1F * hT1090 ;
TH2F * hT1090Energy ;
TH2F * hT3090Energy ;
TH2F * hTDeriMax90Energy ;
TH2F * hTDeriMax90EnergyRel ;
TH2F * hT10DeriMaxEnergy ;
TH2F * hT10DeriMaxEnergyRel ;

TH1D * hEnergyCorr1_T10DeriMaxEnergyNorm_2 = new TH1D("hEnergyCorr1_T10DeriMaxEnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
TH2D * hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hT3090EnergyCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT3090EnergyCorr1_T10DeriMaxEnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyRelCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
TH2D * hT10DeriMaxEnergyCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hT10DeriMaxEnergyRelCorr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH1D("hEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH1D("hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);


TH1D * hEnergyCorr1_T1090EnergyNorm_2 = new TH1D("hEnergyCorr1_T1090EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
TH2D * hT1090EnergyCorr1_T1090EnergyNorm_2 = new TH2D("hT1090EnergyCorr1_T1090EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hT3090EnergyCorr1_T1090EnergyNorm_2 = new TH2D("hT3090EnergyCorr1_T1090EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyCorr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr1_T1090EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyRelCorr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr1_T1090EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
TH2D * hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hT10DeriMaxEnergyRelCorr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr1_T1090EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH1D("hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH1D("hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);


TH1D * hEnergyCorr1_TDeriMax90EnergyNorm_2 = new TH1D("hEnergyCorr1_TDeriMax90EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
TH2D * hT1090EnergyCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hT1090EnergyCorr1_TDeriMax90EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hT3090EnergyCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hT3090EnergyCorr1_TDeriMax90EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr1_TDeriMax90EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hTDeriMax90EnergyRelCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr1_TDeriMax90EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
TH2D * hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
TH2D * hT10DeriMaxEnergyRelCorr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr1_TDeriMax90EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH1D("hEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);

	TH1D * hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH1D("hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";Energy [ADC value];Counts",32000,1,16000);
	TH2D * hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{1090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{3090} [ns];Energy [ADC value]",100,10,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{Derimax90} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{Derimax90,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{10Derimax} [ns];Energy [ADC value]",130,-290,1000,4000,1,4000);
	TH2D * hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 = new TH2D("hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2",";t_{10Derimax,Rel} [ns];Energy [ADC value]",130,-.290,1.000,4000,1,4000);



//get the Tree to fill the histograms

TTree *tDataTree;
TFile *TreeInputFile= new TFile(TreeInputFilename.Data());
	tDataTree=(TTree*)TreeInputFile->Get("DataTree");
	//create energy spectrum
	tDataTree->Draw("Energy>>hEnergy(32000,0,16000)");
		hEnergy =(TH1F*) gDirectory->Get("hEnergy");
		hEnergy->SetTitle(";Energy [ADC value];Counts");
		hEnergy->SetDirectory(0);
	tDataTree->Draw("T1090>>hT1090(100,0,1000)");
		hT1090 =(TH1F*) gDirectory->Get("hT1090");
		hT1090->SetTitle(";t_{1090} [ns];Counts");
		hT1090->SetDirectory(0);
	tDataTree->Draw("Energy:T1090>>hT1090Energy(100,0,1000,4000,0,4000)","","colz");
		hT1090Energy =(TH2F*) gDirectory->Get("hT1090Energy");
		hT1090Energy->SetTitle(";t_{1090} [ns];Energy [ADC value]");
		hT1090Energy->SetDirectory(0);
	tDataTree->Draw("Energy:T3090>>hT3090Energy(100,0,1000,4000,0,4000)","","colz");
		hT3090Energy =(TH2F*) gDirectory->Get("hT3090Energy");
		hT3090Energy->SetTitle(";t_{3090} [ns];Energy [ADC value]");
		hT3090Energy->SetDirectory(0);
	tDataTree->Draw("Energy:TDerimax90>>hTDeriMax90Energy(130,-300,1000,4000,0,4000)","","colz");
		hTDeriMax90Energy =(TH2F*) gDirectory->Get("hTDeriMax90Energy");
		hTDeriMax90Energy->SetTitle(";t_{Derimax90} [ns];Energy [ADC value]");
		hTDeriMax90Energy->SetDirectory(0);
	tDataTree->Draw("Energy:TDerimax90Rel>>hTDeriMax90EnergyRel(130,-.30,1.000,4000,0,4000)","","colz");
		hTDeriMax90EnergyRel =(TH2F*) gDirectory->Get("hTDeriMax90EnergyRel");
		hTDeriMax90EnergyRel->SetTitle(";t_{Derimax90,Rel} [ns];Energy [ADC value]");
		hTDeriMax90EnergyRel->SetDirectory(0);
	tDataTree->Draw("Energy:T10Derimax>>hT10DeriMaxEnergy(130,-300,1000,4000,0,4000)","","colz");
		hT10DeriMaxEnergy =(TH2F*) gDirectory->Get("hT10DeriMaxEnergy");
		hT10DeriMaxEnergy->SetTitle(";t_{10Derimax} [ns];Energy [ADC value]");
		hT10DeriMaxEnergy->SetDirectory(0);
	tDataTree->Draw("Energy:T10DerimaxRel>>hT10DeriMaxEnergyRel(130,-.30,1.000,4000,0,4000)","","colz");
		hT10DeriMaxEnergyRel =(TH2F*) gDirectory->Get("hT10DeriMaxEnergyRel");
		hT10DeriMaxEnergyRel->SetTitle(";t_{10Derimax,Rel} [ns];Energy [ADC value]");
		hT10DeriMaxEnergyRel->SetDirectory(0);
		


//cycle the tree and do corrections
if(ParameterFilename.Contains(".root"))
{
	

	//tDataTree->Print();
	//cout<<tDataTree<<endl;
	double Energy=0;
	double T1090=0;
	double T3090=0;
	double TDeriMax90=0;
	double TDeriMax90Rel=0;
	double T10DeriMax=0;
	double T10DeriMaxRel=0;
	tDataTree->SetBranchAddress("Energy",&Energy);
	tDataTree->SetBranchAddress("T1090",&T1090);
	tDataTree->SetBranchAddress("T3090",&T3090);
	tDataTree->SetBranchAddress("TDerimax90",&TDeriMax90);
	tDataTree->SetBranchAddress("TDerimax90Rel",&TDeriMax90Rel);
	tDataTree->SetBranchAddress("T10Derimax",&T10DeriMax);
	tDataTree->SetBranchAddress("T10DerimaxRel",&T10DeriMaxRel);
	
	
	double ECorr1_T10DeriMax_2;
		double ECorr2_T10DeriMax_2_T1090_2;
		double ECorr2_T10DeriMax_2_TDeriMax90_2;
	double ECorr1_T1090_2;
		double ECorr2_T1090_2_T10DeriMax_2;
		double ECorr2_T1090_2_TDeriMax90_2;
	double ECorr1_TDeriMax90_2;
		double ECorr2_TDeriMax90_2_T10DeriMax_2;
		double ECorr2_TDeriMax90_2_T1090_2;
	
	for(int i = 0;i<tDataTree->GetEntriesFast();i++)
	//for(int i = 0;i<100;i++)
	
	{
		 ECorr1_T10DeriMax_2=-1;
		 ECorr2_T10DeriMax_2_T1090_2=-1;
		 ECorr2_T10DeriMax_2_TDeriMax90_2=-1;
		 ECorr1_T1090_2=-1;
		 ECorr2_T1090_2_T10DeriMax_2=-1;
		 ECorr2_T1090_2_TDeriMax90_2=-1;
		 ECorr1_TDeriMax90_2=-1;
		 ECorr2_TDeriMax90_2_T10DeriMax_2=-1;
		 ECorr2_TDeriMax90_2_T1090_2=-1;
		
		
		tDataTree->GetEntry(i);
		if(funcCorr1_T10DeriMaxEnergyNorm_2)
		{
			ECorr1_T10DeriMax_2=Energy/ funcCorr1_T10DeriMaxEnergyNorm_2->Eval(T10DeriMax);
			if(funcCorr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2)
				ECorr2_T10DeriMax_2_T1090_2=ECorr1_T10DeriMax_2/funcCorr1_T10DeriMaxEnergyNorm_2Corr2_T1090EnergyNorm_2->Eval(T1090);
			if(funcCorr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2)
				ECorr2_T10DeriMax_2_TDeriMax90_2=ECorr1_T10DeriMax_2/funcCorr1_T10DeriMaxEnergyNorm_2Corr2_TDeriMax90EnergyNorm_2->Eval(TDeriMax90);
		}
		if(funcCorr1_T1090EnergyNorm_2)
		{
			ECorr1_T1090_2=Energy/ funcCorr1_T1090EnergyNorm_2->Eval(T1090);
			if(funcCorr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2)
				ECorr2_T1090_2_T10DeriMax_2=ECorr1_T1090_2/funcCorr1_T1090EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2->Eval(T10DeriMax);
			if(funcCorr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2)
				ECorr2_T1090_2_TDeriMax90_2=ECorr1_T1090_2/funcCorr1_T1090EnergyNorm_2Corr2_TDeriMax90EnergyNorm_2->Eval(TDeriMax90);
		}
		if(funcCorr1_TDeriMax90EnergyNorm_2)
		{
			ECorr1_TDeriMax90_2=Energy/ funcCorr1_TDeriMax90EnergyNorm_2->Eval(TDeriMax90);
			if(funcCorr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2)
				ECorr2_TDeriMax90_2_T10DeriMax_2=ECorr1_TDeriMax90_2/funcCorr1_TDeriMax90EnergyNorm_2Corr2_T10DeriMaxEnergyNorm_2->Eval(T10DeriMax);
			if(funcCorr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2)
				ECorr2_TDeriMax90_2_T1090_2=ECorr1_TDeriMax90_2/funcCorr1_TDeriMax90EnergyNorm_2Corr2_T1090EnergyNorm_2->Eval(T1090);
			
		}
		//cout << Energy <<"  " << ECorr1_T10DeriMaxNorm_2<< "   " << T10DeriMax<<  endl;
		
		hEnergyCorr1_T10DeriMaxEnergyNorm_2->Fill(ECorr1_T10DeriMax_2);
		hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2->Fill(T1090,ECorr1_T10DeriMax_2);
		hT3090EnergyCorr1_T10DeriMaxEnergyNorm_2->Fill(T3090,ECorr1_T10DeriMax_2);
		hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2->Fill(TDeriMax90,ECorr1_T10DeriMax_2);
		hTDeriMax90EnergyRelCorr1_T10DeriMaxEnergyNorm_2->Fill(TDeriMax90Rel,ECorr1_T10DeriMax_2);
		hT10DeriMaxEnergyCorr1_T10DeriMaxEnergyNorm_2->Fill(T10DeriMax,ECorr1_T10DeriMax_2);
		hT10DeriMaxEnergyRelCorr1_T10DeriMaxEnergyNorm_2->Fill(T10DeriMaxRel,ECorr1_T10DeriMax_2);
		
			hEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(ECorr2_T10DeriMax_2_T1090_2);
			hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(T1090,ECorr2_T10DeriMax_2_T1090_2);
			hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(T3090,ECorr2_T10DeriMax_2_T1090_2);
			hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2->Fill(TDeriMax90,ECorr2_T10DeriMax_2_T1090_2);
			hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(TDeriMax90Rel,ECorr2_T10DeriMax_2_T1090_2);
			hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(T10DeriMax,ECorr2_T10DeriMax_2_T1090_2);
			hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(T10DeriMaxRel,ECorr2_T10DeriMax_2_T1090_2);
		
			hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2->Fill( ECorr2_T10DeriMax_2_TDeriMax90_2);
			hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2->Fill( T1090,ECorr2_T10DeriMax_2_TDeriMax90_2);
			hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(T3090,ECorr2_T10DeriMax_2_TDeriMax90_2);
			hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill(TDeriMax90,ECorr2_T10DeriMax_2_TDeriMax90_2);
			hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2->Fill(TDeriMax90Rel ,ECorr2_T10DeriMax_2_TDeriMax90_2);
			hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill( T10DeriMax,ECorr2_T10DeriMax_2_TDeriMax90_2);
			hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Fill( T10DeriMaxRel,ECorr2_T10DeriMax_2_TDeriMax90_2);
		
		hEnergyCorr1_T1090EnergyNorm_2->Fill(ECorr1_T1090_2);
		hT1090EnergyCorr1_T1090EnergyNorm_2->Fill(T1090,ECorr1_T1090_2);
		hT3090EnergyCorr1_T1090EnergyNorm_2->Fill(T3090,ECorr1_T1090_2);
		hTDeriMax90EnergyCorr1_T1090EnergyNorm_2->Fill(TDeriMax90,ECorr1_T1090_2);
		hTDeriMax90EnergyRelCorr1_T1090EnergyNorm_2->Fill(TDeriMax90Rel,ECorr1_T1090_2);
		hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2->Fill(T10DeriMax,ECorr1_T1090_2);
		hT10DeriMaxEnergyRelCorr1_T1090EnergyNorm_2->Fill(T10DeriMaxRel,ECorr1_T1090_2);

			hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(ECorr2_T1090_2_TDeriMax90_2);
			hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T1090,ECorr2_T1090_2_TDeriMax90_2);
			hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T3090,ECorr2_T1090_2_TDeriMax90_2);
			hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(TDeriMax90,ECorr2_T1090_2_TDeriMax90_2);
			hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(TDeriMax90Rel, ECorr2_T1090_2_TDeriMax90_2);
			hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T10DeriMax,ECorr2_T1090_2_TDeriMax90_2);
			hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T10DeriMaxRel,ECorr2_T1090_2_TDeriMax90_2);

			hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(ECorr2_T1090_2_T10DeriMax_2);
			hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T1090,ECorr2_T1090_2_T10DeriMax_2);
			hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T3090,ECorr2_T1090_2_T10DeriMax_2);
			hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(TDeriMax90,ECorr2_T1090_2_T10DeriMax_2);
			hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(TDeriMax90Rel,ECorr2_T1090_2_T10DeriMax_2);
			hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T10DeriMax,ECorr2_T1090_2_T10DeriMax_2);
			hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Fill(T10DeriMaxRel,ECorr2_T1090_2_T10DeriMax_2);


		hEnergyCorr1_TDeriMax90EnergyNorm_2->Fill(ECorr1_TDeriMax90_2);
		hT1090EnergyCorr1_TDeriMax90EnergyNorm_2->Fill(T1090,ECorr1_TDeriMax90_2);
		hT3090EnergyCorr1_TDeriMax90EnergyNorm_2->Fill(T3090,ECorr1_TDeriMax90_2);
		hTDeriMax90EnergyCorr1_TDeriMax90EnergyNorm_2->Fill(TDeriMax90,ECorr1_TDeriMax90_2);
		hTDeriMax90EnergyRelCorr1_TDeriMax90EnergyNorm_2->Fill(TDeriMax90Rel,ECorr1_TDeriMax90_2);
		hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2->Fill(T10DeriMax,ECorr1_TDeriMax90_2);
		hT10DeriMaxEnergyRelCorr1_TDeriMax90EnergyNorm_2->Fill(T10DeriMaxRel,ECorr1_TDeriMax90_2);
		
			hEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(ECorr2_TDeriMax90_2_T1090_2);
			hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T1090,ECorr2_TDeriMax90_2_T1090_2);
			hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T3090,ECorr2_TDeriMax90_2_T1090_2);
			hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(TDeriMax90,ECorr2_TDeriMax90_2_T1090_2);
			hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2->Fill(TDeriMax90Rel,ECorr2_TDeriMax90_2_T1090_2); 
			hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T10DeriMax,ECorr2_TDeriMax90_2_T1090_2);
			hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T10DeriMaxRel,ECorr2_TDeriMax90_2_T1090_2);

			hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(ECorr2_TDeriMax90_2_T10DeriMax_2);
			hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T1090,ECorr2_TDeriMax90_2_T10DeriMax_2);
			hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T3090,ECorr2_TDeriMax90_2_T10DeriMax_2);
			hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(TDeriMax90,ECorr2_TDeriMax90_2_T10DeriMax_2);
			hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(TDeriMax90Rel,ECorr2_TDeriMax90_2_T10DeriMax_2);
			hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T10DeriMax,ECorr2_TDeriMax90_2_T10DeriMax_2);
			hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Fill(T10DeriMaxRel,ECorr2_TDeriMax90_2_T10DeriMax_2);
		
		
	}
}
TreeInputFile->Close();

//cout << hEnergy<< endl;
//cout << hT1090<< endl;
//write histograms to file
if(OutputFilename.Contains(".root"))
{
	TFile *OutputFile= new TFile(OutputFilename.Data(),"RECREATE");
		hEnergy->Write(0,TObject::kOverwrite);
		hT1090->Write(0,TObject::kOverwrite);
		hT1090Energy->Write(0,TObject::kOverwrite);
		hT3090Energy->Write(0,TObject::kOverwrite);
		hTDeriMax90Energy->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyRel->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergy->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyRel->Write(0,TObject::kOverwrite);
		
		TDirectory *SubDir= OutputFile->mkdir("Corr1_T10DeriMaxEnergyNorm_2");
		SubDir->cd();
		
		hEnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hT1090EnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hT3090EnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyRelCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyRelCorr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
		
		SubSubDir= SubDir->mkdir("Corr2_T1090EnergyNorm_2");
		SubSubDir->cd();
		
			hEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
		
		SubDir->cd();
		SubSubDir= SubDir->mkdir("Corr2_TDeriMax90EnergyNorm_2");
		SubSubDir->cd();
		
			hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T10DeriMaxEnergyNorm_2 ->Write(0,TObject::kOverwrite);
		
		SubDir= OutputFile->mkdir("Corr1_T1090EnergyNorm_2");
		SubDir->cd();
		hEnergyCorr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
		hT1090EnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT3090EnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyRelCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyRelCorr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		SubSubDir= SubDir->mkdir("Corr2_T10DeriMaxEnergyNorm_2");		
		SubSubDir->cd();
		
			 hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
		
		SubDir->cd();
		SubSubDir= SubDir->mkdir("Corr2_TDeriMax90EnergyNorm_2");
		SubSubDir->cd();
		
			hEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_TDeriMax90EnergyNorm_2Corr1_T1090EnergyNorm_2->Write(0,TObject::kOverwrite);
		
		SubDir= OutputFile->mkdir("Corr1_TDeriMax90EnergyNorm_2");
		SubDir->cd();
		hEnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT1090EnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT3090EnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hTDeriMax90EnergyRelCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		hT10DeriMaxEnergyRelCorr1_TDeriMax90EnergyNorm_2->Write(0,TObject::kOverwrite);
		SubSubDir= SubDir->mkdir("Corr2_T10DeriMaxEnergyNorm_2");
		SubSubDir->cd();
		
			hEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_T10DeriMaxEnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
		
		SubDir->cd();
		SubSubDir= SubDir->mkdir("Corr2_T1090EnergyNorm_2");
		SubSubDir->cd();
		
			hEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT1090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT3090EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hTDeriMax90EnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
			hT10DeriMaxEnergyRelCorr2_T1090EnergyNorm_2Corr1_TDeriMax90EnergyNorm_2 ->Write(0,TObject::kOverwrite);
		
	OutputFile->Close();
}

}
