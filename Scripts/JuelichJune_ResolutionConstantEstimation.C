

void JuelichJune_ResolutionConstantEstimation()
{
	double x[3]={1.175,1.332,1.778};
	double y[3]={4.07,4.19,4.96};
	double y2[3]={pow(4.07,2),pow(4.19,2),pow(4.96,2)};
	
	TGraph *Res = new TGraph(3,x,y);
	//TGraph *Res = new TGraph(2,x+1,y+1);
	TGraph *Res2 = new TGraph(3,x,y2);
	
	TF1 *func = new TF1("func","TMath::Sqrt([0]+x*[1])*2.3548*1000",0,2);
	TF1 *func2 = new TF1("func2","([0]+x * [1])*2.3548*2.3548* 1000000",0,2);
	func->SetParameter(0,3.69E-6);
	//func->SetParLimits(0,2E-6,4E-6);
	//func->FixParameter(0,2.611E-6);
	func->FixParameter(1,4.16193E-7);
	Res->Draw("al");
	
	Res->Fit(func);
	cout <<func->Eval(0)<< endl;
	new TCanvas();
//func->Draw("");

	func2->SetParameter(0,3E-6);
	func2->FixParameter(1,4.16193E-7);
	Res2->Draw("al");
	Res2->Fit(func2);
}
