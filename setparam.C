
void setparam(Bool_t print)
{
   // set here values -----------
   Float_t P1   = 100;
   Float_t P2   = 200;
   Bool_t  Histo = kTRUE; // enable/disable histogramming
   //-----------------------------
	THypGeParameter *fp;
	fp = (THypGeParameter *) go4->GetObject("Parameters/HypGeParameter","Go4");
	if(fp==0) return;
	fp->frP1    = P1;
	fp->frP2    = P2;
	fp->fbHisto = Histo;
	if(print) fp->Print();
}
