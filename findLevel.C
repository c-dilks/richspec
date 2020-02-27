void findLevel(TString infileN="run_000001.bin.hist.root") {

  TFile * infile = new TFile(infileN,"READ");
  TH1I * spec = (TH1I*) infile->Get("hspe100");
  Int_t pedBin;
  Double_t stencil[5];
  Double_t pedPeak,pedAdc;
  Double_t adc,deriv;
  const Int_t h=5;
  Int_t s,cnt;

  TGraph * specDeriv = new TGraph();

  pedBin = spec->GetMaximumBin();
  pedPeak = spec->GetBinContent(pedBin);
  pedAdc = spec->GetBinCenter(pedBin);
  printf("pedestal: bin=%d ADC=%f peak=%f\n",pedBin,pedAdc,pedPeak);
  cnt=0;
  for(int b=pedBin+2*h; b<pedBin+200; b++) {
    adc = spec->GetBinCenter(b);
    // stencil points (s->sb): 0->b-2h  1->b-h  2->b  3->b+h  4->b+2h
    s = 0;
    for(int sb=b-2*h; sb<=b+2*h; sb+=h) stencil[s++] = spec->GetBinContent(sb);
    deriv = (stencil[0] - 8*stencil[1] + 8*stencil[3] - stencil[4]) / (12*h);
    if(deriv<200) specDeriv->SetPoint(cnt++,adc,deriv);
  };


  Double_t ave;
  Double_t adcLev;
  const Int_t np=5;
  for(int p=0; p<specDeriv->GetN(); p++) {
    ave = 0;
    for(int q=0; q<np; q++) {
      specDeriv->GetPoint(p+q,adc,deriv);
      ave += deriv/np;
    };
    if(fabs(ave)<5) {
      specDeriv->GetPoint(p,adcLev,deriv);
      break;
    };
  };

  printf("adcLev=%f\n",adcLev);

  TLine * levLine = new TLine(adcLev,0,adcLev,pedPeak);

  TCanvas * canvDeriv = new TCanvas("canvDeriv","canvDeriv",800,800);
  canvDeriv->Divide(1,2);
  canvDeriv->GetPad(1)->SetLogy();
  canvDeriv->cd(1);
  spec->Draw();
  levLine->Draw();
  canvDeriv->cd(2);
  specDeriv->SetMarkerStyle(kFullCircle);
  specDeriv->SetMarkerColor(kRed);
  specDeriv->Draw("AP");
};



    
  
