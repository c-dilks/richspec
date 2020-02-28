Int_t pedBin;
Double_t pedPeak;
Double_t pedADC;

Double_t findThreshold(TH1I * spec);

void findLevel(TString infileN="run_000001.bin.hist.root") {

  TFile * infile = new TFile(infileN,"READ");
  Double_t threshold;
  Int_t thresholdBin,maxBin;
  Double_t numEvents;
  TString specName;

  TCanvas * canv = new TCanvas("canv","canv",800,800);
  canv->SetLogy();
  TLine * thresholdLine = new TLine();
  thresholdLine->SetLineColor(kBlue);
  thresholdLine->SetLineWidth(3);

  TLatex * numEventsTex = new TLatex();
  numEventsTex->SetNDC(true);
  TString numEventsStr;
  TGraph * numEventsGr = new TGraph();
  numEventsGr->SetMarkerColor(kRed);
  numEventsGr->SetMarkerStyle(kFullCircle);
  numEventsGr->SetTitle("number of events above threshold vs. channel;channel;numEvents");
  Int_t numEventsGrCnt = 0;

  TH1I * spec;
  Int_t chan;
  TKey * key;
  TIter next(gDirectory->GetListOfKeys());
  Bool_t first = true;
  Int_t limit = 0;
  TString suffix;
  while((key=(TKey*)next())) {
    //if(limit++>10) break;
    specName = TString(key->GetName());
    if(specName.Contains(TRegexp("^hspe")) && !strcmp(key->GetClassName(),"TH1I")) {

      // find threshold
      spec = (TH1I*) key->ReadObj();
      threshold = findThreshold(spec);
      thresholdBin = spec->FindBin(threshold);

      // get number of events above threshold
      maxBin = spec->GetNbinsX();
      numEvents = spec->Integral(thresholdBin,maxBin);
      numEventsStr = Form("numEv = %.0f",numEvents);
      numEventsTex->SetText(0.5,0.2,numEventsStr);
      sscanf(spec->GetName(),"hspe%d",&chan);
      numEventsGr->SetPoint(numEventsGrCnt++,chan,numEvents);
      

      // draw
      spec->Draw();
      spec->GetXaxis()->SetRangeUser(pedADC-100,pedADC+500);
      thresholdLine->DrawLine(threshold,0,threshold,pedPeak);
      numEventsTex->Draw();
      suffix = first ? "(":"";
      canv->Print(TString("canv.pdf"+suffix),"pdf");
      first = false;

    };
  };
  canv->SetLogy(0);
  numEventsGr->Draw("AP");
  canv->Print("canv.pdf)","pdf");

};


Double_t findThreshold(TH1I * spec) {

  // locate pedestal
  pedBin = spec->GetMaximumBin();
  pedPeak = spec->GetBinContent(pedBin);
  pedADC = spec->GetBinCenter(pedBin);
  //printf("pedestal: bin=%d ADC=%f peak=%f\n",pedBin,pedADC,pedPeak);

  // evaluate numerical derivative with 5-stencil, with point spacing `h`
  const Int_t h=5;
  Double_t stencil[5];
  Int_t s;
  Double_t adc,deriv;
  TGraph * specDeriv = new TGraph();
  Int_t cnt = 0;
  for(int b=pedBin+2*h; b<pedBin+200; b++) {
    adc = spec->GetBinCenter(b);
    // stencil points (s->sb): 0->b-2h  1->b-h  2->b  3->b+h  4->b+2h
    s = 0;
    for(int sb=b-2*h; sb<=b+2*h; sb+=h) stencil[s++] = spec->GetBinContent(sb);
    deriv = (stencil[0] - 8*stencil[1] + 8*stencil[3] - stencil[4]) / (12*h);
    if(deriv<200) specDeriv->SetPoint(cnt++,adc,deriv);
  };

  // search for threshold, using a moving average of `np` points
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
  //printf("adcLev=%f\n",adcLev);
  if(specDeriv) delete specDeriv;
  return adcLev;
};
