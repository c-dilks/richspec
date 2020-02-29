// options:
Bool_t viewPdf = 1; // view pdf after execution
Bool_t printSpectra = 1; // print ADC spectra to pdf
/////

Int_t pedBin;
Double_t pedPeak;
Double_t pedADC;
Bool_t first = 1;

Double_t findThreshold(TH1I * spec);
void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage=0);

void findLevel(
  TString infileN="../data/726_728_729_2020_02_25_06_33/run_000001.bin.hist.root",
  Bool_t loopMode = 0
  ) {

  TFile * infile = new TFile(infileN,"READ");
  Double_t threshold;
  Int_t thresholdBin,maxBin;
  TString specName;

  TCanvas * canv = new TCanvas("canv","canv",800,800);
  canv->SetLogy();
  TString pdfName = infileN;
  pdfName(TRegexp(".bin.hist.root$")) = ".pdf";
  TLine * thresholdLine = new TLine();
  thresholdLine->SetLineColor(kBlue);
  thresholdLine->SetLineWidth(3);

  Double_t numEvents,mu;
  Double_t muMax = 0;
  TLatex * numEventsTex = new TLatex();
  numEventsTex->SetNDC(true);
  TString numEventsStr;
  TGraph * numEventsGr = new TGraph();
  TGraph * muGr = new TGraph();
  numEventsGr->SetMarkerColor(kBlue); numEventsGr->SetMarkerStyle(kFullCircle);
  numEventsGr->SetTitle("number of events above threshold vs. channel;channel;numEvents");
  muGr->SetMarkerColor(kRed); muGr->SetMarkerStyle(kFullCircle);
  muGr->SetTitle("#mu vs. channel;channel;#mu");
  Int_t grCnt = 0;

  TH1I * spec;
  Int_t chan;
  TKey * key;
  TIter next(gDirectory->GetListOfKeys());
  Int_t limit = 0;
  while((key=(TKey*)next())) {
    //if(limit++>10) break;
    specName = TString(key->GetName());
    if(specName.Contains(TRegexp("^hspe")) && !strcmp(key->GetClassName(),"TH1I")) {
      spec = (TH1I*) key->ReadObj();

      // get channel
      sscanf(spec->GetName(),"hspe%d",&chan);

      // find threshold
      threshold = findThreshold(spec);
      thresholdBin = spec->FindBin(threshold);

      // get number of events above threshold, and calculate mu
      maxBin = spec->GetNbinsX();
      numEvents = spec->Integral(thresholdBin,maxBin);
      numEventsStr = Form("numEv = %.0f",numEvents);
      numEventsTex->SetText(0.5,0.2,numEventsStr);
      numEventsGr->SetPoint(grCnt,chan,numEvents);
      mu = -TMath::Log(1-numEvents/spec->GetEntries());
      muMax = mu>muMax ? mu:muMax;
      muGr->SetPoint(grCnt,chan,mu);
      grCnt++;
      

      // draw spectrum
      if(printSpectra) {
        spec->Draw();
        spec->GetXaxis()->SetRangeUser(pedADC-100,pedADC+500);
        thresholdLine->DrawLine(threshold,0,threshold,pedPeak);
        numEventsTex->Draw();
        printCanv(canv,pdfName);
      };

    };
  };

  // draw other graphs
  canv->SetLogy(0);
  canv->SetGrid(1,1);
  //numEventsGr->GetXaxis()->SetRangeUser(-1,193);
  numEventsGr->Draw("AP");
  printCanv(canv,pdfName);
  muGr->Draw("AP");
  muGr->GetYaxis()->SetRangeUser(0,muMax*1.1);
  printCanv(canv,pdfName);

  // print config
  canv->Clear();
  TString infileLog = infileN;
  infileLog(TRegexp(".bin.hist.root$")) = ".log";
  TString parseCmd = ".! ./parseConfig.sh " + infileLog + " > tempo";
  gROOT->ProcessLine(parseCmd);
  TPaveText * configText = new TPaveText(0.05,0.05,0.95,0.95,"NDC");
  configText->ReadFile("tempo");
  gROOT->ProcessLine(".! rm tempo");
  configText->Draw();
  printCanv(canv,pdfName,1);

  if(viewPdf && !loopMode) gROOT->ProcessLine(TString(".! zathura "+pdfName));
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
  //printf("adcLev=%.0f  numEntries=%.0f\n",adcLev,spec->GetEntries());
  if(specDeriv) delete specDeriv;
  return adcLev;
};

void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage) {
  TString suffix = first ? "(":"";
  if(lastPage) suffix = ")";
  canvas->Print(TString(pdf+suffix),"pdf");
  first = 0;
};
