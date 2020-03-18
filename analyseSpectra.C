// analyse MAPMT spectra

#include "tools.C"

// OPTIONS:
Bool_t viewPdf = 1; // view pdf after execution
Double_t muMaxPlot = 0.18; // if nonzero, override "muMax" for setting
                           // the plot scales below
///////////

// global vars
Int_t pedBin;
Double_t pedPeak;
Double_t pedADC;
Bool_t first = 1;
Double_t findThreshold(TH1I * spec);
void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage=0);
void formatGraphs(TGraph ** gr);

// MAIN
void analyseSpectra(
  TString infileN="datadir/run_000223.bin.hist.root",
  Bool_t loopMode = 0
  ) {

  // open file
  TFile * infile = new TFile(infileN,"READ");

  // define table file
  TString tableFile = infileN;
  tableFile(TRegexp(".bin.hist.root$")) = ".table.dat";
  gSystem->RedirectOutput(tableFile,"w");
  gSystem->RedirectOutput(0);

  // define canvas and pdf name
  TCanvas * canvSpec = new TCanvas("canvSpec","canvSpec",800,800);
  canvSpec->SetLogy();
  TString rootN = infileN;
  rootN(TRegexp(".bin.hist.root$")) = "";
  TString pdfSpecN = rootN + ".spectra.pdf";
  TString pdfPlotN = rootN + ".plots.pdf";

  // define threshold vars
  Double_t threshold;
  Int_t thresholdBin,maxBin;
  TLine * thresholdLine = new TLine();
  thresholdLine->SetLineColor(kBlue);
  thresholdLine->SetLineWidth(3);


  // define numEvents vars
  Double_t numEvents,mu,delta;
  Double_t muMax = 0;
  TLatex * numEventsTex = new TLatex();
  numEventsTex->SetNDC(true);
  TString numEventsStr;

  TMultiGraph * numEventsMgr = new TMultiGraph();
  TMultiGraph * muMgr = new TMultiGraph();
  TMultiGraph * deltaMgr = new TMultiGraph();
  TGraph * numEventsGr[3];
  TGraph * muGr[3];
  TGraph * deltaGr[3];
  Int_t grCnt[3];
  for(int p=0; p<3; p++) {
    numEventsGr[p] = new TGraph();
    muGr[p] = new TGraph();
    deltaGr[p] = new TGraph();
    numEventsMgr->Add(numEventsGr[p]);
    muMgr->Add(muGr[p]);
    deltaMgr->Add(deltaGr[p]);
    grCnt[p] = 0;
  };
  numEventsMgr->SetTitle(
    "number of events above threshold vs. MAROC channel;MAROC channel;numEvents");
  muMgr->SetTitle("#mu vs. MAROC channel;MAROC channel;#mu");
  deltaMgr->SetTitle("#delta vs. MAROC channel;MAROC channel;#delta");
  formatGraphs(numEventsGr);
  formatGraphs(muGr);
  formatGraphs(deltaGr);

  TH2D * muPix[3]; // [pmt]
  TString muPixN,muPixT;
  for(int p=0; p<3; p++) {
    muPixN = Form("muPix%d",p+1);
    muPixT = Form("#mu for each pixel of PMT %d",p+1);
    muPix[p] = new TH2D(muPixN,muPixT,8,0,8,8,0,8);
  };
  Int_t pmt,pix;


  // define spectrum vars
  TH1I * spec;
  TString specN;
  Int_t chan;
  TKey * key;
  TIter next(gDirectory->GetListOfKeys());

  // loop over spectra
  while((key=(TKey*)next())) {
    specN = TString(key->GetName());
    if(specN.Contains(TRegexp("^hspe")) && !strcmp(key->GetClassName(),"TH1I")) {
      spec = (TH1I*) key->ReadObj();

      // get channel, pmt, and pixel numbers
      sscanf(specN,"hspe%d",&chan);
      pmt = chan2pmt(chan);
      pix = chan2pix(chan);

      // find threshold
      threshold = findThreshold(spec);
      thresholdBin = spec->FindBin(threshold);

      // get number of events above threshold
      maxBin = spec->GetNbinsX();
      numEvents = spec->Integral(thresholdBin,maxBin);
      numEventsStr = Form("numEv = %.0f",numEvents);
      numEventsTex->SetText(0.5,0.2,numEventsStr);
      numEventsGr[pmt]->SetPoint(grCnt[pmt],chan,numEvents);

      // calculate delta
      delta = threshold - pedADC;
      deltaGr[pmt]->SetPoint(grCnt[pmt],chan,delta);

      // calculate mu
      mu = numEvents<spec->GetEntries() ? -TMath::Log(1-numEvents/spec->GetEntries()) : 0;
      muMax = mu>muMax ? mu:muMax;
      muGr[pmt]->SetPoint(grCnt[pmt],chan,mu);
      muPix[pmt]->Fill(xPix(pix),yPix(pix),mu);

      // increment graph points counter
      grCnt[pmt]++;

      // output to data table
      gSystem->RedirectOutput(tableFile,"a");
      printf("%d %f %f\n",chan,mu,delta);
      gSystem->RedirectOutput(0);

      
      // draw spectrum
      spec->Draw();
      spec->GetXaxis()->SetRangeUser(pedADC-100,pedADC+500); // zoom range
      thresholdLine->DrawLine(threshold,0,threshold,pedPeak);
      numEventsTex->Draw();
      printCanv(canvSpec,pdfSpecN);

    };
  };
  canvSpec->Clear(); printCanv(canvSpec,pdfSpecN,1);


  // start plots canvas
  TCanvas * canvPlot = new TCanvas("canvPlot","canvPlot",1000,1200);
  first = true;
  canvPlot->Divide(2,3);
  Int_t pad;


  // generate textbox from config log file
  pad = 1; canvPlot->cd(pad);
  TString infileLog = infileN;
  TString tmpFile = infileN;
  infileLog(TRegexp(".bin.hist.root$")) = ".log";
  tmpFile(TRegexp(".bin.hist.root$")) = ".tmp";
  TString parseCmd = ".! ./parseConfig.sh " + infileLog + " > " + tmpFile;
  gROOT->ProcessLine(parseCmd);
  TPaveText * configText = new TPaveText(0.05,0.05,0.95,0.95,"NDC");
  configText->ReadFile(tmpFile);
  gROOT->ProcessLine(TString(".! rm "+tmpFile));
  configText->Draw();

  // draw plots vs. MAROC channel
  /*
  pad = 3; canvPlot->cd(pad);
  canvPlot->GetPad(pad)->SetGrid(1,1);
  numEventsMgr->Draw("AP");
  */
  pad = 3; canvPlot->cd(pad);
  canvPlot->GetPad(pad)->SetGrid(1,1);
  muMgr->Draw("AP");
  muMgr->GetYaxis()->SetRangeUser(0,muMaxPlot>0?muMaxPlot:muMax*1.1);

  pad = 5; canvPlot->cd(pad);
  canvPlot->GetPad(pad)->SetGrid(1,1);
  deltaMgr->Draw("AP");
  //deltaMgr->GetYaxis()->SetRangeUser(0,150);

  // draw muPix
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2g");
  for(int p=0; p<3; p++) { 
    pad = (p+1)*2;
    canvPlot->cd(pad);
    canvPlot->GetPad(pad)->SetGrid(0,0);
    canvPlot->GetPad(pad)->SetLogx(0);
    canvPlot->GetPad(pad)->SetLogy(0);
    canvPlot->GetPad(pad)->SetLogz(0);
    muPix[p]->SetMarkerSize(0.7);
    muPix[p]->SetMinimum(0);
    muPix[p]->SetMaximum(muMaxPlot>0?muMaxPlot:muMax);
    muPix[p]->Draw("colztext"); 
  };
  printCanv(canvPlot,pdfPlotN,1);


  if(viewPdf && !loopMode) {
    gROOT->ProcessLine(TString(".! zathura "+pdfSpecN+" "+pdfPlotN));
  };
};


// method to find threshold
Double_t findThreshold(TH1I * spec) {

  // locate pedestal
  pedBin = spec->GetMaximumBin();
  pedPeak = spec->GetBinContent(pedBin);
  pedADC = spec->GetBinCenter(pedBin);
  //printf("pedestal: bin=%d ADC=%f peak=%f\n",pedBin,pedADC,pedPeak);
  

  ///////////////////////////////
  // ALGORITHM TUNE PARAMS:
  // -- stencil spacing for 5-point numerical derivative
  const Int_t h=3;
  // -- number of points to average in the moving average
  const Int_t np=10; 
  // -- the average must be greater than slopeCut
  const Double_t slopeCut=-30;
  // -- height = difference between first and last point over average;
  //    this height must be smaller than heightCut
  //const Double_t heightCut=300; 
  // -- push the threshold upward by this many ADC counts, to ensure we are away from
  //    crosstalk region
  const Double_t buffer=10; 
  ///////////////////////////////


  // evaluate numerical derivative with 5-stencil, with point spacing `h`
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

  // search for threshold, using a moving average of `np` points of the derivative
  Double_t ave;
  Double_t height;
  Double_t adcLev;

  for(int p=0; p<specDeriv->GetN(); p++) {
    ave = 0;
    for(int q=0; q<np; q++) {
      specDeriv->GetPoint(p+q,adc,deriv);
      ave += deriv/np;
      if(q==0) height=deriv;
      if(q+1==np) height-=deriv;
    };
    if(ave>slopeCut /*&& height<heightCut*/) {
      specDeriv->GetPoint(p,adcLev,deriv);
      break;
    };
  };
  //printf("adcLev=%.0f  numEntries=%.0f\n",adcLev,spec->GetEntries());
  if(specDeriv) delete specDeriv;
  return adcLev + buffer;
};

// print canvas to a pdf page
void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage) {
  TString suffix = first ? "(":"";
  if(lastPage) suffix = ")";
  if(first && lastPage) suffix = "";
  canvas->Print(TString(pdf+suffix),"pdf");
  first = 0;
};

// format graph colors and styles
void formatGraphs(TGraph ** gr) {
  for(int g=0; g<3; g++) {
    gr[g]->SetMarkerStyle(kFullCircle);
  };
  gr[0]->SetMarkerColor(kRed);
  gr[1]->SetMarkerColor(kGreen+1);
  gr[2]->SetMarkerColor(kBlue);
};
