// analyse MAPMT spectra

#include "tools.C"


// OPTIONS:
Bool_t viewPdf = 0; // view pdf after execution
///////////

// global vars
Int_t pedBin;
Double_t pedPeak;
Double_t pedADC;
Bool_t first = 1;
Double_t findThreshold(TH1I * spec);
void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage=0);

// MAIN
void analyseSpectra(
  TString infileN="../data/726_728_729_2020_02_25_06_33/run_000001.bin.hist.root",
  Bool_t loopMode = 0
  ) {

  // open file
  TFile * infile = new TFile(infileN,"READ");

  // load tools
  //gSystem->Load("tools.C");
  //gROOT->LoadMacro("tools.C");

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
  Double_t numEvents,mu;
  Double_t muMax = 0;
  TLatex * numEventsTex = new TLatex();
  numEventsTex->SetNDC(true);
  TString numEventsStr;

  TMultiGraph * numEventsMgr = new TMultiGraph();
  TMultiGraph * muMgr = new TMultiGraph();
  TGraph * numEventsGr[3];
  TGraph * muGr[3];
  Int_t grCnt[3];
  for(int p=0; p<3; p++) {
    numEventsGr[p] = new TGraph();
    muGr[p] = new TGraph();
    numEventsGr[p]->SetMarkerStyle(kFullCircle);
    muGr[p]->SetMarkerStyle(kFullCircle);
    numEventsMgr->Add(numEventsGr[p]);
    muMgr->Add(muGr[p]);
    grCnt[p] = 0;
  };
  numEventsMgr->SetTitle(
    "number of events above threshold vs. MAROC channel;MAROC channel;numEvents");
  muMgr->SetTitle("#mu vs. MAROC channel;MAROC channel;#mu");
  numEventsGr[0]->SetMarkerColor(kRed);
  numEventsGr[1]->SetMarkerColor(kGreen+1);
  numEventsGr[2]->SetMarkerColor(kBlue);
  muGr[0]->SetMarkerColor(kRed);
  muGr[1]->SetMarkerColor(kGreen+1);
  muGr[2]->SetMarkerColor(kBlue);

  TH2D * numEventsPix[3]; // [pmt]
  TString numEventsPixN,numEventsPixT;
  for(int p=0; p<3; p++) {
    numEventsPixN = Form("numEventsPix%d",p+1);
    numEventsPixT = Form("number of events per pixel - PMT %d",p+1);
    numEventsPix[p] = new TH2D(numEventsPixN,numEventsPixT,8,0,8,8,0,8);
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

      // get number of events above threshold, and calculate mu
      maxBin = spec->GetNbinsX();
      numEvents = spec->Integral(thresholdBin,maxBin);
      numEventsStr = Form("numEv = %.0f",numEvents);
      numEventsTex->SetText(0.5,0.2,numEventsStr);
      numEventsGr[pmt]->SetPoint(grCnt[pmt],chan,numEvents);
      mu = numEvents<spec->GetEntries() ? -TMath::Log(1-numEvents/spec->GetEntries()) : 0;
      muMax = mu>muMax ? mu:muMax;
      muGr[pmt]->SetPoint(grCnt[pmt],chan,mu);
      grCnt[pmt]++;
      numEventsPix[pmt]->Fill(xPix(pix),yPix(pix),numEvents);
      
      // draw spectrum
      spec->Draw();
      spec->GetXaxis()->SetRangeUser(pedADC-100,pedADC+500);
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
  printf("parse config file...\n");
  pad = 1; canvPlot->cd(pad);
  TString infileLog = infileN;
  infileLog(TRegexp(".bin.hist.root$")) = ".log";
  TString parseCmd = ".! ./parseConfig.sh " + infileLog + " > tempo";
  gROOT->ProcessLine(parseCmd);
  TPaveText * configText = new TPaveText(0.05,0.05,0.95,0.95,"NDC");
  configText->ReadFile("tempo");
  gROOT->ProcessLine(".! rm tempo");
  configText->Draw();

  // draw numEvents plots
  printf("draw NumEvents...\n");
  pad = 3; canvPlot->cd(pad);
  canvPlot->GetPad(pad)->SetGrid(1,1);
  numEventsMgr->Draw("AP");
  printf("draw mu...\n");
  pad = 5; canvPlot->cd(pad);
  canvPlot->GetPad(pad)->SetGrid(1,1);
  muGr[1]->Print();
  numEventsGr[1]->Print();
  muMgr->Draw("AP");
  muMgr->GetYaxis()->SetRangeUser(0,muMax*1.1);


  // draw numEventsPix
  printf("draw pixels...\n");
  gStyle->SetOptStat(0);
  for(int p=0; p<3; p++) { 
    pad = (p+1)*2;
    canvPlot->cd(pad);
    canvPlot->GetPad(pad)->SetGrid(0,0);
    canvPlot->GetPad(pad)->SetLogx(0);
    canvPlot->GetPad(pad)->SetLogy(0);
    canvPlot->GetPad(pad)->SetLogz(0);
    numEventsPix[p]->Draw("colz"); 
  };
  printCanv(canvPlot,pdfPlotN,1);


  if(viewPdf && !loopMode) gROOT->ProcessLine(TString(".! zathura "+pdfPlotN));
};


// method to find threshold
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
  const Int_t np=3;
  for(int p=0; p<specDeriv->GetN(); p++) {
    ave = 0;
    for(int q=0; q<np; q++) {
      specDeriv->GetPoint(p+q,adc,deriv);
      ave += deriv/np;
    };
    if(ave>-50) {
      specDeriv->GetPoint(p,adcLev,deriv);
      break;
    };
  };
  //printf("adcLev=%.0f  numEntries=%.0f\n",adcLev,spec->GetEntries());
  if(specDeriv) delete specDeriv;
  return adcLev;
};

// print canvas to a pdf page
void printCanv(TCanvas * canvas, TString pdf, Bool_t lastPage) {
  TString suffix = first ? "(":"";
  if(lastPage) suffix = ")";
  if(first && lastPage) suffix = "";
  canvas->Print(TString(pdf+suffix),"pdf");
  first = 0;
};
