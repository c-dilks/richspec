// make plots for alignment analysis

#include "tools.C"

// maximum filterNum
const Int_t filterMax = 8;
enum xy {kX,kY};

TTree * tr;
Int_t runnum,chan;
Int_t laserPos[2];
Float_t mu,thresh;
const Int_t nPMT = 3;
Float_t muMin,muMax;
Int_t muBins;
Float_t lb[2];
Float_t ub[2];
Int_t nb[2];
Int_t pmt,pix;
int c;
  
TH2D * filterPix;
TH3D * dataMu[nPMT];
TH2D * aveMu[nPMT];
TH2D * devMu[nPMT];
TH1D * muHist;
TCanvas * canv;
TString datadir="datadir";

void analyse(Int_t filterNum=0);
Bool_t checkFilter(Int_t filter, Int_t px);


//////////////////////////////////////////////////////


void alignmentAnalysis() {
  // build table of laser positions and mu values
  gROOT->ProcessLine(".! buildAlignmentTable.sh");
  TString table = datadir+"/alignment.dat";
  tr = new TTree("tr","tr");
  tr->ReadFile(table,"runnum/I:x/I:y/I:chan/I:mu/F:thresh/F");
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("x",&laserPos[kX]);
  tr->SetBranchAddress("y",&laserPos[kY]);
  tr->SetBranchAddress("chan",&chan);
  tr->SetBranchAddress("mu",&mu);
  tr->SetBranchAddress("thresh",&thresh);

  // determine laser position matrix; assumes uniform step size
  Int_t minPos[2];
  Int_t maxPos[2];
  Int_t stepSize[2] = {0,0};
  minPos[kX] = tr->GetMinimum("x");
  minPos[kY] = tr->GetMinimum("y");
  maxPos[kX] = tr->GetMaximum("x");
  maxPos[kY] = tr->GetMaximum("y");
  Int_t laserPosTmp[2];
  for(int i=0; i<tr->GetEntries(); i++) {
    tr->GetEntry(i);
    if(i==0) { for(c=0; c<2; c++) laserPosTmp[c]=laserPos[c]; };
    for(c=0; c<2; c++) {
      if(laserPos[c]!=laserPosTmp[c]) {
        stepSize[c] = abs(laserPos[c]-laserPosTmp[c]);
        laserPosTmp[c] = laserPos[c];
      };
    };
  };
  for(c=0; c<2; c++) {
    printf("--- laser %s range: %d - %d, step=%d\n",
    c==kX?"x":"y",minPos[c],maxPos[c],stepSize[c]);
    // binning:
    if(stepSize[c]>0) {
      lb[c] = minPos[c]-stepSize[c]/2;
      ub[c] = maxPos[c]+stepSize[c]/2;
      nb[c] = (maxPos[c]-minPos[c])/stepSize[c]+1;
    } else {
      lb[c] = minPos[c]-0.5;
      ub[c] = maxPos[c]+0.5;
      nb[c] = 1;
    };
  };

  // initialize histos
  gStyle->SetPaintTextFormat(".2g");
  gStyle->SetOptStat(0);
  TString hN,hT;
  muMax = tr->GetMaximum("mu")*1.1;
  muMin = 0;
  muBins = 100;
  filterPix = new TH2D("filterPix","PIXEL FILTER",8,0,8,8,0,8);
  for(int p=0; p<nPMT; p++) {
    hN = Form("dataMu%d",p+1);
    dataMu[p] = new TH3D(
      hN,hN,nb[kX],lb[kX],ub[kX],nb[kY],lb[kY],ub[kY],muBins,muMin,muMax);
    hN = Form("aveMu%d",p+1);
    hT = Form("average #mu for PMT %d",p+1);
    aveMu[p] = new TH2D(hN,hT,nb[kX],lb[kX],ub[kX],nb[kY],lb[kY],ub[kY]);
    hN = Form("devMu%d",p+1);
    hT = Form("stddev #mu for PMT %d",p+1);
    devMu[p] = new TH2D(hN,hT,nb[kX],lb[kX],ub[kX],nb[kY],lb[kY],ub[kY]);
    //aveMu[p]->SetMinimum(0);
    //aveMu[p]->SetMaximum(muMax*0.7);
    //devMu[p]->SetMinimum(0);
    //devMu[p]->SetMaximum(0.015);
    aveMu[p]->SetMarkerSize(0.6);
    devMu[p]->SetMarkerSize(0.6);
  };
  muHist = new TH1D("muHist","muHist",muBins,muMin,muMax);
  canv = new TCanvas("canv","canv",1800,700);
  canv->Divide(nPMT+1,2);

  // call analyse
  for(int f=0; f<=filterMax; f++) analyse(f);
};


//////////////////////////////////////////////////////


void analyse(Int_t filterNum=0) {

  // clear histograms
  filterPix->Reset();
  for(int p=0; p<nPMT; p++) {
    dataMu[p]->Reset();
    aveMu[p]->Reset();
    devMu[p]->Reset();
  };
  muHist->Reset();

  // fill filter pix diagram
  for(int p=1; p<=64; p++) {
    if(checkFilter(filterNum,p)) filterPix->Fill(xPix(p),yPix(p));
  };

  // fill dataMu
  for(int i=0; i<tr->GetEntries(); i++) {
    tr->GetEntry(i);
    pmt = chan2pmt(chan);
    pix = chan2pix(chan);
    if(checkFilter(filterNum,pix)) dataMu[pmt]->Fill(laserPos[kX],laserPos[kY],mu);
  };

  // get average and stddev mu
  Double_t bc;
  for(int p=0; p<nPMT; p++) {
    for(int bx=1; bx<=dataMu[p]->GetNbinsX(); bx++) {
      for(int by=1; by<=dataMu[p]->GetNbinsY(); by++) {
        muHist->Reset();
        for(int bz=1; bz<=dataMu[p]->GetNbinsZ(); bz++) {
          bc = dataMu[p]->GetBinContent(bx,by,bz);
          muHist->SetBinContent(bz,bc);
        };
        aveMu[p]->SetBinContent(bx,by,muHist->GetMean());
        devMu[p]->SetBinContent(bx,by,muHist->GetRMS());
      };
    };
  };

  // draw
  for(int p=0; p<nPMT; p++) {
    canv->cd(p+1);
    aveMu[p]->Draw("colztext");
    canv->cd(nPMT+p+2);
    devMu[p]->Draw("colztext");
  };
  canv->cd(nPMT+1); 
  canv->GetPad(nPMT+1)->SetGrid(1,1);
  filterPix->Draw("col");
  TString pdfName = datadir+"/alignmentPlots.pdf";
  TString suffix = "";
  if(filterNum==0) suffix="(";
  if(filterNum==filterMax) suffix=")";
  canv->Print(TString(pdfName+suffix),"pdf");

  if(filterNum==filterMax) gROOT->ProcessLine(TString(".! ./renameFile.sh "+pdfName));
};


///////////////////////////////////////////


Bool_t checkFilter(Int_t filter, Int_t px) {
  switch(filter) {
    case 0: // all pixels
      return true;
      break;
    case 1: // upper half
      return px <= 32;
      break;
    case 2: // lower half
      return px > 32;
      break;
    case 3: // left half
      return (px-1)%8 < 4;
      break;
    case 4: // right half
      return (px-1)%8 >= 4;
      break;
    case 5: // upper-left quadrant
      return checkFilter(1,px) && checkFilter(3,px);
      break;
    case 6: // upper-right quadrant
      return checkFilter(1,px) && checkFilter(4,px);
      break;
    case 7: // lower-left quadrant
      return checkFilter(2,px) && checkFilter(3,px);
      break;
    case 8: // lower-right quadrant
      return checkFilter(2,px) && checkFilter(4,px);
      break;
    default: return false;
  };
};
