// draws PMT pixels, with channel numbers printed

#include "tools.C"

void drawPixelMap() {

  gStyle->SetOptStat(0);
  gStyle->SetHistMinimumZero();

  TH2D * map[3]; // [pmt]
  TString mapN,mapT;
  int p;
  for(p=0; p<3; p++) {
    mapN = Form("map%d",p+1);
    mapT = Form("channel numbers of PMT %d",p+1);
    map[p] = new TH2D(mapN,mapT,8,0,8,8,0,8);
    map[p]->SetMinimum(-1);
  };

  Int_t pmt,pix,bin;
  for(int chan=0; chan<192; chan++) {
    pmt = chan2pmt(chan);
    pix = chan2pix(chan);
    bin = map[pmt]->FindBin(xPix(pix),yPix(pix));
    map[pmt]->SetBinContent(bin,chan);
  };

  TCanvas * canv = new TCanvas("canv","canv",1200,400);
  canv->Divide(3,1);
  for(p=0; p<3; p++) {
    canv->GetPad(p+1)->SetGrid(1,1);
    canv->cd(p+1);
    map[p]->Draw("text");
  };
};
