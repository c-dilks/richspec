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
TH1D * aveMuProj[nPMT][2];
TH1D * devMuProj[nPMT][2];
enum xy_enum {eX,eY};
TH1D * muHist;
TCanvas * canv;
TString datadir="datadir";

void analyse(Int_t filterNum=0);
Bool_t checkFilter(Int_t filter, Int_t px);
TCanvas * drawMuVsXY(Int_t chanReq=16, Int_t slice=2, Int_t sliceVal=0);


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
  canv = new TCanvas("canv","canv",2000,2000);
  canv->Divide(nPMT+1,6);

  // call analyse
  for(int f=0; f<=filterMax; f++) analyse(f);

  // print extra plots for mu vs. (x,y), or slices in x or y
  TString mupdfName = datadir+"/alignmentMuVsXY.pdf";
  TString mupdfNameL = mupdfName+"(";
  TString mupdfNameR = mupdfName+")";
  Int_t chanRequested = 16;
  drawMuVsXY(chanRequested)->Print(mupdfNameL,"pdf");
  drawMuVsXY(chanRequested,kY,30)->Print(mupdfName,"pdf"); // for 726...
  //drawMuVsXY(chanRequested,kY,110)->Print(mupdfName,"pdf"); // for 753...
  drawMuVsXY(chanRequested,kX,155)->Print(mupdfName,"pdf");
  drawMuVsXY(chanRequested,kX,210)->Print(mupdfName,"pdf");
  drawMuVsXY(chanRequested,kX,270)->Print(mupdfNameR,"pdf");
  gROOT->ProcessLine(TString(".! ./renameFile.sh "+mupdfName));
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

  // projections
  TString aveMuT,devMuT;
  for(int p=0; p<nPMT; p++) {
    aveMuT = aveMu[p]->GetTitle();
    devMuT = devMu[p]->GetTitle();
    aveMuProj[p][eX] = aveMu[p]->ProjectionX();
    aveMuProj[p][eX]->SetTitle(TString("X-projection of "+aveMuT));
    aveMuProj[p][eY] = aveMu[p]->ProjectionY();
    aveMuProj[p][eY]->SetTitle(TString("Y-projection of "+aveMuT));
    devMuProj[p][eX] = devMu[p]->ProjectionX();
    devMuProj[p][eX]->SetTitle(TString("X-projection of "+devMuT));
    devMuProj[p][eY] = devMu[p]->ProjectionY();
    devMuProj[p][eY]->SetTitle(TString("Y-projection of "+devMuT));
    for(int c=0; c<2; c++) {
      aveMuProj[p][c]->SetMinimum(0);
      devMuProj[p][c]->SetMinimum(0);
    };
  };
  

  // draw
  for(int p=0; p<nPMT; p++) {
    canv->cd(0*nPMT+p+1); aveMu[p]->Draw("colztext");
    canv->cd(1*nPMT+p+2); aveMuProj[p][eX]->Draw();
    canv->cd(2*nPMT+p+3); aveMuProj[p][eY]->Draw();
    canv->cd(3*nPMT+p+4); devMu[p]->Draw("colztext");
    canv->cd(4*nPMT+p+5); devMuProj[p][eX]->Draw();
    canv->cd(5*nPMT+p+6); devMuProj[p][eY]->Draw();
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



// draws mu vs. laser (x,y) for a particular channel `chanReq`
// - chanReq is assumed to be in [0,63]; all three PMTs corresponding pixels will be drawn
// - if slice=2, draw 2d plot
// - if slice=kX, draw slice x=sliceVal, a 1d plot binned in y
// - if slice=kY, draw slice y=sliceVal, a 1d plot binned in x
TCanvas * drawMuVsXY(Int_t chanReq=16, Int_t slice=2, Int_t sliceVal=0) {

  TString sfx = Form("%d_%d_%d",chanReq,slice,sliceVal);
  TString canvN = "muXY_"+sfx;
  TCanvas * canvMuXY = new TCanvas(canvN,canvN,2000,500);

  chanReq = chanReq % 64;
  TH1D * d1[3];
  TH2D * d2[3];
  TString dN,dT;
  for(int p=0; p<3; p++) {
    dN = Form("d%d_%s",p,sfx.Data());
    switch(slice) {
      case 2:
        dT = Form("#mu versus laser (x,y) -- PMT%d chan%d;x;y;#mu",p+1,chanReq);
        d2[p] = new TH2D(dN,dT,nb[kX],lb[kX],ub[kX],nb[kY],lb[kY],ub[kY]);
        d2[p]->SetMarkerSize(0.3);
        break;
      case kX:
        dT = Form("#mu vs. laser y for slice x=%d -- PMT%d chan%d;y;#mu",sliceVal,p+1,chanReq); 
        d1[p] = new TH1D(dN,dT,nb[kY],lb[kY],ub[kY]);
        break;
      case kY:
        dT = Form("#mu vs. laser x for slice y=%d -- PMT%d chan%d;x;#mu",sliceVal,p+1,chanReq); 
        d1[p] = new TH1D(dN,dT,nb[kX],lb[kX],ub[kX]);
        break;
      default:
        fprintf(stderr,"ERROR: bad slice value");
        return canvMuXY;
    };
  };

  Int_t bin;
  for(int i=0; i<tr->GetEntries(); i++) {
    tr->GetEntry(i);
    if(chan%64==chanReq) {
      pmt = chan2pmt(chan);
      if(slice==2) {
        bin = d2[pmt]->FindBin(laserPos[kX],laserPos[kY]);
        d2[pmt]->SetBinContent(bin,mu);
      }
      else if(slice==kX) {
        if(laserPos[kX]==sliceVal) {
          bin = d1[pmt]->FindBin(laserPos[kY]);
          d1[pmt]->SetBinContent(bin,mu);
        }
      }
      else if(slice==kY) {
        if(laserPos[kY]==sliceVal) {
          bin = d1[pmt]->FindBin(laserPos[kX]);
          d1[pmt]->SetBinContent(bin,mu);
        }
      }
    };
  };

  canvMuXY->Divide(3,1);
  for(int p=0; p<3; p++) {
    canvMuXY->cd(p+1);
    if(slice==2) {
      d2[p]->SetMinimum(0);
      d2[p]->SetMaximum(muMax);
      d2[p]->Draw("colztext");
    } else {
      canvMuXY->GetPad(p+1)->SetGrid(1,1);
      d1[p]->GetYaxis()->SetRangeUser(0,muMax);
      d1[p]->Draw();
    };
  };

  return canvMuXY;
};
