#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom2.h"
#include <cmath>
#include "TMath.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TText.h"

// large gaussian with a constant background and a small gaussian dip
Double_t myfunction2(Double_t *xin, Double_t *par) {
  Double_t x=xin[0];
  Double_t background = par[0];
  // overall gaussian
  Double_t A=par[1];
  Double_t mu=par[2];
  Double_t sig=par[3];
  Double_t peak=A*TMath::Exp(-0.5*(x-mu)*(x-mu)/sig/sig);
  // small gaussian dip
  Double_t A1=par[4];
  Double_t mu1=par[5];
  Double_t sig1=par[6];
  Double_t dip=A1*TMath::Exp(-0.5*(x-mu1)*(x-mu1)/sig1/sig1);
  return background+peak+dip;
}

// large gaussian with a constant background
Double_t myfunction(Double_t *xin, Double_t *par) {
  Double_t x=xin[0];
  Double_t background = par[0];
  // overall gaussian
  Double_t A=par[1];
  Double_t mu=par[2];
  Double_t sig=par[3];
  Double_t peak=A*TMath::Exp(-0.5*(x-mu)*(x-mu)/sig/sig);
  return background+peak;
}


// large gaussian with a linear background
Double_t myfunction3(Double_t *xin, Double_t *par) {
  Double_t x=xin[0];
  Double_t background = par[0];
  Double_t slope = par[1];
  // overall gaussian
  Double_t A=par[2];
  Double_t mu=par[3];
  Double_t sig=par[4];
  Double_t peak=A*TMath::Exp(-0.5*(x-mu)*(x-mu)/sig/sig);
  return background+peak+slope*x;
}




void User_fit1(int entries=100000) {

     auto *f=new TFile("datadist.root");
     auto h = (TH1F*) f->Get("h");
     TCanvas *tc1 = new TCanvas("tc1","",1200,600);
     tc1->Divide(3,1);
     tc1->cd(1);

     // define graph and histo for residuals
     TGraph *gresiduals = new TGraph();
     gresiduals->SetTitle("residuals;;");
     TH1F *hresiduals = new TH1F("hresiduals", "Pull Distribution", 15, -4, 4);

     /////////////////////////////////////////////////////
     // Here we simulate some physics process
     /////////////////////////////////////////////////////
     // define a TF1 using 4 parameters
     // int Npar = 4;
     // TF1 *f1 = new TF1("f1",myfunction,0,12,4);
     // 5 param fit
     int Npar = 5;
     TF1 *f1 = new TF1("f1",myfunction3,0,12,Npar);
     
     // set the parameter values before fitting
     // f1->SetParameters(30,50,8,3); // initial pars for constant bkg + gaussian
     f1->SetParameters(30,1,50,8,3); // initial pars for linear bkg + gaussian
     // f1->SetParameters(30,50,8,3,-10,6.5,0.2);
     // f1->Draw("same");

     // perform the fit here
     h->Fit("f1");
     // h->Fit("f1","EV");

     // retrieve central parameter values and errors
     TF1 *myfunc = h->GetFunction("f1");
     // myfunc->Print();
     
     cout << "fit parameters" << endl;
     for (int i=0;i<Npar;i++) {
     cout << i << "\t" << myfunc->GetParameter(i) 
          << " +- " << myfunc->GetParError(i) << endl;
     }

     
     // Loop over **all** visible bins and calculate residuals
     for (int i = 1; i <= h->GetNbinsX(); ++i) {
          double binCenter = h->GetBinCenter(i);
          double content   = h->GetBinContent(i);
          double error     = h->GetBinError(i);
          double_t residual = (myfunc->Eval(binCenter)-content)/error;
          gresiduals->SetPoint(i-1, binCenter, residual); // add residual to g
          hresiduals->Fill(residual);
     }

     // fit pull distribution with gaussian
     hresiduals->Fit("gaus","Q");
     TF1 *resfitfunc = hresiduals->GetFunction("gaus");


     // get chi^2 pvalue
     Double_t chi2 = myfunc->GetChisquare();  
     cout << "chi2 / ndf: " << chi2 << " / " << myfunc->GetNDF()
          << " reduced chi2: " << chi2/myfunc->GetNDF() << endl;
     cout << " Pvalue: " 
          << TMath::Prob(chi2,myfunc->GetNDF()) << endl;


     // plot
     TText text; // text for graphs
     text.SetTextSize(0.03);
     text.SetTextColor(kRed);

     h->Draw();
     myfunc->Draw("same");

     tc1->cd(2);
     gresiduals->Draw("AL");

     tc1->cd(3);
     hresiduals->Draw();
     resfitfunc->Draw("same");
     text.DrawTextNDC(0.15, 0.6, Form("mean = %.2f",resfitfunc->GetParameter(1)));
     text.DrawTextNDC(0.15, 0.55, Form("sigma = %.2f",resfitfunc->GetParameter(2)));
 

     tc1->Update();
     tc1->SaveAs("result2.pdf");

     // delete gresiduals;
     // delete hresiduals;
     // delete tc1;
     // delete f1;
     // delete myfunc;
     // delete resfitfunc;
}

