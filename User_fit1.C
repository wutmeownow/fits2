#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom2.h"
#include <cmath>
#include "TMath.h"
#include "TROOT.h"

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


void User_fit1(int entries=100000) {

     auto *f=new TFile("datadist.root");
     auto h = (TH1F*) f->Get("h");
     TCanvas *tc1 = new TCanvas("tc1","",800,600);
     tc1->Divide(1,1);
     h->Draw();

     /////////////////////////////////////////////////////
     // Here we simulate some physics process
     /////////////////////////////////////////////////////
     // define a TF1 using 4 parameters
     TF1 *f1 = new TF1("f1",myfunction,0,12,7);
     
     // set the parameter values before fitting
     f1->SetParameters(30,50,8,3);
     // f1->SetParameters(30,50,8,3,-10,6.5,0.2);
     f1->Draw("same");

     // perform the fit here
     // h->Fit("f1");
     h->Fit("f1","EV");

     // retrieve central parameter values and errors
     TF1 *myfunc = h->GetFunction("f1");
     // myfunc->Print();
     myfunc->Draw("same");
     cout << "fit parameters" << endl;
     for (int i=0;i<4;i++) {
     cout << i << "\t" << myfunc->GetParameter(i) 
          << " +- " << myfunc->GetParError(i) << endl;
     }


     // get chi^2 pvalue
     Double_t chi2 = myfunc->GetChisquare();  
     cout << "chi2 / ndf: " << chi2 << " / " << myfunc->GetNDF()
          << " reduced chi2: " << chi2/myfunc->GetNDF() << endl;
     cout << " Pvalue: " 
          << TMath::Prob(chi2,myfunc->GetNDF()) << endl;

     tc1->Update();
     tc1->SaveAs("result.pdf");

     // delete h;
     // delete tc1;
     // delete f1;
     // delete myfunc;
}

