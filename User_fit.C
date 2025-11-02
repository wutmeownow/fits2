#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom2.h"
#include "TROOT.h"

// User_fit.C
// example of fitting with a user defined function

// A function producing two peaks on top of an exponentialy 
// falling background.  Depends on several parameters.
// Note: this need not be a 1D function
// Generic interface for fcn of n input-values and m parameters
// Functions with this interface may be used to construct a "TFunction"
Double_t myfunction(Double_t *xin, Double_t *par) {
  Double_t x=xin[0];
  Double_t bkgScale=par[0];
  Double_t alpha=par[1];
  Double_t beta=par[2];
  Double_t background = pow(x/beta,-1.0*alpha);
  // gaussian bump 1
  Double_t A1=par[3];
  Double_t mu1=par[4];
  Double_t sig1=par[5];
  Double_t peak1=A1*TMath::Exp(-0.5*(x-mu1)*(x-mu1)/sig1/sig1);
  // gaussian bump 2
  Double_t A2=par[6];
  Double_t mu2=par[7];
  Double_t sig2=par[8];
  Double_t peak2=A2*exp(-0.5*(x-mu2)*(x-mu2)/sig2/sig2);
  return bkgScale*background+peak1+peak2;
}


void User_fit(int entries=100000) {
     //Simple histogram fitting examples

     TFile *f=new TFile("data1.root","recreate");
     TCanvas *tc1 = new TCanvas("tc1","Experiments Results",800,600);
     tc1->Divide(1,1);

     /////////////////////////////////////////////////////
     // Here we simulate some physics process
     /////////////////////////////////////////////////////
     // define a TF1 using 9 parameters
     TF1 *f1 = new TF1("f1",myfunction,300,1000,9);
     // set the parameter values
     f1->SetParameters(1e9,4.7,40,5000,500,2,1200,800,25);
     // fill a histogram with random data using f1 as a pdf
     TH1F *ranHist1 = new TH1F("ranHist1", "Random Histogram",500,300,1000);
     ranHist1->FillRandom("f1",entries);
     ranHist1->Draw("e");

     /////////////////////////////////////////////////////
     
     // Now "pretend" that we don't know the paramaters used to generate
     // the data.

     // all fits begin with initial guesses at the best parameter values
     f1->SetParameters(1.5e6,5,80,200,500,3,50,800,30);
     f1->Draw("same");

     // uncomment this section once you have good fit parameters for examples 
     // of accessing details on the fit results

     //   /*
     
     f1->SetParameters(1e9,4.7,40,5000,500,2,1200,800,25);

     // perform the fit here
     ranHist1->Fit("f1");
     ranHist1->Fit("f1","EV");

     // retrieve central parameter values and errors
     TF1 *myfunc = ranHist1->GetFunction("f1");
     myfunc->Print();
     cout << "fit parameters" << endl;
     for (int i=0;i<9;i++) {
     cout << i << "\t" << myfunc->GetParameter(i) 
          << " +- " << myfunc->GetParError(i) << endl;
     }


     // get chi^2 pvalue
     Double_t chi2 = myfunc->GetChisquare();  
     cout << "chi2 / ndf: " << chi2 << " / " << myfunc->GetNDF()
          << " reduced chi2: " << chi2/myfunc->GetNDF() << endl;
     cout << " Pvalue: " 
          << TMath::Prob(chi2,myfunc->GetNDF()) << endl;
     
     //   */
     tc1->Update();
     tc1->SaveAs("result.pdf");
     f->Write();
}

