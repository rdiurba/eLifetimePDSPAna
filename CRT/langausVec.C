//-----------------------------------------------------------------------
//
//	Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//  to execute this example, do:
//  root > .x langaus.C
// or
//  root > .x langaus.C++
//
//-----------------------------------------------------------------------

#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TRandom.h"
#include "/exp/dune/app/users/rdiurba/protoduneana_v09_66_calib/srcs/protoduneana/protoduneana/singlephase/michelremoving/scripts/LanGausFit.C"
using namespace std;
TFile *ef=new TFile("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
TH3F    *dzMeasNeg=(TH3F*)ef->Get("RecoBkwd_Displacement_X_Neg");
TH3F    *dzMeas=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
/*
TFile *ex=new TFile("/pnfs/dune/persistent/users/rdiurba/efield_5runs_graph.root");
TH1D *xLeft=(TH1D*)ex->Get("efield_fit_trunc_mean");
TH1D *xRight=(TH1D*)ex->Get("efield_fit_trunc_mean_neg");
TFile *efACA=new TFile("/dune/app/users/apaudel/MC_sample_jan26/data_newfiles/Efield_correctedz.root");
TH1D *eXACA=(TH1D*)efACA->Get("med_Ef");
*/
float LAr_density=1.39;
float alp=0.93;
float bet=0.212;
float dedx=1.9;
float recom_correction(float totEf){ 
  float xsi=bet*dedx/(LAr_density*totEf);
  float xsi0=bet*dedx/(LAr_density*0.4867);
  float rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}
float calc_Ef(float xval, float Eo, float Ea){
    float alpha=3.6*0.18;
    float term=(alpha*Eo*xval)/(Ea*3.6);
    return Ea*(1+0.5*term*term);
}
float diff_corr(double t, double weight){
    float nominal=TMath::Exp(-t/52);
    float diffOff=TMath::Exp(-t/35);
    float ratio=nominal/diffOff;
    if (weight<0) ratio=1.f/ratio;
    return ratio;

}
float tot_Ef(float xval,float yval,float zval){
    //float E0value=0.5;
    float E0value=0.4867;
    //E0value=0.4993;
   if(xval>=0){
    float ex=E0value+E0value*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
    float ey=E0value*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
    float ez=E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
    float ex=E0value+E0value*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
    float ey=E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
    float ez=E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
    }
 return E0value;
}
/*
TFile *efACA=new TFile("/dune/app/users/apaudel/MC_sample_jan26/data_newfiles/Efield_correctedz.root");
TH1D *eXACA=(TH1D*)efACA->Get("med_Ef");

float tot_EfACAData(float xval,float yval,float zval){
  float E0value=0.4867;

   if(xval>=0){
    float ex=eXACA->Interpolate(xval);
    float ey=E0value*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
    float ez=E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
    float ex=eXACA->Interpolate(xval);
    float ey=E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
    float ez=E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
    }
 return E0value;
}
*/

void langausVec(std::string runNumber) { 
//std::cout<<dzMeas->GetBinContent(dzMeas->FindBin(-500,300,1))<<std::endl;
   gROOT->Reset();
   gROOT->LoadMacro("protoDUNEStyle.C");
   gROOT->SetStyle("protoDUNEStyle");
   gROOT->ForceStyle();
   
   gStyle->SetTitleX(0.3);
   gStyle->SetOptFit(111);

  //TFile my_file0("/pnfs/dune/persistent/users/apaudel/prod2_calibfactors/YZcalo_r5809.root");
 //TFile my_file0("/pnfs/dune/persistent/users/apaudel/prod2_calibfactors/YZcalo_sce.root");
// TFile my_file0("/dune/app/users/apaudel/prod3_MC/YZcalo_r6256030.root");

std::string fileNameYZ=Form("/exp/dune/app/users/rdiurba/eLifetimeCathode/YZcalo_mich1_r%s.root",runNumber.c_str());
if (runNumber=="MC") fileNameYZ="/exp/dune/app/users/rdiurba/eLifetimeCathode/YZcalo_mich1_mc_sceOn.root";
if (runNumber=="MCNoDiff") fileNameYZ="/exp/dune/app/users/rdiurba/eLifetimeCathode/YZcalo_mich1_mc_sceOn_noDiff.root";
 TFile my_file0(Form("%s",fileNameYZ.c_str()));
  TH2F *YZ_correction_neg=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_negativeX_hist_2");
   TH2F *YZ_correction_pos=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_positiveX_hist_2");

 //TFile my_file1("/dune/app/users/rdiurba/eLifetime/Xcalo_MC6GeV.root");
//  TFile my_file0("/dune/app/users/apaudel/prod3_SCEOFF/YZcalo_r31186687.root");
  //TH1F *x_correction_factor=(TH1F*)my_file1.Get("dqdx_X_correction_hist_2");
   
   //TFile *f = new TFile("/dune/data/users/rdiurba/TwoMCCE.root");
   // TFile *f=new TFile("/dune/data/users/rdiurba/TwoCRTE5796v2.root");
//   TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTProd4ELifetime.root");
  // TFile *f = new TFile("/dune/app/users/rdiurba/dunetpc_v08_50/TwoCRTSingleGen_noSCE_0deg.root");
 //  TFile *f = new TFile("/dune/data/users/rdiurba/TwoProd2NoSCE.root");
   //TFile *f = new TFile("/dune/app/users/rdiurba/TwoCRTMatching.root");   
//    TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTE5841v2.root");
   //TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTNoDiff.root");
   //TFile *f = new TFile("/dune/data/users/rdiurba/TwoProd3SCE_10ms.root");
    //TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTNoDiff.root");
  // TFile *f = new TFile("/dune/data/users/rdiurba/TwoProd3CryoFix.root");
   //TFile *f= new TFile("/dune/data/users/rdiurba/TwoProd3NoDiff_v2.root");
   //TFile *f = new TFile("/dune/data/users/rdiurba/TwoProd3NoDiffNoSCE.root");   
    //TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTCryoFix.root");
   // TFile *f = new TFile("/dune/data/users/rdiurba/TwoProd2NoSCE.root");
    //TFile *f = new TFile("/dune/data/users/rdiurba/TwoCRTE11142.root");
   //TFile *f = new TFile("/dune/app/users/rdiurba/TwoCRTMatchingLegacy_CRTest.root");
  
// TFile *f = new TFile("/dune/app/users/rdiurba/TwoCRT5770_v09_09.root");
// TFile *f=new TFile("/dune/data/users/rdiurba/TwoCRTMatchRun5833.root");

//char set[]="Date 11/11/2018";
   //char s[]="5841Redo";
   gStyle->SetOptFit(1); 
   //char s[]="Prod4MCSCE";
   //char set[]="Prod. 4 MC w/ SCE";

 std::string fileName="/dune/data/users/rdiurba/TwoCRTMatchProd4SCE_1-2.root";
 std::string set="Run "+runNumber;
  std::string s="run"+runNumber;
if(runNumber.find(Form("MCNoDiff"))!=std::string::npos){
std::string set="Simulation No Diff";
std::string s="Prod4aMCNoDiff";
fileName="/pnfs/dune/persistent/users/rdiurba/TwoCRTFilesNov2023/TwoCRTProd4aNoDiff.root";


}
else if (runNumber.find(Form("MC"))!=std::string::npos){
 std::string set="Simulation";
  std::string s="Prod4MCSCE";
 fileName="/pnfs/dune/persistent/users/rdiurba/TwoCRTFilesNov2023/TwoCRTProd4a_0.root";

}
else{
 
fileName=Form("/pnfs/dune/persistent/users/rdiurba/TwoCRTFilesNov2023/TwoCRTELifetime%s.root",runNumber.c_str());


}
//fileName=Form("/dune/data/users/rdiurba/test5759_100.root");
  TFile *f = new TFile(Form("%s",fileName.c_str())); 
   TTree *t1 = (TTree*)f->Get("ana");
   if (runNumber.find(Form("MC"))!=std::string::npos) t1=(TTree*)f->Get("TwoCRTMatching/CRTdQ");
   std::cout<<t1->GetEntries()<<std::endl;
   std::vector<int>* hit=0x0;
   std::vector<int>* decon=0x0;
   std::vector<int>* wireID=0x0;
   Double_t x1,x2,y1,y2,z1,z2,t,crtt;
   std::vector<double>* trkt=0x0;
   std::vector<double>* trkx=0x0;
   std::vector<double>* trky=0x0;
   std::vector<double>* trkz=0x0;
   Double_t energyFrac;
   Int_t hourMinSec, yearMonthDay;
   std::string stampDayString, year, month, day;
   t1->SetBranchAddress("trkhitIntegral",&hit);
   //t1->SetBranchAddress("hsumADC",&hit);
   t1->SetBranchAddress("crtt0",&crtt);
   t1->SetBranchAddress("trkhitt0",&trkt);
   t1->SetBranchAddress("hX_F",&x1);
   t1->SetBranchAddress("hX_B",&x2);
   t1->SetBranchAddress("hY_F",&y1);
   t1->SetBranchAddress("hY_B",&y2);
   t1->SetBranchAddress("hZ_F",&z1);
   t1->SetBranchAddress("hZ_B",&z2);
   t1->SetBranchAddress("hWireID",&wireID);
   t1->SetBranchAddress("trkhity",&trky);
   t1->SetBranchAddress("trkhitz",&trkz);
   t1->SetBranchAddress("trkhitx",&trkx);
   t1->SetBranchAddress("totEnergyFrac",&energyFrac); 
   t1->SetBranchAddress("hdeconChrg",&decon); 
   t1->SetBranchAddress("hourMinSec",&hourMinSec);
   t1->SetBranchAddress("yearMonthDay",&yearMonthDay);
 
  TVector3 vec;
  // Fill Histogram
   //TFile file("hist_all.root");
  // const int tbinsize=150000;
  const int tbinsize = 100000;
  //const int tbinsize = 50000;
	  //const int ntbins = 2200000/tbinsize;  
	  const int ntbins = 2000000/tbinsize;
	  TCanvas *c1 = new TCanvas("c1","c1", 200, 10, 700, 500 );
	  TCanvas *c2 = new TCanvas("c2","c2", 200, 10, 700, 500 );
	  c1->Print("fit.ps[");
	  c2->Print("fitLand.ps[");

	  vector<double> time;
	  vector<double> etime;
	  vector<double> charge;
	  vector<double> echarge;
	  

	  double tt[ntbins];
	  double ett[ntbins];
	  double imp[ntbins];
	  double eimp[ntbins];
	  double dqdx0[ntbins];
	  double edqdx0[ntbins];
	  std::vector<double> deltaQ;
	  std::vector<double> startX;
	  std::vector<double> startY;
	  std::vector<double> startZ;
	  std::vector<double> endX;
	  std::vector<double> endY;
	  std::vector<double> endZ;
	  std::vector<double> crtXPos;
	  std::vector<double> crtYPos;
	  std::vector<double> hitz;
	  std::vector<double> wireVecID;
	  Int_t stampDay, stampTime;
	  
          int n = 0;
	  double trkzCorr=0,dZOffset=0,dZOffset2=0;
	    time.clear();
	    etime.clear();
	    charge.clear();
	    echarge.clear();
	    //TH1D *hSNR = new TH1D("dQ","dQ",80,0,200);
	   std::vector<TH1D*> hSNRVec;
	    for (int i =0; i<ntbins; ++i){
	    hSNRVec.push_back(new TH1D(Form("dQ%d",i),Form("dQ%d",i),100,0,200));
	   }
	   for (Long64_t j=0; j<Long64_t(t1->GetEntries()); j++){
	   
	   t1->GetEntry(j);
	   TVector3 vec(x2-x1,y2-y1,z2-z1);
	   TVector3 vhat=vec.Unit();
	   if(j==0){ stampDay=yearMonthDay; stampTime=hourMinSec;
                stampDayString=std::to_string(yearMonthDay);
                year=stampDayString.substr(0,4);
                month=stampDayString.substr(4,2);
                day=stampDayString.substr(6,2);
                //std::cout<<year<<","<<month<<","<<day<<std::endl;
		}
	   
           /*if (stampDay<=yearMonthDay){
	   if(stampDay<yearMonthDay) stampTime=hourMinSec;

	   if(stampTime<hourMinSec) stampTime=hourMinSec;
	   stampDay=yearMonthDay;

	   }*/
	   if(j%1000==0) std::cout<<j<<"/"<<t1->GetEntries()<<std::endl;
	   if (vhat(2)<.99) continue;
	   for (size_t hitIndex=0; hitIndex<trky->size(); hitIndex++){

	   if (trky->at(hitIndex)>375 || trky->at(hitIndex)<225 || trkz->at(hitIndex)>450 || trkz->at(hitIndex)<240) continue;
	   if (trkx->at(hitIndex)<0) continue;
	   int binIndex=int((trkt->at(hitIndex)-crtt)/tbinsize);

	   if (binIndex>19) continue;
	   trkzCorr=dZOffset+trkz->at(hitIndex);
	   double crty=((y2-y1)/(z2-z1))*(trkzCorr-z1)+y1;
	   double crtx=((x2-x1)/(z2-z1))*(trkzCorr-z1)+x1;
	   float eMeas=tot_Ef(crtx,crty,trkz->at(hitIndex));
	   float coeff=recom_correction(eMeas);
	   double YZCorrectionFactor=1.f;
	   if(crtx>0) YZCorrectionFactor=YZ_correction_pos->GetBinContent(YZ_correction_pos->FindBin(crty,trkz->at(hitIndex)));
	   else  YZCorrectionFactor=YZ_correction_neg->GetBinContent(YZ_correction_neg->FindBin(crty,trkz->at(hitIndex)));
	   double deltaQ=((coeff*(YZCorrectionFactor*hit->at(hitIndex)))/(0.479/vhat(2)));
	   if (hitIndex>trky->size()-3) hSNRVec[(binIndex)]->Fill(deltaQ);
	   else{
	   double deltaQ1=((coeff*(YZCorrectionFactor*hit->at(hitIndex+1)))/(0.479/vhat(2)));
	   double deltaQ2=((coeff*(YZCorrectionFactor*hit->at(hitIndex+2)))/(0.479/vhat(2)));
	    
	   if (deltaQ>80 && deltaQ1>80 && deltaQ2>80 ) {    hitIndex=hitIndex+1;}
	   else{ hSNRVec[(binIndex)]->Fill(deltaQ);}
	   }
	   
	   }
	   
	   }
	    for (int i =0; i<ntbins; ++i){	
	      if (hSNRVec[i]->GetEntries()>300){
	      //TF1* fLand= new TF1("f","landau",30,200);
		    // hSNRVec[i]->Fit("f","RS QN");   
               std::cout<<hSNRVec[i]->GetEntries()<<std::endl;	       
	       TF1 *fitsnr = runlangaufit(hSNRVec[i],2,200);
	       if(fitsnr->GetChisquare()/fitsnr->GetNDF()<20){ 

		  c2->cd();
		  //hSNR->GetXaxis()->SetRange(100,700);
		  float startTime=i*tbinsize*0.000001f;
		  float endTime=(i+1)*tbinsize*0.000001f;
		  hSNRVec[i]->SetTitle(Form("CRT-TPC, Date: %s/%s/%s, %.1f<t<%.1f ms",day.c_str(),month.c_str(),year.c_str(),startTime,endTime));
                  if (runNumber=="MC") hSNRVec[i]->SetTitle(Form("CRT-TPC, Sim., %.1f<t<%.1f ms",startTime,endTime));	
                  if (runNumber=="MCNoDiff") hSNRVec[i]->SetTitle(Form("CRT-TPC, Sim. No Diff., %.1f<t<%.1f ms",startTime,endTime));
          	  hSNRVec[i]->SetXTitle("dQ/dx [(ADC count)#timestick/cm]");
		  hSNRVec[i]->SetYTitle("Number of Hits");
		  hSNRVec[i]->GetXaxis()->CenterTitle();
		  hSNRVec[i]->GetYaxis()->CenterTitle();
	      
		  hSNRVec[i]->Draw("P0E");
		  fitsnr->SetLineColor(2);
		  fitsnr->SetLineWidth(2);
		  
                  fitsnr->Draw("l same");
		  
	    gPad->RedrawAxis();

	    gStyle->SetStatX(0.85);
	    gStyle->SetStatY(0.88);
	  
	    TLatex tL;
	    tL.SetNDC();
	    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	     c2->Print(Form("fitLand%s%d.png",s.c_str(),i));    
             c2->Print(Form("fitLand%s%d.pdf",s.c_str(),i));
		   time.push_back((i*tbinsize+tbinsize/2)/1000000.f);
		  etime.push_back(tbinsize/(2000000.f));
		  
		  printf("%d,%f,%f",(i*tbinsize+tbinsize/2),fitsnr->GetParameter(1),fitsnr->GetParError(1));
		  //printf("\n");
		  //printf("%f,%f,%f",fitsnr->GetParError(1),fitsnr->GetChisquare()/fitsnr->GetNDF(),fLand->GetParError(1));
		  //printf("\n");
		   charge.push_back(fitsnr->GetParameter(1));  

		  echarge.push_back(fitsnr->GetParError(1));

		}
	      }
	    }
	      


	   

	   
	   c1->cd();
            TH1D* hist= new TH1D("lanfithist","lanfithist",20,0.0,2.0);
            int iter2=0;
            for (int i=0; i<time.size();i++){
            int position=hist->FindBin(time.at(i));
            if (position<2 || position>20) continue;
            hist->SetBinContent(position,charge.at(i));
            hist->SetBinError(position, echarge.at(i));
            iter2++;


            }
           if (iter2<10) return;
	    TGraphErrors *gr = new TGraphErrors(time.size(),&time[0],&charge[0],&etime[0],&echarge[0]);
	    double min=time[0];
	    double max=time[time.size()-1];
	    TF1 *fit= new TF1("exp","[0]*exp(x/(-[1]))",0.1,2.0);
	    //fit->SetChisquare(1);
	    fit->SetParName(0,"Constant");
	    fit->SetParName(1,"e^{-} Lifetime [ms]");
	    //fit->SetParName(1,"1/tau");i
	    double initial=charge.at(0);
	    double lifetimeEst=(time.at(time.size()-1)-time.at(0))/(TMath::Log(charge.at(0)/charge.at(charge.size()-1)));
	    if (lifetimeEst<0) lifetimeEst=100;
	    std::cout<<"Testing "<<lifetimeEst<<std::endl; 
	    fit->SetParameter(0,initial);
	    fit->SetParameter(1,lifetimeEst);
	    hist->SetTitle(Form("%s",set.c_str()));

	    //gr->SetTitle(Form("%s - %s",t0.AsString("lc"),t1.AsString("lc")));
	    hist->Draw("E P");
	    //gr->SetTitle(Form("dQ/dx and Hit Time %s",set));
	    hist->GetXaxis()->SetTitle("Drift Time [ms]");
	    hist->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
            hist->SetTitle(Form("CRT-TPC, Date: %s/%s/%s",day.c_str(),month.c_str(),year.c_str()));
	    if (runNumber=="MC") hist->SetTitle("Simulation");
            if (runNumber=="MCNoDiff") hist->SetTitle("Simulation with No Diffusion");
            hist->GetYaxis()->SetRangeUser(0,100);
	    //gr->GetYaxis()->SetRangeUser(150,500);
	    hist->GetXaxis()->CenterTitle();
	    hist->GetYaxis()->CenterTitle();
	    

	    gPad->RedrawAxis();
	    gStyle->SetStatX(0.85);
	    gStyle->SetStatY(0.88);
        TH1D* histCopy=(TH1D*)hist->Clone("histCopy");
double chi2_dof=1000000;
int iter=0;
        while ((chi2_dof>1.05 || chi2_dof<0.95) && iter<100000){
        fit->SetParameter(0,initial);
        fit->SetParameter(1,lifetimeEst);

        histCopy->Fit("exp","QB");
        chi2_dof=fit->GetChisquare()/fit->GetNDF();
        //std::cout<<chi2_dof<<std::endl; 
       if (chi2_dof>1.05){
        for (int bin=0; bin<histCopy->GetXaxis()->GetNbins(); bin++){  if (histCopy->GetBinContent(bin+1)>0) histCopy->SetBinError(bin+1,0.001+histCopy->GetBinError(bin+1));}

        }
        if (chi2_dof<0.95){
        for (int bin=0; bin<histCopy->GetXaxis()->GetNbins(); bin++){ if (histCopy->GetBinContent(bin+1)>0)  histCopy->SetBinError(bin+1,abs(histCopy->GetBinError(bin+1)-0.001));}        }
        iter++;
        }

	    //int status=fit->CovMatrixStatus();
	    std::cout<<fit->GetParameter(1)<<std::endl; 
	    TLatex tL;
	    tL.SetNDC();
	    tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	  c1->Print(Form("fit%s.png",s.c_str()));
	  c1->Print(Form("fit%s.pdf",s.c_str()));
	  ofstream outfile;
	  outfile.open(Form("crt_lifetimeresults_Nov2018.txt"),std::ios_base::app);
	outfile<<runNumber<<','<<stampTime<<','<<stampDay<<','<<fit->GetParameter(0)<<','<<fit->GetParError(0)<<','<<fit->GetParameter(1)<<','<<fit->GetParError(1)<<std::endl;

	double x_pos=-340;
	for(int i=0; i<46; i++){
	x_pos=-360+15.5*i;
	double efield=tot_Ef(x_pos,300,300);
//	std::cout<<recom_correction(efield)<<std::endl;


	}



	}


