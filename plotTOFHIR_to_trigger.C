// g++ -Wall -o plotTOFHIR_to_trigger plotTOFHIR_to_trigger.C functions.hh `root-config --cflags --glibs` -lSpectrum

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"

#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"


#include "TSpline.h"
#include "TCanvas.h"
#include "TSpectrum.h"
// #include "TDirectory.h." 
#include "TObject.h"
#include <algorithm>
#include "functions.hh"

int main(int argc, char** argv)
{
    
    TApplication* theApp = new TApplication("App", &argc, argv);

    TLegend *leg;

    //read input files
    int firstRun = 1;    
    
    if (argc > 1) 
    {
        firstRun= atoi(argv[1]);
    }
    std::cout << "plotting from run = " << firstRun << std::endl;
    
    int lastRun = firstRun;
    if (argc > 2) 
    {
        lastRun= atoi(argv[2]);
    }
    
    bool withTRACKS = 1;
    if (argc > 3) 
    {
        withTRACKS = bool(argv[3]);
    }
    
    std::string data_path;
    if (withTRACKS) data_path = "../data/RecoData/RecoWithTracks/v2/";
    else            data_path = "../data/RecoData/RecoWithoutTracks/v1/";
    
    if (argc > 4) 
    {
        std::cout << " using data path from input line" << std::endl;
        data_path= std::string(argv[4]);
    }
    std::cout << "data_path = " << data_path << std::endl;
    
    
    
    //define tree
    float	step1;
    float	step2;
    float	x_dut;
    float	y_dut;
    int 	ntracks;        
    long long	chTime[400];    
    float	chtot[400];
            
    TChain * tree = new TChain("data", "data");
    
    for (int iRun = firstRun; iRun <= lastRun; iRun++)
    {
//        if ( std::find(bad_runs.begin(), bad_runs.end(), iRun) != bad_runs.end()) continue;       
       if (data_path == "../data/RecoData/RecoWithTracks/v1/")    tree->Add( Form("%s/run%.5d_events.root", data_path.c_str(), iRun) ); 
       
       else if (withTRACKS)    tree->Add( Form("%s/run%.5d_events_withTrack.root", data_path.c_str(), iRun) );       
       else               tree->Add( Form("%s/run%.4d_singles.root", data_path.c_str(), iRun) );
       
       std::cout << "adding run: " << iRun << std::endl;
    }
    
    tree->SetBranchStatus("*",0);
    tree->SetBranchAddress("step1", &step1);
    tree->SetBranchAddress("step2", &step2);

    if (withTRACKS)
    {
        tree->SetBranchAddress("x_dut", &x_dut);
        tree->SetBranchAddress("y_dut", &y_dut);
        tree->SetBranchAddress("ntracks", &ntracks);    
        tree->SetBranchAddress("chtot", &chtot);
        tree->SetBranchAddress("chTime", &chTime);
    }
    else    // old RECO
    {
        tree->SetBranchAddress("tot", &chtot);
        tree->SetBranchAddress("time", &chTime);        
    }
    
    
    
    //count how many steps are contained in the file
    std::cout << "loop to count steps in file" << std::endl;
    int NSTEP1 = 0;
    int NSTEP2 = 0;        
    std::vector<float> step1_vct;
    std::vector<float> step2_vct;
    step1_vct.clear();
    step2_vct.clear();
        
    Long64_t NEVENTS = tree->GetEntries();
//     NEVENTS = 100000;
    std::cout << "nEvents = " << NEVENTS << std::endl;
    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
//         std::cout << "step1 = " << step1 << " :: step2 = " << step2 <<std::endl;
        if (std::find(step1_vct.begin(), step1_vct.end(), step1) == step1_vct.end() )
        {             
            if (step1 > 1.e-2) 
            {
                step1_vct.push_back(step1);
                NSTEP1++;
            }
        }                
        if (std::find(step2_vct.begin(), step2_vct.end(), step2) == step2_vct.end() )
        {            
//             if (step2 > 1.e-2) 
            {
                step2_vct.push_back(step2);
                NSTEP2++;
            }
        }
    }
     
    for (int iStep1 = 0; iStep1<NSTEP1; iStep1++)
    {
        std::cout << "step1 (" << iStep1+1 << "/" << NSTEP1 << "): " << step1_vct.at(iStep1)  << std::endl;        
    }
    for (int iStep2 = 0; iStep2<NSTEP2; iStep2++)
    {
        std::cout << "step2 (" << iStep2+1 << "/" << NSTEP2 << "): " << step2_vct.at(iStep2)  << std::endl;        
    }
    
    


    // define histos    
    double minTot = 0;
    double maxTot = 700;
    
    double minTime = -200000;
    double maxTime = 200000;
    
    int NBINSPOS = 200;
    float minX = 0, maxX = 40;
    
    int REBIN_COEFF = 32;
    
    int NCH = 400;
    int myChList[] = {128, 130, 132, 134, 136, 138, 140, 142, 129, 131, 133, 135, 137, 139, 141, 143};
    int NBARS = 8;
//     int myChList[] = {0, 14, 1, 15};
//     int NBARS = 2;
    
    
    TH1F * hPosX = new TH1F ("hPosX", "hPosX", NBINSPOS, minX, maxX );
    TH1F * hPosY = new TH1F ("hPosY", "hPosY", NBINSPOS, minX, maxX );
    
    
    std::map<float, std::map<float, std::map<int, TH1F * > > > hEffX;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hEffY;
    
    std::map<float, std::map<float, std::map<int, TProfile2D * > > > pXY_Edep;
    
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTot;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTotSel;
    std::map<float, std::map<float, std::map<int, TH1F * > > > hTime;
    
    std::map<float, std::map<float, std::map<int, TH1F * > > > hCTR_UD;
    
        
//     std::map<float, std::map<float, TProfile2D * > > hTot_XY[NCH];
    
//     std::map<float, std::map<float, float > > tot_mean;    
//     std::map<float, std::map<float, float > > time_mean;
        
    
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
    {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
        {
            for (int iCh = 0; iCh < NCH; iCh++)
            {            
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hTot_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
                hTotSel[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTotSel_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hTotSel_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTot, maxTot );
                
                hTime[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTime_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hTime_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), 4000, minTime, maxTime );
                
                hEffX[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hEffX_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hEffX_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINSPOS, minX, maxX );
                
                hEffY[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hEffY_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hEffY_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINSPOS, minX, maxX );
                
                pXY_Edep[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile2D (Form("pXY_Edep_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("pXY_Edep_ch%.3d_step1_%.1f_step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINSPOS, minX, maxX, NBINSPOS, minX, maxX );
            }
            for (int iBar = 0; iBar < NBARS; iBar++)
            {
                hCTR_UD[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iBar] = new TH1F (Form("hCTR_UD_ch%.3d_step1_%.1f_step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), Form("hCTR_UD_ch%.3d_step1_%.1f_step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), 20000, -4000, 4000 );
            }                       
        }
    }
    
    
    //define more histos...               
//     TH2F * hTotTimeWalk  = new TH2F ("hAmpTimeWalk", "hAmpTimeWalk", 1000, -5, 5, 4000, -6000, 6000);
//     TProfile * pTotTimeWalk  = new TProfile ("pAmpTimeWalk", "pAmpTimeWalk", 1000, -5, 5);

    float pos0 = 4.8; 
    float bar_w = 3.05;
    float cut = 2.;
    
    
    
    int myOV = 0, myTH = 0;
//     float cutUp = 135, cutDown = 140;
    float cutUp = 125, cutDown = 120;

    //************************************************************************************//
    //              loop 0
    //************************************************************************************//
    
    std::cout << "(0) looping over events to get MIP peak position" << std::endl;
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
        tree->GetEntry(iEvt);
        if ( myOV != -1 && step1 != step1_vct.at(myOV)) continue;
        if ( myTH != -1 && step2 != step2_vct.at(myTH)) continue;
        
//         if ( myOV != -1 && step1 != myOV) continue;
//         if ( myTH != -1 && step2 != myTH) continue;
        
        if (iEvt%1000 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
//         if (step1!=0) continue;
//         hCTR[step1][step2]->Fill(time2-time1);
        long long time_ref = chTime[384];
        
        hPosX->Fill(x_dut);
        hPosY->Fill(y_dut);
            
        
        for (int iCh = 0; iCh<NCH; iCh++)
        {
            
            if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
            
//             if (chtot[iCh]>0) std::cout << "filling ch[" << iCh << "] with tot = " << chtot[iCh]/1.e3 << " ns :: and time-t_ref = " << chTime[iCh] << "  - " << time_ref << " = " << (chTime[iCh]-time_ref) << " posX = " << x_dut << " :: posY = " << y_dut << std::endl;
            hTot[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
            int barID = (iCh-128)/2;
            
            
            if ( fabs(x_dut-(pos0 + bar_w*barID)) < cut/2) 
            {
                hTotSel[step1][step2][iCh]->Fill(chtot[iCh]/1.e3);
                std::cout << "iCh = " << iCh << " :: barID = " << barID << " posBar = " << pos0 + bar_w*barID << " :: beam_hit = " << fabs(x_dut-pos0 + bar_w*barID) << std::endl;
            }
            hTime[step1][step2][iCh]->Fill( (chTime[iCh] - time_ref));                        
            
            
            if (x_dut != 0 && y_dut != 0)
            {
                pXY_Edep[step1][step2][iCh]->Fill(x_dut, y_dut, chtot[iCh]/1.e3);
            
                if (chtot[iCh]/1.e3>130)
                {
                    hEffX[step1][step2][iCh]->Fill(x_dut, 1);
                    hEffY[step1][step2][iCh]->Fill(y_dut, 1);
                }
            }
//             hTime[step1][step2][iCh]->Fill(chTime[iCh] - time_ref);                        
        }
        
        for (int iBar = 0; iBar<NBARS; iBar++)
        {
            int chBarUp = iBar*2+128;
            int chBarDown = iBar*2+129;
            
//             int chBarUp = iBar;
//             int chBarDown = iBar + 14;
            
            if (chtot[chBarUp]>cutUp && chtot[chBarDown]>cutDown ) 
            {
                hCTR_UD[step1][step2][iBar]->Fill(chTime[chBarUp] - chTime[chBarDown]);                        
//                 std::cout << "CT_UD [" << iBar << "] = " << chTime[chBarUp] << " - " << chTime[chBarDown] << " = " << (chTime[chBarUp] - chTime[chBarDown]) << std::endl;
            }
        }
            
    }
    
    
    
    
    
    float minTotForPeakSearch = 40;
    
    TCanvas *cTots_scan[NSTEP1][NSTEP2][NCH];
    for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
    {
        for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
        {
            for (int iCh = 0; iCh< NCH; iCh++)
            {
                if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
                if ( myOV != -1 && step1 != myOV) continue;
                if ( myTH != -1 && step2 != myTH) continue;
                
                cTots_scan[iStep1][iStep2][iCh] = new TCanvas (Form("cTots_ch%.3d_step1_%.1f_step2_.%1f", iCh, step1_vct.at(iStep1), step1_vct.at(iStep1)), Form("cTots_ch%.3d_step1_%.1f_step2_.%1f", iCh, step1_vct.at(iStep1), step1_vct.at(iStep1)), 800, 400);
//                 cTots_scan[iStep1][iStep2][iCh]->Divide(2,1);
                
                cTots_scan[iStep1][iStep2][iCh]->cd();            
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Rebin(REBIN_COEFF);
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw();
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetXaxis()->SetTitle("tot [ns]");
                hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->GetYaxis()->SetTitle("Counts");
                
                
                hTotSel[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->SetLineColor(kGreen+2);                
                hTotSel[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]->Draw("same");
                
                gPad->SetLogy();
            }
            
            /*
            hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->SetRange(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->FindBin(minTotForPeakSearch) , hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetXaxis()->GetNbins());
//             float x_max = hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetBinCenter(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->GetMaximumBin());
            
            int npeaks = 3;
            TSpectrum *s = new TSpectrum(npeaks);
            Int_t nfound = s->Search(hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)], 20, "nobackground", 0.1);    //sigma, option, threshold
            printf("Found %d candidate peaks in tot1 to fit\n",nfound);
            Double_t *xpeaks = s->GetPositionX();               
            Double_t xp = xpeaks[0];
            for (int ip = 0; ip<npeaks; ip++)
            {
               if (xp<xpeaks[ip]) xp = xpeaks[ip];
            }
            float x_max = xp;
        
            TF1 * fitTots = new TF1 ("fitTots", "gaus", x_max*0.9, x_max*1.15);
            hTot1[step1_vct.at(iStep1)][step2_vct.at(iStep2)]->Fit(fitTots, "QR");
            tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = fitTots->GetParameter(1);
            tot1_sigma[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = fitTots->GetParameter(2);
            
            std::cout << "tot1[" << step1_vct.at(iStep1)<<"]["<<step2_vct.at(iStep2)<<"] = " << tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)] << " +/- " << tot1_sigma[step1_vct.at(iStep1)][step2_vct.at(iStep2)] << std::endl; // << " :: res = " << tot1_sigma/tot1_mean << std::endl;    
            gStep1Tot1[iStep1]->SetPoint(iStep2, step2_vct.at(iStep2), tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)]);
            gStep2Tot1[iStep2]->SetPoint(iStep1, step1_vct.at(iStep1), tot1_mean[step1_vct.at(iStep1)][step2_vct.at(iStep2)]);
            */
        }
    }
    
    
    int selStep1 = myOV;
    int selStep2 = myTH;

    TCanvas * cBeamPosX = new TCanvas ("cBeamPosX", "cBeamPosX", 500, 500);
    cBeamPosX->cd();
    hPosX->Draw();
    hPosX->GetXaxis()->SetTitle("beam X pos [mm]");
    hPosX->GetYaxis()->SetTitle("Counts");        
    gPad->SetLogy();

    TCanvas * cBeamPosY = new TCanvas ("cBeamPosY", "cBeamPosY", 500, 500);
    cBeamPosY->cd();
    hPosY->Draw();
    hPosY->GetXaxis()->SetTitle("beam Y pos [mm]");
    hPosY->GetYaxis()->SetTitle("Counts");        
    gPad->SetLogy();
    
    
    TCanvas * cArrayEffX = new TCanvas ("cArrayEffX", "cArrayEffX", 1600, 300);
    cArrayEffX->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        //normalize efficiency plots
    //
        for (int iBin = 0; iBin < NBINSPOS; iBin++)
        {
            if (hPosX->GetBinContent(iBin+1)>0) hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetBinContent(iBin+1, hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinContent(iBin+1)/hPosX->GetBinContent(iBin+1) );
        }
        cArrayEffX->cd(chId+1);
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("beam X pos [mm]");
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Efficiency");        
        gPad->SetLogy();
    }
    
    TCanvas * cArrayEffY = new TCanvas ("cArrayEffY", "cArrayEffY", 1600, 300);
    cArrayEffY->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEffY->cd(chId+1);
        hEffY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hEffY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("beam Y pos [mm]");
        hEffY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Efficiency");        
        gPad->SetLogy();
    }
    
    /*
    TF1 * fitBarPos = new TF1 ("fitBarPos", fitBarEfficiency, minX, maxX, 4);
    fitBarPos->SetParameters(1., 3., 0.1, 0.8);
    fitBarPos->SetParLimits(0, minX, maxX);
    fitBarPos->SetParLimits(1, 2, 4);
    fitBarPos->SetParLimits(2, 0.01, 0.4);
    fitBarPos->SetParLimits(3, 0.2, 1.);
    */
    
    TF1 * fitBarPos = new TF1 ("fitBarPos", fitBarEffErr, 2, 32, 5);
    fitBarPos->SetParameters(5., 3., 0.01, 0.8, 0.2);
    fitBarPos->SetNpx(5000);
//     fitBarPos->FixParameter(1, 3.);
//     fitBarPos->FixParameter(4, 1.);
//     fitBarPos->FixParameter(2, 0.01);
    
    
    fitBarPos->SetParLimits(0, minX, maxX);
    fitBarPos->SetParLimits(1, 2, 4);
//     fitBarPos->SetParLimits(2, 0.001, 0.03);
    fitBarPos->SetParLimits(3, 0.2, 1.);
    
    
    
//     TF1 * fitBarPos = new TF1 ("fitBarPos", "[0]*pow((x-[1]),4) + [2]", minX, maxX);
    
    
    
    TCanvas * cCompareEffX_top = new TCanvas ("cCompareEffX_top", "cCompareEffX_top", 600, 400);
    cCompareEffX_top->cd();
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->Draw();
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->GetXaxis()->SetTitle("beam X pos [mm]");
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->GetYaxis()->SetTitle("Efficiency");        
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[0]]->GetYaxis()->SetRangeUser(0.01, 2);
    
    for (int chId = 0; chId<NBARS; chId++)
    {
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(chId+1);
        fitBarPos->SetLineColor(chId+1);
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
        fitBarPos->SetParameter(0, 4 + chId*3.);
        for (int i = 0; i< 3; i++) hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitBarPos, "QR");
        float posBar   = fitBarPos->GetParameter(0);
        float widthBar = fitBarPos->GetParameter(1);
        float effBar   = fitBarPos->GetParameter(3);
        float trackRes = fitBarPos->GetParameter(4);
        std::cout << "posBar[" << chId << "] = " << posBar << " :: width = " << widthBar << " :: effBar = " << effBar << " :: trackRes = " << trackRes << std::endl;
        
    }
    gPad->SetLogy();
    
    
    TCanvas * cCompareEffX_bot = new TCanvas ("cCompareEffX_bot", "cCompareEffX_bot", 600, 400);
    cCompareEffX_bot->cd();
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[NBARS]]->Draw();
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[NBARS]]->GetXaxis()->SetTitle("beam X pos [mm]");
    hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[NBARS]]->GetYaxis()->SetTitle("Efficiency");        
    
    for (int chId = NBARS; chId<NBARS*2; chId++)
    {
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(chId - NBARS +1);
        hEffX[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
    }
    gPad->SetLogy();
    
    
    TCanvas * cArrayEdepXY = new TCanvas ("cArrayEdepXY", "cArrayEdepXY", 1600, 300);
    cArrayEdepXY->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEdepXY->cd(chId+1);
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("COLZ");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("beam X pos [mm]");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("beam Y pos [mm]");
        pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetZaxis()->SetTitle("Mean Tot [ns]");
//         gPad->SetLogy();
    }
    
    
    TLine *cutTotUp     = new TLine (cutUp, 0, cutUp, 1000);
    TLine *cutTotDown   = new TLine (cutDown, 0, cutDown, 1000);
    
    cutTotUp   ->SetLineColor(kRed);
    cutTotDown ->SetLineColor(kBlue);
    
    TCanvas * cArrayTots = new TCanvas ("cArrayTots", "cArrayTots", 1600, 300);
    cArrayTots->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayTots->cd(chId+1);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("tot [ns]");
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
        cutTotUp->Draw("same");
        
        
        hTotSel[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(kGreen+2);                
        hTotSel[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
                
    }
    
    
    TCanvas * cArrayTimes = new TCanvas ("cArrayTimes", "cArrayTimes", 1600, 300);
    cArrayTimes->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayTimes->cd(chId+1);
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("Time [ps]");
        hTime[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
    }
    
    
    TCanvas * cArrayCTR_UD = new TCanvas ("cArrayCTR_UD", "cArrayCTR_UD", 1600, 300);
    cArrayCTR_UD->Divide(NBARS, 1);
    for (int barId = 0; barId<NBARS; barId++)
    {
        
        cArrayCTR_UD->cd(barId+1);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Rebin(REBIN_COEFF);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Draw();
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetXaxis()->SetTitle("t_{up} - t_{down} [ps]");
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetYaxis()->SetTitle("Counts");
        
        float mean = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetMean();
        float maxBin = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetMaximumBin();
        float maxX = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetBinCenter(maxBin);
        float rms  = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->GetRMS();
        
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", maxX-rms*0.5 , maxX+rms*0.5);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][barId]->Fit(fitGaus, "QRL");
        
        std::cout << "CTR_UD [" << barId << "] = " << fitGaus->GetParameter(2) << " ps --> sigma_bar = " << fitGaus->GetParameter(2)/2 << " ps " << std::endl;
    }
    
    
    int selBar = 5;
    int selBarCh = selBar*2+128;
//     int chBarDown = iBar*2+129;
    
    TCanvas * cSingleScatterEdep = new TCanvas ("cSingleScatterEdep", "cSingleScatterEdep", 500, 500);
    cSingleScatterEdep->cd();
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->SetStats(0);
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->Draw("COLZ");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetXaxis()->SetTitle("beam X pos [mm]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetYaxis()->SetTitle("beam Y pos [mm]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetZaxis()->SetTitle("Mean Tot [ns]");
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetXaxis()->SetRangeUser(0, 30);
    pXY_Edep[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarCh]->GetYaxis()->SetRangeUser(10, 40);
    
    
    
    TCanvas * cTestFunction = new TCanvas ("cTestFunction", "cTestFunction", 500, 500);
    
//     TF1 * fitTest = new TF1 ("fitTest", fitBarEffErr, minX, maxX, 4);
//     fitTest->SetParameters(1., 3., 0.01, 0.8);
    
    TF1 * fitTest = new TF1 ("fitTest", "[0]/2*(erf((x+[1])/(sqrt(2)*[2])) + 1)", -50, 50);    
    fitTest->SetParameters(0.8, 0., 2);
    fitTest->SetNpx(5000);
    fitTest->Draw();
    
    TF1 * fitTest2 = new TF1 ("fitTest2", "[0]/2*(erf((x+[1])/(sqrt(2)*[2])) + 1)", -50, 50);    
    fitTest2->SetParameters(0.8, 8.5, 2);
    fitTest2->SetNpx(5000);
    fitTest2->SetLineColor(kGreen+1);
    fitTest2->Draw("same");
    
    TF1 * fitTest3 = new TF1 ("fitTest3", "[0]/2*(erf((-x+[1])/(sqrt(2)*[2])) + 1)", -50, 50);    
    fitTest3->SetParameters(0.8, 5.5, 2);
    fitTest3->SetNpx(5000);
    fitTest3->SetLineColor(kBlue+1);
    fitTest3->Draw("same");
    
    TF1 * fitTestMyFunc = new TF1 ("fitTestMyFunc", fitBarEffErr, -20, 20, 5);
    fitTestMyFunc->SetParameters(10., 3., 0.0, 0.8, 0.2);
    fitTestMyFunc->SetLineColor(kBlack);
    fitTestMyFunc->SetNpx(10000);
    fitTestMyFunc->Draw("same");
    
    
//     fitBarPos->FixParameter(1, 3.);
//     fitBarPos->FixParameter(2, 0.01);
    
    
//     fitBarPos->SetParLimits(0, minX, maxX);
//     fitBarPos->SetParLimits(1, 2, 4);
//     fitBarPos->SetParLimits(2, 0.01, 0.4);
//     fitBarPos->SetParLimits(3, 0.2, 1.);
    
    
    
    /*    
    //save plots or histos/graphs to output file
    TFile * outputFile = new TFile(Form("./output/output_run_%.3d.root", firstRun), "RECREATE");
    outputFile->cd();    
    
    outputFile->Write();
    outputFile->Close();    
    */
    
    
    std::cout << "end of program" << std::endl;
    theApp->Run();

    
    

}
