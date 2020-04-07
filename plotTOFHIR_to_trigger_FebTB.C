// to compile: g++ -Wall -o plotTOFHIR_to_triggerFebTB plotTOFHIR_to_trigger_FebTB.C functions.cc `root-config --cflags --glibs` -lSpectrum

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
  // open the file that important things are written to
  std::ofstream myfile;
  myfile.open ("example.txt");
  myfile << "Saving bar position and MPV from the lambda fit." << std::endl;
  
  gStyle->SetTitleXOffset (1.00) ;                                                                                        
  gStyle->SetTitleYOffset (1.2) ;                                                                                                                                                                                                                 
  gStyle->SetPadLeftMargin (0.13) ;                                                                                       
  gStyle->SetPadBottomMargin (0.13) ;                                                                                                                                                                                                              
  gStyle->SetTitleSize (0.05, "xyz") ;                                                                                    
  gStyle->SetLabelSize (0.035,"xyz") ; 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.035);
  TLegend *leg;

//     gROOT->SetBatch(true);
    TApplication* theApp = new TApplication("App", &argc, argv);
    //TLegend *leg;

    //read input files - based on first and last run number listed in the command line arguments
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
    
//     std::string data_path = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1
    std::string data_path = "../data_TB/";  
    if (argc > 3) 
    {
      std::cout << "using data path from input line" << std::endl;
      data_path = std::string(argv[3]);
    }
    std::cout << "data_path = " << data_path << std::endl;
    
    std::string output_plot_folder = "../plots";
    int selStep1     = 0;
    int selStep2     = 0;
    int selBarNumber = 8;
    
    
    // define tree
    float	step1;    
    float	step2;
    float	x_dut;
    float	y_dut;
    float	xIntercept;
    float	yIntercept;
    int 	ntracks;        
    long long   time[400];    
    float	tot[256];
    float       qfine[256];
    float       energy[256];      
    
    
    TChain * tree = new TChain("data", "data");

    
    // going over each of the runs listed as command line arguments
    for (int iRun = firstRun; iRun <= lastRun; iRun++)
      {
	//if ( std::find(bad_runs.begin(), bad_runs.end(), iRun) != bad_runs.end()) continue;       
	tree->Add( Form("%s/run%.5d_events.root", data_path.c_str(), iRun) );  
	std::cout << "adding run: " << iRun << std::endl;
      }

    // extract information from the ROOT tree
    tree->SetBranchStatus("*",0);
    tree->SetBranchAddress("step1", &step1); // overvoltage
    tree->SetBranchAddress("step2", &step2); // (threshold 1 * 10000 + 1) + (threshold 2 * 100 + 1) + (threshold E + 1), first 2 digits are vth1, second 2 are vth2, last 2 are energy 
    tree->SetBranchAddress("x_dut", &x_dut);
    tree->SetBranchAddress("y_dut", &y_dut);
    tree->SetBranchAddress("xIntercept", &xIntercept);
    tree->SetBranchAddress("yIntercept", &yIntercept);
    tree->SetBranchAddress("ntracks", &ntracks);    
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("qfine", &qfine); // qfine is integral of amplified signal in range defined by tot
    tree->SetBranchAddress("energy", &energy); // energy found from calibration plot of qfine vs tot, and then difference between calibration tot and measured
    
    //count how many steps are contained in the file
    std::cout << "loop to count steps in file" << std::endl;
    int NSTEP1 = 0;
    int NSTEP2 = 0;        
    std::vector<double> step1_vct;
    std::vector<double> step2_vct;
    step1_vct.clear();
    step2_vct.clear();
        
    // starting the EVENT loop now
    Long64_t NEVENTS = tree->GetEntries();
    std::cout << "nEvents = " << NEVENTS << std::endl;    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
      {
        tree->GetEntry(iEvt);
        //std::cout << "step1 = " << step1 << " :: step2 = " << step2 <<std::endl;
	// getting vth1 vth2 and vthe from step2 for each event
 	double vth1 = double(int(step2/10000)-1);;
 	double vth2 = double(int((step2-10000*vth1)/100)-1);
 	double vthe = double(int((step2-10000*vth1-step2-100*vth2)/1)-1);

	if (ntracks != 1) continue;
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
            step2_vct.push_back(step2);
            NSTEP2++;
	  }
      }
     
    // step 1 is overvoltage                                                                                                                                           
    for (int iStep1 = 0; iStep1<NSTEP1; iStep1++)
      {
        std::cout << "step1 (" << iStep1+1 << "/" << NSTEP1 << "): " << step1_vct.at(iStep1)  << std::endl;        
      }
    // step 2 is vth1, vth2, vthe in each two digits
    for (int iStep2 = 0; iStep2<NSTEP2; iStep2++)
      {
        std::cout << "step2 (" << iStep2+1 << "/" << NSTEP2 << "): " << step2_vct.at(iStep2)  << std::endl;        
      }    

    // define histos    
    double minTot = 0;
    double maxTot = 700;
    double minEnergy = 0;
    double maxEnergy = 50;
    double minTime = -200000;
    double maxTime = 200000;

    int REBIN_COEFF = 16;
    
    #define NCH 400
    

    // this is for configuration 2, 3 (connector was rotated wrong), 4, or 5
    // conf 6 and 8 have the horizontal array as the low channel numbers. Config 9 has both arrays vertical - mostly focusing on this one
    //    int myChList[] = {0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}; // from channelMapping1 in mtd drawMatrices.cfg, VERTICAL
    //    int myChTop[] = {3,10,0,1,7,14,5,12,21,20,23,22,16,18,17,19}; // one side from channelMapping1
    //    int myChBot[] = {4,6,15,8,13,2,11,27,32,31,30,29,28,26,24,25}; // one side from channelMapping1
    //    int myChList[] = {57,63,50,60,59,55,61,56,58,53,62,54,9,51,38,52,34,43,33,42,36,44,35,46,37,45,39,47,41,48,40,49}; // channelMapping2 HORIZONTAL
    //    int myChTop[] = {57,50,59,61,58,62,9,38,34,33,36,35,37,39,41,40}; // one side from channelMapping2, totalEnergy[0]
    //    int myChBot[] = {63,60,55,56,53,54,51,52,43,42,44,46,45,47,48,49}; // one side from channelMapping2, totalEnergy[1]
    
    
    // mapping for the pin connector array, caltech array
    int myChList[32] = {63,57,60,50,55,59,56,61,53,58,54,62,51,9,52,38,
                        40,49,41,48,39,47,37,45,35,46,36,44,33,42,34,43}; // VERTICAL
                        
    int myChTop[16] = {63,57,60,50,55,59,56,61,53,58,54,62,51,9,52,38}; // one side from channelMapping2, totalEnergy[0]
    int myChBot[16] = {40,49,41,48,39,47,37,45,35,46,36,44,33,42,34,43}; // one side from channelMapping2, totalEnergy[1]
    bool HORIZONTAL = true;
    
    int NBARS = 16; // full array has 16 bars, 32 SiPM readouts




 
 // declare the histograms, these will be filled in the channel loop    
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTot_cut;
 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hEnergy;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hEnergy_cut;
 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hTime;
 std::map<float, std::map<float, std::map<int, TProfile * > > > pBarResponse_vs_X;
 std::map<float, std::map<float, std::map<int, TProfile * > > > pBarResponse_vs_Y;
 std::map<float, std::map<float, std::map<int, TProfile2D * > > > pBarResponse_vs_XY;
 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hEff_vs_X; 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hEff_vs_Y;
 std::map<float, std::map<float, std::map<int, TH2F * > > > hEff_vs_XY;
 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hCTR_UD;
 std::map<float, std::map<float, TProfile * > > pTot_vs_Xpos_overlay;
 std::map<float, std::map<float, TProfile * > > pTot_vs_Ypos_overlay;
 
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT_cut;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT_cut_ln;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT_cut_lln;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT_cut_rn;
 std::map<float, std::map<float, std::map<int, TH1F * > > > hXT_cut_rrn;
 

 int NBINPOS = 500;
 double minXpos = 0;
 double maxXpos = 50;
 double minYpos = 0;
 double maxYpos = 50;
 
 std::vector<double> * trk_blacklist = new std::vector<double>();
 
 int NBL = 14;
 double bl[NBL] = {0., 5.35123,  26.4772, 31.3572,  18.2681, 22.6292,   5.76614,  5.31838,  23.264,  18.1366,  27.6117,  21.0186,  19.3968};
 for (int i = 0; i<NBL ; i++)
 {
    trk_blacklist->push_back(bl[i]);
 }
 
 TH1F * hPosX = new TH1F ("hPosX", "hPosX", NBINPOS, minXpos, maxXpos);
 TH1F * hPosY = new TH1F ("hPosY", "hPosY", NBINPOS, minYpos, maxYpos);
 
 int NBINS = 2000;
 
 // step 1, step 2, and channel listing loops. Histograms are defined inside the loops
 for (int iStep1 = 0; iStep1< NSTEP1; iStep1++)
   {
     for (int iStep2 = 0; iStep2< NSTEP2; iStep2++)
       {
	 for (int iCh = 0; iCh < NCH; iCh++)
	   {       
               
             if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
             
             std::string this_caption = Form("ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2);
             std::string this_title   = Form(", ch%.3d, step1_%.1f, step2_%.1f", iCh, step1_vct.at(iStep1), step2_vct.at(iStep2));
             
             
	     hTot[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_%s", this_caption.c_str()), Form("Time over threshold, %s", this_title.c_str()), NBINS, minTot, maxTot );             	     
	     hTot_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTot_cut_%s", this_caption.c_str()), Form("Time over threshold with x cut %s", this_title.c_str()), NBINS, minTot, maxTot );
             
             hEnergy[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hEnergy_%s", this_caption.c_str()), Form("Energy %s", this_title.c_str()), NBINS, minEnergy, maxEnergy );
             hEnergy_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hEnergy_cut_%s", this_caption.c_str()), Form("Energy (cut) %s", this_title.c_str()), NBINS, minEnergy, maxEnergy );
             
	     hTime[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hTime_%s", this_caption.c_str()), Form("Time Hist %s", this_title.c_str()), NBINS, minTime, maxTime );
             

             pBarResponse_vs_X[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]  = new TProfile (Form("pBarResponse_vs_Xpos_%s", this_caption.c_str()), Form("Energy vs. X pos %s", this_title.c_str()), NBINPOS, minXpos, maxXpos );
	     pBarResponse_vs_Y[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]  = new TProfile (Form("pBarResponse_vs_Ypos_%s", this_caption.c_str()), Form("Energy vs. Y pos %s", this_title.c_str()), NBINPOS, minYpos, maxYpos );
             pBarResponse_vs_XY[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TProfile2D (Form("pBarResponse_vs_XY_%s", this_caption.c_str()), Form("XY Energy dep %s", this_title.c_str()), NBINPOS, minXpos, maxXpos, NBINPOS, minYpos, maxYpos );
             
	     hEff_vs_X[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]  = new TH1F (Form("hEff_vs_Xpos_%s", this_caption.c_str()), Form("Efficiency vs. X pos %s", this_title.c_str()), NBINPOS, minXpos, maxXpos);
	     hEff_vs_Y[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh]  = new TH1F (Form("hEff_vs_Ypos_%s", this_caption.c_str()), Form("Efficiency vs. Y pos %s", this_title.c_str()), NBINPOS, minYpos, maxYpos);            
 	     hEff_vs_XY[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH2F (Form("hEff_vs_XY_%s", this_caption.c_str()), Form("Efficiency vs X-Y %s", this_title.c_str()), NBINPOS, minXpos, maxXpos, NBINPOS, minYpos, maxYpos );
	     
	     
             
	     hXT[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_%s", this_caption.c_str()), Form("Cross Talk  %s", this_title.c_str()), 400, 0, 1);
             hXT_cut[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_cut_%s", this_caption.c_str()), Form("Cross Talk  %s", this_title.c_str()), 400, 0, 1);
             hXT_cut_ln[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_cut_ln_%s", this_caption.c_str()), Form("Cross Talk ln bar %s", this_title.c_str()), 400, 0, 1);
             hXT_cut_lln[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_cut_lln_%s", this_caption.c_str()), Form("Cross Talk lln bar %s", this_title.c_str()), 400, 0, 1);
             hXT_cut_rn[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_cut_rn_%s", this_caption.c_str()), Form("Cross Talk rn bar %s", this_title.c_str()), 400, 0, 1);
             hXT_cut_rrn[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iCh] = new TH1F (Form("hXT_cut_rrn_%s", this_caption.c_str()), Form("Cross Talk rrn bar %s", this_title.c_str()), 400, 0, 1);
             
	   }

	 pTot_vs_Xpos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Xpos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. X pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINS, minXpos, maxXpos );
	 pTot_vs_Ypos_overlay[step1_vct.at(iStep1)][step2_vct.at(iStep2)] = new TProfile (Form("pTot_vs_Ypos_over_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("ToT vs. Y pos overlay, step1_%.1f, step2_%.1f", step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINS, minYpos, maxYpos );

	 // bar loop (outside of channel loop) to define time resolution for a single bar
	 for (int iBar = 0; iBar < NBARS; iBar++)
	   {
	     hCTR_UD[step1_vct.at(iStep1)][step2_vct.at(iStep2)][iBar] = new TH1F (Form("hCTR_UD_ch%.3d_step1vct_%.1f_step2vct_%.1f_step1_%i_step2_%i", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2), iStep1, iStep2), Form("Bar Time Res, UD_ch%.3d, step1_%.1f, step2_%.1f", iBar, step1_vct.at(iStep1), step2_vct.at(iStep2)), NBINS, -2000, 2000 );
	   }                       
       }
   }
    
 //************************************************************************************//
 //              loop 0
 //************************************************************************************//
 
 std::cout << "(0) looping over events to get MIP peak position" << std::endl;
 for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
 {
     //std::cout << "processing event before 1 MIP req "<< iEvt<< std::endl;
     tree->GetEntry(iEvt);
     if (ntracks != 1) continue; // require 1 MIP per event
     
     if (iEvt%1000 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
     //	if (step1!=6) continue;
     
     

     if (std::find(trk_blacklist->begin(), trk_blacklist->end(), double(x_dut)) != trk_blacklist->end() ) continue;
     
     long long time_ref = time[384];
     
     hPosX->Fill(xIntercept);
     hPosY->Fill(yIntercept);

     for (int iCh = 0; iCh<NCH; iCh++) // channel loop in the event loop
     {
	 if (std::find(std::begin(myChList), std::end(myChList), iCh) == std::end(myChList) ) continue;
         
         // skip any thing where the energy isnt a normal value (should always be postive, default value is -9999) for both energy and tot
	 
         if (!passBugEventCuts(iCh, energy, qfine, tot)) continue;
         
	 //no cuts
	 hTot[step1][step2][iCh]->Fill(tot[iCh]/1.e3); 
         hEnergy[step1][step2][iCh]->Fill(energy[iCh]);
         hTime[step1][step2][iCh]->Fill(time[iCh] - time_ref);
                  
         
         double my_signal = energy[iCh];	 
	 if (xIntercept != 0 && yIntercept != 0)
         {
             pBarResponse_vs_X[step1][step2][iCh]->Fill(xIntercept, my_signal);
             pBarResponse_vs_Y[step1][step2][iCh]->Fill(yIntercept, my_signal);
	     pBarResponse_vs_XY[step1][step2][iCh]->Fill(xIntercept, yIntercept, my_signal);            
             
         }

         // make efficiency plots
//          if (my_signal >= 0.9 * MIP[iCh] && my_signal <= 3 * MIP[iCh] && xIntercept > 0)
         if (my_signal > 10 && my_signal < 100
             && yIntercept > 0 && xIntercept > 0
         )
         {
            hEff_vs_X[step1][step2][iCh]->Fill(xIntercept, 1);
            hEff_vs_Y[step1][step2][iCh]->Fill(yIntercept, 1);
            hEff_vs_XY[step1][step2][iCh]->Fill(xIntercept, yIntercept, 1);
            
            hEnergy_cut[step1][step2][iCh]->Fill(my_signal); // cut out low region where energy and tot are non-linear with each other
            
         }
    }
  }
 
 
 
 
 
    // plotting first round of plots
    std::cout << std::endl;
    std::cout << "first set of plots" << std::endl;
    
    TCanvas * cBeamPosX = new TCanvas ("cBeamPosX", "cBeamPosX", 500, 500);
    cBeamPosX->cd();
    hPosX->Rebin(2);
    hPosX->Draw();
    hPosX->GetXaxis()->SetTitle("beam X pos [mm]");
    hPosX->GetYaxis()->SetTitle("Counts");        
    gPad->SetLogy();

    TCanvas * cBeamPosY = new TCanvas ("cBeamPosY", "cBeamPosY", 500, 500);
    cBeamPosY->cd();
    hPosY->Rebin(2);
    hPosY->Draw();
    hPosY->GetXaxis()->SetTitle("beam Y pos [mm]");
    hPosY->GetYaxis()->SetTitle("Counts");
    gPad->SetLogy();

    
    
    TCanvas * cArrayTots = new TCanvas ("cArrayTots", "cArrayTots", 1600, 300);
    cArrayTots->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayTots->cd(chId+1);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
	hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRangeUser(200, 750);
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("Tot [ns]");
        hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
    }
    
    TCanvas * cTotsZoom = new TCanvas ("cTotsZoom", "cTotsZoom", 500, 800);
    cTotsZoom->Divide(1, 2);
    cTotsZoom->cd(1);
    hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw();
    hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->GetXaxis()->SetRangeUser(200,750);

    
    cTotsZoom->cd(2);
    hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw();
    hTot[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->GetXaxis()->SetRangeUser(200,750);

      
 
    std::cout << "ICs" << std::endl;
 
    double MIP_peak[NBARS*2] = {0};
    double MIP_peak_err[NBARS*2] = {0};
    double IC[NBARS*2] = {0};
    double avgMIP = 0; 
    int peak_count = 0;
    
    TH1F * hIC_tot = new TH1F ("hIC_tot", "hIC_tot", 20, 0.7, 1.3);
    TH1F * hIC_top = new TH1F ("hIC_top", "hIC_top", 20, 0.7, 1.3);
    TH1F * hIC_bot = new TH1F ("hIC_bot", "hIC_bot", 20, 0.7, 1.3);
    
    TCanvas * cArrayEnergies = new TCanvas ("cArrayEnergies", "cArrayEnergies", 1600, 300);
    cArrayEnergies->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEnergies->cd(chId+1);
        hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
 	hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRangeUser(0,100);
        hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("Energy [a.u.]");
        hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
        
        hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(kGreen+2);
        hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetMarkerColor(kGreen+2);
        hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
        
        TF1 * fitLandau = new TF1 ("fitLandau","landau", 0, 100);
	hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitLandau, "QRL");
        MIP_peak[chId] = fitLandau->GetParameter(1);        
        MIP_peak_err[chId] = fitLandau->GetParError(1);        
        if (MIP_peak[chId]>0)
        {
            avgMIP += MIP_peak[chId];
            peak_count++;
        }
    }
    
        
    
    TCanvas * cEnergiesZoom = new TCanvas ("cEnergiesZoom", "cEnergiesZoom", 500, 800);
    cEnergiesZoom->Divide(1, 2);
    cEnergiesZoom->cd(1);
    hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw();
    hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("same");
    cEnergiesZoom->cd(2);
    hEnergy[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw();
    hEnergy_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("same");
 
    
    
//     avgMIP /= NBARS*2;
    avgMIP /= peak_count; // we have 29 active channels!!    
    TGraphErrors * gMIP_top = new TGraphErrors();
    TGraphErrors * gMIP_bot = new TGraphErrors();
    TGraphErrors * gIC_top = new TGraphErrors();
    TGraphErrors * gIC_bot = new TGraphErrors();
    
    for (int chId = 0; chId<NBARS*2; chId++) 
    {
        IC[chId] = MIP_peak[chId] / avgMIP;
        if (chId<NBARS) 
        {
            gMIP_top->SetPoint(chId, chId, MIP_peak[chId]);
            gMIP_top->SetPointError(chId, 0, MIP_peak_err[chId]);
            
            gIC_top->SetPoint(chId, chId, IC[chId]);
            gIC_top->SetPointError(chId, 0, MIP_peak_err[chId]/avgMIP);
            hIC_top->Fill(IC[chId]);
        }
        else
        {
            gMIP_bot->SetPoint(chId-NBARS, chId-NBARS, MIP_peak[chId]);
            gMIP_bot->SetPointError(chId-NBARS, 0, MIP_peak_err[chId]);
            
            gIC_bot->SetPoint(chId-NBARS, chId-NBARS, IC[chId]);
            gIC_bot->SetPointError(chId-NBARS, 0, MIP_peak_err[chId]/avgMIP);
            hIC_bot->Fill(IC[chId]);
        }
        hIC_tot->Fill(IC[chId]);
    }
    
    TCanvas * cMIP_peak = new TCanvas ("cMIP_peak", "cMIP_peak", 600, 500);
    cMIP_peak->cd();
    gMIP_top->Draw("APE");
    gMIP_top->GetYaxis()->SetRangeUser(0, 20);
    gMIP_top->GetYaxis()->SetTitle("MIP peak [a.u.]");
    gMIP_top->SetMarkerStyle(20);
    gMIP_top->SetMarkerColor(kBlue+1);
    gMIP_top->SetLineColor(kBlue+1);
    
    gMIP_bot->SetMarkerStyle(21);
    gMIP_bot->SetMarkerColor(kRed+1);
    gMIP_bot->SetLineColor(kRed+1);
    gMIP_bot->Draw("same PE");
    
    leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
    leg->AddEntry(gMIP_top, "top SiPMs", "lp" );
    leg->AddEntry(gMIP_bot, "bottom SiPMs", "lp" );
    leg->Draw();
    gPad->SetGridy();
    
    
    TCanvas * cIC = new TCanvas ("cIC", "cIC", 600, 500);
    cIC->cd();
    gIC_top->Draw("APE");
    gIC_top->GetYaxis()->SetRangeUser(0.7, 1.3);
    gIC_top->GetYaxis()->SetTitle("MIP peak / avg MIP ");
    gIC_top->SetMarkerStyle(20);
    gIC_top->SetMarkerColor(kBlue+1);
    gIC_top->SetLineColor(kBlue+1);
    
    gIC_bot->SetMarkerStyle(21);
    gIC_bot->SetMarkerColor(kRed+1);
    gIC_bot->SetLineColor(kRed+1);
    gIC_bot->Draw("same PE");
    
    leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
    leg->AddEntry(gIC_top, "top SiPMs", "lp" );
    leg->AddEntry(gIC_bot, "bottom SiPMs", "lp" );
    leg->Draw();
    gPad->SetGridy();
    
    
    TCanvas * cIC_RMS = new TCanvas ("cIC_RMS", "cIC_RMS", 1000, 500);
    cIC_RMS->Divide(3,1);
    cIC_RMS->cd(1);
    hIC_tot->Draw();    
    hIC_tot->GetXaxis()->SetTitle("IC (all)");
    hIC_tot->SetLineColor(kBlack);
    hIC_tot->SetFillColor(kBlack);
    hIC_tot->SetFillStyle(3001);
    gPad->SetGridx();
    
    cIC_RMS->cd(2);
    hIC_top->Draw();    
    hIC_top->GetXaxis()->SetTitle("IC (top)");
    hIC_top->SetLineColor(kBlue+1);
    hIC_top->SetFillColor(kBlue+1);
    hIC_top->SetFillStyle(3001);
    gPad->SetGridx();
    
    cIC_RMS->cd(3);
    hIC_bot->Draw();    
    hIC_bot->GetXaxis()->SetTitle("IC (bottom)");
    hIC_bot->SetLineColor(kRed+1);
    hIC_bot->SetFillColor(kRed+1);
    hIC_bot->SetFillStyle(3001);
    gPad->SetGridx();
    

    
    //2D energy vs position plots
    TCanvas * cArrayEnergy_vsXY = new TCanvas ("cArrayEnergy_vsXY", "cArrayEnergy_vsXY", 1600, 300);
    cArrayEnergy_vsXY->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEnergy_vsXY->cd(chId+1);
        pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("COLZ");
        pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("X [mm]");
        pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Y [mm]");
        pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetZaxis()->SetTitle("Energy [a.u.]");
    }
    
    TCanvas * cEnergy_vsXYZoom = new TCanvas ("cEnergy_vsXYZoom", "cEnergy_vsXYZoom", 500, 800);
    cEnergy_vsXYZoom->Divide(1, 2);
    cEnergy_vsXYZoom->cd(1);
    pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("COLZ");
    cEnergy_vsXYZoom->cd(2);
    pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("COLZ");
    
    
    //1D energy vs position plots
    TCanvas * cArrayEnergy_vsPos = new TCanvas ("cArrayEnergy_vsPos", "cArrayEnergy_vsPos", 1600, 300);
    cArrayEnergy_vsPos->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEnergy_vsPos->cd(chId+1);
        TProfile * this_profile;
        std::string axis_label;
        if (HORIZONTAL)
        {
            axis_label = "X [mm]";
            this_profile = pBarResponse_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]];
        }
        else 
        {
            axis_label = "Y [mm]";
            this_profile = pBarResponse_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]];
        }        
        this_profile->Draw();
        this_profile->GetXaxis()->SetRangeUser(0,100);
        this_profile->GetXaxis()->SetTitle(axis_label.c_str());        
        this_profile->GetYaxis()->SetTitle("Energy [a.u.]");            
    }
    
    TCanvas * cEnergy_vsPosZoom = new TCanvas ("cEnergy_vsPosZoom", "cEnergy_vsPosZoom", 500, 800);
    cEnergy_vsPosZoom->Divide(1, 2);
    
    TProfile * this_profile_bot;
    TProfile * this_profile_top;
    if (HORIZONTAL)
    {
        this_profile_top = pBarResponse_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]];
        this_profile_bot = pBarResponse_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]];
    }
    else
    {
        this_profile_top = pBarResponse_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]];
        this_profile_bot = pBarResponse_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]];
    }
    cEnergy_vsPosZoom->cd(1);
    this_profile_top->Draw();    
    cEnergy_vsPosZoom->cd(2);
    this_profile_bot->Draw();
    
    
    //2D efficiency vs position plots
    TCanvas * cArrayEfficiency_vsXY = new TCanvas ("cArrayEfficiency_vsXY", "cArrayEfficiency_vsXY", 1600, 300);
    cArrayEfficiency_vsXY->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEfficiency_vsXY->cd(chId+1);
        hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("COLZ");
//  	pBarResponse_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRangeUser(0,100);
        hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("X [mm]");
        hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Y [mm]");
        hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetZaxis()->SetTitle("Efficiency [a.u.]");
        hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetZaxis()->SetRangeUser(0,1);
    }
    
    TCanvas * cEfficiency_vsXYZoom = new TCanvas ("cEfficiency_vsXYZoom", "cEfficiency_vsXYZoom", 500, 800);
    cEfficiency_vsXYZoom->Divide(1, 2);
    cEfficiency_vsXYZoom->cd(1);
    hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("COLZ");    
    cEfficiency_vsXYZoom->cd(2);
    hEff_vs_XY[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("COLZ");

    
    
    //1D efficiency vs position plots
    TCanvas * cArrayEfficiency_vsPos = new TCanvas ("cArrayEfficiency_vsPos", "cArrayEfficiency_vsPos", 1600, 300);
    cArrayEfficiency_vsPos->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayEfficiency_vsPos->cd(chId+1);
        TH1F * this_histo;
        std::string axis_label;
        if (HORIZONTAL)
        {
            axis_label = "X [mm]";
            this_histo = hEff_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]];
        }
        else 
        {
            axis_label = "Y [mm]";
            this_histo = hEff_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]];
        }        
        this_histo->Draw();
        this_histo->GetXaxis()->SetRangeUser(0,100);
        this_histo->GetXaxis()->SetTitle(axis_label.c_str());        
        this_histo->GetYaxis()->SetTitle("Efficiency [a.u.]");            
    }
    
    TCanvas * cEfficiency_vsPosZoom = new TCanvas ("cEfficiency_vsPosZoom", "cEfficiency_vsPosZoom", 500, 800);
    cEfficiency_vsPosZoom->Divide(1, 2);    
    TH1F * this_histo_bot;
    TH1F * this_histo_top;    
    if (HORIZONTAL)
    {
        this_histo_top = hEff_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]];
        this_histo_bot = hEff_vs_X[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]];
    }
    else
    {
        this_histo_top = hEff_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]];
        this_histo_bot = hEff_vs_Y[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]];
    }
    cEfficiency_vsPosZoom->cd(1);
    this_histo_top->Draw();    
    cEfficiency_vsPosZoom->cd(2);
    this_histo_bot->Draw();

    
    
//     double center[NCH] = {0};

    
    
 //************************************************************************************//
 //              loop 1
 //************************************************************************************//
    
    
    
    
//	 calculate the cross talk, and plot this
// 	 plot Tot, normalized by MIP peak energy, and then as a fraction of the total energy deposited in all channels
 
    std::cout << "(1) looping over events to calculate cross talk and time resolution" << std::endl;
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {
     
            
        tree->GetEntry(iEvt);
        if (ntracks != 1) continue; // require 1 MIP per event
            
        if (iEvt%1000 == 0) std::cout << "processing event: " << iEvt << "\r" << std::flush;
        //	if (step1!=6) continue;
        
        
        if (std::find(trk_blacklist->begin(), trk_blacklist->end(), double(x_dut)) != trk_blacklist->end() ) continue;
        
//         long long time_ref = time[384];
               
        //loop over the bars within the array to calculate cross-talk
        for (int chId = 0; chId<NBARS*2; chId++)
        {
            
            if (!passBugEventCuts(myChList[chId], energy, qfine, tot)) continue;
            if (!neighborsPassBugEventCuts(chId, false, energy, qfine, tot, myChList)) continue;
         
             //calculate cross-talk
             double in_me     = getMyLightFraction(chId,  0, false, energy, IC, myChList);
             double in_my_ln  = getMyLightFraction(chId, -1, false, energy, IC, myChList);
             double in_my_rn  = getMyLightFraction(chId,  1, false, energy, IC, myChList);
             
//              std::cout << "in_me[" << chId << "] = " << in_me << std::endl;
             hXT[step1][step2][myChList[chId]]->Fill(in_me);
             
             if (energy[myChList[chId]] > MIP_peak[chId]*0.8 && energy[myChList[chId]] < MIP_peak[chId]*5)
             {                
                hXT_cut[step1][step2][myChList[chId]]->Fill(in_me);
                hXT_cut_ln[step1][step2][myChList[chId]]->Fill(in_my_ln);
                hXT_cut_rn[step1][step2][myChList[chId]]->Fill(in_my_rn);
             }                          
        }
        
        
        //loop over the bars within the array to calculate cross-talk
        for (int iBar = 0; iBar<NBARS; iBar++)
        {
            
            if (!passBugEventCuts(myChList[iBar], energy, qfine, tot)) continue;    // check that top SiPM is ok
            if (!passBugEventCuts(myChList[iBar+16], energy, qfine, tot)) continue; // check that bottom SiPM is ok
            
            //select mip peak in both sipms
            if (   energy[myChList[iBar]] > MIP_peak[iBar]*0.8 && energy[myChList[iBar]] < MIP_peak[iBar]*5
                && energy[myChList[iBar+16]] > MIP_peak[iBar+16]*0.8 && energy[myChList[iBar+16]] < MIP_peak[iBar+16]*5
            )
	    {
                hCTR_UD[step1][step2][iBar]->Fill(time[myChList[iBar]] - time[myChList[iBar+16]]);                        
            }                         
        }
    }
 
 
 
    double XT[NBARS*2] = {0};
    double XT_err[NBARS*2] = {0};
    double XT_ln[NBARS*2] = {0};
    double XT_err_ln[NBARS*2] = {0};
    double XT_rn[NBARS*2] = {0};
    double XT_err_rn[NBARS*2] = {0};
    
    double XT_norm[NBARS*2] = {0};
    double avgXT = 0;        
    int XT_count = 0;
    
    TH1F * hXT_norm_tot = new TH1F ("hXT_norm_tot", "hXT_norm_tot", 20, 0.7, 1.3);
    TH1F * hXT_norm_top = new TH1F ("hXT_norm_top", "hXT_norm_top", 20, 0.7, 1.3);
    TH1F * hXT_norm_bot = new TH1F ("hXT_norm_bot", "hXT_norm_bot", 20, 0.7, 1.3);
    
    TCanvas * cArrayXT = new TCanvas ("cArrayXT", "cArrayXT", 1600, 300);
    cArrayXT->Divide(NBARS, 2);
    for (int chId = 0; chId<NBARS*2; chId++)
    {
        cArrayXT->cd(chId+1);
        hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw();
        hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetRangeUser(0,100);
        hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetXaxis()->SetTitle("XT");
        hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetYaxis()->SetTitle("Counts");
        
        hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(kGreen+2);
        hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetMarkerColor(kGreen+2);
        hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
        
        hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(kRed+1);
        hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetMarkerColor(kRed+1);
        hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");
        
        hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetLineColor(kYellow+2);
        hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->SetMarkerColor(kYellow+2);
        hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Draw("same");

        //fitting histos
        //main bar
        float maxBin    = hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetMaximumBin();
        float maxCenter = hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinCenter(maxBin);        
        TF1 * fitGaus = new TF1 ("fitGaus","gaus", 0.9*maxCenter, 1);
        fitGaus->SetLineColor(kGreen+1);
        hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitGaus, "QRL");
        XT[chId] = fitGaus->GetParameter(1);        
        XT_err[chId] = fitGaus->GetParError(1);
        
        //left bar
        maxBin    = hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetMaximumBin();
        maxCenter = hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinCenter(maxBin);
        fitGaus = new TF1 ("fitGaus","gaus", 0.7*maxCenter, 1.5*maxCenter);
        fitGaus->SetLineColor(kRed+2);
        hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitGaus, "QRL");
        XT_ln[chId] = fitGaus->GetParameter(1);        
        XT_err_ln[chId] = fitGaus->GetParError(1);

        //right bar
        maxBin    = hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetMaximumBin();
        maxCenter = hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->GetBinCenter(maxBin);
        fitGaus = new TF1 ("fitGaus","gaus", 0.7*maxCenter, 1.5*maxCenter);
        fitGaus->SetLineColor(kYellow+1);
        hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChList[chId]]->Fit(fitGaus, "QRL");
        XT_rn[chId] = fitGaus->GetParameter(1);        
        XT_err_rn[chId] = fitGaus->GetParError(1);
        
        //counting good fits
        if (XT[chId]>0) 
        {
            avgXT += XT[chId];
            XT_count++;           
        }
    }
    
    TCanvas * cXTZoom = new TCanvas ("cXTZoom", "cXTZoom", 500, 800);
    cXTZoom->Divide(1, 2);
    cXTZoom->cd(1);
    hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw();
    hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("same");
    hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("same");
    hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChTop[selBarNumber]]->Draw("same");
    cXTZoom->cd(2);
    hXT[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw();
    hXT_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("same");
    hXT_cut_ln[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("same");
    hXT_cut_rn[step1_vct.at(selStep1)][step2_vct.at(selStep2)][myChBot[selBarNumber]]->Draw("same");
            
    
    
    avgXT /= XT_count; 
    TGraphErrors * gXT_top = new TGraphErrors();
    TGraphErrors * gXT_bot = new TGraphErrors();
    TGraphErrors * gXT_top_ln = new TGraphErrors();
    TGraphErrors * gXT_bot_ln = new TGraphErrors();
    TGraphErrors * gXT_top_rn = new TGraphErrors();
    TGraphErrors * gXT_bot_rn = new TGraphErrors();
    
    TGraphErrors * gXT_norm_top = new TGraphErrors();
    TGraphErrors * gXT_norm_bot = new TGraphErrors();
        
    for (int chId = 0; chId<NBARS*2; chId++) 
    {
        XT_norm[chId] = XT[chId] / avgXT;
        if (chId<NBARS) 
        {
            gXT_top->SetPoint(chId, chId, XT[chId]);
            gXT_top->SetPointError(chId, 0, XT_err[chId]);
            if(XT_ln[chId]>0 && XT_err_ln[chId] < 0.2)
            {
                gXT_top_ln->SetPoint(chId, chId, XT_ln[chId]);
                gXT_top_ln->SetPointError(chId, 0, XT_err_ln[chId]);
            }
            if(XT_rn[chId]>0 && XT_err_rn[chId] < 0.2)
            {
                gXT_top_rn->SetPoint(chId, chId, XT_rn[chId]);
                gXT_top_rn->SetPointError(chId, 0, XT_err_rn[chId]);
            }            
            gXT_norm_top->SetPoint(chId, chId, XT_norm[chId]);
            gXT_norm_top->SetPointError(chId, 0, XT_err[chId]/avgXT);
            hXT_norm_top->Fill(XT_norm[chId]);
            
        }
        else
        {
            gXT_bot->SetPoint(chId-NBARS, chId-NBARS, XT[chId]);
            gXT_bot->SetPointError(chId-NBARS, 0, XT_err[chId]);
            if(XT_ln[chId]>0 && XT_err_ln[chId] < 0.2)
            {
                gXT_bot_ln->SetPoint(chId-NBARS, chId-NBARS, XT_ln[chId]);
                gXT_bot_ln->SetPointError(chId-NBARS, 0, XT_err_ln[chId]);
            }
            if(XT_rn[chId]>0 && XT_err_rn[chId] < 0.2)
            {
                gXT_bot_rn->SetPoint(chId-NBARS, chId-NBARS, XT_rn[chId]);
                gXT_bot_rn->SetPointError(chId-NBARS, 0, XT_err_rn[chId]);
            }            
            gXT_norm_bot->SetPoint(chId-NBARS, chId-NBARS, XT_norm[chId]);
            gXT_norm_bot->SetPointError(chId-NBARS, 0, XT_err[chId]/avgXT);
            hXT_norm_bot->Fill(XT_norm[chId]);
        }
        hXT_norm_tot->Fill(XT_norm[chId]);
    }
        
    
    TCanvas * cXT_abs = new TCanvas ("cXT_abs", "cXT_abs", 600, 500);
    cXT_abs->cd();
    gXT_top->Draw("APE");
    gXT_top->GetYaxis()->SetRangeUser(0.5, 1.);
    gXT_top->GetYaxis()->SetTitle("Fraction of light within bar (XT)");
    gXT_top->SetMarkerStyle(20);
    gXT_top->SetMarkerColor(kBlue+1);
    gXT_top->SetLineColor(kBlue+1);
    
    gXT_bot->SetMarkerStyle(21);
    gXT_bot->SetMarkerColor(kRed+1);
    gXT_bot->SetLineColor(kRed+1);
    gXT_bot->Draw("same PE");
    
    leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
    leg->AddEntry(gXT_top, "top SiPMs", "lp" );
    leg->AddEntry(gXT_bot, "bottom SiPMs", "lp" );
    leg->Draw();
    gPad->SetGridy();
    
    TCanvas * cXT_abs_nb = new TCanvas ("cXT_abs_nb", "cXT_abs_nb", 600, 500);
    cXT_abs_nb->cd();
    gXT_top_ln->Draw("APE");
    gXT_top_ln->GetYaxis()->SetRangeUser(0., 0.2);
    gXT_top_ln->GetYaxis()->SetTitle("Fraction of light inside each neighboring bar");    
    gXT_top_ln->SetMarkerStyle(20);
    gXT_top_ln->SetMarkerColor(kBlue+1);
    gXT_top_ln->SetLineColor(kBlue+1);
    
    gXT_bot_ln->SetMarkerStyle(20);
    gXT_bot_ln->SetMarkerColor(kRed+1);
    gXT_bot_ln->SetLineColor(kRed+1);
    gXT_bot_ln->Draw("same PE");
    
    gXT_top_rn->SetMarkerStyle(24);
    gXT_top_rn->SetMarkerColor(kBlue+1);
    gXT_top_rn->SetLineColor(kBlue+1);
    gXT_top_rn->Draw("same PE");
    
    gXT_bot_rn->SetMarkerStyle(24);
    gXT_bot_rn->SetMarkerColor(kRed+1);
    gXT_bot_rn->SetLineColor(kRed+1);
    gXT_bot_rn->Draw("same PE");
    
    leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
    leg->AddEntry(gXT_top_ln, "top SiPMs (left bar)", "lp" );
    leg->AddEntry(gXT_top_rn, "top SiPMs (right bar)", "lp" );
    leg->AddEntry(gXT_bot_ln, "bot SiPMs (left bar)", "lp" );
    leg->AddEntry(gXT_bot_rn, "bot SiPMs (right bar)", "lp" );
    
    leg->Draw();
    gPad->SetGridy();
    
    
    
    
    TCanvas * cXT_norm = new TCanvas ("cXT_norm", "cXT_norm", 600, 500);
    cXT_norm->cd();
    gXT_norm_top->Draw("APE");
    gXT_norm_top->GetYaxis()->SetRangeUser(0.7, 1.3);
    gXT_norm_top->GetYaxis()->SetTitle("Fraction of light within bar / avg XT");
    gXT_norm_top->SetMarkerStyle(20);
    gXT_norm_top->SetMarkerColor(kBlue+1);
    gXT_norm_top->SetLineColor(kBlue+1);
    
    gXT_norm_bot->SetMarkerStyle(21);
    gXT_norm_bot->SetMarkerColor(kRed+1);
    gXT_norm_bot->SetLineColor(kRed+1);
    gXT_norm_bot->Draw("same PE");
    
    leg = new TLegend(0.15,0.75,0.5,0.88,NULL,"brNDC");
    leg->AddEntry(gXT_norm_top, "top SiPMs", "lp" );
    leg->AddEntry(gXT_norm_bot, "bottom SiPMs", "lp" );
    leg->Draw();
    gPad->SetGridy();
    
    
    TCanvas * cXT_norm_RMS = new TCanvas ("cXT_norm_RMS", "cXT_norm_RMS", 1000, 500);
    cXT_norm_RMS->Divide(3,1);
    cXT_norm_RMS->cd(1);
    hXT_norm_tot->Draw();    
    hXT_norm_tot->GetXaxis()->SetTitle("XT_norm (all)");
    hXT_norm_tot->SetLineColor(kBlack);
    hXT_norm_tot->SetFillColor(kBlack);
    hXT_norm_tot->SetFillStyle(3001);
    gPad->SetGridx();    
    cXT_norm_RMS->cd(2);
    hXT_norm_top->Draw();    
    hXT_norm_top->GetXaxis()->SetTitle("XT_norm (top)");
    hXT_norm_top->SetLineColor(kBlue+1);
    hXT_norm_top->SetFillColor(kBlue+1);
    hXT_norm_top->SetFillStyle(3001);
    gPad->SetGridx();
    
    cXT_norm_RMS->cd(3);
    hXT_norm_bot->Draw();    
    hXT_norm_bot->GetXaxis()->SetTitle("XT_norm (bottom)");
    hXT_norm_bot->SetLineColor(kRed+1);
    hXT_norm_bot->SetFillColor(kRed+1);
    hXT_norm_bot->SetFillStyle(3001);
    gPad->SetGridx();
    
    
    

    //do timing studies

    std::cout << "time resolution" << std::endl;
 
    double CTR[NBARS] = {0};
    double CTR_err[NBARS] = {0};
    double CTR_norm[NBARS] = {0};
    double avgCTR = 0; 
    int CTR_count = 0;
    
    TH1F * hCTR_norm = new TH1F ("hCTR_norm", "hCTR_norm", 20, 0.7, 1.3);
    
    TCanvas * cArrayCTR_UD = new TCanvas ("cArrayCTR_UD", "cArrayCTR_UD", 1600, 300);
    cArrayCTR_UD->Divide(NBARS);
    for (int iBar = 0; iBar<NBARS; iBar++)
    {
        cArrayCTR_UD->cd(iBar+1);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->Draw();
 	hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->GetXaxis()->SetRangeUser(-1000,1000);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->GetXaxis()->SetTitle("Energy [a.u.]");
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->GetYaxis()->SetTitle("Counts");
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->Rebin(REBIN_COEFF);
        
        //fitting histo
        float maxBin    = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->GetMaximumBin();
        float maxCenter = hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->GetBinCenter(maxBin);        
        TF1 * fitGaus = new TF1 ("fitGaus","gaus", maxCenter-200, maxCenter+200);
//         fitGaus->SetLineColor(kGreen+1);
        hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][iBar]->Fit(fitGaus, "QRL");
        CTR[iBar] = fitGaus->GetParameter(2);        
        CTR_err[iBar] = fitGaus->GetParError(2);   
        
        if (CTR[iBar]>0)
        {
            avgCTR += CTR[iBar];
            CTR_count++;
        }
    }
    
        
    
    TCanvas * cCTR_UDZoom = new TCanvas ("cCTR_UDZoom", "cCTR_UDZoom", 500, 500);    
    cCTR_UDZoom->cd();
    hCTR_UD[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarNumber]->Draw();
//     hCTR_UD_cut[step1_vct.at(selStep1)][step2_vct.at(selStep2)][selBarNumber]->Draw("same");

 
    avgCTR /= CTR_count; 
    TGraphErrors * gCTR_UD      = new TGraphErrors();    
    TGraphErrors * gCTR_UD_norm = new TGraphErrors();    
    
    for (int iBar = 0; iBar<NBARS; iBar++) 
    {
        CTR_norm[iBar] = CTR[iBar] / avgCTR;
        
        if (CTR[iBar]>0)
        {
            gCTR_UD->SetPoint(iBar, iBar, CTR[iBar]);
            gCTR_UD->SetPointError(iBar, 0, CTR_err[iBar]);
        
            gCTR_UD_norm->SetPoint(iBar, iBar, CTR_norm[iBar]);
            gCTR_UD_norm->SetPointError(iBar, 0, CTR_err[iBar]/avgCTR);
            hCTR_norm->Fill(IC[iBar]);
        }
                
    }
    
    TCanvas * cCTR_UD_abs = new TCanvas ("cCTR_UD_abs", "cCTR_UD_abs", 600, 500);
    cCTR_UD_abs->cd();
    gCTR_UD->Draw("APE");
    gCTR_UD->GetYaxis()->SetRangeUser(0, 300);
    gCTR_UD->GetYaxis()->SetTitle("CTR [ps]");
    gCTR_UD->SetMarkerStyle(20);
    gCTR_UD->SetMarkerColor(kBlue+1);
    gCTR_UD->SetLineColor(kBlue+1);    
    gPad->SetGridy();
     
    TCanvas * cCTR_UD_norm = new TCanvas ("cCTR_UD_norm", "cCTR_UD_norm", 600, 500);
    cCTR_UD_norm->cd();
    gCTR_UD_norm->Draw("APE");
    gCTR_UD_norm->GetYaxis()->SetRangeUser(0, 2);
    gCTR_UD_norm->GetYaxis()->SetTitle("CTR_{i} / CTR_{ave}");
    gCTR_UD_norm->SetMarkerStyle(20);
    gCTR_UD_norm->SetMarkerColor(kBlue+1);
    gCTR_UD_norm->SetLineColor(kBlue+1);    
    gPad->SetGridy();

    TCanvas * cCTR_norm_RMS = new TCanvas ("cCTR_norm_RMS", "cCTR_norm_RMS", 1000, 500);    
    cCTR_norm_RMS->cd();
    hCTR_norm->Draw();    
    hCTR_norm->GetXaxis()->SetTitle("CTR_{i} / CTR_{ave}");
    hCTR_norm->SetLineColor(kGreen+1);
    hCTR_norm->SetFillColor(kGreen+1);
    hCTR_norm->SetFillStyle(3001);
    gPad->SetGridx(); 
    
//     //save plots or histos/graphs to output file
//     TFile * outputFile = new TFile(Form("./output/output_run_%.3d.root", firstRun), "RECREATE");
//     outputFile->cd();    
//     
//     outputFile->Write();
//     outputFile->Close();    
//     
    
    std::cout << "end of program" << std::endl;
    theApp->Run();
    myfile.close();

}
