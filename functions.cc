#include "functions.hh"
#include "TMath.h"
#include <iostream>


bool passBugEventCuts (int iCh, float energy[256], float qfine[256], float tot[256])
{
    if (   tot[iCh]    < -100 || energy[iCh] < -100
        || qfine[iCh]  < 15   || qfine[iCh]  > 500 
        || energy[iCh] == 0   || energy[iCh] == 20 
        || energy[iCh] == 40  || energy[iCh] == 60
        )      return false; 
    
    else       return true;
    
}

bool neighborsPassBugEventCuts (int chId, bool extended, float energy[256], float qfine[256], float tot[256], int myChList[32])
{
    
    bool passedExtended = true;
    bool passedCompact  = true;
//     std::cout << "testing neighbors... [" << chId << "]" << std::endl;
//     
    int myBar;
    if (chId < 16) myBar = chId;
    else           myBar = chId-16;
                                
    if (myBar>0 && !passBugEventCuts(myChList[chId-1], energy, qfine, tot) ) 
    {
        passedExtended = false;
        passedCompact  = false;
    }
//     std::cout << "passed ln test..." << std::endl;
    if (myBar>1 && !passBugEventCuts(myChList[chId-2], energy, qfine, tot) ) 
    {
        passedExtended = false;                
    }
//     std::cout << "passed lln test..." << std::endl;
    if (myBar<15 && !passBugEventCuts(myChList[chId+1], energy, qfine, tot) ) 
    {
        passedExtended = false;
        passedCompact  = false;
    }
//     std::cout << "passed rn test..." << std::endl;
    if (myBar<14 && !passBugEventCuts(myChList[chId+2], energy, qfine, tot) ) 
    {
        passedExtended = false;                
    }
//     std::cout << "passed rnn test..." << std::endl;
    
    if (extended) {return passedExtended;}
    else          {return passedCompact;}
    
}

double getMyLightFraction(int chId, int neighbor, bool extended, float energy[256], double IC[32], int myChList[32])
{
    
    double this_bar_signal = 0.;
    double lnn = 0., ln = 0., rn = 0., rnn = 0.;
//     float IC_cut = 0.1;
//     std::cout << "entering function" << std::endl;
    
    int myBar;
    if (chId < 16) myBar = chId;     //top
    else           myBar = chId-16;  //bot
    
         
    if (IC[chId]>0 && energy[myChList[chId]] > 0)  this_bar_signal    = energy[myChList[chId]]   / IC[chId];
    if (myBar>0 && IC[chId-1]>0  && energy[myChList[chId-1]] > 0) ln  = energy[myChList[chId-1]] / IC[chId-1];
    if (myBar>1 && IC[chId-2]>0  && energy[myChList[chId-2]] > 0) lnn = energy[myChList[chId-2]] / IC[chId-2];
    if (myBar<15 && IC[chId+1]>0 && energy[myChList[chId+1]] > 0) rn  = energy[myChList[chId+1]] / IC[chId+1];
    if (myBar<14 && IC[chId+2]>0 && energy[myChList[chId+2]] > 0) rnn = energy[myChList[chId+2]] / IC[chId+2];
//         std::cout << "chId [" << chId << "] --> this_bar = " << this_bar_signal << " :: ln = " << ln << " :: lnn = " << lnn << " :: rn = " << rn << " :: rnn = " << rnn << std::endl;
    
    
    
    double tot;
    if (extended) tot = this_bar_signal+ln+lnn+rn+rnn;  //use sum of 5 bars as denominator
    else          tot = this_bar_signal+ln+rn;          //use sum of 3 bars as denominator
    
    double frac;
    if (tot!= 0) 
    {
        if (neighbor == 0)  frac = this_bar_signal/tot; // fraction of light in central bar
        if (neighbor == -1) frac = ln/tot;              // fraction of light in -1 bar
        if (neighbor == -2) frac = lnn/tot;             // fraction of light in -2 bar
        if (neighbor == 1)  frac = rn/tot;              // fraction of light in +1 bar
        if (neighbor == 2)  frac = rnn/tot;             // fraction of light in +2 bar
    }
    else         frac = 0;
    

    if (frac>0.3 && frac<0.45)
    {
//         std::cout << "IC = [" << IC[chId-2] << ", " << IC[chId-1] << ", " << IC[chId] << ", " << IC[chId+1] << ", " << IC[chId+2] << "]" <<  std::endl;
        
        
    }
    
//     std::cout << "frac=" << frac << " :: energy["<<chId <<"] = [" << energy[myChList[chId-2]] << ", " << energy[myChList[chId-1]] << ", " << energy[myChList[chId]] << ", " << energy[myChList[chId+1]] << ", " << energy[myChList[chId+2]] << "]" <<  std::endl;

    
    return frac;    
    
    
}


double fitBarEffErr(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double t = x[0];
  
  double center		= par[0];
  double width 		= par[1];
  double baseline	= par[2];
  double max     	= par[3];
  double sigma_x   	= par[4];
  
  double f;
  double edge1 = center-width/2;
  double edge2 = center+width/2;
  
  if ( t <= center)
  {
      f = baseline + max/2*(erf( (t-edge1) / (sqrt(2)*sigma_x)) + 1);
  }
  if (t > center)
  {
      f = baseline + max/2*(erf((-t+edge2) / (sqrt(2)*sigma_x)) + 1);
  }
  
  
  return f;
}





// 
/*
        float in_lnb  = getL_NBFraction(chId);
        float in_rnb  = getR_NBFraction(chId);
        
        float in_lnnb = getL_NNBFraction(chId);
        float in_rnnb = getR_NNBFraction(chId);*/


double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  std::sort(vals->begin(),vals->end());

  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);

//   unsigned int minPoint = 0;
//   unsigned int maxPoint = 0;
  double delta = 999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
  {
    double tmpMin = vals->at(point);
    double tmpMax = vals->at(point+maxPoints-1);
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
      min = tmpMin;
      max = tmpMax;
    }
  }

  
  return delta;
  
}

double FindSmallestInterval(std::vector<double>* vals, const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  std::sort(vals->begin(),vals->end());

  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);

  double delta = 9999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
  {
    double tmpMin = vals->at(point);
    double tmpMax = vals->at(point+maxPoints-1);
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
    }
  }

  
  return delta;
  
}



double poissonf(double*x, double*par)                                         
{                                                    
//   Double_t res=0.;
  Double_t xx=x[0];
  if (xx<=0) return  0;

  // Poisson distribution
  // par[1] - distribution parameter
  return par[0]*TMath::Power(par[1],xx)/TMath::Gamma(xx+1)/TMath::Exp(par[1]);

//   return par[0]*TMath::Poisson(x[0],par[1]);
}
