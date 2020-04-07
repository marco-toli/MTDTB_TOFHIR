#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>
#include <algorithm>
#include "TMath.h"


bool passBugEventCuts (int iCh, float energy[256], float qfine[256], float tot[256]);

bool neighborsPassBugEventCuts (int chId, bool extended, float energy[256], float qfine[256], float tot[256], int myChList[32]);

double getMyLightFraction(int chId, int neighbor, bool extended, float energy[256], double IC[32], int myChList[32]);


double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity);
double FindSmallestInterval(std::vector<double>* vals, const double fraction, const bool verbosity);


double poissonf(double*x,double*par);




double fitBarEffErr(double* x, double* par);
