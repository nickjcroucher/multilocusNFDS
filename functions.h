//
//  functions.h
//  frequencyDependentSimulation
//
//  Created by Nicholas Croucher on 28/09/2015.
//  Copyright (c) 2015 Imperial College. All rights reserved.
//

#ifndef __frequencyDependentSimulation__functions__
#define __frequencyDependentSimulation__functions__

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include <vector>
#include "parms.h"

///////////////////////
// Generic functions //
///////////////////////

// random seed
unsigned long int random_seed();
extern gsl_rng * rgen;

// usage message
void usage(char *);

////////////////////////
// Input file parsing //
////////////////////////

// parse input file
int parseInputFile(std::vector<isolate*> *pop, std::vector<cog*> *accessoryLoci, double lower, double upper, std::vector<int> &samplingList, std::vector<std::string> *serotypeList, std::vector<int> *scList, std::vector<std::string> *cogList, char *inputFilename, char * vtCogName,int&);

// check input values
bool checkInputValues(struct parms *sp, char* inputFilename, char* vtCogName, char* propFile, char *weightFile);

// parse marker file
int parseMarkerFile(std::vector<isolate*> *pop,char *markerFilename,std::vector<std::string> *markerList);

// parse frequency file
int parseFrequencyFile(char *frequencyFilename,std::vector<cog*> *accessoryLoci);

// parse weighting file
int parseWeightingFile(char* weightingFilename,std::vector<cog*> *accessoryLoci);

// parse ordering file
int parseOrderingFile(char* orderingFilename,std::vector<cog*> *accessoryLoci,struct parms *sp);

////////////////////////////////
// Pre-processing information //
////////////////////////////////

// diving population for immigration
int dividePopulationForImmigration(std::vector<isolate*> *pop,std::vector<int> *scList,std::vector<std::vector<isolate*> > *popBySc,int maxScNum);

// get isolates for starting population
int getStartingIsolates(std::vector<isolate*> *pop,std::vector<isolate*> *first,std::vector<cog*> *accessoryLoci, int,std::vector<double> &eqFreq,std::vector<double> &cogWeights,std::vector<double> &cogDeviations,std::vector<int> &startingVtScFrequencies,std::vector<int> &startingNvtScFrequencies,std::vector<int> &scList);

// sample selected from first generation
int firstSample(std::vector<isolate*> *currentIsolates,int firstSample,std::ofstream& sampleOutFile,int minGen);

// summarise statistics per generation
int summariseGeneration(std::vector<isolate*> *pop,int sampleSize,std::vector<int> *scs,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<std::string> *seros,std::vector<int> *serotypeF);

// calculate COG deviations
int getCogDeviations(std::vector<double> * ef,std::vector<isolate*> *currentStrains,std::vector<double> *cogWeights,std::vector<double> *cogDeviations,std::vector<double> *startingCogFrequencies,std::vector<double> &startingScVtFrequencies,std::vector<double> &startingScNvtFrequencies,int);

// alter vaccine formulation
int alterVaccineFormulation(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<std::vector<isolate*> > *populationBySc);

//////////////////////////////
// Per-generation functions //
//////////////////////////////

// select next generation
int reproduction(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop,std::vector<std::vector<isolate*> > *popBySc,std::vector<double> *cogWeights,std::vector<double> *cogDeviations,struct parms *sp, std::vector<double> * ef, std::vector<int> * vtScFreq,std::vector<int> * nvtScFreq,std::vector<double> * piGen,std::vector<int> *scList,int gen);

// recombination
int recombination(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop,char* markerFilename,double transformationProportion,double transformationRate,double transformationAsymmetryLoci, double transformationAsymmetryMarker, std::vector<double> *cogWeights, std::vector<double> *cogDeviations, std::vector<double> * ef);

// move isolate into next generation
int nextGeneration(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop);

// compare samples
int compareSamples(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<cog*> *accessoryLoci,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<int> &sampledSeroFreq,std::vector<std::string> &serotypeList,std::vector<double> &vtCogFittingStatsList,std::vector<double> &nvtCogFittingStatsList,std::vector<double> &strainFittingStatsList,std::ofstream& sampleOutFile);

// just record sample statistics in pure simulation mode
int justRecordStats(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<cog*> *accessoryLoci);

// calculate Pearson correlation
double pearson(std::vector<double> *x,std::vector<double> *y);

//////////////////////
// Output functions //
//////////////////////

// calculate reproductive fitness differences
int rFitMetricCalculation(int minGen,int maxScNum,int numGen,std::vector<int> &samplingList,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<isolate*> *population,std::vector<double> &rFitVector);

// print output
int printOutput(char* outputFilename,std::vector<std::string> *seroList,std::vector<std::vector<int> > &sampledSeroFreq,std::vector<int> *scList,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,std::vector<std::string> *cogList,std::vector<double> *cogWeights,int gen,int minGen,std::vector<cog*> *accessoryLoci,std::vector<int> samplingList,std::vector<std::vector<double> > &piGen,struct parms *sp);

// print populations
int printPop(char* outputFilename,std::string suffix,std::vector<isolate*> *currentIsolates, char *markerFilename, std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList);

// print population genotype sample
int printPopSample (char* prefixStar,std::string suffix,std::vector<isolate*> *currentIsolates,char* markerFilename,std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList, int sampleSize);

/////////
// OLD //
/////////

// calculate SC deviations
//int getStrainStats(std::vector<isolate*> *pop,std::vector<isolate*> *currentPop, int year, int numSc, double &popDev);


// calculate final COG frequencies
//int calculateFinalCogFrequency(std::vector<isolate*> *pop,std::vector<double> *endingCogFrequencies);







// calculate strain change statistics
//int getStrainChangeStats(std::vector<isolate*> *pop,int firstYear,int midYear,int finalYear,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,int mGen,int fGen,int numSc,double &changeDev);

// get COG frequencies
//int getCogFreq(std::vector<isolate*> *pop,int,std::vector<int> *freq,int &strains);

// get summary statistics
//int getStats(std::vector<double> *startingFreq,std::vector<double> *endingFreq,int vtCog,double &frequencyCorrelation,double &vtCogDecrease);

// find change in COG proportions
//int getProp(std::vector<isolate*> *pop,struct parms *sp,std::vector<int> *starting,std::vector<int> *middling,std::vector<int> *ending,std::vector<double> *ef,std::vector<double> *mf,std::vector<double> *ff,std::vector<std::string> *cogList,int &vc, int ns, int &ng,int ss, int ms, int es);





#endif /* defined(__frequencyDependentSimulation__functions__) */
