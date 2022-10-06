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
int parseInputFile(std::vector<isolate*> *pop, std::vector<cog*> *accessoryLoci, double lower, double upper, std::vector<int> *samplingList, std::vector<std::string> *serotypeList, std::vector<int> *scList, std::vector<std::string> *cogList, char *inputFilename, char * vtCogName,int&, bool useCogList);

// check input values
bool checkInputValues(struct parms *sp, char* inputFilename, char* vtCogName, char* propFile, char *weightFile);

// parse marker file
int parseMarkerFile(std::vector<isolate*> *pop,char *markerFilename,std::vector<std::string> *markerList, bool useMarkerList);

// validate input files
int compareInputPopulations(std::vector<isolate*> *popA, std::vector<isolate*> *popB, bool check_markers);

// parse frequency file
int parseFrequencyFile(char *frequencyFilename,std::vector<cog*> *accessoryLoci);

// parse weighting file
int parseWeightingFile(char* weightingFilename,std::vector<cog*> *accessoryLoci);

// parse ordering file
int parseOrderingFile(char* orderingFilename,std::vector<cog*> *accessoryLoci,struct parms *sp, char* outputFilename);

////////////////////////////////
// Pre-processing information //
////////////////////////////////

// picking out valid strains for migration
std::vector<int> getValidStrains(std::vector<std::vector<isolate*> > migrantInput);

// overall population division function
int generateMigrantPool(std::vector<std::vector<std::vector<isolate*> > > *migrantPool, std::vector<isolate*> *population, std::vector<isolate*> *migrant_population, char* migrantFilename, std::vector<int> *scList, int maxScNum, int minGen,struct parms *p);

// diving population for immigration
int dividePopulationForImmigration(std::vector<isolate*> *pop,std::vector<int> *scList,std::vector<std::vector<isolate*> > *popBySc,int maxScNum);

// dividing population for immigration by time
int dividePopulationForImmigrationByTime(std::vector<isolate*> *pop, int minGen, int numgen,std::vector<std::vector<isolate*> > *popByTime);

// get isolates for starting population
int getStartingIsolates(std::vector<isolate*> *pop,struct parms *sp,std::vector<isolate*> *first,std::vector<cog*> *accessoryLoci, int,std::vector<double> &eqFreq,std::vector<double> &cogWeights,std::vector<double> &cogDeviations,std::vector<int> &startingVtScFrequencies,std::vector<int> &startingNvtScFrequencies,std::vector<int> *scList, int minGen, float seedStartingPopulation, char* migrantFilename, std::vector<isolate*> *migrant_population, int maxScNum);

// sample selected from first generation
int firstSample(std::vector<isolate*> *currentIsolates,int firstSample,std::ofstream& sampleOutFile,int minGen);

// summarise statistics per generation
int summariseGeneration(std::vector<isolate*> *pop,int sampleSize,std::vector<int> *scs,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<std::string> *seros,std::vector<int> *serotypeF);

// calculate COG deviations
int getCogDeviations(std::vector<double> * ef,std::vector<isolate*> *currentStrains,std::vector<double> *cogWeights,std::vector<double> *cogDeviations,std::vector<double> *startingCogFrequencies,std::vector<double> &startingScVtFrequencies,std::vector<double> &startingScNvtFrequencies,int);

// alter vaccine formulation
int alterVaccineFormulation(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<std::vector<std::vector<isolate*> > > *populationBySc);

// parse disease data
int parse_disease_data(char* epiFilename,
                       std::vector<int> *diseaseTime,
                       std::vector<std::string> *diseaseSeroList,
                       std::vector<int> *diseaseScList,
                       std::vector<int> *diseaseVt,
                       std::vector<double> *diseaseInvasiveness,
                       std::vector<int> *diseasePopulation,
                       std::vector<int> *diseaseCount);

//////////////////////////////
// Per-generation functions //
//////////////////////////////

// select next generation
int reproduction(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<std::vector<std::vector<isolate*> > > *migrantPool,std::vector<double> *cogWeights,std::vector<double> *cogDeviations,struct parms *sp, std::vector<double> * ef, std::vector<int> * vtScFreq,std::vector<int> * nvtScFreq,std::vector<double> * piGen,std::vector<int> *scList, int gen,std::vector<double> * timeGen,std::vector<double> * fitGen,std::vector<std::string> * isolateGen,std::vector<int> * countGen, double popLimitFactor, int minGen, int secondVaccinationGeneration, float partialVaccine);

// recombination
int recombination(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,char* markerFilename,double transformationProportion,double transformationRate,double transformationAsymmetryLoci, double transformationAsymmetryMarker);

// locus frequency updating after recombination
int update_locus_freq(std::vector<isolate*> *futureIsolates, std::vector<double> *cogWeights, std::vector<double> *cogDeviations, std::vector<double> *ef);

// move isolate into next generation
int nextGeneration(std::vector<isolate*> *pop,std::vector<isolate*> *new_pop,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates, std::vector<std::vector<std::vector<isolate*> > > *migrantPool);

// update population
int updatePopulation(std::vector<isolate*> *pop,std::vector<isolate*> *new_pop);

// compare samples
int compareSamples(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<cog*> *accessoryLoci,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<int> &sampledSeroFreq,std::vector<std::string> &serotypeList,std::vector<double> &vtCogFittingStatsList,std::vector<double> &nvtCogFittingStatsList,std::vector<double> &strainFittingStatsList,std::ofstream& sampleOutFile,struct parms *sp);

// compare to disease samples
int compare_to_disease_data(std::vector<double> &diseaseDivergence,
                            int simulation_time,
                            std::vector<isolate*> *currentIsolates,
                            std::vector<int> *diseaseTime,
                            std::vector<std::string> *diseaseSeroList,
                            std::vector<int> *diseaseScList,
                            std::vector<int> *diseaseVt,
                            std::vector<double> *diseaseInvasiveness,
                            std::vector<int> *diseasePopulation,
                            std::vector<int> *diseaseCount,
                            std::ofstream& diseaseOutFile);

// just record sample statistics in pure simulation mode
int justRecordStats(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<cog*> *accessoryLoci);

// calculate Pearson correlation
double pearson(std::vector<double> *x,std::vector<double> *y);

//////////////////////
// Output functions //
//////////////////////

// calculate reproductive fitness differences
int rFitMetricCalculation(int minGen,std::vector<int> *samplingList,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<isolate*> *population,std::vector<double> &rFitVector);

// print output
int printOutput(char* outputFilename,std::vector<std::string> *seroList,std::vector<std::vector<int> > &sampledSeroFreq,std::vector<int> *scList,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,int gen,int minGen,std::vector<cog*> *accessoryLoci,std::vector<int> *samplingList,std::vector<std::vector<double> > &piGen,struct parms *sp,std::vector<double> * timeGen,std::vector<double> * fitGen,std::vector<std::string> * isolateGen,std::vector<int> * countGen);

// print populations
int printPop(char* outputFilename,std::string suffix,std::vector<isolate*> *currentIsolates, char *markerFilename, std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList);

// print population genotype sample
int printPopSample (char* prefixStar,std::string suffix,std::vector<isolate*> *currentIsolates,char* markerFilename,std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList, int sampleSize);

// tidy up isolates
void tidyUpIsolates(std::vector<isolate*> *a_list, std::vector<isolate*> *b_list, std::vector<isolate*> *c_list, std::vector<isolate*> *d_list);

// tidy up isolates
void tidyUpLoci(std::vector<cog*> *cog_list);

#endif /* defined(__frequencyDependentSimulation__functions__) */
