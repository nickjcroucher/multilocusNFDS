//
//  main.cpp
//  frequencyDependentSimulation
//
//  Created by Nicholas Croucher on 27/09/2015.
//  Copyright (c) 2015 Imperial College. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "functions.h"
#include "parms.h"

// define external function for random number generation
gsl_rng * rgen;

int main(int argc, char * argv[]) {

    ////////////////////////////////
    // Initialise data structures //
    ////////////////////////////////
    
    // Seeds the rng
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rgen = gsl_rng_alloc (T);
    gsl_rng_set(rgen, random_seed());
    
    // initiate parameterisation
    parms p;
    p.selectedProp = 1.0;
    p.lowerSelection = 0;
    p.higherSelection = 1;
    
    ////////////////////////
    // Parse command line //
    ////////////////////////
    
    // parse command line
    char* inputFilename = 0;
    char* outputFilename = 0;
    char* frequencyFilename = 0;
    char* weightingFilename = NULL;
    char* orderingFilename = NULL;
    char* vtCogName = 0;

    if (argc == 1) {
        usage(argv[0]);
        return 1;
    } else {
        int opt = 0;
        while ((opt = getopt(argc,argv,"hc:p:s:v:i:t:n:g:u:l:y:j:k:f:x:w:r:o:")) != EOF) {
            switch (opt) {
                case 'h':
                    usage(argv[0]);
                    return 0;
                case 'p':
                    p.programme = (*optarg);
                    break;
                case 's':
                    p.fSelection = atof(optarg);
                    break;
                case 'v':
                    p.vSelection = atof(optarg);
                    break;
                case 'c':
                    vtCogName = optarg;
                    break;
                case 'i':
                    p.immigrationRate = atof(optarg);
                    break;
                case 't':
                    p.immigrationType = atoi(optarg);
                    break;
                case 'n':
                    p.popSize = atoi(optarg);
                    break;
                case 'g':
                    p.numGen = atoi(optarg);
                    break;
                case 'u':
                    p.upperLimit = atof(optarg);
                    p.upperLimit = p.upperLimit + 0.000001;
                    break;
                case 'l':
                    p.lowerLimit = atof(optarg);
                    p.lowerLimit = p.lowerLimit - 0.0000001;
                    break;
                case 'y':
                    p.selectedProp = atof(optarg);
                    break;
                case 'j':
                    p.lowerSelection = atof(optarg);
                    break;
                case 'k':
                    p.higherSelection = atof(optarg);
                    break;
                case 'f':
                    inputFilename = optarg;
                    break;
                case 'x':
                    frequencyFilename = optarg;
                    break;
                case 'w':
                    weightingFilename = optarg;
                    break;
                case 'r':
                    orderingFilename = optarg;
                    break;
                case 'o':
                    outputFilename = optarg;
                    break;
            }
        }
    }
    
    // validate command line input
    bool inputValid = checkInputValues(&p,inputFilename,vtCogName,orderingFilename,weightingFilename);
    if (inputValid == 0) {
        std::cerr << "Input is not valid!" << std::endl;
        usage(argv[0]);
        return 1;
    }
    
    //////////////////////
    // Parse input file //
    //////////////////////
    
    // parse input file
    std::vector<isolate*> *population = new std::vector<isolate*>;
    std::vector<cog*> *accessoryLoci = new std::vector<cog*>;
    std::vector<int> samplingList(p.numGen+1,0);
    std::vector<std::string> serotypeList;
    std::vector<int> scList;
    std::vector<std::string> cogList;
    int minGen = 0;
    
    int fileRead = parseInputFile(population,accessoryLoci,p.lowerLimit,p.upperLimit,samplingList,&serotypeList,&scList,&cogList,inputFilename,vtCogName,minGen);
    if (fileRead != 0) {
        usage(argv[0]);
        return 1;
    }
    int genLimit = p.numGen+minGen;
    int maxScNum = 1+(*std::max_element(std::begin(scList),std::end(scList)));
    
    //////////////////////
    // Pre-process data //
    //////////////////////
    
    // split population for immigration by SC
    std::vector<std::vector<isolate*> > *populationBySc = new std::vector<std::vector<isolate*> >;
    if (p.immigrationRate > 0 && p.immigrationType == 1) {
        int divCheck = dividePopulationForImmigration(population,&scList,populationBySc,maxScNum);
        if (divCheck != 0) {
            std::cerr << "Unable to split population into sequence clusters" << std::endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    // parse actual statistics
    if (p.programme == "s" && !(strcmp(frequencyFilename,"0"))) {
        int parseCheck = parseFrequencyFile(frequencyFilename,accessoryLoci);
        if (parseCheck != 0) {
            std::cerr << "Unable to parse frequency file" << std::endl;
            usage(argv[0]);
            return 1;
        }
    }

    // Differentially weight COGs by fixed input file or parameterisation
    if (weightingFilename != NULL) {
        int weightCheck = parseWeightingFile(weightingFilename,accessoryLoci);
        if (weightCheck != 0) {
            std::cerr << "Unable to parse weighting file" << std::endl;
            usage(argv[0]);
            return 1;
        }
    } else if (orderingFilename != NULL) {
        int orderCheck = parseOrderingFile(orderingFilename,accessoryLoci,&p);
        if (orderCheck != 0) {
            std::cerr << "Unable to parse ordering file" << std::endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    // extract COG information for simulation calculations
    std::vector<double> eqFreq;
    std::vector<double> cogWeights;
    std::vector<cog*>::iterator cit;
    for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
        eqFreq.push_back((*cit)->eqFreq);
        cogWeights.push_back((*cit)->weight);
    }
    
    /////////////////////////////
    // Select first generation //
    /////////////////////////////
    
    // run initial simulation - data structures
    // vectors for recording the current and next generation
    std::vector<isolate*> *currentIsolates = new std::vector<isolate*>;
    std::vector<isolate*> *futureIsolates = new std::vector<isolate*>;
    // 2D vectors for recording actual population history (integers)
    std::vector<std::vector<int> > vtScFreq(p.numGen+1,std::vector<int>(scList.size()));
    std::vector<std::vector<int> > nvtScFreq(p.numGen+1,std::vector<int>(scList.size()));
    std::vector<std::vector<int> > cogFreq(p.numGen+1,std::vector<int>(accessoryLoci->size()));
    // 2D vectors for recording cog deviations
    std::vector<std::vector<double> > piGen(p.numGen+1,std::vector<double>(accessoryLoci->size(),0.0));
    // 2D vectors for recording sampled population history (doubles)
    std::vector<std::vector<int> > sampledSeroFreq(p.numGen+1,std::vector<int>(serotypeList.size()));
    std::vector< std::vector<double> > sampledVtScFreq(p.numGen+1,std::vector<double>(scList.size(),0.0));
    std::vector< std::vector<double> > sampledNvtScFreq(p.numGen+1,std::vector<double>(scList.size(),0.0));
    // data structures for COG frequency measurements
    std::vector<double> cogDeviations(eqFreq.size());
    
    // initialise population in first generation, record simulated population statistics
    int gen = 0;
    int initialiseCheck = getStartingIsolates(population,currentIsolates,accessoryLoci,p.popSize,eqFreq,cogWeights,cogDeviations, vtScFreq[0],nvtScFreq[0],scList);
    if (initialiseCheck != 0) {
        std::cerr << "Unable to initialise population" << std::endl;
        usage(argv[0]);
        return 1;
    }
    
    // get sample from first generation
    int numberComparisons = 0;
    
    int summaryCheck = summariseGeneration(currentIsolates,samplingList[0],&scList,sampledVtScFreq,sampledNvtScFreq,&serotypeList,&sampledSeroFreq[gen]);
    if (summaryCheck != 0) {
        std::cerr << "Unable to summarise output of first generation" << std::endl;
        usage(argv[0]);
        return 1;
    }
    
    // get initial COG deviations
    std::vector<std::vector<double> > simulatedCogFrequencies;
    std::vector<double> vtCogFittingStatsList;
    std::vector<double> nvtCogFittingStatsList;
    std::vector<double> strainFittingStatsList;
    
    /////////////////////////////////
    // Iterate through generations //
    /////////////////////////////////
    
    // iterate through generations
    for (gen = minGen+1; gen <= genLimit; gen++) {
        
        // allow cells to reproduce and update COG deviations array
        int reproCheck = reproduction(currentIsolates,futureIsolates,population,populationBySc,&cogWeights,&cogDeviations,&p,&eqFreq,&vtScFreq[gen-minGen],&nvtScFreq[gen-minGen],&piGen[gen-minGen],&scList,gen);
        if (reproCheck != 0) {
            std::cerr << "Population failed to reproduce at generation " << gen << std::endl;
            usage(argv[0]);
            return 1;
        }
        
        // move on to next generation
        currentIsolates->clear();
        currentIsolates->insert(currentIsolates->begin(),futureIsolates->begin(),futureIsolates->end());
        futureIsolates->clear();
        
        if ((gen-minGen) < samplingList.size() && samplingList[gen-minGen] > 0 && p.programme != "s" && p.programme != "x") {
            int compareSamplesCheck = compareSamples(gen,minGen,samplingList[gen-minGen],currentIsolates,population,accessoryLoci,scList,sampledVtScFreq,sampledNvtScFreq,sampledSeroFreq[gen-minGen],serotypeList,vtCogFittingStatsList,nvtCogFittingStatsList,strainFittingStatsList);
            if (compareSamplesCheck != 0) {
                std::cerr << "Unable to compare simulated and actual frequencies" << std::endl;
                usage(argv[0]);
                return 1;
            } else {
                numberComparisons++;
            }
        }
    }
    
    ///////////////////////////////
    // Print summary information //
    ///////////////////////////////
    
    if (p.programme != "s" && p.programme != "x") {
        
        // calculate reproductive fitness metric
        std::vector<double> rFitVector(scList.size(),0.0);
        int rFitMetricCheck = rFitMetricCalculation(minGen,maxScNum,p.numGen,samplingList,scList,sampledVtScFreq,sampledNvtScFreq,population,rFitVector);
        if (rFitMetricCheck != 0) {
            std::cerr << "Unable to calculate reproductive fitness metric!" << std::endl;
            usage(argv[0]);
            return 1;
        }
        
        // sum up deviations
        double totalVtCogDeviation = 0.0;
        double totalNvtCogDeviation = 0.0;
        double totalStrainDeviation = 0.0;
        for (int i = 0; i < vtCogFittingStatsList.size(); i++) {
            totalVtCogDeviation+=vtCogFittingStatsList[i];
            totalNvtCogDeviation+=nvtCogFittingStatsList[i];
            totalStrainDeviation+=strainFittingStatsList[i];
        }
        
        // now add reproductive fitness deviations
        double totalRFitnessDeviation = 0.0;
        for (int i = 0; i < scList.size(); i++) {
            totalRFitnessDeviation+=rFitVector[i];
        }
        
        // print summaries
        std::cout << totalVtCogDeviation << "\t" << totalNvtCogDeviation << "\t" << totalStrainDeviation << "\t" << totalRFitnessDeviation << "\t" << numberComparisons << std::endl;
    }
    
    /////////////////////////////
    // Print full output files //
    /////////////////////////////
    
    // print output files if simulating
    if (p.programme != "f") {
        int printCheck = printOutput(outputFilename,&serotypeList,sampledSeroFreq,&scList,vtScFreq,nvtScFreq,&cogList,&cogWeights,p.numGen,minGen,accessoryLoci,samplingList,piGen,&p);
        if (printCheck != 0) {
            std::cerr << "Could not write to output files" << std::endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    return 0;
}

