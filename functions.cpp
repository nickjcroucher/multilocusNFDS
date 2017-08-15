//
//  functions.cpp
//  frequencyDependentSimulation
//
//  Created by Nicholas Croucher on 28/09/2015.
//  Copyright (c) 2015 Imperial College. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <sys/time.h>
#include <cmath>
#include "functions.h"
#include "parms.h"

///////////////////////
// Generic functions //
///////////////////////

/////////////////////////////
// random number generator //
/////////////////////////////

unsigned long int random_seed() {
    
    unsigned int seed;
    struct timeval tv;
    FILE *devrandom;
    
    if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
        gettimeofday(&tv,0);
        seed = tv.tv_sec + tv.tv_usec;
    } else {
        size_t nread = 0;
        nread = fread(&seed,sizeof(seed),1,devrandom);
        fclose(devrandom);
    }
    
    return(seed);
    
}

///////////////////
// usage message //
///////////////////

void usage (char* fn) {
    
    std::cerr << "Programme:" << std::endl << "\tp\tprogramme type:" << std::endl << "\t\t's' - simulation without fitting" << std::endl << "\t\t'x' - extended simulation output" << std::endl << "\t\t'f' - just return fitting metrics" << std::endl << "\t\t'b' - both fit to data and print simulation" << std::endl << std::endl << "Simulation parameters:" << std::endl << "\ts\tfrequency dependent selection pressure [double between 0 and 1]" << std::endl << "\tv\tvaccine selection pressure  [double between 0 and 1]" << std::endl << "\ti\timmigration rate [double between 0 and 1]" << std::endl << "\tt\timmigration type [0 - by strain, 1 - by SC]" << std::endl << "\tn\tpopulation carrying capacity  [integer]" << std::endl << "\tg\tnumber of generations [integer]" << std::endl << "\tu\tupper gene frequency limit [double between 0 and 1]" << std::endl << "\tl\tlower gene frequency limit [double between 0 and 1]" << std::endl << "\tq\tgeneration in which vaccine formulation is changed" << std::endl << std::endl << "Model fitting parameters:" << std::endl << "\tc\tvaccine target COG name [string]" << std::endl << "\tb\tbeginning year [integer]" << std::endl << "\tm\tmid year [integer]" << std::endl << "\te\tending year [integer]" << std::endl << std::endl << "Filenames:" << std::endl << "\tf\tinputFilename" << std::endl << "\tx\tfrequency file name [only for simulating]" << std::endl << "\to\toutput file prefix" << std::endl << "\tw\tweighting file" << std::endl  << "\tr\tcog reordering file" << std::endl;
    
}

////////////////////////
// Input file parsing //
////////////////////////

//////////////////////
// parse input file //
//////////////////////

int parseInputFile(std::vector<isolate*> *pop, std::vector<cog*> *accessoryLoci, double lower, double upper,  std::vector<int> &samplingList,std::vector<std::string> *st, std::vector<int> *sc, std::vector<std::string> *cogList, char *inputFilename, char * vtCogName, int &minGen) {
    
    // indices
    int s = 0;
    int v = -1;
    
    // COG sampling data structure
    std::vector<std::string> tmpCogList;
    std::vector<int> samplingTimes;
    int eqPop = 0;
    
    // parse file
    std::ifstream infile;
    infile.open(inputFilename, std::ifstream::in);
    if (infile) {
        std::string line;
        // parse lines
        while (std::getline(infile, line)) {
            // temporary information stores
            std::string sample_id;
            int sample_time = -1;
            int sample_sc = -1;
            std::string sample_serotype = "noSero";
            double sample_vt;
            double sample_latent_vt;
            std::vector<bool> sample_genotype;
            std::string iname;
            std::istringstream iss(line);
            int sIndex = 0;
            while (iss) {
                std::string temp;
                while (getline(iss, temp, '\t')) {
                    if (sIndex == 0) {
                        iname = temp;
                        sample_id = iname;
                    } else if (sIndex == 1) {
                        sample_time = atoi(temp.c_str());
                        if (sample_time < minGen) {
                            minGen = sample_time;
                        }
                    } else if (sIndex == 2) {
                        sample_serotype = temp;
                    } else if (sIndex == 3) {
                        double vt_double = atof(temp.c_str());
                        if (vt_double == 0) {
                            sample_vt = 0.0;
                            sample_latent_vt = 0.0;
                        } else if (vt_double <= 1) {
                            sample_vt = vt_double;
                            sample_latent_vt = 0.0;
                        } else if (vt_double > 1) {
                            sample_vt = 0.0;
                            sample_latent_vt = vt_double - 1.0;
                        } else {
                            std::cerr << "Unknown VT status for " << sample_id << std::endl;
                        }
                    } else if (sIndex == 4) {
                        sample_sc = atoi(temp.c_str());
                    } else if (sIndex > 4) {
                        if (iname == "Taxon") {
                            // record COG names
                            tmpCogList.push_back(temp);
                            // identify the key VT COG
                            if (vtCogName != 0 && temp == vtCogName) {
                                v = sIndex - 5;
                            }
                        } else {
                            sample_genotype.push_back(atoi(temp.c_str()));
                        }
                    }
                    sIndex++;
                }
            }
            if (iname != "Taxon" && sample_time != -1 && sample_sc != -1 && sample_serotype.compare("noSero") != 0) {
                std::vector<bool> sample_markers(0);
                isolate* tmp = new isolate(sample_id,sample_time,sample_sc,sample_serotype,sample_vt,sample_latent_vt,&sample_genotype,&sample_markers);
                pop->push_back(tmp);
                // record isolate information
                samplingTimes.push_back(sample_time);
                sc->push_back(sample_sc);
                st->push_back(sample_serotype);
                // calculate pre- and peri-vaccine population size
                if (sample_time <= 0) {
                    eqPop++;
                }
                // increment index
                s++;
            }
        }
        infile.close();
        
        // get unique serotypes and sequence clusters
        sort(st->begin(),st->end());
        st->erase(unique(st->begin(),st->end()),st->end());
        sort(sc->begin(),sc->end());
        sc->erase(unique(sc->begin(),sc->end()),sc->end());
        
        // calculate sample timings and sizes
        int maxTime = 0;
        std::vector<int>::iterator iter;
        for (iter = samplingTimes.begin(), samplingTimes.end() ; iter != samplingTimes.end(); ++iter) {
            if ((*iter) > maxTime) {
                maxTime = (*iter);
            }
        }

        for (iter = samplingTimes.begin(), samplingTimes.end() ; iter != samplingTimes.end(); ++iter) {
            samplingList[(*iter)-minGen]++;
        }
        
        // create vector of accessory locus COG objects
        std::vector<int> includeLocus(tmpCogList.size(),0);
        std::vector<isolate*>::iterator iiter;
        std::vector<std::string> intCogList;
        for (int i = 0; i < tmpCogList.size(); i++) {
            // data structures for recording frequencies
            int overallFreq = 0;
            std::vector<double> cogFrequencies(samplingList.size(),0);
            double eqFreq = 0;
            // calculate gene frequencies from isolates
            for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                overallFreq+=(*iiter)->genotype[i];
                cogFrequencies[((*iiter)->year)-minGen]+=(double((*iiter)->genotype[i])/double(samplingList[((*iiter)->year)-minGen]));
                if ((*iiter)->year <= 0) {
                    eqFreq+=(double((*iiter)->genotype[i])/double(eqPop));
                }
            }
            // retain if present at intermediate frequency OR vt-defining COG
//            if ((double(overallFreq)/double(pop->size()) >= lower && double(overallFreq)/double(pop->size()) <= upper) || i == v) {
            if ((eqFreq >= (lower-1e-07) && eqFreq <= (upper+1e-07)) || i == v) {
                includeLocus[i] = 1;
                int vtType = 0;
                if (i == v) {
                    vtType = 1;
                }
                intCogList.push_back(tmpCogList[i]);
                cog* tmpCog = new cog(tmpCogList[i],vtType,1.0,eqFreq,&cogFrequencies);
                if (tmpCog->id.length() > 0) {
                    accessoryLoci->push_back(tmpCog);
                } else {
                    
                    std::cerr << "Undefined COG ID at line " << i << std::endl;
                    
                }
            } else {
                includeLocus[i] = 0;
            }
        }
        
        // recalculate isolate genotypes to only include COGs at intermediate frequency
        for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
            std::vector<bool> tmpGenotype;
            for (int i = 0; i < includeLocus.size(); i++) {
                if (includeLocus[i] == 1) {
                    tmpGenotype.push_back((*iiter)->genotype[i]);
                }
            }
            (*iiter)->genotype = tmpGenotype;
        }
        
        // return cog list
        *cogList = intCogList;
        
    } else {
        std::cerr << "Problem with input file: " << strerror(errno) << std::endl;
        return 1;
    }
    
    return 0;
}

//////////////////////////////
// parse marker information //
//////////////////////////////

int parseMarkerFile(std::vector<isolate*> *pop,char *markerFilename,std::vector<std::string> *markerList) {
    
    // marker sampling data structure
    int markerLength = 0;
    
    // parse file
    std::ifstream infile;
    infile.open(markerFilename, std::ifstream::in);
    if (infile) {
        std::string line;
        // parse lines
        while (std::getline(infile, line)) {
            // temporary information stores
            std::string sample_id;
            std::vector<bool> sample_markers;
            std::string iname;
            std::istringstream iss(line);
            int sIndex = 0;
            while (iss) {
                std::string temp;
                while (getline(iss, temp, '\t')) {
                    if (sIndex == 0) {
                        sample_id = temp;

                    } else if (sIndex > 4) {
                        if (sample_id == "Taxon") {
                            // record COG names
                            markerList->push_back(temp);
                        } else {
                            sample_markers.push_back(atoi(temp.c_str()));
                        }
                    }
                    sIndex++;
                }
            }
            if (iname != "Taxon") {
                std::vector<isolate*>::iterator iiter;
                for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                    if ((*iiter)->id == sample_id) {
                        (*iiter)->markers = sample_markers;
                        if (sample_markers.size() > 0) {
                            markerLength = sample_markers.size();
                        }
                    }
                }
            }
        }
        infile.close();
        
        // check all marker lengths are the same
        std::vector<isolate*>::iterator iiter;
        for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
            if (markerLength != (*iiter)->markers.size()) {
                std::cerr << "Isolate " << (*iiter)->id << " has incorrect marker information; expecting " << markerLength << " but found " << (*iiter)->markers.size() << std::endl;
                return 1;
            }
        }
        
    } else {
        std::cerr << "Problem with marker file: " << strerror(errno) << std::endl;
        return 1;
    }
    
    return 0;
    
}

////////////////////////////
// parse input parameters //
////////////////////////////

bool checkInputValues(struct parms *sp,char * inputFilename,char * vtCogName, char* propFile, char* weightFile) {
    bool tmpvalid = 1;
    // check mode
    if (sp->programme != "f" && sp->programme != "s" && sp->programme != "b" && sp->programme != "x") {
        std::cerr << "Must select a programme: 'f' for fitting, 's' for simulating ('x' for extended output), 'b' for both" << std::endl;
        tmpvalid = 0;
    }
    if (propFile != 0 && weightFile != 0) {
        std::cerr << "You can't use a proportion file and a weight file. You just can't, stop it." << std::endl;
        tmpvalid = 0;
    }
    // check values for parameters
    if (sp->fSelection < 0 || sp->fSelection > 1000) {
        std::cerr << "Invalid frequency dependent selection value: " << sp->fSelection << "; should be between 0 and 1" << std::endl;
        tmpvalid = 0;
    }
    if (sp->vSelection < 0 || sp->vSelection > 1) {
        std::cerr << "Invalid vaccine selection value: " << sp->vSelection << "; should be between 0 and 1" << std::endl;
        tmpvalid = 0;
    }
    if (sp->immigrationRate < 0 || sp->immigrationRate > 1) {
        std::cerr << "Invalid immigration rate value: " << sp->immigrationRate << "; should be between 0 and 1" << std::endl;
        tmpvalid = 0;
    }
    if (sp->immigrationRate > 0) {
        if (sp->immigrationType != 0 && sp->immigrationType != 1) {
            std::cerr << "Invalid immigration type: " << sp->immigrationType << "; should be either '0' (by isolate) or '1' (by sc)" << std::endl;
            tmpvalid = 0;
        }
    }
    if (sp->upperLimit < 0 || sp->upperLimit > 1) {
        std::cerr << "Invalid upper limit value: " << sp->upperLimit << "; should be between 0 and 1" << std::endl;
        tmpvalid = 0;
    }
    if (sp->lowerLimit < 0 || sp->lowerLimit > 1 || sp->lowerLimit >= sp->upperLimit) {
        std::cerr << "Invalid lower limit value: " << sp->lowerLimit << "; should be between 0 and 1 and lower than the upper limit of " << sp->upperLimit << std::endl;
        tmpvalid = 0;
    }
    if (sp->selectedProp < 0 || sp->selectedProp > 1) {
        std::cerr << "Proportion of intermediate-frequency COGs under frequency dependent selection: " << sp->selectedProp << "; needs to be between zero and one" << std::endl;
        tmpvalid = 0;
    } else if (sp->selectedProp < 1 && propFile == 0) {
        std::cerr << "You are limiting the proportion of COGs under frequency-dependent selection in a nonsensical manner - please stop" << std::endl;
    }
    if (sp->popSize < 0) {
        std::cerr << "Invalid population size value: " << sp->popSize << "; needs to be greater than zero" << std::endl;
        tmpvalid = 0;
    }
    if (sp->numGen < 0) {
        std::cerr << "Invalid number of generations: " << sp->numGen << "; needs to be greater than zero" << std::endl;
        tmpvalid = 0;
    }
    if (vtCogName == 0 && (sp->programme != "s" && sp->programme != "x")) {
        std::cerr << "Need a COG name: currently not defined" << std::endl;
        tmpvalid = 0;
    }
    // check recombination parameters
    if (sp->transformationProportion*sp->transformationRate == 0 && sp->transformationProportion+sp->transformationRate > 0) {
        std::cerr << "Warning! Need to set transformation proportion (z) and transformation rate (e) both greater than zero for recombination to occur" << std::endl;
        tmpvalid = 0;
    } else if (sp->transformationProportion+sp->transformationRate > 0) {
        if (!(sp->transformationAsymmetryLoci >= 0 && sp->transformationAsymmetryLoci <= 1)) {
            std::cerr << "Transformation asymmetry needs to be between 0 and 1" << std::endl;
            tmpvalid = 0;
        }
        if (!(sp->transformationAsymmetryMarker >= 0 && sp->transformationAsymmetryMarker <= 1)) {
            std::cerr << "Transformation asymmetry needs to be between 0 and 1" << std::endl;
            tmpvalid = 0;
        }
    }
    // check input file exists
    if (inputFilename != 0) {
        std::ifstream infile(inputFilename);
        if (!infile) {
            tmpvalid = 0;
            std::cerr << "Problem with input file: " << strerror(errno) << std::endl;
        }
        infile.close();
    } else {
        std::cerr << "No input file name provided!" << std::endl;
        tmpvalid = 0;
    }
    // return value
    return tmpvalid;
}

///////////////////////////////////////////////////
// parse COG frequency file for pure simulations //
///////////////////////////////////////////////////

int parseFrequencyFile(char *frequencyFilename,std::vector<cog*> *accessoryLoci) {
    
    // sc output
    std::ifstream fFile;
    fFile.open(frequencyFilename,std::ifstream::in);
    
    // parse frequency file line-by-line
    if (fFile.is_open()) {
        std::string line;
        while (std::getline(fFile, line)) {
            std::istringstream iss(line);
            int sIndex = 0;
            std::string cname;
            double frac = -1;
            while (iss) {
                std::string temp;
                while (getline(iss, temp, '\t')) {
                    if (sIndex == 0) {
                        cname = temp;
                    } else if (sIndex == 1) {
                        frac = std::stod(temp);
                    }
                    sIndex++;
                }
                // replace calculated equilibrium frequency with that specified by separate input file
                // only retain loci included by previous criteria, and included in this second file
                if (frac >= 0 && frac <= 1) {
                    std::vector<cog*>::iterator cit;
                    for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                        if ((*cit)->id.compare(cname) == 0) {
//                        if ((*cit)->id == cname) {
//                            cog* newLocus = (*cit);
//                            std::vector<double> newFreq(newLocus->actualFreq.size(),frac);
//                            (*cit)->actualFreq = newFreq;
                            (*cit)->eqFreq = frac;
                        }
                    }
                } else {
                    std::cerr << "Fraction needs to be between 0 and 1: " << frac << std::endl;
                    return 1;
                }
            }
        }
        fFile.close();
    } else {
        std::cerr << "Unable to read file " << frequencyFilename << std::endl;
        return 1;
    }
    
    
    return 0;
}

//////////////////////////////
// parse COG weighting file //
//////////////////////////////

int parseWeightingFile(char* weightingFilename,std::vector<cog*> *accessoryLoci) {
    
    // open weighting file
    std::ifstream wFile;
    wFile.open(weightingFilename,std::ifstream::in);
    
    // parse weighting file
    if (wFile.is_open()) {
        std::string line;
        while (std::getline(wFile, line)) {
            // parse values from line
            std::istringstream iss(line);
            int index = 0;
            std::string cname;
            double weight = 1.0;
            while (iss) {
                std::string temp;
                while (getline(iss, temp, '\t')) {
                    if (index == 0) {
                        cname = temp;
                    } else {
                        weight = std::stof(temp);
                    }
                    std::vector<cog*>::iterator cit;
                    for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                        if ((*cit)->id.compare(cname) == 0) {
                            (*cit)->weight = weight;
                        }
                    }
                    index++;
                }
            }
        }
        wFile.close();
    } else {
        std::cerr << "Unable to read file " << weightingFilename << std::endl;
        return 1;
    }
    
    return 0;
}

/////////////////////////////
// parse COG ordering file //
/////////////////////////////


int parseOrderingFile(char* orderingFilename,std::vector<cog*> *accessoryLoci,struct parms *sp) {
    
    // reordered list
    std::vector<std::string> orderedAccessoryLoci;
    std::vector<cog*>::iterator cit;
    
    // open ordering file
    std::ifstream oFile;
    oFile.open(orderingFilename,std::ifstream::in);
    
    // parse weighting file
    if (oFile.is_open()) {
        std::string line;
        while (std::getline(oFile, line)) {
            // parse values from line
            std::istringstream iss(line);
            while (iss) {
                std::string cogName;
                while (getline(iss, cogName)) {
                    orderedAccessoryLoci.push_back(cogName);
                }
            }
        }
        oFile.close();
    } else {
        std::cerr << "Unable to read file " << orderingFilename << std::endl;
        return 1;
    }
    
    // check no loci have been duplicated, or lost, during reordering
    if (accessoryLoci->size() != orderedAccessoryLoci.size()) {
        std::cerr << "Duplicate or missing COGs found in COG reordering file: expecting " << accessoryLoci->size() << ", found " << orderedAccessoryLoci.size() << "; beneath the lists are compared:" << std::endl;
        for (int j = 0; j < orderedAccessoryLoci.size(); j++) {
            std::cerr << orderedAccessoryLoci[j] << "\t";
            for (int i = 0; i < accessoryLoci->size(); i++) {
                if ((*accessoryLoci)[i]->id == orderedAccessoryLoci[j]) {
                    std::cerr << (*accessoryLoci)[i]->id;
                }
            }
            std::cerr << std::endl;
        }
        std::cerr << "Alternative ordering:" << std::endl;
        for (int i = 0; i < accessoryLoci->size(); i++) {
            std::cerr << (*accessoryLoci)[i]->id << "\t";
            for (int j = 0; j < orderedAccessoryLoci.size(); j++) {
                if ((*accessoryLoci)[i]->id == orderedAccessoryLoci[j]) {
                    std::cerr << orderedAccessoryLoci[j] << std::endl;
                }
            }
        }
        return 1;
    }
    
    // assign weights appropriately
    for (int cindex = 0; cindex < orderedAccessoryLoci.size(); cindex++) {
        for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
            if ((*cit)->id.compare(orderedAccessoryLoci[cindex]) == 0) {
                if ((float(cindex)/float(orderedAccessoryLoci.size())) <= (1.0-sp->selectedProp)) {
                    (*cit)->weight = sp->lowerSelection;
                } else {
                    (*cit)->weight = sp->higherSelection;
                }
            }
        }
    }
    
    return 0;
}

////////////////////////////////
// Pre-processing information //
////////////////////////////////


///////////////////////////////////////////
// divide isolates by SC for immigration //
///////////////////////////////////////////

int dividePopulationForImmigration(std::vector<isolate*> *pop,std::vector <int> *scList,std::vector<std::vector<isolate*> > *popBySc, int maxScNum) {
    
//    int maxScNum = *std::max_element(std::begin(*scList),std::end(*scList));
    
    // check there are > 0 sequence clusters
    if (maxScNum == 0) {
        std::cerr << "Cannot find any sequence clusters" << std::endl;
    }
    
//    maxScNum++;
//    std::vector<std::vector<isolate*> > tmpStrains(maxScNum);
    std::vector<std::vector<isolate*> > tmpStrains(scList->size());
    
    for (int s = 0; s < scList->size(); s++) {
        std::vector<isolate*>::iterator cit;
        for (cit = pop->begin(), pop->end(); cit != pop->end(); ++cit) {
            if ((*cit)->sc == (*scList)[s]) {
//                tmpStrains[(*scList)[s]].push_back((*cit));
                tmpStrains[s].push_back((*cit));
            }
        }
    }
    
    (*popBySc) = tmpStrains;
    
    return 0;
}

///////////////////////////
// get first year sample //
///////////////////////////

int getStartingIsolates(std::vector<isolate*> *pop,std::vector<isolate*> *first,std::vector<cog*> *accessoryLoci,int psize,std::vector<double> &eqFreq,std::vector<double> &cogWeights,std::vector<double> &cogDeviations,std::vector<int> &startingVtScFrequencies,std::vector<int> &startingNvtScFrequencies,std::vector<int> &scList) {
    
    // get all isolates observed in the pre- or peri-vaccine samples
    std::vector<isolate*> possibleFirst;
    std::vector<isolate*>::iterator iter;
    for (iter = pop->begin(), pop->end() ; iter != pop->end(); ++iter) {
        if ((*iter)->year <= 0) {
            possibleFirst.push_back(*iter);
        }
    }
    
    // fill first timepoint with random sample of isolates from pre-/peri-vaccination samples
    // record starting COG frequencies
    std::vector<int> observedVtSc;
    std::vector<int> observedNvtSc;
    while (first->size() < psize) {
        //int selection = rand()%possibleFirst.size();
        int selection = int(double(gsl_rng_uniform(rgen))*int(possibleFirst.size()));
        first->push_back(possibleFirst[selection]);
        // record sequence clusters
        if (possibleFirst[selection]->vt > 0) {
            observedVtSc.push_back(possibleFirst[selection]->sc);
        } else {
            observedNvtSc.push_back(possibleFirst[selection]->sc);
        }
        // calculate gene frequencies
        for (int i = 0; i < possibleFirst[selection]->genotype.size();i++) {
            (*accessoryLoci)[i]->simFreq[0]+=(double(possibleFirst[selection]->genotype[i])/double(psize));
        }
    }
    
    // record sequence cluster statistics
    for (int i = 0; i < scList.size(); ++i) {
        startingVtScFrequencies[i] = std::count(observedVtSc.begin(),observedVtSc.end(),scList[i]);
        startingNvtScFrequencies[i] = std::count(observedNvtSc.begin(),observedNvtSc.end(),scList[i]);
    }
    
    // record gene frequency statistics
    std::vector<double> startingCogFrequencies(accessoryLoci->size(),0.0);
    for (int i = 0; i < accessoryLoci->size(); i++) {
        startingCogFrequencies[i] = (*accessoryLoci)[i]->simFreq[0];
    }
    std::transform(eqFreq.begin(), eqFreq.end(), startingCogFrequencies.begin(), cogDeviations.begin(), std::minus<double>());
    std::transform(cogWeights.begin(), cogWeights.end(), cogDeviations.begin(), cogDeviations.begin(), std::multiplies<double>());
    
    return 0;
}

////////////////////////////////////////////////////
// Select first sample for input file replication //
////////////////////////////////////////////////////

int firstSample(std::vector<isolate*> *currentIsolates,int firstSample,std::ofstream& sampleOutFile,int minGen) {
    
    // data structures for sample
    std::vector<isolate*> isolateSample;
    
    // get appropriately sized random sample from simulation
    while (isolateSample.size() < firstSample) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample.push_back(selectedIsolate);
        // print record of sample to file
        int vtInt = 0;
        if (selectedIsolate->latent_vt > 0) {
            vtInt = 2;
        } else if (selectedIsolate->vt > 0) {
            vtInt = 1;
        }
        sampleOutFile << selectedIsolate->id << "\t" << minGen << "\t" << selectedIsolate->serotype << "\t" << vtInt << "\t" << selectedIsolate->sc << std::endl;
    }
    
    return 0;
    
}

///////////////////////////////////////////
// get summary statistics per generation //
///////////////////////////////////////////

int summariseGeneration(std::vector<isolate*> *pop,int sampleSize,std::vector<int> *scs,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<std::string> *seros,std::vector<int> *serotypeF) {

    // record sequence clusters and serotypes from random sample
    std::vector<std::string> genSerotypes;
    std::vector<int> sampledVtSequenceClusters;
    std::vector<int> sampledNvtSequenceClusters;
    
    // select random sample and record sequence clusters and serotypes
    for (int i = 0; i <= sampleSize; ++i) {
        //int selection = rand()%pop->size();
        int selection = int(double(gsl_rng_uniform(rgen))*pop->size());
        isolate selectedIsolate = *(*pop)[selection];
        if (selectedIsolate.vt) {
            sampledVtSequenceClusters.push_back(selectedIsolate.sc);
        } else {
            sampledNvtSequenceClusters.push_back(selectedIsolate.sc);
        }
        // record serotypes
        genSerotypes.push_back(selectedIsolate.serotype);
    }
    
    // now summarise serotypes in sample
    for (int i = 0; i < seros->size(); i++) {
        int tmpSeroNum = std::count(genSerotypes.begin(),genSerotypes.end(),(*seros)[i]);
        (*serotypeF)[i] = tmpSeroNum;
    }
    
    // summarise sequence clusters in sample
    for (int i = 0; i < scs->size(); i++) {
        sampledVtScFreq[0][i]+=(double(std::count(sampledVtSequenceClusters.begin(),sampledVtSequenceClusters.end(),(*scs)[i]))/double(sampleSize));
        sampledNvtScFreq[0][i]+=(double(std::count(sampledNvtSequenceClusters.begin(),sampledNvtSequenceClusters.end(),(*scs)[i]))/double(sampleSize));
    }
    
    return 0;
}

/////////////////////////////////////////////////
// alter vaccine formulation during simulation //
/////////////////////////////////////////////////

int alterVaccineFormulation(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<std::vector<isolate*> > *popBySc) {

    std::vector<isolate*>::iterator iter;
    
    // change VT of current population
    for (iter = currentIsolates->begin(), currentIsolates->end(); iter != currentIsolates->end(); ++iter) {
        if ((*iter)->vt == 0) {
            (*iter)->vt = (*iter)->latent_vt;
        }
    }
    
    // change VT of underlying population
    for (iter = pop->begin(), pop->end(); iter != pop->end(); ++iter) {
        if ((*iter)->vt == 0) {
            (*iter)->vt = (*iter)->latent_vt;
        }
    }
    
    // change VT of potential immigrant isolates
    for (int i = 0; i < popBySc->size(); i++) {
        for (iter = (*popBySc)[i].begin(), (*popBySc)[i].end(); iter != (*popBySc)[i].end(); ++iter) {
            if ((*iter)->vt == 0.0) {
                (*iter)->vt = (*iter)->latent_vt;
            }
        }
    }
    
    return 0;
}

//////////////////////////////
// Per-generation functions //
//////////////////////////////

//////////////////////////////////////////////
// select next generation in the population //
//////////////////////////////////////////////

int reproduction(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop,std::vector<std::vector<isolate*> > *popBySc, std::vector<double> *cogWeights, std::vector<double> *cogDeviations,struct parms *sp, std::vector<double> * ef, std::vector<int> * vtScFreq,std::vector<int> * nvtScFreq,std::vector<double> * piGen,std::vector<int> *scList,int gen) {
    
    // new COG deviations array
    std::vector<int> futureCogCount(ef->size());
    std::fill (futureCogCount.begin(),futureCogCount.end(),0);
    
    // record population statistics
    std::vector<int> futureVtScs;
    std::vector<int> futureNvtScs;
    std::vector<std::string> futureSerotypes;
    
    // basic reproduction number based on immigration and population size
    double baseR = (1-sp->immigrationRate)*(double(sp->popSize)/double(currentIsolates->size()));
    
    // sort current population
    std::vector<isolate*>::iterator iter;
    std::sort(currentIsolates->begin(),currentIsolates->end());
    
    // allow old generation to reproduce
    std::string oldId = "";
    double oldFitness = 0.0;
    for (iter = currentIsolates->begin(), currentIsolates->end(); iter != currentIsolates->end(); ++iter) {
        // calculate fitness of each new genotype in population
        if ((*iter)->id != oldId) {
            
            std::vector<double> fitnesses(cogDeviations->size());
            
            std::transform(cogDeviations->begin(), cogDeviations->end(), (*iter)->genotype.begin(), fitnesses.begin(), std::multiplies<double>());
            double freqDepFitSum = std::accumulate(fitnesses.begin(),fitnesses.end(),0.0);
            double freqDepFit = pow((1+sp->fSelection),freqDepFitSum);
            double vaccineFit = 1.0;
            
            // only switch on vaccine selection pressure after vaccine is introduced
            if (gen >= 0 && (*iter)->vt > 0) {
                vaccineFit = 1.0 - (double(sp->vSelection)*(*iter)->vt);
            }
            // calculate overall fitness
            double overallFitness = baseR*vaccineFit*freqDepFit;
            oldFitness = overallFitness;
            oldId = (*iter)->id;
            
            // debug
//            if (gen == 200) {
//                std::cerr << (*iter)->id << std::endl;
//            }
//                if ((*iter)->sc == 25) {
//                    std::cerr << gen << "\t" << freqDepFit << std::endl;
//                }
            // debug
            
        }
        // select offspring by Poisson distribution
        int progeny = gsl_ran_poisson(rgen,oldFitness);
        for (int p = 0; p < progeny; p++) {
            futureIsolates->push_back((*iter));
            // record statistics
            if ((*iter)->vt > 0) {
                futureVtScs.push_back((*iter)->sc);
            } else {
                futureNvtScs.push_back((*iter)->sc);
            }
            futureSerotypes.push_back((*iter)->serotype);
        }
        // tally their COGs in new matrix
        if (progeny > 0) {
            for (int c = 0; c < (*iter)->genotype.size(); c++) {
                futureCogCount[c]+=(*iter)->genotype[c]*progeny;
            }
        }
    }
    
    // allow for immigration
    double adjustedImmigrationRate = (sp->immigrationRate)*(double(sp->popSize)/double(currentIsolates->size()));
    int immigrants = gsl_ran_binomial(rgen,adjustedImmigrationRate,sp->popSize);
    
    for (int f = 0; f < immigrants; f++) {
        isolate* selectedIsolate;
        if (sp->immigrationType == 0) {
//            int selection = int(rand()%pop->size());
            int selection = int(double(gsl_rng_uniform(rgen))*pop->size());
            selectedIsolate = (*pop)[selection];
        } else if (sp->immigrationType == 1) {
//            int selectedScIndex = int(rand()%scList->size());
//            int selectedSc = (*scList)[selectedScIndex];
            int selectedScIndex = int(double(gsl_rng_uniform(rgen))*scList->size());
            std::vector<isolate*> candidates = (*popBySc)[selectedScIndex];
//            int selection = int(rand()%candidates.size());
            int selection = int(double(gsl_rng_uniform(rgen))*candidates.size());
            selectedIsolate = candidates[selection];
            
        } else {
            return 1;
        }
        futureIsolates->push_back(selectedIsolate);
        // record statistics
        if (selectedIsolate->vt > 0) {
            futureVtScs.push_back(selectedIsolate->sc);
        } else {
            futureNvtScs.push_back(selectedIsolate->sc);
        }
        futureSerotypes.push_back(selectedIsolate->serotype);
        // tally COGs in new matrix
        for (int c = 0; c < selectedIsolate->genotype.size(); c++) {
            futureCogCount[c]+=selectedIsolate->genotype[c];
        }
    }
    
    // calculate COG deviations in next generation
    std::vector<int>::iterator iit;
    std::vector<double> cogFractions;
    for (iit = futureCogCount.begin(), futureCogCount.end(); iit != futureCogCount.end(); ++iit) {
        cogFractions.push_back(double(*iit)/double(futureIsolates->size()));
    }
    std::fill(cogDeviations->begin(),cogDeviations->end(),0);
    std::transform(ef->begin(), ef->end(), cogFractions.begin(), cogDeviations->begin(), std::minus<double>());
    // now summarise COG deviations before they are weighted in the extended output mode
    std::transform(cogDeviations->begin(), cogDeviations->end(), piGen->begin(), piGen->begin(), std::plus<double>());
    // now weight the COG deviations
    std::transform(cogWeights->begin(), cogWeights->end(), cogDeviations->begin(), cogDeviations->begin(), std::multiplies<double>());
    
    // now summarise sequence clusters - vaccine types
    for (int i = 0; i < scList->size(); i++) {
        (*vtScFreq)[i] = std::count(futureVtScs.begin(),futureVtScs.end(),(*scList)[i]);
    }
    
    // now summarise sequence clusters - non-vaccine types
    for (int i = 0; i < scList->size(); i++) {
        (*nvtScFreq)[i] = std::count(futureNvtScs.begin(),futureNvtScs.end(),(*scList)[i]);
    }
    
    return 0;
}

///////////////////
// Recombination //
///////////////////

bool alleleExchange (bool r, bool d, double a) {
    
    if (r) {
        if (d) {
            return true;
        } else {
            double p =  (double)rand() / RAND_MAX;
            if (a <= 1) {
                return false;
            } else if (p <= double(1.0/a)) {
                return false;
            } else {
                return true;
            }
        }
    } else {
        if (d) {
            double p =  (double)rand() / RAND_MAX;
            if (a >= 1) {
                return true;
            } else if (p <= a) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
    
    std::cerr << "Failed to read genotype" << std::endl;
    exit(1);
    return false;
}

const std::vector<std::string> explode(const std::string& s, const char& c) {
    std::string buff = "";
    std::vector<std::string> v;
    
    for (std::string::size_type i = 0; i < s.size(); ++i) {
        char n = s[i];
        if (n != c) {
            buff+=n;
        } else if (n == c && buff != "") {
            v.push_back(buff);
            buff = "";
        }
    }
    
    if (buff != "") {
        v.push_back(buff);
    }
    
    return v;
}

int recombination(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop,char* markerFilename,double transformationProportion,double transformationRate,double transformationAsymmetryLoci, double transformationAsymmetryMarker,std::vector<double> *cogWeights, std::vector<double> *cogDeviations, std::vector<double> *ef) {

    // recombination
    std::vector<isolate*>::iterator iit;
    for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
        isolate recipient = *(*iit);
        bool recHappened = false;
        // test if recombination occurs in this isolate at this timestep
        double rTrans = (double)rand() / RAND_MAX;
        if (rTrans <= transformationProportion) {
            // if so, how much of the genome is replaced
            for (int i = 0; i < recipient.genotype.size(); i++) {
                // recombination at selected loci
                double pTrans = (double)rand() / RAND_MAX;
                if (pTrans <= transformationRate) {
                    int j = int(rand() % pop->size());
                    bool donorAllele = (*pop)[j]->genotype[i];
                    bool recipientAllele = recipient.genotype[i];
                    recipient.genotype[i] = alleleExchange(recipientAllele,donorAllele,transformationAsymmetryLoci);
                    if (recipient.genotype[i] != recipientAllele) {
                        recHappened = true;
                    }
                    
                }
            }
            // recombination at unselected markers
            if (markerFilename != NULL) {
                for (int i = 0; i < recipient.markers.size(); i++) {
                    double mTrans = (double)rand() / RAND_MAX;
                    if (mTrans <= transformationRate) {
                        recHappened = true;
                        int j = int(rand() % pop->size());
                        bool donorAllele = (*pop)[j]->markers[i];
                        bool recipientAllele = recipient.markers[i];
                        recipient.markers[i] = alleleExchange(recipientAllele,donorAllele,transformationAsymmetryMarker);
                        if (recipient.markers[i] != recipientAllele) {
                            recHappened = true;
                        }
                    }
                }
            }
            // update ID if recombination has happened
            if (recHappened) {
                
                std::vector<std::string> seglist;
                seglist = explode(recipient.id,'_');
                
                // generate new hexadecimal suffix
                char newSuffix[7];
                for(int i = 0; i < 7; i++) {
                    sprintf(newSuffix + i, "%x", rand() % 16);
                }
                
                // assign new ID
                recipient.id = seglist[0]+"_"+newSuffix;
                
                // store in new population
                
                isolate* tmp = new isolate(recipient.id,recipient.year,recipient.sc,recipient.serotype,recipient.vt,recipient.latent_vt,&recipient.genotype,&recipient.markers);
    //            pop->push_back(tmp);
                futureIsolates->push_back(&(*tmp));
                
            } else {
                futureIsolates->push_back((*iit));
            }
        } else {
            futureIsolates->push_back((*iit));
        }

    }

    // calculate COG deviations in next generation
    std::vector<double> cogFractions(ef->size(),0.0);
    for (iit = futureIsolates->begin(), futureIsolates->end(); iit != futureIsolates->end(); ++iit) {
        for (int c = 0; c < (*iit)->genotype.size(); c++) {
            cogFractions[c]+=(double((*iit)->genotype[c])/double(futureIsolates->size()));
        }
    }
    // alter COG deviation vector according to recombinant phenotypes
    std::fill(cogDeviations->begin(),cogDeviations->end(),0.0);
    std::transform(ef->begin(), ef->end(), cogFractions.begin(), cogDeviations->begin(), std::minus<double>());
    // now weight the COG deviations
    std::transform(cogWeights->begin(), cogWeights->end(), cogDeviations->begin(), cogDeviations->begin(), std::multiplies<double>());
    
//    std::vector<isolate>::iterator isoi;
//    for (isoi = recombinantGenotypes.begin(), recombinantGenotypes.end() ; isoi != recombinantGenotypes.end(); ++isoi) {
//        futureIsolates->push_back(&(*isoi));
//    }
    
    return 0;
    
}

////////////////////////////////////////
// Move isolates into next generation //
////////////////////////////////////////

int nextGeneration(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<isolate*> *pop) {
    
    sort(currentIsolates->begin(), currentIsolates->end());
    std::vector<isolate*>::iterator iit;
    iit = unique(currentIsolates->begin(), currentIsolates->end());
    currentIsolates->resize(distance(currentIsolates->begin(),iit));

    for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
        if (!(std::find(futureIsolates->begin(), futureIsolates->end(),(*iit))!=futureIsolates->end())) {
            if (!(std::find(pop->begin(), pop->end(),(*iit))!=pop->end())) {
                delete (*iit);
            }
        }
    }
    
//    std::vector <isolate*> isolateList = *currentIsolates;
//    
//    for (int i = 0; i < isolateList.size(); i++) {
//        if (i == 0 || isolateList[i] != isolateList[i-1]) {
//            if (!(std::find(futureIsolates->begin(), futureIsolates->end(),isolateList[i])!=futureIsolates->end())) {
//                delete isolateList[i];
//            }
//        }
//    }
    
//    std::vector<isolate*>::iterator iit;
//    for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
//        if (!(std::find(futureIsolates->begin(), futureIsolates->end(),(*iit))!=futureIsolates->end())) {
//            delete (*iit);
//        }
//    }
    
    currentIsolates->clear();
    currentIsolates->insert(currentIsolates->begin(),futureIsolates->begin(),futureIsolates->end());
    futureIsolates->clear();
    
    return 0;
    
}

////////////////////////////////////////////
// Compare simulation and genomic samples //
////////////////////////////////////////////

int compareSamples(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<cog*> *accessoryLoci,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<int> &sampledSeroFreq,std::vector<std::string> &serotypeList,std::vector<double> &vtCogFittingStatsList,std::vector<double> &nvtCogFittingStatsList,std::vector<double> &strainFittingStatsList,std::ofstream& sampleOutFile) {
    
    // data structures for sample
    std::vector<isolate*> isolateSample;
    std::vector<std::string> currentSerotypeObservations;
    std::vector<int> currentVtScObservations;
    std::vector<int> currentNvtScObservations;
    
    // get appropriately sized random sample from simulation
    while (isolateSample.size() < sampleSize) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample.push_back(selectedIsolate);
        // record serotype frequencies
        currentSerotypeObservations.push_back(selectedIsolate->serotype);
        // record sequence cluster frequencies
        if (selectedIsolate->vt > 0) {
            currentVtScObservations.push_back(selectedIsolate->sc);
        } else {
            currentNvtScObservations.push_back(selectedIsolate->sc);
        }
        // print record of sample to file
        int vtInt = 0;
        if (selectedIsolate->latent_vt > 0) {
            vtInt = 2;
        } else if (selectedIsolate->vt > 0) {
            vtInt = 1;
        }
        sampleOutFile << selectedIsolate->id << "\t" << gen << "\t" << selectedIsolate->serotype << "\t" << vtInt << "\t" << selectedIsolate->sc << std::endl;
        // calculate gene frequencies
        for (int i = 0; i < selectedIsolate->genotype.size();i++) {
            (*accessoryLoci)[i]->simFreq[gen-minGen]+=(double(selectedIsolate->genotype[i])/double(sampleSize));
        }
    }
    
    // summarise serotype information
    for (int i = 0; i < serotypeList.size(); i++) {
        sampledSeroFreq[i] = std::count(currentSerotypeObservations.begin(),currentSerotypeObservations.end(),serotypeList[i]);
    }

    // summarise sequence cluster information from genomic data
    std::vector<int> currentGenomicVtScObservations;
    std::vector<int> currentGenomicNvtScObservations;
    
    std::vector<isolate*>::iterator iiter;
    for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
        if ((*iiter)->year == gen && (*iiter)->vt > 0) {
            currentGenomicVtScObservations.push_back((*iiter)->sc);
        } else if ((*iiter)->year == gen && (*iiter)->vt == 0.0) {
            currentGenomicNvtScObservations.push_back((*iiter)->sc);
        }
    }
    
    // summarise sequence cluster observations into counts
    std::vector<int> currentVtScCounts(scList.size(),0);
    std::vector<int> currentNvtScCounts(scList.size(),0);
    std::vector<int> currentGenomicVtScCount(scList.size(),0);
    std::vector<int> currentGenomicNvtScCount(scList.size(),0);
    
    for (int i = 0; i < scList.size(); ++i) {
        // observations from simulated data
        currentVtScCounts[i] = std::count(currentVtScObservations.begin(),currentVtScObservations.end(),scList[i]);
        currentNvtScCounts[i] = std::count(currentNvtScObservations.begin(),currentNvtScObservations.end(),scList[i]);
        // general storage of frequencies for printing
        sampledVtScFreq[gen-minGen][i] = double(currentVtScCounts[i])/double(sampleSize);
        sampledNvtScFreq[gen-minGen][i] = double(currentNvtScCounts[i])/double(sampleSize);
        // observations from genomic data
        currentGenomicVtScCount[i] = std::count(currentGenomicVtScObservations.begin(),currentGenomicVtScObservations.end(),scList[i]);
        currentGenomicNvtScCount[i] = std::count(currentGenomicNvtScObservations.begin(),currentGenomicNvtScObservations.end(),scList[i]);
    }
    
    // calculate divergence between real and simulated samples
    double cumulativeScDifference = 0.0;
//    double cumulativeJSD = 0.0;
    for (int i = 0; i < scList.size(); i++) {
//        double dkl_vt = std::abs(((double(currentGenomicVtScCount[i])+0.5)/double(sampleSize))*log((double(currentGenomicVtScCount[i])+0.5)/(double(currentVtScCounts[i])+0.5)));
//        double dkl_nvt = std::abs(((double(currentGenomicNvtScCount[i])+0.5)/double(sampleSize))*log((double(currentGenomicNvtScCount[i])+0.5)/(double(currentNvtScCounts[i])+0.5)));
        
        // calculate Jensen-Shannon divergence
        double jsd_vt = 0;
        double jsd_nvt = 0;
        
        double currentVtSimFreq = double(currentVtScCounts[i])/double(sampleSize);
        double currentNvtSimFreq = double(currentNvtScCounts[i])/double(sampleSize);
        double currentVtGenomicFreq = double(currentGenomicVtScCount[i])/double(sampleSize);
        double currentNvtGenomicFreq = double(currentGenomicNvtScCount[i])/double(sampleSize);
        double m_vt = (0.5/double(sampleSize))*(double(currentGenomicVtScCount[i])+double(currentVtScCounts[i]));
        double m_nvt = (0.5/double(sampleSize))*(double(currentNvtScCounts[i])+double(currentGenomicNvtScCount[i]));
        if (m_vt > 0.0) {
            if (currentVtSimFreq > 0) {
                jsd_vt += 0.5*currentVtSimFreq*log(currentVtSimFreq/m_vt);
            }
            if (currentVtGenomicFreq > 0) {
                jsd_vt += 0.5*currentVtGenomicFreq*log(currentVtGenomicFreq/m_vt);
            }
        }
        if (m_nvt > 0.0) {
            if (currentNvtSimFreq > 0) {
                jsd_nvt += 0.5*currentNvtSimFreq*log(currentNvtSimFreq/m_nvt);
            }
            if (currentNvtGenomicFreq > 0) {
                jsd_nvt += 0.5*currentNvtGenomicFreq*log(currentNvtGenomicFreq/m_nvt);
            }
        }
        
        cumulativeScDifference+=(jsd_vt+jsd_nvt);
//        std::cerr << scList[i] << "\t" << dkl_vt << "\t" << dkl_nvt << "\t" << jsd_vt << "\t" << jsd_nvt << std::endl;
//        std::cerr << scList[i] << "\t" << dkl_vt << "\t" << dkl_nvt << "\t" << jsd_vt << "\t" << jsd_nvt << "\t" << currentVtSimFreq << "\t" << currentNvtSimFreq << "\t" << currentVtGenomicFreq << "\t" << currentNvtGenomicFreq << std::endl;
//        std::cerr << "sc " << scList[i] << " vt KLD " << dkl_vt << " nvt KLD " << dkl_nvt << " total KLD " << cumulativeScDifference << std::endl;
//        std::cerr << "\tsc " << scList[i] << " vt JSD " << jsd_vt << " nvt JSD " << jsd_nvt << " total JSD " << cumulativeJSD << std::endl;
    }
    
    // record COG frequencies from both simulation sample and genomic sample
    std::vector<double> currentSampleVtFreq;
    std::vector<double> currentSampleNvtFreq;
    std::vector<double> startingSampleVtFreq;
    std::vector<double> startingSampleNvtFreq;
    std::vector<double> currentGenomicVtFreq;
    std::vector<double> currentGenomicNvtFreq;
    std::vector<double> startingGenomicVtFreq;
    std::vector<double> startingGenomicNvtFreq;
    
    // recalculated each time - allows flexibility for changing VT definition over time for multiple
    // vaccine introductions
    for (int i = 0; i < accessoryLoci->size(); i++) {
        // separate by VT for fitting statistic calculation
        if ((*accessoryLoci)[i]->vt > 0) {
//            currentSampleVtFreq.push_back(generationCogFrequencies[i]);
            startingSampleVtFreq.push_back((*accessoryLoci)[i]->simFreq[0]);
            currentSampleVtFreq.push_back((*accessoryLoci)[i]->simFreq[gen-minGen]);
            startingGenomicVtFreq.push_back((*accessoryLoci)[i]->actualFreq[0]);
            currentGenomicVtFreq.push_back((*accessoryLoci)[i]->actualFreq[gen-minGen]);
        } else {
//            currentSampleNvtFreq.push_back(generationCogFrequencies[i]);
            startingSampleNvtFreq.push_back((*accessoryLoci)[i]->simFreq[0]);
            currentSampleNvtFreq.push_back((*accessoryLoci)[i]->simFreq[gen-minGen]);
            startingGenomicNvtFreq.push_back((*accessoryLoci)[i]->actualFreq[0]);
            currentGenomicNvtFreq.push_back((*accessoryLoci)[i]->actualFreq[gen-minGen]);
        }
    }
    
    // calculate vaccine type COGs - absolute, not squared, value
    double vtCogStat = 0;
    for (int i = 0; i < currentSampleVtFreq.size(); i++) {
        if (startingSampleVtFreq[i] > 0 && startingGenomicVtFreq[i] > 0) {
            vtCogStat+=std::abs((currentSampleVtFreq[i]/startingSampleVtFreq[i])-(currentGenomicVtFreq[i]/startingGenomicVtFreq[i]));
        }
    }
    
    // calculate correlations - absolute, not squared, value
    double simulationNvtCorrelation = pearson(&currentSampleNvtFreq,&startingSampleNvtFreq);
    double genomicNvtCorrelation = pearson(&currentGenomicNvtFreq,&startingGenomicNvtFreq);
//    double nvtCogStat = std::abs((simulationNvtCorrelation-genomicNvtCorrelation)/(1-genomicNvtCorrelation));
    double nvtCogStat = std::abs(simulationNvtCorrelation-genomicNvtCorrelation);
    
    // record final statistics
    vtCogFittingStatsList.push_back(vtCogStat);
    nvtCogFittingStatsList.push_back(nvtCogStat);
    strainFittingStatsList.push_back(cumulativeScDifference);
    
    return 0;
}

//////////////////////////////////////////////
// just record simple simulation statistics //
//////////////////////////////////////////////

int justRecordStats(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<cog*> *accessoryLoci) {
    
    // data structures for sample
    std::vector<isolate*> isolateSample;
    std::vector<std::string> currentSerotypeObservations;
    std::vector<int> currentVtScObservations;
    std::vector<int> currentNvtScObservations;
    
    // get appropriately sized random sample from simulation
    while (isolateSample.size() < sampleSize) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample.push_back(selectedIsolate);
        // record serotype frequencies
        currentSerotypeObservations.push_back(selectedIsolate->serotype);
        // record sequence cluster frequencies
        if (selectedIsolate->vt > 0) {
            currentVtScObservations.push_back(selectedIsolate->sc);
        } else {
            currentNvtScObservations.push_back(selectedIsolate->sc);
        }
        // calculate gene frequencies
        for (int i = 0; i < selectedIsolate->genotype.size();i++) {
            (*accessoryLoci)[i]->simFreq[gen-minGen]+=(double(selectedIsolate->genotype[i])/double(sampleSize));
        }
    }
    
    return 0;
}

//////////////////////////////////
// calculate summary statistics //
//////////////////////////////////

double pearson(std::vector<double> *x,std::vector<double> *y) {
    
    double xi = 0;
    double yi = 0;
    double xs = 0;
    double ys = 0;
    double xy = 0;
    
    // check length of arrays first
    if (x->size() != y->size()) {
        std::cerr << "Cannot calculate correlation between arrays of different lengths" << std::endl;
        return -1;
    }
    
    int n = x->size()-1;
    
    for (int i = 0; i <= n; i++) {
        xi+=(*x)[i];
        xs+=pow((*x)[i],2);
        yi+=(*y)[i];
        ys+=pow((*y)[i],2);
        xy+=((*x)[i]*(*y)[i]);
    }
    
    double frequencyCorrelation = (n*xy-xi*yi)/(sqrt(n*xs-pow(xi,2))*sqrt(n*ys-pow(yi,2)));
    
    return frequencyCorrelation;
    
}

//////////////////////
// Output functions //
//////////////////////

////////////////////////////////////////////////
// calculate reproductive fitness comparisons //
////////////////////////////////////////////////

int rFitMetricCalculation(int minGen,int maxScNum,int numGen,std::vector<int> &samplingList,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<isolate*> *population,std::vector<double> &rFitVector) {
    
    // data structures
    std::vector<double> actualFoldChanges(samplingList.size(),0.0);
    std::vector<double> simulatedFoldChanges(samplingList.size(),0.0);
    
    // calculate actual statistics
    std::vector< std::vector<double> > realVtScFreq(numGen+1,std::vector<double>(scList.size(),0.0));
//    realVtScFreq[genIndex][scIndex]+=(1/double(samplingList[genIndex]));
    std::vector< std::vector<double> > realNvtScFreq(numGen+1,std::vector<double>(scList.size(),0.0));
//    realNvtScFreq[genIndex][scIndex]+=(1/double(samplingList[genIndex]));
    
    // record VT and NVT observations
    std::vector< std::vector<int> > genomicVtObservations(numGen+1);
    std::vector< std::vector<int> > genomicNvtObservations(numGen+1);
    
    std::vector<isolate*>::iterator iiter;
    for (iiter = population->begin(), population->end(); iiter != population->end(); ++iiter) {
        int genIndex = (*iiter)->year-minGen;
        if (samplingList[genIndex] > 0) {
            if ((*iiter)->vt > 0) {
                genomicVtObservations[genIndex].push_back((*iiter)->sc);
            } else {
                genomicNvtObservations[genIndex].push_back((*iiter)->sc);
            }
        } else {
            std::cerr << "Misalignment between sampling generations!" << std::endl;
            return 1;
        }
    }
    
    // summarise SC information
    for (int genIndex = 0; genIndex < samplingList.size(); genIndex++) {
        if (samplingList[genIndex] > 0) {
            for (int scIndex = 0; scIndex < scList.size(); scIndex++) {
                realVtScFreq[genIndex][scIndex] = double(std::count(genomicVtObservations[genIndex].begin(),genomicVtObservations[genIndex].end(),scList[scIndex]))/double(samplingList[genIndex]);
                realNvtScFreq[genIndex][scIndex] = double(std::count(genomicNvtObservations[genIndex].begin(),genomicNvtObservations[genIndex].end(),scList[scIndex]))/double(samplingList[genIndex]);
            }
        }
    }
    
    // count number of timepoints at which samples are taken
    // for calculating mean frequencies
    double numberOfSamples = 0.0;
    for (int i = 0; i < samplingList.size(); ++i) {
        if (samplingList[i] > 0) {
            numberOfSamples+=1.0;
        }
    }
    
    // return comparison values
    for (int scIndex = 0; scIndex < scList.size(); scIndex++) {
        
        // record observed timesteps
        std::vector<int> realVtScObservationTimepoints;
        std::vector<int> realNvtScObservationTimepoints;
        std::vector<int> simVtScObservationTimepoints;
        std::vector<int> simNvtScObservationTimepoints;
        
        // record frequencies
        std::vector<double> realVtScObservations;
        std::vector<double> realNvtScObservations;
        std::vector<double> simVtScObservations;
        std::vector<double> simNvtScObservations;
        
        // iterate through generations
        for (int genIndex = 0; genIndex <= numGen; ++genIndex) {
            if (samplingList[genIndex] > 0) {
                // real VT observations
                if (realVtScFreq[genIndex][scIndex] > 0 || realVtScObservationTimepoints.size() > 0) {
                    realVtScObservations.push_back(realVtScFreq[genIndex][scIndex]);
                    realVtScObservationTimepoints.push_back(genIndex);
                }
                // real NVT observations
                if (realNvtScFreq[genIndex][scIndex] > 0 || realNvtScObservationTimepoints.size() > 0) {
                    realNvtScObservations.push_back(realNvtScFreq[genIndex][scIndex]);
                    realNvtScObservationTimepoints.push_back(genIndex);
                }
                // simulated & sampled VT observations
                if (sampledVtScFreq[genIndex][scIndex] > 0 || simVtScObservationTimepoints.size() > 0) {
                    simVtScObservations.push_back(sampledVtScFreq[genIndex][scIndex]);
                    simVtScObservationTimepoints.push_back(genIndex);
                }
                // simulated & sampled NVT observations
                if (sampledNvtScFreq[genIndex][scIndex] > 0 || simNvtScObservationTimepoints.size() > 0) {
                    simNvtScObservations.push_back(sampledNvtScFreq[genIndex][scIndex]);
                    simNvtScObservationTimepoints.push_back(genIndex);
                }
            }
        }
        // distance metric calculation
        double RmetricDeviation = 0.0;
        // calculate metric comparison if >1 timepoint for VT
        if (realVtScObservations.size() > 1 && simVtScObservations.size() > 1) {
            // calculate metric for real data
            std::vector<double> realTmp;
            double cumulativeFrequency = realVtScObservations[0];
            for (int j = 1; j < realVtScObservations.size(); ++j) {
                cumulativeFrequency+=realVtScObservations[j];
                double invExponent = double(realVtScObservationTimepoints[j]) - double(realVtScObservationTimepoints[0]);
                realTmp.push_back(pow(realVtScObservations[j]/realVtScObservations[0],(1/invExponent)));
            }
            double tmp_t = 0.0;
            double tmp_c = 0.0;
            for (int j = 0; j < realTmp.size(); ++j) {
                tmp_t+=realTmp[j];
                tmp_c++;
            }
            double realRmetric = tmp_t/tmp_c;
            // calculate metric for simulated data
            std::vector<double> simTmp;
            for (int j = 1; j < simVtScObservations.size(); ++j) {
                double invExponent = double(simVtScObservationTimepoints[j]) - double(simVtScObservationTimepoints[0]);
                simTmp.push_back(pow(simVtScObservations[j]/simVtScObservations[0],(1/invExponent)));
            }
            tmp_t = 0.0;
            tmp_c = 0.0;
            for (int j = 0; j < simTmp.size(); ++j) {
                tmp_t+=simTmp[j];
                tmp_c++;
            }
            double simRmetric = tmp_t/tmp_c;
            
            // summarise difference between metric estimates weighted by SC frequency in actual data
            RmetricDeviation+=((cumulativeFrequency/numberOfSamples)*fabs(realRmetric-simRmetric));

        }
        // calculate metric comparison if >1 timepoint for NVT
        if (realNvtScObservations.size() > 1 && simNvtScObservations.size() > 1) {
            
            // calculate metric for real data
            std::vector<double> realTmp;
            double cumulativeFrequency = realNvtScObservations[0];
            for (int j = 1; j < realNvtScObservations.size(); ++j) {
                cumulativeFrequency+=realNvtScObservations[j];
                double invExponent = double(realNvtScObservationTimepoints[j]) - double(realNvtScObservationTimepoints[0]);
                realTmp.push_back(pow(realNvtScObservations[j]/realNvtScObservations[0],(1/invExponent)));
            }
            double tmp_t = 0.0;
            double tmp_c = 0.0;
            for (int j = 0; j < realTmp.size(); ++j) {
                tmp_t+=realTmp[j];
                tmp_c++;
            }
            double realRmetric = tmp_t/tmp_c;
            
            // calculate metric for simulated data
            std::vector<double> simTmp;
            for (int j = 1; j < simNvtScObservations.size(); ++j) {
                double invExponent = double(simNvtScObservationTimepoints[j]) - double(simNvtScObservationTimepoints[0]);
                simTmp.push_back(pow(simNvtScObservations[j]/simNvtScObservations[0],(1/invExponent)));
            }
            tmp_t = 0.0;
            tmp_c = 0.0;
            for (int j = 0; j < simTmp.size(); ++j) {
                tmp_t+=simTmp[j];
                tmp_c++;
            }
            double simRmetric = tmp_t/tmp_c;
            
            // summarise difference between metric estimates weighted by SC frequency in actual data
            RmetricDeviation+=((cumulativeFrequency/numberOfSamples)*fabs(realRmetric-simRmetric));
        }
        // record deviation
        rFitVector[scIndex] = RmetricDeviation;
        
    }

    return 0;
    
}

///////////////////////////
// write output to files //
///////////////////////////

int printOutput(char* outputFilename,std::vector<std::string> *seroList,std::vector<std::vector<int> > &sampledSeroFreq,std::vector<int> *scList,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,std::vector<std::string> *cogList,std::vector<double> *cogWeights,int gen,int minGen,std::vector<cog*> *accessoryLoci,std::vector<int> samplingList,std::vector<std::vector<double> > &piGen,struct parms *sp) {
    
    // parse file names
    std::stringstream prefixStream;
    std::string prefix(outputFilename);
    // serotype output
    std::string seroOutFilename = prefix + ".sero.out";
    std::ofstream seroOutFile;
    seroOutFile.open(seroOutFilename,std::ios::out);
    
    // write serotype output
    if (seroOutFile.is_open()) {
        // write header
        seroOutFile << "Serotype";
        std::vector<std::string>::iterator siter;
        for (siter = seroList->begin(), seroList->end(); siter != seroList->end(); ++siter) {
            seroOutFile <<  "\t" << (*siter);
        }
        seroOutFile << std::endl;
        // write values
        int pseudoGen = minGen;
        std::vector<std::vector<int> >::iterator iiter;
        for (iiter = sampledSeroFreq.begin(), sampledSeroFreq.end(); iiter != sampledSeroFreq.end(); ++iiter) {
            seroOutFile << pseudoGen;
            std::vector<int>::iterator jiter;
            for (jiter = (*iiter).begin(), (*iiter).end(); jiter != (*iiter).end(); ++jiter) {
                seroOutFile << "\t" << (*jiter);
            }
            seroOutFile << std::endl;
            pseudoGen++;
        }
        seroOutFile.close();
    } else {
        std::cerr << "Unable to write to file " << seroOutFilename << std::endl;
        return 1;
    }
    
    // sc output
    std::string scOutFilename = prefix + ".sc.out";
    std::ofstream scOutFile;
    scOutFile.open(scOutFilename,std::ios::out);
    
    // write serotype output
    if (scOutFile.is_open()) {
        // write header
        scOutFile << "SC";
        std::vector<int>::iterator siter;
        for (siter = scList->begin(), scList->end(); siter != scList->end(); ++siter) {
            scOutFile <<  "\tSC" << (*siter) << "_VT\tSC" << (*siter) << "_NVT";
        }
        scOutFile << std::endl;
        // write values
        for (int pseudoGen = minGen; pseudoGen < (gen+1+minGen); pseudoGen++) {
            scOutFile << pseudoGen;
            for (int i = 0; i != scList->size(); i++) {
                scOutFile << "\t" << vtScFreq[pseudoGen-minGen][i] << "\t" << nvtScFreq[pseudoGen-minGen][i];
            }
            scOutFile << std::endl;
        }
        scOutFile.close();
    } else {
        std::cerr << "Unable to write to file " << scOutFilename << std::endl;
        return 1;
    }
    
    // cog output
    std::string cogOutFilename = prefix + ".cog.out";
    std::ofstream cogOutFile;
    cogOutFile.open(cogOutFilename,std::ios::out);
    
    // write cog output
    if (cogOutFile.is_open()) {
        // write header
        cogOutFile << "Name\tVT\tWeight";
        int t = 0;
        for (int g = 0; g < samplingList.size(); g++) {
            if (samplingList[g] > 0) {
                if (t == 0) {
                    // need to replace with minimum generation information
                    cogOutFile << "\t" << "InitialGenomicFrequency" << "\t" << "InitialSimulatedFrequency";
                } else {
                    // also need to alter based on minimum generation information
                    cogOutFile << "\t" << "Genome_" << g+minGen << "\t" << "Simulation_" << g+minGen;
                }
                t++;
            }
        }
        cogOutFile << std::endl;
        // write content
        for (int c = 0; c < accessoryLoci->size(); c++) {
            cog tmpCog = *(*accessoryLoci)[c];
            cogOutFile << tmpCog.id << "\t" << tmpCog.vt << "\t" << tmpCog.weight;
            for (int g = 0; g < tmpCog.actualFreq.size(); g++) {
                if (samplingList[g] > 0) {
                    cogOutFile << "\t" << tmpCog.actualFreq[g] << "\t" << tmpCog.simFreq[g];
                }
            }
            cogOutFile << std::endl;
        }
        
        // close
        cogOutFile.close();
    } else {
        std::cerr << "Unable to write to file " << cogOutFilename << std::endl;
        return 1;
    }
    
    // COG deviation output
    if (sp->programme == "x") {
        // pi output
        std::string piOutFilename = prefix + ".pi.out";
        std::ofstream piOutFile;
        piOutFile.open(piOutFilename,std::ios::out);
        
        // write pi output
        if (piOutFile.is_open()) {
            // write header
            piOutFile << "Gen";
            for (int g = 0; g < accessoryLoci->size(); g++) {
                piOutFile << "\t" << (*accessoryLoci)[g]->id;
            }
            piOutFile << std::endl;
            // write content
            for (int pseudoGen = minGen; pseudoGen < (gen+1+minGen); pseudoGen++) {
                piOutFile << pseudoGen;
                for (int c = 0; c < accessoryLoci->size(); c++) {
                    piOutFile << "\t" << piGen[pseudoGen][c];
                }
                piOutFile << std::endl;
            }
            // close
            piOutFile.close();
        } else {
            std::cerr << "Unable to write to file " << piOutFilename << std::endl;
            return 1;
        }
        
        
    }
    
    // end
    return 0;
}

//////////////////////
// print population //
//////////////////////

int printPop(char* prefixStar,std::string suffix,std::vector<isolate*> *currentIsolates,char* markerFilename,std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList) {
    
    // population size
    double pSize = double(currentIsolates->size());
    // parse file names
    std::stringstream prefixStream;
    std::string prefix(prefixStar);
    // selected loci output
    std::string popFileName = prefix + "." + suffix;
    std::ofstream popOutFile;
    popOutFile.open(popFileName,std::ios::out);
    // marker loci output
    std::string marFileName = prefix + ".markers." + suffix;
    std::ofstream marOutFile;
    if (markerFilename != NULL) {
        marOutFile.open(marFileName,std::ios::out);
    }
    
    // write accessory output
    if (popOutFile.is_open()) {

        popOutFile << "Locus\tFrequency" << std::endl;
        std::vector<double> currentFreq(accessoryLoci->size(),0.0);
        std::vector<isolate*>::iterator iit;
        for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
            for (int i = 0; i < (*iit)->genotype.size(); i++) {
                currentFreq[i]+=double((*iit)->genotype[i])/pSize;
            }
        }
        
        for (int i = 0; i < accessoryLoci->size(); i++) {
            popOutFile << (*accessoryLoci)[i]->id << "\t" << currentFreq[i] << std::endl;
        }
        
        popOutFile.close();
    } else {
        std::cerr << "Unable to write to file " << popFileName << std::endl;
        return 1;
    }

    // write marker output
    if (markerFilename != NULL) {
        if (marOutFile.is_open()) {

            marOutFile << "Marker\tFrequency" << std::endl;
            std::vector<double> currentFreq(markerList->size(),0.0);
            std::vector<isolate*>::iterator iit;
            for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
                for (int i = 0; i < (*iit)->markers.size(); i++) {
                    currentFreq[i]+=double((*iit)->markers[i])/pSize;
                }
            }
            
            for (int i = 0; i < markerList->size(); i++) {
                marOutFile << (*markerList)[i] << "\t" << currentFreq[i] << std::endl;
            }
            
            marOutFile.close();
        } else {
            std::cerr << "Unable to write to file " << markerFilename << std::endl;
            return 1;
        }
    }
    
    return 0;
}

int printPopSample (char* prefixStar,std::string suffix,std::vector<isolate*> *currentIsolates,char* markerFilename,std::vector<cog*> *accessoryLoci,std::vector<std::string> *markerList, int sampleSize) {

    // parse file names
    std::stringstream prefixStream;
    std::string prefix(prefixStar);
    // selected loci output
    std::string popFileName = prefix + "." + suffix;
    std::ofstream popOutFile;
    popOutFile.open(popFileName,std::ios::out);
    // marker loci output
    std::string marFileName = prefix + ".markers." + suffix;
    std::ofstream marOutFile;
    if (markerFilename != NULL) {
        marOutFile.open(marFileName,std::ios::out);
    }

    // pick random sample of specified size
    std::vector<isolate*> *currentSample = new std::vector<isolate*>;
//    std::vector<isolate*> currentSample;
    while (currentSample->size() < sampleSize) {
//        int selectedIsolateIndex = int(rand()%currentIsolates->size());
        int selectedIsolateIndex = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        currentSample->push_back((*currentIsolates)[selectedIsolateIndex]);
    }
    
    // write accessory output
    if (popOutFile.is_open()) {
    
        // print headers
        popOutFile << "Taxon    Time    Serotype    VT  SC";
        std::vector<cog*>::iterator cit;
        for (cit = accessoryLoci->begin(), accessoryLoci->end() ; cit != accessoryLoci->end(); ++cit) {
            popOutFile << "\t" << (*cit)->id;
        }
        popOutFile << std::endl;

        // print data
        std::vector<isolate*>::iterator iit;
        for (iit = currentSample->begin(), currentSample->end() ; iit != currentSample->end(); ++iit) {
            std::string outLine = (*iit)->id;
            // add metadata
            std::ostringstream metaString;
            metaString << "\t" << (*iit)->year << "\t" << (*iit)->serotype << "\t" << (*iit)->vt << "\t" << (*iit)->sc;
            outLine.append(metaString.str());
            // end add metadata
            for (int i = 0; i < (*iit)->genotype.size(); i++) {
                outLine.append("\t");
                outLine.append(std::string((*iit)->genotype[i] ? "1" : "0"));
            }
            popOutFile << outLine << std::endl;
        }
    } else {
        std::cerr << "Unable to write to population sample locus file " << popFileName << std::endl;
        return 1;
    }
    
    // write marker output
    if (markerFilename != NULL) {
        if (marOutFile.is_open()) {
    
            // print headers
            std::string titleLine = "Marker";
            for (int i = 0; i < markerList->size(); i++) {
                titleLine.append("\t");
                titleLine.append((*markerList)[i]);
            }
            marOutFile << titleLine << std::endl;

            // print data
            std::vector<isolate*>::iterator iit;
            for (iit = currentSample->begin(), currentSample->end() ; iit != currentSample->end(); ++iit) {
                std::string outLine = (*iit)->id;
                for (int i = 0; i < (*iit)->markers.size(); i++) {
                    outLine.append("\t");
                    outLine.append(std::string((*iit)->markers[i] ? "1" : "0"));
                }
                marOutFile << outLine << std::endl;
            }
    
        } else {
            std::cerr << "Unable to write to file " << markerFilename << std::endl;
            return 1;
        }
    }
            
    return 0;
}

//// Unsorted



//int getStats(std::vector<double> *x,std::vector<double> *y,int vtCog,double &frequencyCorrelation,double &vtCogDecrease) {
//    
//    double xi = 0;
//    double yi = 0;
//    double xs = 0;
//    double ys = 0;
//    double xy = 0;
//    
//    // check length of arrays first
//    if (x->size() != y->size()) {
//        std::cerr << "Cannot calculate correlation between arrays of different lengths" << std::endl;
//        return 1;
//    }
//    
//    int n = x->size()-1;
//    
//    for (int i = 0; i <= n; i++) {
//        if (i == vtCog) {
//            vtCogDecrease = (*y)[i]/(*x)[i];
//        } else {
//            xi+=(*x)[i];
//            xs+=pow((*x)[i],2);
//            yi+=(*y)[i];
//            ys+=pow((*y)[i],2);
//            xy+=((*x)[i]*(*y)[i]);
//        }
//    }
//    
//    frequencyCorrelation = (n*xy-xi*yi)/(sqrt(n*xs-pow(xi,2))*sqrt(n*ys-pow(yi,2)));
//    
//    return 0;
//}

/////////////////////////////////////////
// calculate strain summary statistics //
/////////////////////////////////////////

//int getStrainStats(std::vector<isolate*> *pop,std::vector<isolate*> *currentPop, int year, int numSc, double &popDev) {
//    
//    double tmpDev = 0;
//    
//    for (int s = 0; s < numSc; s++) {
//        // store frequencies
//        int totalStrains = 0;
//        int vtScStrains = 0;
//        int nvtScStrains = 0;
//        // calculate frequencies in actual population
//        std::vector<isolate*>::iterator iter;
//        for (iter = pop->begin(), pop->end() ; iter != pop->end(); ++iter) {
//            if ((*iter)->year == year) {
//                if ((*iter)->sc == s) {
//                    if ((*iter)->vt == 1) {
//                        vtScStrains++;
//                    } else {
//                        nvtScStrains++;
//                    }
//                }
//                totalStrains++;
//            }
//        }
//        double vtActualPrev = double(vtScStrains)/double(totalStrains);
//        double nvtActualPrev = double(nvtScStrains)/double(totalStrains);
//        vtScStrains = 0;
//        nvtScStrains = 0;
//        totalStrains = 0;
//        // calculate frequencies in simulated population
//        for (iter = currentPop->begin(), currentPop->end() ; iter != currentPop->end(); ++iter) {
//            if ((*iter)->sc == s) {
//                if ((*iter)->vt == 1) {
//                    vtScStrains++;
//                } else {
//                    nvtScStrains++;
//                }
//            }
//            totalStrains++;
//        }
//        double vtSimPrev = double(vtScStrains)/double(totalStrains);
//        double nvtSimPrev = double(nvtScStrains)/double(totalStrains);
//        // increment deviation stat
//        tmpDev+=pow((vtSimPrev-vtActualPrev),2);
//        tmpDev+=pow((nvtSimPrev-nvtActualPrev),2);
//    }
//    
//    popDev = tmpDev;
//    
//    return 0;
//}







/////////////////////////////////////
// calculate final COG frequencies //
/////////////////////////////////////

//int calculateFinalCogFrequency(std::vector<isolate*> *pop,std::vector<double> *endingCogFrequencies) {
//    
//    std::vector<isolate*>::iterator cit;
//    double increment = 1/double(pop->size());
//    for (cit = pop->begin(), pop->end(); cit != pop->end(); ++cit) {
//        for (int c = 0; c < (*cit)->genotype.size(); c++) {
//            if ((*cit)->genotype[c] == 1) {
//                (*endingCogFrequencies)[c]+=increment;
//            }
//        }
//    }
//    return 0;
//}





////////////////////////////////////////////
// summarise changes in strain population //
////////////////////////////////////////////

//int getStrainChangeStats(std::vector<isolate*> *pop,int firstYear,int midYear,int finalYear,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,int mGen,int fGen,int numSc,double &changeDev) {
//    
//    // compare strain frequencies between each pair of years
//    double tmpDev = 0;
//    
//    for (int s = 1; s <= numSc; s++) {
//        // store frequencies
//        int totalStrainsFirst = 0;
//        int vtScStrainsFirst = 0;
//        int nvtScStrainsFirst = 0;
//        int totalStrainsMid = 0;
//        int vtScStrainsMid = 0;
//        int nvtScStrainsMid = 0;
//        int totalStrainsFinal = 0;
//        int vtScStrainsFinal = 0;
//        int nvtScStrainsFinal = 0;
//        // calculate frequencies in actual population
//        std::vector<isolate*>::iterator iter;
//        for (iter = pop->begin(), pop->end() ; iter != pop->end(); ++iter) {
//            if ((*iter)->year == firstYear) {
//                if ((*iter)->sc == s) {
//                    if ((*iter)->vt == 1) {
//                        vtScStrainsFirst++;
//                    } else {
//                        nvtScStrainsFirst++;
//                    }
//                }
//                totalStrainsFirst++;
//            } else if ((*iter)->year == midYear) {
//                if ((*iter)->sc == s) {
//                    if ((*iter)->vt == 1) {
//                        vtScStrainsMid++;
//                    } else {
//                        nvtScStrainsMid++;
//                    }
//                }
//                totalStrainsMid++;
//            } else if ((*iter)->year == finalYear) {
//                if ((*iter)->sc == s) {
//                    if ((*iter)->vt == 1) {
//                        vtScStrainsFinal++;
//                    } else {
//                        nvtScStrainsFinal++;
//                    }
//                }
//                totalStrainsFinal++;
//            }
//        }
//        double vtFirstChange = (double(vtScStrainsMid)/double(totalStrainsMid))-(double(vtScStrainsFirst)/double(totalStrainsFirst));
//        double vtSecondChange = (double(vtScStrainsFinal)/double(totalStrainsFinal))-(double(vtScStrainsMid)/double(totalStrainsMid));
//        double nvtFirstChange = (double(nvtScStrainsMid)/double(totalStrainsMid))-(double(nvtScStrainsFirst)/double(totalStrainsFirst));
//        double nvtSecondChange = (double(nvtScStrainsFinal)/double(totalStrainsFinal))-(double(vtScStrainsMid)/double(totalStrainsMid));
//        
//        // reset frequencies
//        totalStrainsFirst = 0;
//        vtScStrainsFirst = 0;
//        nvtScStrainsFirst = 0;
//        totalStrainsMid = 0;
//        vtScStrainsMid = 0;
//        nvtScStrainsMid = 0;
//        totalStrainsFinal = 0;
//        vtScStrainsFinal = 0;
//        nvtScStrainsFinal = 0;
//        // calculate frequencies in simulated population
//        for (int i = 0; i != numSc; i++) {
//            totalStrainsFirst+=(vtScFreq[0][i]+nvtScFreq[0][i]);
//            totalStrainsMid+=(vtScFreq[mGen][i]+nvtScFreq[mGen][i]);
//            totalStrainsFinal+=(vtScFreq[fGen][i]+nvtScFreq[fGen][i]);
//            if (i == (s-1)) {
//                vtScStrainsFirst+=vtScFreq[0][i];
//                nvtScStrainsFirst+=nvtScFreq[0][i];
//                vtScStrainsMid+=vtScFreq[mGen][i];
//                nvtScStrainsMid+=nvtScFreq[mGen][i];
//                vtScStrainsFinal+=vtScFreq[fGen][i];
//                nvtScStrainsFinal+=nvtScFreq[fGen][i];
//            }
//        }
//        double vtFirstChangeSim = (double(vtScStrainsMid)/double(totalStrainsMid))-(double(vtScStrainsFirst)/double(totalStrainsFirst));
//        double vtSecondChangeSim = (double(vtScStrainsFinal)/double(totalStrainsFinal))-(double(vtScStrainsMid)/double(totalStrainsMid));
//        double nvtFirstChangeSim = (double(nvtScStrainsMid)/double(totalStrainsMid))-(double(nvtScStrainsFirst)/double(totalStrainsFirst));
//        double nvtSecondChangeSim = (double(nvtScStrainsFinal)/double(totalStrainsFinal))-(double(nvtScStrainsMid)/double(totalStrainsMid));
//        // increment deviation stat
//        tmpDev+=pow((vtFirstChange-vtFirstChangeSim),2);
//        tmpDev+=pow((nvtFirstChange-nvtFirstChangeSim),2);
//        
////        std::cerr << "here sc" << s << " changes by " << vtFirstChange << " for vt/nvt " << nvtFirstChange << "; in Sim, " << vtFirstChangeSim << " for vt/nvt " << nvtFirstChangeSim << "; deviation is " << tmpDev << std::endl;
//        
//        tmpDev+=pow((vtSecondChange-vtSecondChangeSim),2);
//        tmpDev+=pow((nvtSecondChange-nvtSecondChangeSim),2);
//        
////        std::cerr << "here sc" << s << " changes by " << vtSecondChange << " for vt/nvt " << nvtSecondChange << "; in Sim, " << vtSecondChangeSim << " for vt/nvt " << nvtSecondChangeSim << "; deviation is " << tmpDev << std::endl << std::endl;
//    }
//    
//    // return final value of statistic
//    changeDev = tmpDev;
//    
//    return 0;
//}



//int parseOrderingFile(char* orderingFilename,std::vector<double> *cogWeights,std::vector<std::string> *cogList,struct parms *sp,std::vector<double> *eqFreq,std::vector<double> *midFreq,std::vector<double> *finalFreq,int &vtCog,std::vector<isolate*> *pop) {
//    
//    // initialise new cogList
//    std::vector<std::string> newCogList;
//    
//    // initialise new frequency lists
//    std::vector<double> newEqFreq;
//    std::vector<double> newMidFreq;
//    std::vector<double> newFinalFreq;
//    std::vector<int> reorderOrder;
//    
//    // open ordering file
//    std::ifstream oFile;
//    oFile.open(orderingFilename,std::ifstream::in);
//    
//    // keep track of VT COG
//    int newVtCog = -1;
//    
//    // parse weighting file
//    if (oFile.is_open()) {
//        std::string line;
//        while (std::getline(oFile, line)) {
//            // parse values from line
//            std::istringstream iss(line);
//            while (iss) {
//                std::string cogName;
//                while (getline(iss, cogName)) {
//                    int seen = 0;
//                    int newIndex = -1;
//                    for (int cindex = 0; cindex < cogList->size(); cindex++) {
//                        if ((*cogList)[cindex].compare(cogName) == 0) {
//                            seen = 1;
//                            newIndex = cindex;
//                            if (cindex == vtCog) {
//                                newVtCog = newCogList.size();
//                            }
//                        }
//                    }
//                    if (seen == 1 && newIndex != -1) {
//                        newCogList.push_back(cogName);
//                        newEqFreq.push_back((*eqFreq)[newIndex]);
//                        newMidFreq.push_back((*midFreq)[newIndex]);
//                        newFinalFreq.push_back((*finalFreq)[newIndex]);
//                        reorderOrder.push_back(newIndex);
//                    }
//                }
//            }
//        }
//        oFile.close();
//    } else {
//        std::cerr << "Unable to read file " << orderingFilename << std::endl;
//        return 1;
//    }
//    
//    // use new order
//    cogList->assign(newCogList.begin(),newCogList.end());
//    eqFreq->assign(newEqFreq.begin(),newEqFreq.end());
//    midFreq->assign(newMidFreq.begin(),newMidFreq.end());
//    finalFreq->assign(newFinalFreq.begin(),newFinalFreq.end());
//
//    // check VT COG is included
//    if (newVtCog != -1) {
//        vtCog = newVtCog;
//    } else {
//        std::cerr << "Reordering file does not include vaccine-type COG; please restart" << std::endl;
//        return 1;
//    }
//    
//    // reorder genotypes
//    std::vector<isolate*>::iterator iit;
//    for (iit = pop->begin(), pop->end() ; iit != pop->end(); ++iit) {
//        std::vector<bool> newGenotype;
//        for (int y = 0; y < reorderOrder.size(); y++) {
//            newGenotype.push_back((*iit)->genotype[reorderOrder[y]]);
//        }
//        (*iit)->genotype.assign(newGenotype.begin(),newGenotype.end());
//    }
//    
//    // assign weights appropriately
//    for (int cindex = 0; cindex < cogList->size(); cindex++) {
//        if ((float(cindex)/float(cogList->size())) <= sp->selectedProp) {
//            (*cogWeights)[cindex] = sp->higherSelection;
//        } else {
//            (*cogWeights)[cindex] = sp->lowerSelection;
//        }
//    }
//    
//    return 0;
//}



///////////////
// Graveyard //
///////////////


//////////////////////////////////
// calculate summary statistics //
//////////////////////////////////

//int getCogFreq(std::vector<isolate*> *pop,int year,std::vector<int> *freq,int& strains) {
//    // iterate through isolates in population
//    std::vector<isolate*>::iterator iter;
//    for (iter = pop->begin(), pop->end() ; iter != pop->end(); ++iter) {
//        if ((*iter)->year == year) {
//            std::transform((*iter)->genotype.begin(), (*iter)->genotype.end(), freq->begin(), freq->begin(), std::plus<int>());
//            strains++;
//        }
//    }
//    return 0;
//}

///////////////////////////////
// calculate cog proportions //
///////////////////////////////

//int getProp(std::vector<isolate*> *pop,struct parms *sp,std::vector<int> *starting,std::vector<int> *middling,std::vector<int> *ending,std::vector<double> *ef,std::vector<double> *mf,std::vector<double> *ff,std::vector<std::string> *cogList,int &vc, int ns, int &ng,int ss, int ms, int es) {
//    std::vector<int>::iterator iter;
//    // store names of informative genes
//    std::vector<std::string> tmpCogNames;
//    // debug
//    const int numberOfGenes = ng;
//    bool includeGenes[numberOfGenes];
//    int index = 0;
//    for (iter = starting->begin(), starting->end() ; iter != starting->end(); ++iter) {
//        // check on starting criteria
//        if ((double)(*iter)/(double)ss > sp->lowerLimit && (double)(*iter)/(double)ss < sp->upperLimit) {
//            includeGenes[index] = 1;
//        } else {
//            includeGenes[index] = 0;
//        }
//        index++;
//    }
//    // ensure VT COG always included
//    includeGenes[vc] = 1;
//    // now eliminate all unnecessary loci
//    int newvc = vc;
//    int rIndex = 0;
//    // get new list of COG names
//    index = 0;
//    std::vector<std::string>::iterator siter;
//    for (siter = cogList->begin(), cogList->end() ; siter != cogList->end(); ++siter) {
//        if (includeGenes[index] == 1) {
//            tmpCogNames.push_back(*siter);
//        }
//        index++;
//    }
//    (*cogList) = tmpCogNames;
//    tmpCogNames.clear();
//    // reduce the starting COGs array and find equilibrium frequencies
//    index = 0;
//    std::vector<int> newStarting;
//    for (iter = starting->begin(), starting->end() ; iter != starting->end(); ++iter) {
//        if (index == vc) {
//            newvc = rIndex;
//        }
//        if (includeGenes[index] == 1) {
//            newStarting.push_back((*iter));
//            ef->push_back((double(*iter)/double(ss)));
//
//            rIndex++;
//        }
//        index++;
//    }
//    vc = newvc;
//    starting = &newStarting;
//    ng = newStarting.size();
//    // record midpoint COG proportions and reduce the middle COGs array
//    index = 0;
//    std::vector<int> newMiddle;
//    for (iter = middling->begin(), middling->end() ; iter != middling->end(); ++iter) {
//        if (includeGenes[index] == 1) {
//            newMiddle.push_back((*iter));
//            mf->push_back((double(*iter)/double(ms)));
//        }
//        index++;
//    }
//    middling = &newMiddle;
//    // record final COG proportions and reduce the ending COGs array
//    index = 0;
//    std::vector<int> newEnding;
//    for (iter = ending->begin(), ending->end() ; iter != ending->end(); ++iter) {
//        if (includeGenes[index] == 1) {
//            newEnding.push_back((*iter));
//            ff->push_back((double(*iter)/double(es)));
//        }
//        index++;
//    }
//    ending = &newEnding;
//    // reduce the genotypes for each isolate
//    std::vector<isolate*>::iterator iit;
//    for (iit = pop->begin(), pop->end() ; iit != pop->end(); ++iit) {
//        std::vector<bool> newGenotype;
//        std::vector<bool>::iterator git;
//        index = 0;
//        for (git = (*iit)->genotype.begin(), (*iit)->genotype.end() ; git != (*iit)->genotype.end(); ++git) {
//            if (includeGenes[index] == 1) {
//                newGenotype.push_back(*git);
//            }
//            index++;
//        }
//        (*iit)->genotype.clear();
//        for (git = newGenotype.begin(), newGenotype.end() ; git != newGenotype.end(); ++git) {
//            (*iit)->genotype.push_back(*git);
//        }
//    }
//    // completed
//    return 0;
//}

//////////////////////////////////////////////
// calculate COG deviation from equilibrium //
//////////////////////////////////////////////
//
//int getCogDeviations(std::vector<double> * ef,std::vector<isolate*> *currentStrains,std::vector<double> *cogWeights,std::vector<double> *cogDeviations,std::vector<double> *startingCogFrequencies,std::vector<double> &startingVtScFrequencies,std::vector<double> &startingNvtScFrequencies,int initialSampleSize) {
//
//    std::vector<double> currentPrev((*currentStrains)[0]->genotype.size(),0.0);
//    std::vector<int>::iterator cogIter;
//    std::vector<isolate*> initialPopulation = *currentStrains;
//    std::vector<isolate*> initialPopulationSample;
//    if (initialSampleSize <= 0) {
//        std::cerr << "Zero sample size at first timestep; need to alter input timescale" << std::endl;
//        return 1;
//    }
//
//    while (initialPopulationSample.size() < initialSampleSize) {
//        // sample isolates at random
//        int selection = rand()%initialPopulation.size();
//        isolate selectedIsolate = (*initialPopulation[selection]);
//        initialPopulationSample.push_back(&selectedIsolate);
//        // record sequence cluster frequencies
//        if (selectedIsolate.vt) {
//            startingVtScFrequencies[selectedIsolate.sc-1]+=(1.0/double(initialSampleSize));
//        } else {
//            startingNvtScFrequencies[selectedIsolate.sc-1]+=(1.0/double(initialSampleSize));
//        }
//        // calculate gene frequencies
//        for (int i = 0; i < selectedIsolate.genotype.size();i++) {
//            currentPrev[i]+=(double(selectedIsolate.genotype[i])/double(initialSampleSize));
//        }
//    }
//
//    (*startingCogFrequencies) = currentPrev;
//    std::transform(ef->begin(), ef->end(), currentPrev.begin(), cogDeviations->begin(), std::minus<double>());
//    std::transform(cogWeights->begin(), cogWeights->end(), cogDeviations->begin(), cogDeviations->begin(), std::multiplies<double>());
//
//    return 0;
//}