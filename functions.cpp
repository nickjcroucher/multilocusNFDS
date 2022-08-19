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
    
    std::cerr << "Programme:" << fn << std::endl << "\tp\tprogramme type:" << std::endl << "\t\t's' - simulation without fitting" << std::endl << "\t\t'x' - extended simulation output" << std::endl << "\t\t'f' - just return fitting metrics" << std::endl << "\t\t'b' - both fit to data and print simulation" << std::endl << std::endl << "Simulation parameters:" << std::endl << "\ts\tfrequency dependent selection pressure [double between 0 and 1]" << std::endl << "\tv\tvaccine selection pressure  [double between 0 and 1]" << std::endl << "\ti\timmigration rate [double between 0 and 1]" << std::endl << "\tt\timmigration type [0 - by strain, 1 - by SC]" << std::endl << "\tn\tpopulation carrying capacity  [integer]" << std::endl << "\tg\tnumber of generations [integer]" << std::endl << "\tu\tupper gene frequency limit [double between 0 and 1]" << std::endl << "\tl\tlower gene frequency limit [double between 0 and 1]" << std::endl << "\tq\tgeneration in which vaccine formulation is changed" << std::endl << std::endl << "Model fitting parameters:" << std::endl << "\tc\tvaccine target COG name [string]" << std::endl << "\tb\tbeginning year [integer]" << std::endl << "\tm\tmid year [integer]" << std::endl << "\te\tending year [integer]" << std::endl << std::endl << "Filenames:" << std::endl << "\tf\tinputFilename" << std::endl << "\tx\tfrequency file name [only for simulating]" << std::endl << "\to\toutput file prefix" << std::endl << "\tw\tweighting file" << std::endl  << "\tr\tcog reordering file" << std::endl;
    
}

////////////////////////
// Input file parsing //
////////////////////////

//////////////////////
// parse input file //
//////////////////////

int parseInputFile(std::vector<isolate*> *pop, std::vector<cog*> *accessoryLoci, double lower, double upper,  std::vector<int> *samplingList,std::vector<std::string> *st, std::vector<int> *sc, std::vector<std::string> *cogList, char *inputFilename, char * vtCogName, int &minGen, bool useCogList) {
    
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
            bool sample_vt = 0;
            bool sample_latent_vt = 0;
            bool partial_vt = 0;
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
                        int vt_int = atoi(temp.c_str());
                        if (vt_int == 0) {
                            sample_vt = false;
                            sample_latent_vt = false;
                            partial_vt = false;
                        } else  if (vt_int == 1) {
                            sample_vt = true;
                            sample_latent_vt = false;
                            partial_vt = false;
                        } else if (vt_int == 2) {
                            sample_vt = false;
                            sample_latent_vt = true;
                            partial_vt = false;
                        } else  if (vt_int == 3) {
                            sample_vt = true;
                            sample_latent_vt = false;
                            partial_vt = true;
                        } else if (vt_int == 4) {
                            sample_vt = false;
                            sample_latent_vt = true;
                            partial_vt = true;
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
                isolate* tmp = new isolate(sample_id,
                                           sample_time,
                                           sample_sc,
                                           sample_serotype,
                                           sample_vt,
                                           sample_latent_vt,
                                           partial_vt,
                                           &sample_genotype,
                                           &sample_markers,
                                           1.0);
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
            if ((*iter)-minGen > maxTime) {
                maxTime = (*iter)-minGen;
            }
        }
        
        samplingList->resize(maxTime+1,0);
        for (iter = samplingTimes.begin(), samplingTimes.end() ; iter != samplingTimes.end(); ++iter) {
            (*samplingList)[(*iter)-minGen]++;
        }

        // create vector of accessory locus COG objects
        if (useCogList) {
            
            // find matching positions between pre-calculated coglist
            // and new input file
            std::vector<int> cogMatches;
            if (cogMatches.size() < cogList->size()) {
                for (unsigned int i = 0; i < cogList->size(); i++) {
                    int match_pos = -1;
                    for (unsigned int j = 0; j < tmpCogList.size(); j++) {
                        if ((*cogList)[i] == tmpCogList[j]) {
                            match_pos = j;
                            j = tmpCogList.size();
                        }
                    }
                    if (match_pos == -1) {
                        std::cerr << "Unable to find COG " << (*cogList)[i] << " in migrant strain input file" << std::endl;
                        return 1;
                    } else {
                        cogMatches.push_back(match_pos);
                    }
                }
            }
            
            // add ordered cogs to input file
            std::vector<isolate*>::iterator iiter;
            for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                std::vector<bool> tmpGenotype;
                for (unsigned int j = 0; j < cogMatches.size(); j++) {
                    tmpGenotype.push_back((*iiter)->genotype[j]);
                }
                (*iiter)->genotype = tmpGenotype;
            }
            
        } else {
            std::vector<int> includeLocus(tmpCogList.size(),0);
            std::vector<isolate*>::iterator iiter;
            std::vector<std::string> intCogList;
            for (unsigned int i = 0; i < tmpCogList.size(); i++) {
                // data structures for recording frequencies
                int overallFreq = 0;
                std::vector<double> cogFrequencies(samplingTimes.size(),0);
                double eqFreq = 0;
                // calculate gene frequencies from isolates
                for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                    overallFreq+=(*iiter)->genotype[i];
                    cogFrequencies[((*iiter)->year)-minGen]+=(double((*iiter)->genotype[i])/double((*samplingList)[((*iiter)->year)-minGen]));
                    if ((*iiter)->year <= 0) {
                        eqFreq+=(double((*iiter)->genotype[i])/double(eqPop));
                    }
                }
                // retain if present at intermediate frequency OR vt-defining COG
                if ((eqFreq >= (lower-1e-07) && eqFreq <= (upper+1e-07)) || i == unsigned(v)) {
                    includeLocus[i] = 1;
                    int vtType = 0;
                    if (i == unsigned(v)) {
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
                for (unsigned int i = 0; i < includeLocus.size(); i++) {
                    if (includeLocus[i] == 1) {
                        tmpGenotype.push_back((*iiter)->genotype[i]);
                    }
                }
                (*iiter)->genotype = tmpGenotype;
            }
        
            // return cog list
            *cogList = intCogList;
        }
        
    } else {
        std::cerr << "Problem with input file: " << strerror(errno) << std::endl;
        return 1;
    }
    
    return 0;
}

//////////////////////////////
// parse marker information //
//////////////////////////////

int parseMarkerFile(std::vector<isolate*> *pop,char *markerFilename,std::vector<std::string> *markerList, bool useMarkerList) {
    
    // data structures
    std::vector<std::string> tmpMarkerList;
    
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
                            tmpMarkerList.push_back(temp);
                        } else {
                            sample_markers.push_back(atoi(temp.c_str()));
                        }
                    }
                    sIndex++;
                }
            }
            // record marker genotypes and attach to isolate objects
            std::vector<isolate*>::iterator iiter;
            for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                if ((*iiter)->id == sample_id) {
                    (*iiter)->markers = sample_markers;
                }
            }

        }
        infile.close();
        
        
        // check whether markers need reordering
        if (useMarkerList) {
            
            // find matching positions between pre-calculated coglist
            // and new input file
            
            // data structure
            std::vector<int> markerMatches;
            if (markerMatches.size() < markerList->size()) {
                for (unsigned int i = 0; i < markerList->size(); i++) {
                    int match_pos = -1;
                    for (unsigned int j = 0; j < tmpMarkerList.size(); j++) {
                        if ((*markerList)[i] == tmpMarkerList[j]) {
                            match_pos = j;
                            j = tmpMarkerList.size();
                        }
                    }
                    if (match_pos == -1) {
                        std::cerr << "Unable to find marker " << (*markerList)[i] << " in migrant strain input file" << std::endl;
                        return 1;
                    } else {
                        markerMatches.push_back(match_pos);
                    }
                }
            }
            
            // add ordered markers to correct isolate
            std::vector<isolate*>::iterator iiter;
            for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                std::vector<bool> tmpGenotype;
                std::vector<bool> sample_markers = (*iiter)->markers;
                for (unsigned int j = 0; j < markerMatches.size(); j++) {
                    tmpGenotype.push_back(sample_markers[j]);
                }
                (*iiter)->markers = tmpGenotype;
            }
            
        } else {
            // modify marker list
            (*markerList) = tmpMarkerList;
            // modify individual isolates
        }

        // check all marker lengths are the same
        std::vector<isolate*>::iterator iiter;
        long markerLength = tmpMarkerList.size();
        for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
            if (unsigned(markerLength) != (*iiter)->markers.size()) {
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

//////////////////////////
// validate input files //
//////////////////////////

int compareInputPopulations(std::vector<isolate*> *popA, std::vector<isolate*> *popB, bool check_markers) {

    long genotype_length = -1;
    std::vector<isolate*>::iterator iter;
    
    // check populationA genotype lengths are consistent
    for (iter = popA->begin(), popA->end() ; iter != popA->end(); ++iter) {
        long current_genotype_size = (*iter)->genotype.size(); // explicit cast for precision reasons
        if (genotype_length != -1 && genotype_length != current_genotype_size) {
            std::cerr << "Inconsistent genotype length for isolate " << (*iter)->id << " in first population" << std::endl;
            return 1;
        } else {
            genotype_length = (*iter)->genotype.size();
        }
    }

    // check populationB genotype lengths are consistent
    for (iter = popB->begin(), popB->end() ; iter != popB->end(); ++iter) {
        long current_genotype_size = (*iter)->genotype.size();
        if (genotype_length != current_genotype_size) {
            std::cerr << "Inconsistent genotype length for isolate " << (*iter)->id << " in second population - may not match that in the first population" << std::endl;
            return 1;
        }
    }
    
    if (check_markers) {
        
        genotype_length = -1;
        
        // check populationA genotype lengths are consistent
        for (iter = popA->begin(), popA->end() ; iter != popA->end(); ++iter) {
            long current_marker_size = (*iter)->markers.size();
            if (genotype_length != -1 && genotype_length != current_marker_size) {
                std::cerr << "Inconsistent marker genotype length for isolate " << (*iter)->id << " in first population" << std::endl;
                return 1;
            } else {
                genotype_length = (*iter)->markers.size();
            }
        }
        
        // check populationB genotype lengths are consistent
        for (iter = popB->begin(), popB->end() ; iter != popB->end(); ++iter) {
            long current_marker_size = (*iter)->markers.size();
            if (genotype_length != current_marker_size) {
                std::cerr << "Inconsistent marker genotype length for isolate " << (*iter)->id << " in second population - may not match that in the first population" << std::endl;
                return 1;
            }
        }
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
        if (sp->immigrationType != 0 && sp->immigrationType != 1 && sp->immigrationType != 2 && sp->immigrationType != 3) {
            std::cerr << "Invalid immigration type: " << sp->immigrationType << "; should be either '0' (by isolate),  '1' (by sc), or '2' (by time)" << std::endl;
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
    if (sp->vaccineLag < 0) {
        std::cerr << "Vaccine lag (in months) must be equal to or greater than zero" << std::endl;
        tmpvalid = 0;
    }
    if (sp->nfdsLag < 1) {
        std::cerr << "NFDS lag (in months) must be equal to or greater than one" << std::endl;
        tmpvalid = 0;
    }
    // check recombination parameters
    if (sp->transformationProportion*sp->transformationRate == 0 && sp->transformationProportion+sp->transformationRate > 0) {
        std::cerr << "Warning! Need to set transformation proportion (z) and transformation rate (e) both greater than zero for recombination to occur" << std::endl;
        tmpvalid = 0;
    } else if (sp->transformationProportion+sp->transformationRate > 0) {
        if (sp->transformationAsymmetryLoci < 0) {
            std::cerr << "Transformation asymmetry needs to be be greater than 0" << std::endl;
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
        for (unsigned int j = 0; j < orderedAccessoryLoci.size(); j++) {
            std::cerr << orderedAccessoryLoci[j] << "\t";
            for (unsigned int i = 0; i < accessoryLoci->size(); i++) {
                if ((*accessoryLoci)[i]->id == orderedAccessoryLoci[j]) {
                    std::cerr << (*accessoryLoci)[i]->id;
                }
            }
            std::cerr << std::endl;
        }
        std::cerr << "Alternative ordering:" << std::endl;
        for (unsigned int i = 0; i < accessoryLoci->size(); i++) {
            std::cerr << (*accessoryLoci)[i]->id << "\t";
            for (unsigned int j = 0; j < orderedAccessoryLoci.size(); j++) {
                if ((*accessoryLoci)[i]->id == orderedAccessoryLoci[j]) {
                    std::cerr << orderedAccessoryLoci[j] << std::endl;
                }
            }
        }
        return 1;
    }
    
    // assign weights appropriately - most conserved frequencies are at the END of the list
    unsigned int cindex = 0;
    if (sp->het_mode == "l") { // logistic
        double max_value = sp->lowerSelection + 1.0 / (1.0 + exp(-1.0 * sp->decayRate * (float(accessoryLoci->size()) - (1 - sp->selectedProp) * float(accessoryLoci->size()))));
        for (unsigned int cindex = 0; cindex < orderedAccessoryLoci.size(); cindex++) {
            for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                if ((*cit)->id.compare(orderedAccessoryLoci[cindex]) == 0) {
                    double relative_weight_from_order = (sp->lowerSelection + 1.0 / (1.0 + exp(-1.0 * sp->decayRate * (float(cindex) + 1.0 - (1 - sp->selectedProp) * float(accessoryLoci->size()))))) * (1 / max_value);
                    (*cit)->weight = sp->higherSelection * relative_weight_from_order;
                    break;
                }
            }
        }
    } else if (sp->het_mode == "e") { // exponential
        double max_value = sp->lowerSelection + exp(-1 * sp->decayRate * 0.0);
        for (unsigned int cindex = 0; cindex < orderedAccessoryLoci.size(); cindex++) {
            for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                if ((*cit)->id.compare(orderedAccessoryLoci[cindex]) == 0) {
                    double relative_weight_from_order = (sp->lowerSelection + exp(-1 * sp->decayRate * float(accessoryLoci->size() - cindex - 1))) * (1 / max_value);
                    (*cit)->weight = sp->higherSelection * relative_weight_from_order;
                    break;
                }
            }
        }
    } else if (sp->het_mode == "r") { // linear
        double max_value = sp->lowerSelection + 1;
        for (unsigned int cindex = 0; cindex < orderedAccessoryLoci.size(); cindex++) {
            for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                if ((*cit)->id.compare(orderedAccessoryLoci[cindex]) == 0) {
                    double relative_weight_from_order = sp->lowerSelection + 1 - (float(accessoryLoci->size() - cindex - 1) * sp->decayRate);
                    relative_weight_from_order = relative_weight_from_order * (1 / max_value);
                    if (relative_weight_from_order < 0.0) {
                        relative_weight_from_order = 0.0;
                    }
                    (*cit)->weight = sp->higherSelection * relative_weight_from_order;
                    break;
                }
            }
        }
    } else if (sp->het_mode == "s") { // step
        for (unsigned int cindex = 0; cindex < orderedAccessoryLoci.size(); cindex++) {
            for (cit = accessoryLoci->begin(), accessoryLoci->end(); cit != accessoryLoci->end(); ++cit) {
                if ((*cit)->id.compare(orderedAccessoryLoci[cindex]) == 0) {
                    if ((float(cindex)/float(orderedAccessoryLoci.size())) <= (1.0-sp->selectedProp)) {
                        (*cit)->weight = sp->lowerSelection;
                    } else {
                        (*cit)->weight = sp->higherSelection;
                    }
                    break;
                }
            }
        }
    } else {
        std::cerr << "Allowed heterogeneous functions are (l)ogistic, (e)xponential, linea(r) and (s)tep" << std::endl;
        exit(1);
    }
    return 0;
}

////////////////////////////////
// Pre-processing information //
////////////////////////////////

///////////////////////////
// Generate migrant pool //
///////////////////////////

int generateMigrantPool(std::vector<std::vector<std::vector<isolate*> > > *migrantPool, std::vector<isolate*> *population, std::vector<isolate*> *migrant_population, char* migrantFilename, std::vector<int> *scList, int maxScNum, int minGen,struct parms *p) {
    
    // Split population for immigration by SC
    if (p->immigrationType == 1) {
        
        int divCheck = 1;
        std::vector<std::vector<isolate*> > *populationBySc = new std::vector<std::vector<isolate*> >;
        if (migrantFilename != NULL) {
            divCheck = dividePopulationForImmigration(migrant_population,scList,populationBySc,maxScNum);
        } else {
            divCheck = dividePopulationForImmigration(population,scList,populationBySc,maxScNum);
        }
        if (divCheck != 0) {
            std::cerr << "Unable to split population into sequence clusters" << std::endl;
            return 1;
        }
        migrantPool->push_back(*populationBySc);
    // Split by time
    } else if (p->immigrationType == 2) {
        // split population for immigration by time
        int divCheck = 1;
        std::vector<std::vector<isolate*> > *populationByTime = new std::vector<std::vector<isolate*> >;
        if (migrantFilename != NULL) {
            divCheck = dividePopulationForImmigrationByTime(migrant_population,minGen,p->numGen,populationByTime);
        } else {
            divCheck = dividePopulationForImmigrationByTime(population,minGen,p->numGen,populationByTime);
        }
        if (divCheck != 0) {
            std::cerr << "Unable to split population by isolation times" << std::endl;
            return 1;
        }
        migrantPool->push_back(*populationByTime);
        
    // Split by time and strain
    } else if (p->immigrationType == 3) {
        // split population for immigration by time
        int divCheck = 1;
        std::vector<std::vector<isolate*> > *populationByTime = new std::vector<std::vector<isolate*> > (p->numGen+1, std::vector<isolate*>());
        if (migrantFilename != NULL) {
            divCheck = dividePopulationForImmigrationByTime(migrant_population,minGen,p->numGen,populationByTime);
        } else {
            divCheck = dividePopulationForImmigrationByTime(population,minGen,p->numGen,populationByTime);
        }
        if (divCheck != 0) {
            std::cerr << "Unable to split population by isolation times" << std::endl;
            return 1;
        }
        // then split each time point by strain
        for (int g = 0; g <= p->numGen; g++) {
            std::vector<std::vector<isolate*> > *populationByTimeAndSc = new std::vector<std::vector<isolate*> >;
            if ((*populationByTime)[g].size() >= 1) {
                std::vector<isolate*> *tmpStrains = new std::vector<isolate*>;
                tmpStrains = &(*populationByTime)[g]; // needs fixing
                divCheck = dividePopulationForImmigration(tmpStrains,scList,populationByTimeAndSc,maxScNum);
                if (divCheck != 0) {
                    std::cerr << "Unable to split population by SC for time " << g << std::endl;
                    return 1;
                } else {
                    migrantPool->push_back(*populationByTimeAndSc);
                }
            } else {
                migrantPool->push_back(*populationByTimeAndSc);
            }
        }
    // Do not split, just sample at random
    } else {
        std::vector<std::vector<isolate*> > *tmpStrains = new std::vector<std::vector<isolate*> >;
        if (migrantFilename != NULL) {
            tmpStrains->push_back(*migrant_population);
        } else {
            tmpStrains->push_back(*population);
        }
        migrantPool->push_back(*tmpStrains);
    }
    
    return 0;
}


///////////////////////////////////////////
// divide isolates by SC for immigration //
///////////////////////////////////////////

int dividePopulationForImmigration(std::vector<isolate*> *pop,std::vector <int> *scList,std::vector<std::vector<isolate*> > *popBySc, int maxScNum) {
    
    // check there are > 0 sequence clusters
    if (maxScNum == 0) {
        std::cerr << "Cannot find any sequence clusters" << std::endl;
    }
    
    std::vector<std::vector<isolate*> > tmpStrains(scList->size());
    
    for (unsigned int s = 0; s < scList->size(); s++) {
        std::vector<isolate*>::iterator cit;
        for (cit = pop->begin(), pop->end(); cit != pop->end(); ++cit) {
            if ((*cit)->sc == (*scList)[s]) {
                tmpStrains[s].push_back((*cit));
            }
        }
    }
    
    (*popBySc) = tmpStrains;
    
    return 0;
}

/////////////////////////////////////////////
// divide isolates by time for immigration //
/////////////////////////////////////////////

int dividePopulationForImmigrationByTime(std::vector<isolate*> *pop, int minGen, int numGen,std::vector<std::vector<isolate*> > *popByTime) {
    

    std::vector<std::vector<isolate*> > tmpStrains(numGen+1);

    std::vector<isolate*>::iterator cit;
    for (cit = pop->begin(), pop->end(); cit != pop->end(); ++cit) {
        for (int t = 0; t <= numGen; t++) {
            if (((*cit)->year)-minGen <= numGen) { // && ((*cit)->year)-minGen >= 0) {
                if (((*cit)->year)-minGen == t) {
                    tmpStrains[t].push_back((*cit));
                }
            }
        }
    }
    
    (*popByTime) = tmpStrains;

    return 0;
    
}

///////////////////////////
// get first year sample //
///////////////////////////

int getStartingIsolates(std::vector<isolate*> *pop,struct parms *sp,std::vector<isolate*> *first,std::vector<cog*> *accessoryLoci,int psize,std::vector<double> &eqFreq,std::vector<double> &cogWeights,std::vector<double> &cogDeviations,std::vector<int> &startingVtScFrequencies,std::vector<int> &startingNvtScFrequencies,std::vector<int> *scList, int minGen, float seedStartingPopulation, char* migrantFilename, std::vector<isolate*> *migrant_population, int maxScNum) {
    
    // get all isolates observed in the pre- or peri-vaccine samples
    std::vector<isolate*> *possibleFirst = new std::vector<isolate*>;
    std::vector<isolate*> *possibleFirst_unsampled = new std::vector<isolate*>;
    std::vector<isolate*>::iterator iter;
    std::vector<int> observedVtSc;
    std::vector<int> observedNvtSc;
    
    int first_sample_size = 0;
    for (iter = pop->begin(), pop->end() ; iter != pop->end(); ++iter) {
        if (minGen < 0) {
            // Use pre-vaccine population if possible
            if ((*iter)->year < 0) {
                possibleFirst->push_back(*iter);
                first_sample_size++;
            } else {
                possibleFirst_unsampled->push_back(*iter);
            }
        } else {
            // if no pre-vaccine population, use the peri-vaccination population
            if ((*iter)->year == 0) {
                possibleFirst->push_back(*iter);
                first_sample_size++;
            } else {
                possibleFirst_unsampled->push_back(*iter);
            }
        }
    }
    
    // add in genotypes not detected in first sample if seeding unsampled genotypes
    if (seedStartingPopulation > 1e-6) {
        
        // data structure for the seeding genotypes
        std::vector< std::vector<isolate*> > *isolates_for_seeding = new std::vector<std::vector<isolate*> >(scList->size(),std::vector<isolate*>());
        
        // calculate the number of unsampled bacteria to add
        // seedStartingPopulation is the probability of not sampling each undetected genotype
        float upper_freq = 1.0 - exp(log(seedStartingPopulation)/float(first_sample_size));
        int number_of_unsampled_bacteria = round(upper_freq*psize);
        std::vector< int > unseen_sc;
        
        // if simplest migration type, select isolates from later generations
        if (sp->immigrationType == 0) {
            
            (*isolates_for_seeding)[0] = *possibleFirst_unsampled;
//            isolates_for_seeding[0] = possibleFirst_unsampled;
//            unseen_sc.push_back(0);
            
        } else {
            
            // migrationType == 1 - divide by strain and seed from any later timeperiod
            // works with and without migration file
            if (sp->immigrationType == 1) {
                // divide population by strain
                int divCheck = 1;
                std::vector<std::vector<isolate*> > *isolatesBySc = new std::vector<std::vector<isolate*> >;
                if (migrantFilename != NULL) {
                    divCheck = dividePopulationForImmigration(migrant_population,scList,isolatesBySc,maxScNum);
                } else {
                    divCheck = dividePopulationForImmigration(pop,scList,isolatesBySc,maxScNum);
                }
                isolates_for_seeding = isolatesBySc;
                if (divCheck != 0) {
                    std::cerr << "Unable to split population into sequence clusters" << std::endl;
                    return 1;
                }

            // Split by time
            } else if (sp->immigrationType == 2) {
                // split population for immigration by time
                std::vector<std::vector<isolate*> > *isolatesByTime = new std::vector<std::vector<isolate*> >;
                int divCheck = 1;
                if (migrantFilename == NULL) {
                    std::cerr << "Need a separate migration file when seeding the initial population with migration mode 2" << std::endl;
                    return 1;
                } else {
                    divCheck = dividePopulationForImmigrationByTime(migrant_population,minGen,sp->numGen,isolatesByTime);
                }
                if (divCheck != 0) {
                    std::cerr << "Unable to split population by isolation times" << std::endl;
                    return 1;
                }
                // add to main data structure here
                if (minGen < 0) {
                    for (int gen = 0; gen < (0-minGen); gen++) {
                        if ((*isolatesByTime)[gen].size() > 0) {
                            for (int x = 0; x < (*isolatesByTime)[gen].size(); x++) {
                                (*isolates_for_seeding)[0].push_back((*isolatesByTime)[gen][x]);
                            }
                        }
                    }
                } else {
                    (*isolates_for_seeding)[0] = (*isolatesByTime)[0];
                }
                // Check possible seed isolates were found
                if ((*isolates_for_seeding)[0].size() == 0) {
                    std::cerr << "No possible seed isolates found" << std::endl;
                    return 1;
                }
            // Split by time and strain
            } else if (sp->immigrationType == 3) {
                // split population for immigration by time
                int divCheck = 1;
                std::vector<std::vector<isolate*> > *populationByTime = new std::vector<std::vector<isolate*> > (sp->numGen+1, std::vector<isolate*>());
                if (migrantFilename == NULL) {
                    std::cerr << "Need a separate migration file when seeding the initial population with migration mode 3" << std::endl;
                    return 1;
                } else {
                    divCheck = dividePopulationForImmigrationByTime(migrant_population,minGen,sp->numGen,populationByTime);
                }
                if (divCheck != 0) {
                    std::cerr << "Unable to split population by isolation times" << std::endl;
                    return 1;
                }
                // Split the initial population
                std::vector<std::vector<std::vector<isolate*> > > *isolatesByTimeAndSc = new std::vector<std::vector<std::vector<isolate*> > >;
                // then split each time point by strain
                for (int g = 0; g <= sp->numGen; g++) {
                    std::vector<std::vector<isolate*> > *isolatesBySc = new std::vector<std::vector<isolate*> >;
                    if ((*populationByTime)[g].size() >= 1) {
                        std::vector<isolate*> *tmpStrains = new std::vector<isolate*>;
                        tmpStrains = &(*populationByTime)[g]; // needs fixing
                        divCheck = dividePopulationForImmigration(tmpStrains,scList,isolatesBySc,maxScNum);
                        if (divCheck != 0) {
                            std::cerr << "Unable to split population->by SC for time " << g << std::endl;
                            return 1;
                        }
                        else {
                            isolatesByTimeAndSc->push_back(*isolatesBySc);
                        }
                    }
                    else {
                        isolatesByTimeAndSc->push_back(*isolatesBySc);
                    }
                }
                if (divCheck != 0) {
                    std::cerr << "Unable to split population by SC for time " << std::endl;
                    return 1;
                }
                // add to main data structure here
                if (minGen < 0) {
                    for (int gen = 0; gen < (0-minGen); gen++) {
                        if ((*isolatesByTimeAndSc)[gen].size() > 0) {
                            for (int strain_index = 0; strain_index < scList->size(); strain_index++) {
                                for (int x = 0; x < (*isolatesByTimeAndSc)[gen][strain_index].size(); x++) {
                                    isolate *isolate_for_seeding = (*isolatesByTimeAndSc)[gen][strain_index][x];
                                    (*isolates_for_seeding)[strain_index].push_back(isolate_for_seeding);
                                    
                                }
                            }
                        }
                    }
                } else {
                    *isolates_for_seeding = (*isolatesByTimeAndSc)[0];
                }
                // Check possible seed isolates were found
                int max_strain_size = 0;
                for (int strain_index = 0; strain_index < (*isolates_for_seeding).size(); strain_index++) {
                    if ((*isolates_for_seeding)[strain_index].size() > max_strain_size) {
                        max_strain_size = (*isolates_for_seeding)[strain_index].size();
                        break;
                    }
                }
                if (max_strain_size == 0) {
                    std::cerr << "No possible seed isolates found" << std::endl;
                    return 1;
                }
            }
            
            // iterate through strains - only add in strains that were not observed in the starting
            // population - in modes 1 and 3
            if (sp->immigrationType == 1 || sp->immigrationType == 3) {
                std::vector<isolate*>::iterator first_iter;
                for (int strain_index = 0; strain_index < scList->size(); strain_index++) {
                    // first check whether the strain has been observed at the starting timepoint
                    int sc = (*scList)[strain_index];
                    bool seen = 0;
                    for (first_iter = possibleFirst->begin(), possibleFirst->end() ; first_iter != possibleFirst->end(); ++first_iter) {
                        if ((*first_iter)->sc == sc) {
                            seen = 1;
                            break;
                        }
                    }
                    
                    // to maintain consistency with definition above, use all pre-vaccine strains
                    // when such samples are available
                    // only consider the cases where there are multiple pre-vaccine generations
    //                if (minGen < -1) {
                        
    //                    for (int pv_gen = 1; pv_gen < (0-minGen); pv_gen++) {
    //                        if (isolates_for_seeding[pv_gen].size() > 0) {
    //                            if (!isolates_for_seeding[pv_gen][strain_index].empty()) {
    //                                for (int isolate_index = 0; isolate_index < isolates_for_seeding[pv_gen][strain_index].size(); isolate_index++) {
    //                                    isolates_for_seeding[0][strain_index].push_back(isolates_for_seeding[pv_gen][strain_index][isolate_index]);
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
                    
                    // then add unseen strains
                    if (!seen) {
                        // check whether strain is observed in later timesteps
                        if ((*isolates_for_seeding)[strain_index].size() > 0) {
                            unseen_sc.push_back(strain_index);
                        }
                    }
                }
            } else {
                unseen_sc.push_back(0);
            }
            
            // iterate up to the determined sample size
            for (int index = 0; index < unseen_sc.size(); index++) {
                int unseen_strain_index = unseen_sc[index];
                int unseen_sc_number = (*scList)[unseen_strain_index];
                for (int unsampled_index = 0; unsampled_index < number_of_unsampled_bacteria; unsampled_index++) {
                    int selection = int(double(gsl_rng_uniform(rgen))*int((*isolates_for_seeding)[unseen_strain_index].size()));
                    isolate *selected_isolate = (*isolates_for_seeding)[unseen_strain_index][selection];
                    first->push_back(selected_isolate);
                    // record sequence clusters
                    if (selected_isolate->vt) {
                        observedVtSc.push_back(selected_isolate->sc);
                    } else {
                        observedNvtSc.push_back(selected_isolate->sc);
                    }
                    // calculate gene frequencies
                    for (unsigned int i = 0; i < selected_isolate->genotype.size();i++) {
                        (*accessoryLoci)[i]->simFreq[0]+=(double(selected_isolate->genotype[i])/double(psize));
                    }
                }
            }
        }
    }
    
    // fill first timepoint with random sample of isolates from pre-/peri-vaccination samples
    // record starting COG frequencies
    while (first->size() < unsigned(psize)) {
        //int selection = rand()%possibleFirst.size();
        int selection = int(double(gsl_rng_uniform(rgen))*int(possibleFirst->size()));
        first->push_back((*possibleFirst)[selection]);
        // record sequence clusters
        if ((*possibleFirst)[selection]->vt) {
            observedVtSc.push_back((*possibleFirst)[selection]->sc);
        } else {
            observedNvtSc.push_back((*possibleFirst)[selection]->sc);
        }
        // calculate gene frequencies
        for (unsigned int i = 0; i < (*possibleFirst)[selection]->genotype.size();i++) {
            (*accessoryLoci)[i]->simFreq[0]+=(double((*possibleFirst)[selection]->genotype[i])/double(psize));
        }
    }
    
    // record sequence cluster statistics
    for (unsigned int i = 0; i < scList->size(); ++i) {
        startingVtScFrequencies[i] = std::count(observedVtSc.begin(),observedVtSc.end(),(*scList)[i]);
        startingNvtScFrequencies[i] = std::count(observedNvtSc.begin(),observedNvtSc.end(),(*scList)[i]);
    }
    
    // record gene frequency statistics
    std::vector<double> startingCogFrequencies(accessoryLoci->size(),0.0);
    for (unsigned int i = 0; i < accessoryLoci->size(); i++) {
        startingCogFrequencies[i] = (*accessoryLoci)[i]->simFreq[0];
    }
    std::transform(eqFreq.begin(), eqFreq.end(), startingCogFrequencies.begin(), cogDeviations.begin(), std::minus<double>());
    std::transform(cogWeights.begin(), cogWeights.end(), cogDeviations.begin(), cogDeviations.begin(), std::multiplies<double>());
    
    // tidy
    delete possibleFirst;
    
    return 0;
}

////////////////////////////////////////////////////
// Select first sample for input file replication //
////////////////////////////////////////////////////

int firstSample(std::vector<isolate*> *currentIsolates,int firstSample,std::ofstream& sampleOutFile,int minGen) {
    
    // data structures for sample
    std::vector<isolate*> *isolateSample = new std::vector<isolate*>;
    
    // get appropriately sized random sample from simulation
    while (isolateSample->size() < unsigned(firstSample)) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample->push_back(selectedIsolate);
        // print record of sample to file
        int vtInt = 0;
        if (selectedIsolate->latent_vt) {
            vtInt = 2;
        } else if (selectedIsolate->vt) {
            vtInt = 1;
        }
        sampleOutFile << selectedIsolate->id << "\t" << minGen << "\t" << selectedIsolate->serotype << "\t" << vtInt << "\t" << selectedIsolate->sc << std::endl;
    }
    
    // tidy
    delete isolateSample;
    
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
    for (unsigned int i = 0; i < seros->size(); i++) {
        int tmpSeroNum = std::count(genSerotypes.begin(),genSerotypes.end(),(*seros)[i]);
        (*serotypeF)[i] = tmpSeroNum;
    }
    
    // summarise sequence clusters in sample
    for (unsigned int i = 0; i < scs->size(); i++) {
        sampledVtScFreq[0][i]+=(double(std::count(sampledVtSequenceClusters.begin(),sampledVtSequenceClusters.end(),(*scs)[i]))/double(sampleSize));
        sampledNvtScFreq[0][i]+=(double(std::count(sampledNvtSequenceClusters.begin(),sampledNvtSequenceClusters.end(),(*scs)[i]))/double(sampleSize));
    }
    
    return 0;
}

/////////////////////////////////////////////////
// alter vaccine formulation during simulation //
/////////////////////////////////////////////////

int alterVaccineFormulation(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<std::vector<std::vector<isolate*> > > *popBySc) {

    std::vector<isolate*>::iterator iter;
    
    // change VT of current population
    for (iter = currentIsolates->begin(), currentIsolates->end(); iter != currentIsolates->end(); ++iter) {
        if (!(*iter)->vt) {
            (*iter)->vt = (*iter)->latent_vt;
        }
    }
    
    // change VT of underlying population
    for (iter = pop->begin(), pop->end(); iter != pop->end(); ++iter) {
        if (!(*iter)->vt) {
            (*iter)->vt = (*iter)->latent_vt;
        }
    }
    
    // change VT of potential immigrant isolates
    for (unsigned int i = 0; i < popBySc->size(); i++) {
        for (unsigned int j = 0; j < (*popBySc)[i].size(); j++) {
            for (iter = (*popBySc)[i][j].begin(), (*popBySc)[i][j].end(); iter != (*popBySc)[i][j].end(); ++iter) {
                if (!(*iter)->vt) {
                    (*iter)->vt = (*iter)->latent_vt;
                }
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

std::vector<int> getValidStrains(std::vector<std::vector<isolate*> > migrantInput) {
    
    std::vector<int> validStrains;
    
    for (unsigned int i = 0; i < migrantInput.size(); i++) {
        if (migrantInput[i].size() >= 1) {
            validStrains.push_back(i);
        }
    }
    
    return validStrains;
    
}

int reproduction(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,std::vector<std::vector<std::vector<isolate*> > > *migrantPool, std::vector<double> *cogWeights, std::vector<double> *cogDeviations,struct parms *sp, std::vector<double> * ef, std::vector<int> * vtScFreq,std::vector<int> * nvtScFreq,std::vector<double> * piGen,std::vector<int> *scList, int gen,std::vector<double> * timeGen,std::vector<double> * fitGen,std::vector<std::string> * isolateGen,std::vector<int> * countGen, double popLimitFactor, int minGen, int secondVaccinationGeneration, float partialVaccine) {
    
    // calculate population limit for memory management
    int popLimit = round(popLimitFactor * double(sp->popSize));
    if (popLimit < 0) {
        popLimit = std::numeric_limits<int>::max();
    }
    
    // calculate vaccine strength
    double firstVaccineSelection = double(sp->vSelection);
    double secondVaccineSelection = double(sp->vSelection);
    if (sp->vaccineLag >= 0 && gen >= 0) {
        // time since first vaccine introduction
        if (gen <= sp->vaccineLag) {
            firstVaccineSelection = double(gen)/double(sp->vaccineLag)*double(sp->vSelection);
        }
        // date of last vaccine introduction
        if (gen >= secondVaccinationGeneration) {
            // is simulation in lag period?
            int timeSinceVaccination = gen - secondVaccinationGeneration;
            if (timeSinceVaccination < sp->vaccineLag) {
                secondVaccineSelection = double(timeSinceVaccination)/double(sp->vaccineLag)*double(sp->vSelection);
            }
        }
    }
    
    // new COG deviations array
    std::vector<int> futureCogCount(ef->size());
    std::fill (futureCogCount.begin(),futureCogCount.end(),0);
    
    // record population statistics
    std::vector<int> futureVtScs;
    std::vector<int> futureNvtScs;
    std::vector<std::string> futureSerotypes;
    int genotypeCount = 0;
    
    // basic reproduction number based on immigration and population size
    double baseR = (1-sp->immigrationRate);

    // sort current population
    std::vector<isolate*>::iterator iter;
    std::sort(currentIsolates->begin(),currentIsolates->end());
    
    // allow old generation to reproduce
    std::string oldId = "";
    double oldFitness = 0.0;
    
    // Store non-standardised fitnesses
    std::vector<double> unstandardised_fitnesses;
    
    // Calculate non-standardised fitnesses
    for (iter = currentIsolates->begin(), currentIsolates->end(); iter != currentIsolates->end(); ++iter) {
        
        // calculate fitness of each new genotype in population
        if ((*iter)->id != oldId) {
            
            // record isolate fitnesses for extended output
            if (sp->programme == "x" && oldFitness != 0.0) {
                isolateGen->push_back((*iter)->id);
                fitGen->push_back(oldFitness);
                timeGen->push_back(gen);
                countGen->push_back(genotypeCount);
                genotypeCount = 0;
            }
            
            std::vector<double> fitnesses(cogDeviations->size());
            
            std::transform(cogDeviations->begin(), cogDeviations->end(), (*iter)->genotype.begin(), fitnesses.begin(), std::multiplies<double>());
            double freqDepFitSum = std::accumulate(fitnesses.begin(),fitnesses.end(),0.0);
            double freqDepFit = pow((1+sp->fSelection),freqDepFitSum);
            double vaccineFit = 1.0;
            
            // only switch on vaccine selection pressure after vaccine is introduced
            if (gen >= 0 && (*iter)->vt) {
                // if vaccine lag, then only apply to latent_vt for the second vaccine
                if (sp->vaccineLag > 0) {
                    if (gen >= secondVaccinationGeneration) {
                        // Here it is assumed anything cross-reactive in first vaccine
                        // is included in the second vaccine - needs to be fixed
                        if ((*iter)->latent_vt || (*iter)->partial_vt) {
                            // reduced selection on serotypes unique to second vaccine
                            if ((*iter)->latent_vt) {
                                if ((*iter)->partial_vt) {
                                    vaccineFit = 1.0 - (secondVaccineSelection*partialVaccine);
                                } else {
                                    vaccineFit = 1.0 - secondVaccineSelection;
                                }
                            } else {
                                vaccineFit = 1.0 - secondVaccineSelection;
                            }
                        } else {
                            // full selection on serotypes in first vaccine
                            // adjusted for time since introduction of first vaccine
                            if ((*iter)->partial_vt) {
                                vaccineFit = 1.0 - firstVaccineSelection*partialVaccine;
                            } else {
                                vaccineFit = 1.0 - firstVaccineSelection;
                            }
                        }
                    } else {
                        // full selection on serotypes in first vaccine
                        // adjusted for time since introduction of first vaccine
                        if ((*iter)->partial_vt) {
                            vaccineFit = 1.0 - firstVaccineSelection*partialVaccine;
                        } else {
                            vaccineFit = 1.0 - firstVaccineSelection;
                        }
                    }
                } else {
                    if ((*iter)->partial_vt) {
                        vaccineFit = 1.0 - sp->vSelection*partialVaccine;
                    } else {
                        vaccineFit = 1.0 - sp->vSelection;
                    }
                }
            }
            
            // calculate overall fitness
            double overallFitness = baseR*vaccineFit*freqDepFit;
            oldFitness = overallFitness;
            oldId = (*iter)->id;
        }
        
        (*iter)->fitness = oldFitness;
        unstandardised_fitnesses.push_back(oldFitness);
        
    }
    
    // Standardise fitnesses for density dependent-regulation
    double standardisation_factor = 1.0;
    if (sp->densdepMode == 0) {
        standardisation_factor = double(sp->popSize)/double(currentIsolates->size());
    } else if (sp->densdepMode == 1) {
        double densdep = double(sp->popSize)/double(currentIsolates->size());
        double mean_unstandardised_fitness = std::accumulate(unstandardised_fitnesses.begin(),
                                                             unstandardised_fitnesses.end(),
                                                             0.0)/unstandardised_fitnesses.size();
        standardisation_factor = densdep/mean_unstandardised_fitness;
    } else {
        std::cerr << "Unrecognsised densdep mode: " << sp->densdepMode << std::endl;
    }
    
    // select offspring by Poisson distribution using standardised fitness
    for (iter = currentIsolates->begin(), currentIsolates->end(); iter != currentIsolates->end(); ++iter) {

        int progeny = gsl_ran_poisson(rgen,standardisation_factor * (*iter)->fitness);
        for (int p = 0; p < progeny; p++) {
            futureIsolates->push_back((*iter));
            // record statistics
            if ((*iter)->vt) {
                futureVtScs.push_back((*iter)->sc);
            } else {
                futureNvtScs.push_back((*iter)->sc);
            }
            futureSerotypes.push_back((*iter)->serotype);
        }
        genotypeCount++;
        
        // halt simulations if next population is above limit
        // add check here for rapidly-expanding populations
        if (futureIsolates->size() > popLimit) {
            // end process
            return 8888;
        }
        
        // tally their COGs in new matrix
        if (progeny > 0) {
            for (unsigned int c = 0; c < (*iter)->genotype.size(); c++) {
                futureCogCount[c]+=(*iter)->genotype[c]*progeny;
            }
        }
    }
    
    // allow for immigration
    double adjustedImmigrationRate = (sp->immigrationRate)*(double(sp->popSize)/double(currentIsolates->size()));
    int immigrants = gsl_ran_binomial(rgen,adjustedImmigrationRate,sp->popSize);
    
    // if mode == 2 immigration, chose the timestep
    // from which to select the migrants
    int migration_gen = -1;
    if (sp->immigrationType == 2) {
        for (int g = 0; g <= gen-minGen; g++) {
            std::vector<isolate*> *tmp_candidates = &(*migrantPool)[0][g];
            if (tmp_candidates->size() >= 1) {
                migration_gen = g;
            }
        }
    } else if (sp->immigrationType == 3) {
        for (int g = 0; g <= gen-minGen; g++) {
            std::vector<std::vector<isolate*> > *tmp_candidates = &(*migrantPool)[g];
            if (tmp_candidates->size() >= 1) {
                migration_gen = g;
            }
        }
    }

    std::vector<isolate*> *candidates = new std::vector<isolate*>;
    int selection = 0;
    
    for (int f = 0; f < immigrants; f++) {
        
        if (sp->immigrationType == 0) {
            candidates = &(*migrantPool)[0][0];
            selection = int(double(gsl_rng_uniform(rgen))*candidates->size());
        } else if (sp->immigrationType == 1) {
            std::vector<int> migrantStrains = getValidStrains((*migrantPool)[0]);
            int selectedScIndex = int(double(gsl_rng_uniform(rgen))*migrantStrains.size());
            candidates = &(*migrantPool)[0][migrantStrains[selectedScIndex]];
            selection = int(double(gsl_rng_uniform(rgen))*candidates->size());
        } else if (sp->immigrationType == 2 || sp->immigrationType == 3) {
            if (migration_gen == -1) {
                std::cerr << "No migrant genotypes available at generation " << gen << std::endl;
                return 1;
            } else {
                if (sp->immigrationType == 2) {
                    candidates = &(*migrantPool)[0][migration_gen];
                    selection = int(double(gsl_rng_uniform(rgen))*candidates->size());
                } else if (sp->immigrationType == 3) {
                    std::vector<int> migrantStrains = getValidStrains((*migrantPool)[migration_gen]);
                    int selectedScIndex = int(double(gsl_rng_uniform(rgen))*migrantStrains.size());
                    candidates = &(*migrantPool)[migration_gen][migrantStrains[selectedScIndex]];
                    selection = int(double(gsl_rng_uniform(rgen))*candidates->size());
                }
            }

        } else {
            return 1;
        }

        // select isolates
        isolate* selectedIsolate = (*candidates)[selection];
        
        futureIsolates->push_back(selectedIsolate);
        // record statistics
        if (selectedIsolate->vt) {
            futureVtScs.push_back(selectedIsolate->sc);
        } else {
            futureNvtScs.push_back(selectedIsolate->sc);
        }
        futureSerotypes.push_back(selectedIsolate->serotype);
        // tally COGs in new matrix
        for (unsigned int c = 0; c < selectedIsolate->genotype.size(); c++) {
            futureCogCount[c]+=selectedIsolate->genotype[c];
        }
    }
    
    // halt simulations if next population is above limit
    if (futureIsolates->size() > popLimit) {
        // end process
        return 8888;
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
    if (sp->programme == "x") {
        std::transform(cogDeviations->begin(), cogDeviations->end(), piGen->begin(), piGen->begin(), std::plus<double>());
    }
    // now weight the COG deviations
    std::transform(cogWeights->begin(), cogWeights->end(), cogDeviations->begin(), cogDeviations->begin(), std::multiplies<double>());
    
    // now summarise sequence clusters - vaccine types
    for (unsigned int i = 0; i < scList->size(); i++) {
        (*vtScFreq)[i] = std::count(futureVtScs.begin(),futureVtScs.end(),(*scList)[i]);
    }
    
    // now summarise sequence clusters - non-vaccine types
    for (unsigned int i = 0; i < scList->size(); i++) {
        (*nvtScFreq)[i] = std::count(futureNvtScs.begin(),futureNvtScs.end(),(*scList)[i]);
    }
    
    return 0;
}

///////////////////
// Recombination //
///////////////////

bool alleleExchange (bool r, bool d, double a) {
    
    // present in recipient
    if (r) {
        // also present in donor
        if (d) {
            return true;
        // absent in donor
        } else {
            double p =  double(gsl_rng_uniform(rgen));
            // if asymmetry favours deletion
            if (a <= 1) {
                return false;
            // if asymmetry favours insertion
            } else if (p <= double(1.0/a)) {
                return false;
            } else {
                return true;
            }
        }
    // absent in recipient
    } else {
        // present in donor
        if (d) {
            double p =  double(gsl_rng_uniform(rgen));
            // asymmetry favours insertions
            if (a >= 1) {
                return true;
            // asymmetry favours deletions
            } else if (p <= a) {
                return true;
            } else {
                return false;
            }
        // also absent in donor
        } else {
            return false;
        }
    }

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


std::string generateRandomString(size_t length) {
    const char* charmap = "ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567";
    const size_t charmapLength = strlen(charmap);
    auto generator = [&](){ return charmap[rand()%charmapLength]; };
    std::string result;
    result.reserve(length);
    generate_n(back_inserter(result), length, generator);
    return result;
}

int recombination(std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates,char* markerFilename,double transformationProportion,double transformationRate,double transformationAsymmetryLoci, double transformationAsymmetryMarker) {

    // recombination
    std::vector<isolate*>::iterator iit;
    
    std::vector<isolate*> *isolate_array = currentIsolates;
    
    for (iit = isolate_array->begin(), isolate_array->end() ; iit != isolate_array->end(); ++iit) {
        bool recHappened = false;
        // test if recombination occurs in this isolate at this timestep
        double rTrans = gsl_rng_uniform(rgen);
        if (rTrans <= transformationRate) {
            // find recipient if recombination happened
            isolate recipient = *(*iit);
            // select the donor isolate
            int j = int(gsl_rng_uniform(rgen) * isolate_array->size());
            isolate donor = *(*isolate_array)[j];
            // if so, how much of the genome is replaced
            for (unsigned int i = 0; i < recipient.genotype.size(); i++) {
                // recombination at selected loci
                double pTrans = gsl_rng_uniform(rgen);
                if (pTrans <= transformationProportion) {
                    bool donorAllele = donor.genotype[i];
                    bool recipientAllele = recipient.genotype[i];
                    recipient.genotype[i] = alleleExchange(recipientAllele,donorAllele,transformationAsymmetryLoci);
                    if (recipient.genotype[i] != recipientAllele) {
                        recHappened = true;
                    }
                    
                }
            }
            // recombination at unselected markers
            if (markerFilename != NULL) {
                for (unsigned int i = 0; i < recipient.markers.size(); i++) {
                    double mTrans = gsl_rng_uniform(rgen);
                    if (mTrans <= transformationProportion) {
                        bool donorAllele = donor.markers[i];
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
                
                // assign ID with new hexadecimal suffix
                std::string new_id = seglist[0]+"_"+generateRandomString(7);
                
                // store in new population
                
                isolate* tmp = new isolate(new_id,
                                           recipient.year,
                                           recipient.sc,
                                           recipient.serotype,
                                           recipient.vt,
                                           recipient.latent_vt,
                                           recipient.partial_vt,
                                           &recipient.genotype,
                                           &recipient.markers,
                                           recipient.fitness
                                           );
                futureIsolates->push_back(tmp);
                
            } else {
                futureIsolates->push_back((*iit));
            }
        } else {
            futureIsolates->push_back((*iit));
        }

    }

    return 0;
}

///////////////////////////////////////////////////////
// Update equilibrium frequencies post-recombination //
///////////////////////////////////////////////////////

int update_locus_freq(std::vector<isolate*> *futureIsolates, std::vector<double> *cogWeights, std::vector<double> *cogDeviations, std::vector<double> *ef) {

//    std::cout << "Recombination count per generation: " << recombination_count << " num. loci transferred is: " << transferred_loci << " of attempted " << try_transfer << " representing - 11: " << one_one << " 10: " << one_zero << " 01: " << zero_one << " 00: " << zero_zero << std::endl; // debug
    
    std::vector<isolate*>::iterator iit;
    
    // calculate COG deviations in next generation
    std::vector<double> cogFractions(ef->size(),0.0);
    for (iit = futureIsolates->begin(), futureIsolates->end(); iit != futureIsolates->end(); ++iit) {
        for (unsigned int c = 0; c < (*iit)->genotype.size(); c++) {
            cogFractions[c]+=(double((*iit)->genotype[c])/double(futureIsolates->size()));
        }
    }
    // alter COG deviation vector according to recombinant phenotypes
    std::fill(cogDeviations->begin(),cogDeviations->end(),0.0);
    std::transform(ef->begin(), ef->end(), cogFractions.begin(), cogDeviations->begin(), std::minus<double>());
    // now weight the COG deviations
    std::transform(cogWeights->begin(), cogWeights->end(), cogDeviations->begin(), cogDeviations->begin(), std::multiplies<double>());
    
    return 0;
    
}

////////////////////////////////////////
// Move isolates into next generation //
////////////////////////////////////////

int updatePopulation(std::vector<isolate*> *pop,std::vector<isolate*> *new_pop) {
    
    // swap populations
    pop->clear();
    pop->shrink_to_fit();
    pop->swap(*new_pop);
    
    // finish
    return 0;
    
}

int nextGeneration(std::vector<isolate*> *pop,std::vector<isolate*> *new_pop,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *futureIsolates, std::vector<std::vector<std::vector<isolate*> > > *migrantPool) {
    
    // list of objects to delete
    std::vector<isolate*> *to_delete = new std::vector<isolate*>;
    
    // iterate through currentIsolates to identify obsolete genotypes
    sort(currentIsolates->begin(), currentIsolates->end());
    std::vector<isolate*>::iterator iit;
    iit = unique(currentIsolates->begin(), currentIsolates->end());
    currentIsolates->resize(distance(currentIsolates->begin(),iit));
    
    for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
        if (!(std::find(futureIsolates->begin(), futureIsolates->end(),(*iit))!=futureIsolates->end())) {
            if (!(std::find(new_pop->begin(), new_pop->end(), (*iit))!=new_pop->end())) {
                bool found = false;
                for (int a = 0; a < migrantPool->size(); a++) {
                    for (int b = 0; b < migrantPool[a].size(); b++) {
                        for (int c = 0; c < migrantPool[a][b].size(); c++) {
                            if (std::find(migrantPool[a][b][c].begin(), migrantPool[a][b][c].end(), (*iit))!=migrantPool[a][b][c].end()) {
                                found = true;
                            }
                        }
                    }
                }
                if (!(found)) {
                    to_delete->push_back(*iit);
                }
            }
        }
    }

    // iterate through population to identify obsolete genotypes
    sort(pop->begin(), pop->end());
    iit = unique(pop->begin(), pop->end());
    pop->resize(distance(pop->begin(),iit));
    
    for (iit = pop->begin(), pop->end() ; iit != pop->end(); ++iit) {
        if (!(std::find(futureIsolates->begin(), futureIsolates->end(),(*iit))!=futureIsolates->end())) {
            if (!(std::find(new_pop->begin(), new_pop->end(), (*iit))!=new_pop->end())) {
                if (!(std::find(to_delete->begin(), to_delete->end(), (*iit))!=to_delete->end())) {
                    to_delete->push_back(*iit);
                }
            }
        }
    }

    // delete obsolete isolates
    sort(to_delete->begin(), to_delete->end());
    iit = unique(to_delete->begin(), to_delete->end());
    pop->resize(distance(to_delete->begin(),iit));

    int i = 0;
    for (iit = to_delete->begin(), to_delete->end() ; iit != to_delete->end(); ++iit) {
        delete (*iit);
        ++i;
    }
    
    // replace current isolates with future isolates
    currentIsolates->clear();
    currentIsolates->shrink_to_fit();
    currentIsolates->swap(*futureIsolates);
    futureIsolates->clear();
    futureIsolates->shrink_to_fit();
    
    // replace old population with new population
    // if old population not in extant population
    // or new migrant pool, then dropped from simulation
    pop->clear();
    pop->shrink_to_fit();
    pop->swap(*new_pop);
    new_pop->clear();
    new_pop->shrink_to_fit();
    
    // tidy up
    delete to_delete;
    
    return 0;
    
}

////////////////////////////////////////////
// Compare simulation and genomic samples //
////////////////////////////////////////////

int compareSamples(int gen,int minGen,int sampleSize,std::vector<isolate*> *currentIsolates,std::vector<isolate*> *pop,std::vector<cog*> *accessoryLoci,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<int> &sampledSeroFreq,std::vector<std::string> &serotypeList,std::vector<double> &vtCogFittingStatsList,std::vector<double> &nvtCogFittingStatsList,std::vector<double> &strainFittingStatsList,std::ofstream& sampleOutFile,struct parms *sp) {
    
    // data structures for sample
    std::vector<isolate*> isolateSample;
    std::vector<std::string> currentSerotypeObservations;
    std::vector<int> currentVtScObservations;
    std::vector<int> currentNvtScObservations;
    
    // get appropriately sized random sample from simulation
    while (isolateSample.size() < unsigned(sampleSize)) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample.push_back(selectedIsolate);
        // record serotype frequencies
        currentSerotypeObservations.push_back(selectedIsolate->serotype);
        // record sequence cluster frequencies
        if (selectedIsolate->vt) {
            currentVtScObservations.push_back(selectedIsolate->sc);
        } else {
            currentNvtScObservations.push_back(selectedIsolate->sc);
        }
        // print record of sample to file
        int vtInt = 0;
        if (selectedIsolate->latent_vt) {
            vtInt = 2;
        } else if (selectedIsolate->vt) {
            vtInt = 1;
        }
        sampleOutFile << selectedIsolate->id << "\t" << gen << "\t" << selectedIsolate->serotype << "\t" << vtInt << "\t" << selectedIsolate->sc << std::endl;
        // calculate gene frequencies
        for (unsigned int i = 0; i < selectedIsolate->genotype.size();i++) {
            (*accessoryLoci)[i]->simFreq[gen-minGen]+=(double(selectedIsolate->genotype[i])/double(sampleSize));
        }
    }
    
    // summarise serotype information
    for (unsigned int i = 0; i < serotypeList.size(); i++) {
        sampledSeroFreq[i] = std::count(currentSerotypeObservations.begin(),currentSerotypeObservations.end(),serotypeList[i]);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::vector<std::string> > serotypes_in_sc(scList.size());
    
    // First identify all serotypes in all SC - in the actual and simulated populations
    int total_sample_size = 0;
    int simulated_sample_size = 0;
    for (int i = 0; i < scList.size(); ++i) {
        std::vector<isolate*>::iterator iiter;
        // Iterate through actual population
        for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
            if ((*iiter)->year == gen) {
                if ((*iiter)->sc == scList[i]) {
                    bool unseen_serotype = true;
                    for (int j = 0; j < serotypes_in_sc[i].size(); ++j) {
                        if ((*iiter)->serotype == serotypes_in_sc[i][j]) {
                            unseen_serotype = false;
                        }
                    }
                    if (unseen_serotype) {
                        serotypes_in_sc[i].push_back((*iiter)->serotype);
                    }
                    total_sample_size++;
                }
            }
        }
        // Iterate through simulated population
        for (iiter = isolateSample.begin(), isolateSample.end() ; iiter != isolateSample.end(); ++iiter) {
            if ((*iiter)->sc == scList[i]) {
                bool unseen_serotype = true;
                for (int j = 0; j < serotypes_in_sc[i].size(); ++j) {
                    if ((*iiter)->serotype == serotypes_in_sc[i][j]) {
                        unseen_serotype = false;
                    }
                }
                if (unseen_serotype) {
                    serotypes_in_sc[i].push_back((*iiter)->serotype);
                }
                simulated_sample_size++;
            }
        }
    }
    
    // Check samples sizes match
    if (simulated_sample_size != total_sample_size) {
        std::cerr << "Observed sample size (" << total_sample_size << ") does not match that of simulation (" << simulated_sample_size << ")" << std::endl;
        exit(1);
    }
    
    // Record frequencies and calculate deviations
    double jsd = 0.0;
    for (int i = 0; i < scList.size(); ++i) {
        for (int j = 0; j < serotypes_in_sc[i].size(); ++j) {
            int actual_count = 0;
            int simulated_count = 0;
            std::vector<isolate*>::iterator iiter;
            // Iterate through actual population
            for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
                if ((*iiter)->year == gen && (*iiter)->sc == scList[i] && (*iiter)->serotype == serotypes_in_sc[i][j]) {
                    actual_count++;
                }
            }
            // Iterate through simulated population
            for (iiter = isolateSample.begin(), isolateSample.end() ; iiter != isolateSample.end(); ++iiter) {
                if ((*iiter)->sc == scList[i] && (*iiter)->serotype == serotypes_in_sc[i][j]) {
                    simulated_count++;
                }
            }
            // Calculate frequencies
            double actual_frequency = actual_count/(1.0*total_sample_size);
            double simulated_frequency = simulated_count/(1.0*simulated_sample_size);
            double M = 0.5*(actual_frequency+simulated_frequency);
            // Calculate JSD
            double individual_jsd = 0.0;
            if (actual_frequency > 0.0) {
                individual_jsd += 0.5*actual_frequency*log(actual_frequency/M);
            }
            if (simulated_frequency > 0.0) {
                individual_jsd += 0.5*simulated_frequency*log(simulated_frequency/M);
            }
            jsd += individual_jsd;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
//    // summarise sequence cluster information from genomic data
//    std::vector<int> currentGenomicVtScObservations;
//    std::vector<int> currentGenomicNvtScObservations;
//    
//    std::vector<isolate*>::iterator iiter;
//    for (iiter = pop->begin(), pop->end() ; iiter != pop->end(); ++iiter) {
//        if ((*iiter)->year == gen && (*iiter)->vt) {
//            currentGenomicVtScObservations.push_back((*iiter)->sc);
//        } else if ((*iiter)->year == gen && !((*iiter)->vt)) {
//            currentGenomicNvtScObservations.push_back((*iiter)->sc);
//        }
//    }
//    
//    // summarise sequence cluster observations into counts
//    std::vector<int> currentVtScCounts(scList.size(),0);
//    std::vector<int> currentNvtScCounts(scList.size(),0);
//    std::vector<int> currentGenomicVtScCount(scList.size(),0);
//    std::vector<int> currentGenomicNvtScCount(scList.size(),0);
//    
//    for (unsigned int i = 0; i < scList.size(); ++i) {
//        // observations from simulated data
//        currentVtScCounts[i] = std::count(currentVtScObservations.begin(),currentVtScObservations.end(),scList[i]);
//        currentNvtScCounts[i] = std::count(currentNvtScObservations.begin(),currentNvtScObservations.end(),scList[i]);
//        // general storage of frequencies for printing
//        sampledVtScFreq[gen-minGen][i] = double(currentVtScCounts[i])/double(sampleSize);
//        sampledNvtScFreq[gen-minGen][i] = double(currentNvtScCounts[i])/double(sampleSize);
//        // observations from genomic data
//        currentGenomicVtScCount[i] = std::count(currentGenomicVtScObservations.begin(),currentGenomicVtScObservations.end(),scList[i]);
//        currentGenomicNvtScCount[i] = std::count(currentGenomicNvtScObservations.begin(),currentGenomicNvtScObservations.end(),scList[i]);
//    }
//    
//    // calculate divergence between real and simulated samples
//    double cumulativeScDifference = 0.0;
////    double cumulativeJSD = 0.0;
//    for (unsigned int i = 0; i < scList.size(); i++) {
////        double dkl_vt = std::abs(((double(currentGenomicVtScCount[i])+0.5)/double(sampleSize))*log((double(currentGenomicVtScCount[i])+0.5)/(double(currentVtScCounts[i])+0.5)));
////        double dkl_nvt = std::abs(((double(currentGenomicNvtScCount[i])+0.5)/double(sampleSize))*log((double(currentGenomicNvtScCount[i])+0.5)/(double(currentNvtScCounts[i])+0.5)));
//        
//        // calculate Jensen-Shannon divergence
//        double jsd_vt = 0;
//        double jsd_nvt = 0;
//        
//        double currentVtSimFreq = double(currentVtScCounts[i])/double(sampleSize);
//        double currentNvtSimFreq = double(currentNvtScCounts[i])/double(sampleSize);
//        double currentVtGenomicFreq = double(currentGenomicVtScCount[i])/double(sampleSize);
//        double currentNvtGenomicFreq = double(currentGenomicNvtScCount[i])/double(sampleSize);
//        double m_vt = (0.5/double(sampleSize))*(double(currentGenomicVtScCount[i])+double(currentVtScCounts[i]));
//        double m_nvt = (0.5/double(sampleSize))*(double(currentNvtScCounts[i])+double(currentGenomicNvtScCount[i]));
//        if (m_vt > 0.0) {
//            if (currentVtSimFreq > 0) {
//                jsd_vt += 0.5*currentVtSimFreq*log(currentVtSimFreq/m_vt);
//            }
//            if (currentVtGenomicFreq > 0) {
//                jsd_vt += 0.5*currentVtGenomicFreq*log(currentVtGenomicFreq/m_vt);
//            }
//        }
//        if (m_nvt > 0.0) {
//            if (currentNvtSimFreq > 0) {
//                jsd_nvt += 0.5*currentNvtSimFreq*log(currentNvtSimFreq/m_nvt);
//            }
//            if (currentNvtGenomicFreq > 0) {
//                jsd_nvt += 0.5*currentNvtGenomicFreq*log(currentNvtGenomicFreq/m_nvt);
//            }
//        }
//        
//        cumulativeScDifference+=(jsd_vt+jsd_nvt);
////        std::cerr << scList[i] << "\t" << dkl_vt << "\t" << dkl_nvt << "\t" << jsd_vt << "\t" << jsd_nvt << std::endl;
////        std::cerr << scList[i] << "\t" << dkl_vt << "\t" << dkl_nvt << "\t" << jsd_vt << "\t" << jsd_nvt << "\t" << currentVtSimFreq << "\t" << currentNvtSimFreq << "\t" << currentVtGenomicFreq << "\t" << currentNvtGenomicFreq << std::endl;
////        std::cerr << "sc " << scList[i] << " vt KLD " << dkl_vt << " nvt KLD " << dkl_nvt << " total KLD " << cumulativeScDifference << std::endl;
////        std::cerr << "\tsc " << scList[i] << " vt JSD " << jsd_vt << " nvt JSD " << jsd_nvt << " total JSD " << cumulativeJSD << std::endl;
//    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
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
    for (unsigned int i = 0; i < accessoryLoci->size(); i++) {
        // separate by VT for fitting statistic calculation
        if ((*accessoryLoci)[i]->vt) {
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
    for (unsigned int i = 0; i < currentSampleVtFreq.size(); i++) {
        if (startingSampleVtFreq[i] > 0 && startingGenomicVtFreq[i] > 0) {
            vtCogStat+=std::abs((currentSampleVtFreq[i]/startingSampleVtFreq[i])-(currentGenomicVtFreq[i]/startingGenomicVtFreq[i]));
        }
    }
    
    // calculate correlations - absolute, not squared, value
    double simulationNvtCorrelation = 0.0;
    double genomicNvtCorrelation = 0.0;
    if (sp->programme == "x") {
        double simulationNvtCorrelation = pearson(&currentSampleNvtFreq,&startingSampleNvtFreq);
        double genomicNvtCorrelation = pearson(&currentGenomicNvtFreq,&startingGenomicNvtFreq);
    }
//    double nvtCogStat = std::abs((simulationNvtCorrelation-genomicNvtCorrelation)/(1-genomicNvtCorrelation));
    double nvtCogStat = std::abs(simulationNvtCorrelation-genomicNvtCorrelation);
    
    // record final statistics
    vtCogFittingStatsList.push_back(vtCogStat);
    nvtCogFittingStatsList.push_back(nvtCogStat);
    strainFittingStatsList.push_back(jsd);
    
    return 0;
}

//////////////////////////////////////////////
// Parse disease information //
//////////////////////////////////////////////

int parse_disease_data(char* epiFilename,
                       std::vector<int> *diseaseTime,
                       std::vector<std::string> *diseaseSeroList,
                       std::vector<int> *diseaseScList,
                       std::vector<int> *diseaseVt,
                       std::vector<double> *diseaseInvasiveness,
                       std::vector<int> *diseasePopulation,
                       std::vector<int> *diseaseCount) {
    
    std::ifstream disease_file;
    disease_file.open(epiFilename, std::ifstream::in);
    if (disease_file) {
        std::string line;
        // parse lines
        while (std::getline(disease_file, line)) {
            // temporary information stores
            int disease_time = -1;
            std::string disease_serotype;
            int disease_vt = -1;
            int disease_sc = -1;
            double disease_invasiveness = 0.0;
            int disease_population = -1;
            int disease_count = -1;
            // read line by line
            std::istringstream iss(line);
            int sIndex = 0;
            while (iss) {
                std::string temp;
                while (getline(iss, temp, '\t')) {
                    if (sIndex == 1) {
                        disease_time = atoi(temp.c_str());
                    } else if (sIndex == 2) {
                        disease_serotype = temp.c_str();
                    } else if (sIndex == 3) {
                        disease_vt = atoi(temp.c_str());
                    } else if (sIndex == 4) {
                        disease_sc = atoi(temp.c_str());
                    } else if (sIndex == 5) {
                        disease_invasiveness = atof(temp.c_str());
                    } else if (sIndex == 6) {
                        disease_population = atoi(temp.c_str());
                    } else if (sIndex == 7) {
                        disease_count = atoi(temp.c_str());
                    }
                    sIndex++;
                }
                if (disease_serotype != "Serotype") {
                    diseaseTime->push_back(disease_time);
                    diseaseSeroList->push_back(disease_serotype);
                    diseaseScList->push_back(disease_sc);
                    diseaseVt->push_back(disease_vt);
                    diseaseInvasiveness->push_back(disease_invasiveness);
                    diseasePopulation->push_back(disease_population);
                    diseaseCount->push_back(disease_count);
                }
            }
        }
    }
    
    return 0;
    
}

//////////////////////////////////////////////
// Compare to disease samples //
//////////////////////////////////////////////

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
                            std::ofstream& diseaseOutFile) {
    
    // first calculate total population size for frequency calculation
    int population_size = currentIsolates->size();
    
    // store the simulated carriage and disease counts
    int total_disease_count = 0;
    int actual_disease_count = 0;
    std::vector<int> disease_counts;
    std::vector<int> carriage_counts;
    
    // then iterate through strain/serotype combinations
    // calculate simulated disease frequencies
    double total_deviation = 0.0;
    for (unsigned int i = 0; i < diseaseTime->size(); i++) {
        int simulated_disease_count = 0;
        int simulated_carriage_count = 0;
        // only test at the appropriate time point
        if ((*diseaseTime)[i] == simulation_time) {
            // get serotype and strain
            std::string serotype = (*diseaseSeroList)[i];
            int sc = (*diseaseScList)[i];
            // calculate the carriage frequency of each
            std::vector<isolate*>::iterator iiter;
            for (iiter = currentIsolates->begin(), currentIsolates->end(); iiter != currentIsolates->end(); ++iiter) {
                if (sc == (*iiter)->sc && serotype == (*iiter)->serotype) {
                    simulated_carriage_count++;
                }
            }
            double carriage_frequency = ((1.0)*simulated_carriage_count)/((1.0)*population_size);
            // calculate the Poisson variable
            double poisson_var = carriage_frequency*(*diseasePopulation)[i]*(*diseaseInvasiveness)[i];
            // draw from the distribution
            simulated_disease_count = gsl_ran_poisson(rgen,poisson_var);
            total_disease_count += simulated_disease_count;
            // record the actual disease total
            actual_disease_count += (*diseaseCount)[i];
        }
        disease_counts.push_back(simulated_disease_count);
        carriage_counts.push_back(simulated_carriage_count);
    }
    
    // Iterate a second time using the total_disease value
    for (unsigned int i = 0; i < diseaseTime->size(); i++) {
        // only test at the appropriate time point
        if ((*diseaseTime)[i] == simulation_time) {
            // isolate characteristics
            std::string serotype = (*diseaseSeroList)[i];
            int sc = (*diseaseScList)[i];
            double simulated_carriage_count = carriage_counts[i];
            double carriage_frequency = ((1.0)*simulated_carriage_count)/((1.0)*population_size);
            // calculate frequencies
            double simulated_disease_frequency = 1.0*disease_counts[i]/total_disease_count;
            double actual_disease_frequency = 1.0*(*diseaseCount)[i]/actual_disease_count;
            // add deviation to the appropriate total
            double JSD = 0.0;
            double M = (simulated_disease_frequency + actual_disease_frequency)/2.0;
            if (M > 0) {
                double D_PM = 0.0;
                if (actual_disease_frequency > 0) {
                    D_PM = actual_disease_frequency*log(actual_disease_frequency/M);
                }
                double D_QM = 0.0;
                if (simulated_disease_frequency > 0) {
                    D_QM = simulated_disease_frequency*log(simulated_disease_frequency/M);
                }
                JSD = 0.5*(D_PM+D_QM);
            }
            total_deviation += JSD;
            // print output
            diseaseOutFile << simulation_time << "\t" << carriage_frequency << "\t" << (*diseasePopulation)[i] << "\t" << (*diseaseInvasiveness)[i] << "\t" << serotype << "\t" << (*diseaseVt)[i] << "\t" << sc << "\t" << (*diseaseCount)[i] << "\t" << disease_counts[i] << "\t" << JSD << std::endl;
        }
        
    }
    
    // Record
    diseaseDivergence.push_back(total_deviation);
    
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
    while (isolateSample.size() < unsigned (sampleSize)) {
//        int selection = rand()%currentIsolates->size();
        int selection = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        isolate *selectedIsolate = (*currentIsolates)[selection];
        isolateSample.push_back(selectedIsolate);
        // record serotype frequencies
        currentSerotypeObservations.push_back(selectedIsolate->serotype);
        // record sequence cluster frequencies
        if (selectedIsolate->vt) {
            currentVtScObservations.push_back(selectedIsolate->sc);
        } else {
            currentNvtScObservations.push_back(selectedIsolate->sc);
        }
        // calculate gene frequencies
        for (unsigned int i = 0; i < selectedIsolate->genotype.size();i++) {
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

int rFitMetricCalculation(int minGen,std::vector<int> *samplingList,std::vector<int> &scList,std::vector< std::vector<double> > &sampledVtScFreq,std::vector< std::vector<double> > &sampledNvtScFreq,std::vector<isolate*> *population,std::vector<double> &rFitVector) {
    
    // data structures
    std::vector<double> actualFoldChanges(samplingList->size(),0.0);
    std::vector<double> simulatedFoldChanges(samplingList->size(),0.0);
    
    // calculate actual statistics
    std::vector< std::vector<double> > realVtScFreq(samplingList->size()+1,std::vector<double>(scList.size(),0.0));
//    realVtScFreq[genIndex][scIndex]+=(1/double(samplingList[genIndex]));
    std::vector< std::vector<double> > realNvtScFreq(samplingList->size()+1,std::vector<double>(scList.size(),0.0));
//    realNvtScFreq[genIndex][scIndex]+=(1/double(samplingList[genIndex]));
    
    // record VT and NVT observations
    std::vector< std::vector<int> > genomicVtObservations(samplingList->size()+1);
    std::vector< std::vector<int> > genomicNvtObservations(samplingList->size()+1);
    
    std::vector<isolate*>::iterator iiter;
    for (iiter = population->begin(), population->end(); iiter != population->end(); ++iiter) {
        int genIndex = (*iiter)->year-minGen;
        if ((*samplingList)[genIndex] > 0) {
            if ((*iiter)->vt) {
                genomicVtObservations[genIndex].push_back((*iiter)->sc);
            } else {
                genomicNvtObservations[genIndex].push_back((*iiter)->sc);
            }
        } else {
            std::cerr << "Misalignment between sampling generations at year " << (*iiter)->year << " with min gen " << minGen << std::endl;
            return 1;
        }
    }
    
    // summarise SC information
    for (unsigned int genIndex = 0; genIndex < samplingList->size(); genIndex++) {
        if ((*samplingList)[genIndex] > 0) {
            for (unsigned int scIndex = 0; scIndex < scList.size(); scIndex++) {
                realVtScFreq[genIndex][scIndex] = double(std::count(genomicVtObservations[genIndex].begin(),genomicVtObservations[genIndex].end(),scList[scIndex]))/double((*samplingList)[genIndex]);
                realNvtScFreq[genIndex][scIndex] = double(std::count(genomicNvtObservations[genIndex].begin(),genomicNvtObservations[genIndex].end(),scList[scIndex]))/double((*samplingList)[genIndex]);
            }
        }
    }
    
    // count number of timepoints at which samples are taken
    // for calculating mean frequencies
    int numberOfSamples = 0;
    for (unsigned int i = 0; i < samplingList->size(); ++i) {
        if ((*samplingList)[i] > 0) {
            numberOfSamples+=1;
        }
    }
    
    // return comparison values
    for (unsigned int scIndex = 0; scIndex < scList.size(); scIndex++) {
        
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
        for (unsigned int genIndex = 0; genIndex <= samplingList->size(); ++genIndex) {
            if ((*samplingList)[genIndex] > 0) {
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
            for (unsigned int j = 1; j < realVtScObservations.size(); ++j) {
                cumulativeFrequency+=realVtScObservations[j];
                double invExponent = double(realVtScObservationTimepoints[j]) - double(realVtScObservationTimepoints[0]);
                realTmp.push_back(pow(realVtScObservations[j]/realVtScObservations[0],(1/invExponent)));
            }
            double tmp_t = 0.0;
            double tmp_c = 0.0;
            for (unsigned int j = 0; j < realTmp.size(); ++j) {
                tmp_t+=realTmp[j];
                tmp_c++;
            }
            double realRmetric = tmp_t/tmp_c;
            // calculate metric for simulated data
            std::vector<double> simTmp;
            for (unsigned int j = 1; j < simVtScObservations.size(); ++j) {
                double invExponent = double(simVtScObservationTimepoints[j]) - double(simVtScObservationTimepoints[0]);
                simTmp.push_back(pow(simVtScObservations[j]/simVtScObservations[0],(1/invExponent)));
            }
            tmp_t = 0.0;
            tmp_c = 0.0;
            for (unsigned int j = 0; j < simTmp.size(); ++j) {
                tmp_t+=simTmp[j];
                tmp_c++;
            }
            double simRmetric = tmp_t/tmp_c;
            
            // summarise difference between metric estimates weighted by SC frequency in actual data
            RmetricDeviation+=((cumulativeFrequency/double(numberOfSamples))*fabs(realRmetric-simRmetric));

        }
        // calculate metric comparison if >1 timepoint for NVT
        if (realNvtScObservations.size() > 1 && simNvtScObservations.size() > 1) {
            
            // calculate metric for real data
            std::vector<double> realTmp;
            double cumulativeFrequency = realNvtScObservations[0];
            for (unsigned int j = 1; j < realNvtScObservations.size(); ++j) {
                cumulativeFrequency+=realNvtScObservations[j];
                double invExponent = double(realNvtScObservationTimepoints[j]) - double(realNvtScObservationTimepoints[0]);
                realTmp.push_back(pow(realNvtScObservations[j]/realNvtScObservations[0],(1/invExponent)));
            }
            double tmp_t = 0.0;
            double tmp_c = 0.0;
            for (unsigned int j = 0; j < realTmp.size(); ++j) {
                tmp_t+=realTmp[j];
                tmp_c++;
            }
            double realRmetric = tmp_t/tmp_c;
            
            // calculate metric for simulated data
            std::vector<double> simTmp;
            for (unsigned int j = 1; j < simNvtScObservations.size(); ++j) {
                double invExponent = double(simNvtScObservationTimepoints[j]) - double(simNvtScObservationTimepoints[0]);
                simTmp.push_back(pow(simNvtScObservations[j]/simNvtScObservations[0],(1/invExponent)));
            }
            tmp_t = 0.0;
            tmp_c = 0.0;
            for (unsigned int j = 0; j < simTmp.size(); ++j) {
                tmp_t+=simTmp[j];
                tmp_c++;
            }
            double simRmetric = tmp_t/tmp_c;
            
            // summarise difference between metric estimates weighted by SC frequency in actual data
            RmetricDeviation+=((cumulativeFrequency/double(numberOfSamples))*fabs(realRmetric-simRmetric));
        }
        // record deviation
        rFitVector[scIndex] = RmetricDeviation;
        
    }

    return 0;
    
}

///////////////////////////
// write output to files //
///////////////////////////

int printOutput(char* outputFilename,std::vector<std::string> *seroList,std::vector<std::vector<int> > &sampledSeroFreq,std::vector<int> *scList,std::vector<std::vector<int> > &vtScFreq,std::vector<std::vector<int> > &nvtScFreq,int gen,int minGen,std::vector<cog*> *accessoryLoci,std::vector<int> *samplingList,std::vector<std::vector<double> > &piGen,struct parms *sp,std::vector<double> * timeGen,std::vector<double> * fitGen,std::vector<std::string> * isolateGen,std::vector<int> * countGen) {
    
    // debug
//    for (unsigned int y = 0; y < samplingList->size(); y++) {
//        std::cerr << y << "\t" << (*samplingList)[y] << std::endl;
//    }
    
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
            for (unsigned int i = 0; i != scList->size(); i++) {
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
        for (unsigned int g = 0; g < samplingList->size(); g++) {
            if ((*samplingList)[g] > 0) {
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
        for (unsigned int c = 0; c < accessoryLoci->size(); c++) {
            cog tmpCog = *(*accessoryLoci)[c];
            cogOutFile << tmpCog.id << "\t" << tmpCog.vt << "\t" << tmpCog.weight;
//            for (unsigned int g = 0; g < tmpCog.actualFreq.size(); g++) {
            for (unsigned int g = 0; g < samplingList->size(); g++) {
                if ((*samplingList)[g] > 0) {
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
            for (unsigned int g = 0; g < accessoryLoci->size(); g++) {
                piOutFile << "\t" << (*accessoryLoci)[g]->id;
            }
            piOutFile << std::endl;
            // write content
            for (int pseudoGen = minGen; pseudoGen < (gen+1+minGen); pseudoGen++) {
                piOutFile << pseudoGen;
                for (unsigned int c = 0; c < accessoryLoci->size(); c++) {
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
        
        // isolate fitness output
        std::string fitOutFilename = prefix + ".fitnesses.out";
        std::ofstream fitOutFile;
        fitOutFile.open(fitOutFilename,std::ios::out);
        
        // write pi output
        if (fitOutFile.is_open()) {
            // write header
            fitOutFile << "Gen\tIsolate\tCount\tFitness" << std::endl;
            // write content
            for (int index = 0; index < timeGen->size(); index++) {
                fitOutFile << (*timeGen)[index] << "\t" << (*isolateGen)[index] << "\t" << (*countGen)[index] << "\t" << (*fitGen)[index] << std::endl;
            }
            // close
            fitOutFile.close();
        } else {
            std::cerr << "Unable to write to file " << fitOutFilename << std::endl;
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
        std::vector<int> currentFreq(accessoryLoci->size(),0.0);
        std::vector<isolate*>::iterator iit;
        for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
            std::transform (currentFreq.begin(), currentFreq.end(), (*iit)->genotype.begin(), currentFreq.begin(), std::plus<int>());
        }
        
        for (unsigned int i = 0; i < accessoryLoci->size(); i++) {
            popOutFile << (*accessoryLoci)[i]->id << "\t" << double(currentFreq[i])/currentIsolates->size() << std::endl;
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
            std::vector<int> currentFreq(markerList->size(),0.0);
            std::vector<isolate*>::iterator iit;
            for (iit = currentIsolates->begin(), currentIsolates->end() ; iit != currentIsolates->end(); ++iit) {
                std::transform (currentFreq.begin(), currentFreq.end(), (*iit)->markers.begin(), currentFreq.begin(), std::plus<int>());
            }
            
            for (unsigned int i = 0; i < markerList->size(); i++) {
                marOutFile << (*markerList)[i] << "\t" << double(currentFreq[i])/currentIsolates->size() << std::endl;
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
    while (currentSample->size() < unsigned(sampleSize)) {
//        int selectedIsolateIndex = int(rand()%currentIsolates->size());
        int selectedIsolateIndex = int(double(gsl_rng_uniform(rgen))*currentIsolates->size());
        currentSample->push_back((*currentIsolates)[selectedIsolateIndex]);
    }
    
    // write accessory output
    if (popOutFile.is_open()) {
    
        // print headers
        popOutFile << "Taxon";
        std::vector<cog*>::iterator cit;
        for (cit = accessoryLoci->begin(), accessoryLoci->end() ; cit != accessoryLoci->end(); ++cit) {
            popOutFile << "\t" << (*cit)->id;
        }
        popOutFile << std::endl;

        // print data
        std::vector<isolate*>::iterator iit;
        for (iit = currentSample->begin(), currentSample->end() ; iit != currentSample->end(); ++iit) {
            std::string outLine = (*iit)->id;
            for (unsigned int i = 0; i < (*iit)->genotype.size(); i++) {
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
            for (unsigned int i = 0; i < markerList->size(); i++) {
                titleLine.append("\t");
                titleLine.append((*markerList)[i]);
            }
            marOutFile << titleLine << std::endl;

            // print data
            std::vector<isolate*>::iterator iit;
            for (iit = currentSample->begin(), currentSample->end() ; iit != currentSample->end(); ++iit) {
                std::string outLine = (*iit)->id;
                for (unsigned int i = 0; i < (*iit)->markers.size(); i++) {
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

    // free memory
    delete currentSample;
    
    return 0;
}

/////////////
// Tidy up //
/////////////

void tidyUpIsolates(std::vector<isolate*> *a_list, std::vector<isolate*> *b_list, std::vector<isolate*> *c_list, std::vector<isolate*> *d_list) {
    
    // merge lists
    std::vector<isolate*> *isolate_list = new std::vector<isolate*>;
    isolate_list->reserve( a_list->size() + b_list->size() + c_list->size() + d_list->size() ); // preallocate memory
    isolate_list->insert( isolate_list->end(), a_list->begin(), a_list->end());
    isolate_list->insert( isolate_list->end(), b_list->begin(), b_list->end());
    isolate_list->insert( isolate_list->end(), c_list->begin(), c_list->end());
    isolate_list->insert( isolate_list->end(), d_list->begin(), d_list->end());
    
    sort(isolate_list->begin(), isolate_list->end());
    std::vector<isolate*>::iterator iit;
    iit = unique(isolate_list->begin(), isolate_list->end());
    isolate_list->resize(distance(isolate_list->begin(),iit));
    
    for (iit = isolate_list->begin(), isolate_list->end() ; iit != isolate_list->end(); ++iit) {
        if ((*iit) != NULL) {
            delete (*iit);
        }
    }
    
    delete isolate_list;
}

void tidyUpLoci(std::vector<cog*> *cog_list) {
    
    sort(cog_list->begin(), cog_list->end());
    std::vector<cog*>::iterator cit;
    
    for (cit = cog_list->begin(), cog_list->end() ; cit != cog_list->end(); ++cit) {
        if ((*cit) != NULL) {
            delete (*cit);
        }
    }

    delete cog_list;
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
