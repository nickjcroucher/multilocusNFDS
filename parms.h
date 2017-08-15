//
//  parms.h
//  frequencyDependentSimulation
//
//  Created by Nicholas Croucher on 28/09/2015.
//  Copyright (c) 2015 Imperial College. All rights reserved.
//

#ifndef frequencyDependentSimulation_parms_h
#define frequencyDependentSimulation_parms_h

#include <iostream>
#include <string>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include <vector>
#include "functions.h"

// structure for storing parameters
struct parms {
    std::string programme;
    double fSelection;
    double vSelection;
    double immigrationRate;
    int immigrationType;
    double upperLimit;
    double lowerLimit;
    double selectedProp;
    double lowerSelection;
    double higherSelection;
    double transformationProportion;
    double transformationRate;
    double transformationAsymmetryLoci;
    double transformationAsymmetryMarker;
    int popSize;
    int numGen;
    int genotypeSampleSize;
};

// structure for isolate objects
struct isolate {
    std::string id;
    int year;
    int sc;
    std::string serotype;
    double vt;
    double latent_vt;
    std::vector<bool> genotype;
    std::vector<bool> markers;
    
    // constructor for struct 'isolate'
    isolate(std::string init_id, int init_year, int init_sc, std::string init_serotype, double init_vt, double second_vt, std::vector<bool> *init_genotype, std::vector<bool> *init_markers) {
        id = init_id;
        year = init_year;
        sc = init_sc;
        serotype = init_serotype;
        vt = init_vt;
        latent_vt = second_vt;
        genotype = (*init_genotype);
        markers = (*init_markers);
    };
    
    // copy constructor for struct 'isolate'
    isolate(const isolate& other): id(other.id),year(other.year),sc(other.sc),serotype(other.serotype),vt(other.vt),latent_vt(other.latent_vt),genotype(other.genotype),markers(other.markers) {}
    
    // destructor for struct 'isolate'
//    ~isolate ();

};

// structure for COG objects
struct cog {
    std::string id;
    bool vt;
    double weight;
    double currentFreq;
    double eqFreq;
    std::vector<double> actualFreq;
    std::vector<double> simFreq;
    
    // constructor for struct 'cog'
    cog(std::string init_id, int init_vt, double init_weight, double init_eqFreq, std::vector<double> *init_actualFreq) {
        id = init_id;
        vt = init_vt;
        weight = init_weight;
        currentFreq = 0.0;
        eqFreq = init_eqFreq;
        actualFreq = (*init_actualFreq);
        std::vector<double> tmpSimFreq(actualFreq.size(),0.0);
        simFreq = tmpSimFreq;
    };
    
    // copy constructor for struct 'isolate'
    cog(const cog& other): id(other.id),vt(other.vt),weight(other.weight),currentFreq(other.currentFreq),actualFreq(other.actualFreq),simFreq(other.simFreq) {}
};

#endif
