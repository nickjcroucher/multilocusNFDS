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
    bool vt;
    bool latent_vt;
    std::vector<bool> genotype;
    std::vector<bool> markers;
    
    // new constructor for struct 'isolate'
//    isolate(std::string init_id, int init_year, int init_sc, std::string init_serotype, bool init_vt, bool second_vt, std::vector<bool> *init_genotype, std::vector<bool> *init_markers) : id(init_id),year(init_year),sc(init_sc),serotype(init_serotype),vt(init_vt),latent_vt(second_vt),genotype(*init_genotype),markers(*init_markers) {}
    
    // constructor for struct 'isolate'
    isolate(std::string init_id, int init_year, int init_sc, std::string init_serotype, bool init_vt, bool second_vt, std::vector<bool> *init_genotype, std::vector<bool> *init_markers) {
        id = init_id;
        year = init_year;
        sc = init_sc;
        serotype = init_serotype;
        vt = init_vt;
        latent_vt = second_vt;
        genotype = (*init_genotype);
        markers = (*init_markers);
    };
    
    // destructor for struct 'isolate'
//    ~isolate() {
//        genotype.clear();
//        genotype.shrink_to_fit();
//        markers.clear();
//        markers.shrink_to_fit();
//    }
//    
//    // copy constructor for struct 'isolate'
//    isolate(const isolate& other): id(other.id),year(other.year),sc(other.sc),serotype(other.serotype),vt(other.vt),latent_vt(other.latent_vt),genotype(other.genotype),markers(other.markers) {}
//    
//    // copy assignment operator
//    isolate& operator=(const isolate& other) {
//        // check for self-assignment
//        if (&other != this) {
//            id = other.id;
//            year = other.year;
//            sc = other.sc;
//            serotype = other.serotype;
//            vt = other.vt;
//            latent_vt = other.latent_vt;
//            genotype = other.genotype;
//            markers = other.markers;
//        }
//        return *this;
//    }
//    
//    // move constructor
//    isolate(isolate&& other) {
//        
//        id = other.id;
//        other.id = nullptr;
//        year = other.year;
//        other.year = nullptr;
//        sc = other.sc;
//        other.sc = nullptr;
//        serotype = other.serotype;
//        other.serotype = nullptr;
//        vt = other.vt;
//        other.vt = nullptr;
//        latent_vt = other.latent_vt;
//        other.latent_vt = nullptr;
//        genotype = other.genotype;
//        other.genotype = nullptr;
//        markers = other.markers;
//        other.markers = nullptr;
//        
//    }
    
    // move assignment operator
//    isolate& operator=(isolate&& other) {
//        if (&other != this) {
//            delete p;
//            p = other.p;
//            other.p = nullptr;
//
//        }
//        return *this;
//    }

};

// constructor
//isolate::isolate(std::string init_id, int init_year, int init_sc, std::string init_serotype, bool init_vt, bool second_vt, std::vector<bool> *init_genotype, std::vector<bool> *init_markers) {
//    id = init_id;
//    year = init_year;
//    sc = init_sc;
//    serotype = init_serotype;
//    vt = init_vt;
//    latent_vt = second_vt;
//    genotype = (*init_genotype);
//    markers = (*init_markers);
//}

// destructor
//isolate::~isolate() {
////    delete[] genotype;
//    genotype = nullptr;
////    delete[] markers;
//    markers = nullptr;
//}

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
    
    // copy constructor for struct 'cog'
    cog(const cog& other): id(other.id),vt(other.vt),weight(other.weight),currentFreq(other.currentFreq),actualFreq(other.actualFreq),simFreq(other.simFreq) {}
};

#endif
