//
//  Strain.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#ifndef Strain_hpp
#define Strain_hpp

class Gene;
class Strain;
class Group;
class Relation;
class Contig;

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "Settings.cpp"
#include "Gene.hpp"
#include "Contig.hpp"
#include "Database.hpp"

class Gene;
//class unordered_map;
using namespace std;

class Strain
{
public:
    string id;
    string id_full;
    string id_alt;
    string species;

    int bp;
    int ubp;
    int CDS;
    int tRNA;
    int tmRNA;
    int rRNA;
    double GC;
    //vector<Gene *> genes;
    vector<Contig *> contigs;
    unordered_map<string, Contig *> contigpool;
    Strain();

    void Parse(string, unordered_map<string, Gene *> *);
    void ExportChartData(string, Database *);
    void PrintContigs();
    void ExportPTT(string);
    /*
     Gene *GetGene(string geneID)
     {
     vector<Gene>::iterator gene;
     
     for (gene = this->genes.begin(); gene != this->genes.end(); gene++)
     {
     if (gene->id == geneID)
     return gene.base();
     }
     return NULL;
     }*/
};
//void ParseStrain(Strain *strain, string INPUT_DIRECTORY);
//void ParseStrain(Strain *strain);

#endif /* Strain_hpp */
