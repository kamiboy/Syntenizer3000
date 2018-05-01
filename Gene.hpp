//
//  Gene.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#ifndef Gene_hpp
#define Gene_hpp

class Gene;
class Strain;
class Group;
class Relation;
class Contig;

#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include "Strain.hpp"
#include "Group.hpp"
#include "Relation.hpp"
#include "Score.hpp"
#include "Dataset.hpp"
#include "Contig.hpp"

using namespace std;

class Gene
{
public:
    string id;
    string id_full;
    string product;
    string sequence;
    string type;
    //string strain;
    Contig *contig;
    Strain *strain;
    //set<Strain*> strains;
    double coverage;
    //double gc;
    double GC3s;
    int start;
    int end;
    bool orientation;
    bool paralog;
    int copies;
    Group *group;
    //vector<Group*> neighbours;
    //vector<Gene*> neighbours;
    //Gene* neighbours[2][SIZE_OF_NEIGHTBOURHOOD];
    Gene* neighbours[SIZE_OF_NEIGHTBOURHOOD];
    //Group *neighbourhood[SIZE_OF_NEIGHTBOURHOOD];
    Gene *left_neighbour;
    Gene *right_neighbour;
    vector<Relation*> relations;
    vector<double> scoresMatrix;
    int length;
    int votes;
    int antiVotes;
    bool excluded;
    bool examined;
    double score;
    int matchHistogram[SIZE_OF_NEIGHTBOURHOOD+1];

    unordered_map<Group*, double> groupscores;
    unordered_map<Gene*, double> scores;
    unordered_map<Gene*, double> mirrorscores;

    void CompareNeighbours(Gene *gene);
    vector<Group *> Match(Dataset *);
    double Match(Group *);
    double CalculateGC3s();
    /*
    ~Gene(void)
    {
        printf("Gene %s is being deconstructed.\n", id.c_str());
    }*/
};

#endif /* Gene_hpp */
