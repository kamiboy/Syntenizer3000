//
//  main.cpp
//  Syntenizer3000
//
//  Created by Camous.M. on 18/07/16.
//
//
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <tgmath.h>
#include <iostream>
#include <string>
#include <functional>
#include <unordered_map>
#include <set>
#include "Settings.hpp"
#include "Gene.hpp"
#include "Relation.hpp"
#include "Strain.hpp"
#include "Group.hpp"
//#include "Needleman-Wunsch.hpp"
#include "Parse.hpp"
#include "Database.hpp"
#include "Timer.hpp"
//hash<string> HashString;
#define GROUP_HISTOGRAM_SIZE    280
#define TEST_GENES              10
#define TEST_GROUPS             SIZE_OF_NEIGHTBOURHOOD

using namespace std;

/*
double ScoreSynteny(Gene* g1, Gene* g2, bool reverse, unordered_map<int64_t, double> *similaritypool)
{
    double score;
    //vector<Gene*> seq_1 = g1->neighbours;
    //vector<Gene*> seq_2 = g2->neighbours;
    vector<Gene*> seq_1_al;
    vector<Gene*> seq_2_al;
    double d = GAP_PENALTY;
    //int  L1 = SIZE_OF_NEIGHTBOURHOOD;
    //int  L2 = SIZE_OF_NEIGHTBOURHOOD;
    // Create alignment
    score = nw_align(  g1->neighbours, g2->neighbours, seq_1_al, seq_2_al, d, similaritypool );
    
#if DEBUG
    //int  L_al = seq_1_al.size();
    //cout << "Length after alignment: " << L_al << endl;
#endif

    //for( int i = 0; i <= L2; i++ )  delete F[ i ];
    //delete[] F;
    //for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
    //delete[] traceback;

    return  score;
}
*/
/*
double SequenceSimilarity(Gene *g1, Gene *g2)
{
    for (auto r = g1->relations.begin(); r != g1->relations.end(); r++)
        if ((*r)->gene == g2)
            return (1.0 - (*r)->score);

    return 0;
}
double ScoreSynteny(Gene *gene1, Gene *gene2, int n, bool mirrored)
{
    int i, j;
    double matches = 0.0f;
    double weight = 1.0;
    
    for (i = 0; i < n; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < n; j++)
        {
            double similarity = SequenceSimilarity(gene1->neighbours[i], gene2->neighbours[j]);
            
            if (similarity > 0)
            {
                int iloc = (i / 2)+1;
                int jloc = (j / 2)+1;
                
                if (i%2)
                    iloc *= -1;
                if (j%2)
                    jloc *= -1;
                
                if (mirrored)
                    jloc *= -1;
                
                partialmatch = similarity * (1.0f / (double)(abs(iloc-jloc)+1.0f));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        matches += maxpartialmatch * weight;
        
        if (i%2)
            weight *= 0.9;
    }
    
    return matches;
}
*/
/*
double Syntenize3000(Group *group, unordered_map<int64_t, double> *similaritypool)
{
    double scoresum = 0;
    
    if (group->genes.size() > 1000)
    {
        printf("Group %s - genes: %Lf is too large, skipping\n", group->id.c_str(), (long double)group->genes.size());
        return 0.0f;
    }
    printf("Syntenizing Group %s - genes: %Lf\n", group->id.c_str(), (long double)group->genes.size());
    
    for (auto g1 = group->genes.begin(); g1 != group->genes.end(); g1++)
        for (auto g2 = g1+1; g2 != group->genes.end(); g2++)
            if (*g1 == *g2)
                continue;
            else
            {
                //double score1 = ScoreSynteny(*g1, *g2, false, similaritypool);
                //double score2 = ScoreSynteny(*g1, *g2, true, similaritypool);
                double score1 = nw_align((*g1)->neighbours, (*g2)->neighbours, GAP_PENALTY, false, similaritypool);
                double score2 = nw_align((*g1)->neighbours, (*g2)->neighbours, GAP_PENALTY, true, similaritypool);
                
                if (score1 > 0 || score2 > 0)
                {
                    printf("none zero score!\n");
                }
                
                scoresum += (score1 > score2) ? score1 : score2;
            }
    
    printf("Group %s - scoresum: %f, genes: %Lf\n", group->id.c_str(), scoresum, (long double)group->genes.size());
    
    return(scoresum / (double)group->genes.size());
    
}

void PopulateNeighbours(vector<Strain*> *strains)
{
    vector<Strain*>::iterator strain;
    vector<Gene*>::iterator gene;
    
    for (strain = strains->begin(); strain != strains->end(); strain++)
    {
        for (gene = (*strain)->genes.begin(); gene != (*strain)->genes.end(); gene++)
        {
            int n;
            Gene *g;
            g = *gene;
            
            //linear neighbourhood
            for (n = 0; n < SIZE_OF_NEIGHTBOURHOOD/2; n++)
                g = g->left_neighbour;
            
            for (n = 0; n < SIZE_OF_NEIGHTBOURHOOD; n++)
            {
                if (g == *gene)
                    g = g->right_neighbour;

                (*gene)->neighbours[n] = g;
                g = g->right_neighbour;
                
            }

            // gene centered neighbourhood
            //g = *gene;
            //for (n = 0; n < SIZE_OF_NEIGHTBOURHOOD; n += 2)
            //{
            //    g = g->left_neighbour;
            //    (*gene)->neighbours[n] = g;
            //}
            
            //g = *gene;
            //for (n = 1; n < SIZE_OF_NEIGHTBOURHOOD; n += 2)
            //{
            //    g = g->right_neighbour;
            //    (*gene)->neighbours[n] = g;
            //}
        }
    }
}
*/
/*
int Assimilate(unordered_map<Gene *, Group *> *grouppool, vector<Group*> *groups, Group *group, Gene *gene, int counter)
{
    if (group == NULL)
    {
        counter++;
        group = new Group;
        group->id = "group"+to_string(counter);
        groups->push_back(group);
    }
    
    try
    {
        grouppool->at(gene);
    }
    catch (...)
    {
        (*grouppool)[gene] = group;
        gene->group = group;
        group->genes.push_back(gene);
        
        for (auto r = gene->relations.begin(); r != gene->relations.end(); r++)
            counter = Assimilate(grouppool, groups, group, (*r)->gene, counter);
    }

    //if (!group->HasGene(gene))
    //{
    //    grouppool[gene] = group;
    //    group->genes.push_back(gene);

    //    for (auto r = gene->relations.begin(); r != gene->relations.end(); r++)
    //        counter = Assimilate(group, (*r)->gene, counter);
    //}

    return counter;
}
int PopulateGroups(unordered_map<Gene *, Group *> *grouppool, vector<Group*> *groups, unordered_map<string, Gene *> *genepool)
{
    int counter = 0;

    for (auto g = genepool->begin(); g != genepool->end(); g++)
    {
        Gene *gene = g->second;
        
        try
        {
            grouppool->at(gene);
        }
        catch (...)
        {
            counter = Assimilate(grouppool, groups, NULL, gene, counter);
        }
    }
    
    return(counter);
}
 */

void Histogram(Dataset *dataset1, Dataset *dataset2, string filename)
{
    int maxmappedGenesSize = 0;
    int maxunmappedGenesSize = 0;
    int maxstrainsSize = 0;
    int maxscoreSize = 0;
    int maxgroupsetSize = 0;
    int maxmappedSizeHistogram = 0;
    int maxunmappedSizeHistogram = 0;
    int groups1to1 = 0;
    int groups1to2 = 0;
    int genes1to1 = 0;
    int genes1to2 = 0;
    int groupMappedGenesHistogram[GROUP_HISTOGRAM_SIZE];
    int groupUnmappedGenesHistogram[GROUP_HISTOGRAM_SIZE];
    int groupStrainsHistogram[GROUP_HISTOGRAM_SIZE];
    int groupStrains1to1Histogram[GROUP_HISTOGRAM_SIZE];
    int groupStrains1to2Histogram[GROUP_HISTOGRAM_SIZE];

    double groupStrainsSyntenySophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSyntenyOldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSyntenySimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSyntenyAdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSyntenyFastHistogram[GROUP_HISTOGRAM_SIZE];

    double groupScoreSyntenySophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSyntenyOldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSyntenySimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSyntenyAdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSyntenyFastHistogram[GROUP_HISTOGRAM_SIZE];

    double groupStrainsSynteny1to1SophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to1OldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to1SimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to1AdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to1FastHistogram[GROUP_HISTOGRAM_SIZE];
    
    double groupScoreSynteny1to1SophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to1OldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to1SimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to1AdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to1FastHistogram[GROUP_HISTOGRAM_SIZE];

    double groupStrainsSynteny1to2SophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to2OldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to2SimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to2AdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupStrainsSynteny1to2FastHistogram[GROUP_HISTOGRAM_SIZE];
    
    double groupScoreSynteny1to2SophisticatedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to2OldHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to2SimpleHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to2AdjustedHistogram[GROUP_HISTOGRAM_SIZE];
    double groupScoreSynteny1to2FastHistogram[GROUP_HISTOGRAM_SIZE];

    int groupToGroupsHistogram[GROUP_HISTOGRAM_SIZE];
    int mappedSizeHistogram[GROUP_HISTOGRAM_SIZE];
    int unmappedSizeHistogram[GROUP_HISTOGRAM_SIZE];
    //set<Strain*> superstrainset;
    printf("Generating %s\n", filename.c_str());
    for (int i = 0; i < GROUP_HISTOGRAM_SIZE; i++)
    {
        groupMappedGenesHistogram[i] = 0;
        groupUnmappedGenesHistogram[i] = 0;
        
        groupStrainsHistogram[i] = 0;
        groupStrains1to1Histogram[i] = 0;
        groupStrains1to2Histogram[i] = 0;
        
        groupStrainsSyntenySophisticatedHistogram[i] = 0;
        groupStrainsSyntenyOldHistogram[i] = 0;
        groupStrainsSyntenySimpleHistogram[i] = 0;
        groupStrainsSyntenyAdjustedHistogram[i] = 0;
        groupStrainsSyntenyFastHistogram[i] = 0;

        groupScoreSyntenySophisticatedHistogram[i] = 0;
        groupScoreSyntenyOldHistogram[i] = 0;
        groupScoreSyntenySimpleHistogram[i] = 0;
        groupScoreSyntenyAdjustedHistogram[i] = 0;
        groupScoreSyntenyFastHistogram[i] = 0;

        groupStrainsSynteny1to1SophisticatedHistogram[i] = 0;
        groupStrainsSynteny1to1OldHistogram[i] = 0;
        groupStrainsSynteny1to1SimpleHistogram[i] = 0;
        groupStrainsSynteny1to1AdjustedHistogram[i] = 0;
        groupStrainsSynteny1to1FastHistogram[i] = 0;

        groupScoreSynteny1to1SophisticatedHistogram[i] = 0;
        groupScoreSynteny1to1OldHistogram[i] = 0;
        groupScoreSynteny1to1SimpleHistogram[i] = 0;
        groupScoreSynteny1to1AdjustedHistogram[i] = 0;
        groupScoreSynteny1to1FastHistogram[i] = 0;

        groupStrainsSynteny1to2SophisticatedHistogram[i] = 0;
        groupStrainsSynteny1to2OldHistogram[i] = 0;
        groupStrainsSynteny1to2SimpleHistogram[i] = 0;
        groupStrainsSynteny1to2AdjustedHistogram[i] = 0;
        groupStrainsSynteny1to2FastHistogram[i] = 0;

        groupScoreSynteny1to2SophisticatedHistogram[i] = 0;
        groupScoreSynteny1to2OldHistogram[i] = 0;
        groupScoreSynteny1to2SimpleHistogram[i] = 0;
        groupScoreSynteny1to2AdjustedHistogram[i] = 0;
        groupScoreSynteny1to2FastHistogram[i] = 0;

        groupToGroupsHistogram[i] = 0;
        mappedSizeHistogram[i] = 0;
        unmappedSizeHistogram[i] = 0;
    }

    for (auto group = dataset1->groups.begin(); group != dataset1->groups.end(); group++)
    {
        set<Group*> groupset;
        set<Strain*> strainset;
        groupset.clear();
        strainset.clear();
        int mappedGenes = 0;
        int unmappedGenes = 0;
        //int strains = 0;
    
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            //mappedGenes += (*gene)->votes;
            //unmappedGenes += (*gene)->votes;

            strainset.insert((*gene)->strain);

            bool mapped = false;
            try
            {
                //Group *g = dataset2->grouppool.at((*gene)->strain->id + '|' + (*gene)->id);
                Group *g = dataset2->grouppool.at(*gene);
                groupset.insert(g);
                mapped = true;
            }
            catch (const out_of_range & e){}
            
            if (mapped)
                mappedSizeHistogram[(*gene)->length / 100]++;
            else
                unmappedSizeHistogram[(*gene)->length / 100]++;
            
            if (mapped && (*gene)->length / 100 > maxmappedSizeHistogram)
                maxmappedSizeHistogram = (*gene)->length/ 100;
            
            if (!mapped && (*gene)->length / 100 > maxunmappedSizeHistogram)
                maxunmappedSizeHistogram = (*gene)->length / 100;

        }
        unmappedGenes = (*group)->genes.size()-mappedGenes;
        //test
        if (mappedGenes >= GROUP_HISTOGRAM_SIZE /*|| unmappedGenes >= GROUP_HISTOGRAM_SIZE*/ || strainset.size() >= GROUP_HISTOGRAM_SIZE || groupset.size() >= GROUP_HISTOGRAM_SIZE)
        {
            printf("GROUP_HISTOGRAM_SIZE of %i is too small!\n", GROUP_HISTOGRAM_SIZE);
            exit(1);
        }
        
        groupMappedGenesHistogram[(mappedGenes*100) / (mappedGenes+unmappedGenes)]++;
        groupUnmappedGenesHistogram[(unmappedGenes *100) / (mappedGenes+unmappedGenes)]++;
        groupStrainsHistogram[strainset.size()]++;
        
        groupStrainsSyntenySophisticatedHistogram[strainset.size()] += (*group)->syntenyScoreSophisticated;
        groupStrainsSyntenyOldHistogram[strainset.size()] += (*group)->syntenyScoreOld;
        groupStrainsSyntenySimpleHistogram[strainset.size()] += (*group)->syntenyScoreSimple;
        groupStrainsSyntenyAdjustedHistogram[strainset.size()] += (*group)->syntenyScoreAdjusted;
        groupStrainsSyntenyFastHistogram[strainset.size()] += (*group)->syntenyScoreFast;
        
        groupScoreSyntenySophisticatedHistogram[(int)(*group)->syntenyScoreSophisticated]++;
        groupScoreSyntenyOldHistogram[(int)(*group)->syntenyScoreOld]++;
        groupScoreSyntenySimpleHistogram[(int)(*group)->syntenyScoreSimple]++;
        groupScoreSyntenyAdjustedHistogram[(int)(*group)->syntenyScoreAdjusted]++;
        groupScoreSyntenyFastHistogram[(int)(*group)->syntenyScoreFast]++;
    
        groupToGroupsHistogram[groupset.size()]++;
        
        if (groupset.size() == 1)
        {
            groupStrains1to1Histogram[strainset.size()]++;
            groups1to1++;
            genes1to1 += (*group)->genes.size();
            groupStrainsSynteny1to1SophisticatedHistogram[strainset.size()] += (*group)->syntenyScoreSophisticated;
            groupStrainsSynteny1to1OldHistogram[strainset.size()] += (*group)->syntenyScoreOld;
            groupStrainsSynteny1to1SimpleHistogram[strainset.size()] += (*group)->syntenyScoreSimple;
            groupStrainsSynteny1to1AdjustedHistogram[strainset.size()] += (*group)->syntenyScoreAdjusted;
            groupStrainsSynteny1to1FastHistogram[strainset.size()] += (*group)->syntenyScoreFast;

            groupScoreSynteny1to1SophisticatedHistogram[(int)(*group)->syntenyScoreSophisticated]++;
            groupScoreSynteny1to1OldHistogram[(int)(*group)->syntenyScoreOld]++;
            groupScoreSynteny1to1SimpleHistogram[(int)(*group)->syntenyScoreSimple]++;
            groupScoreSynteny1to1AdjustedHistogram[(int)(*group)->syntenyScoreAdjusted]++;
            groupScoreSynteny1to1FastHistogram[(int)(*group)->syntenyScoreFast]++;
        }
        else
        {
            groupStrains1to2Histogram[strainset.size()]++;
            groups1to2++;
            genes1to2 += (*group)->genes.size();
    
            groupStrainsSynteny1to2SophisticatedHistogram[strainset.size()] += (*group)->syntenyScoreSophisticated;
            groupStrainsSynteny1to2OldHistogram[strainset.size()] += (*group)->syntenyScoreOld;
            groupStrainsSynteny1to2SimpleHistogram[strainset.size()] += (*group)->syntenyScoreSimple;
            groupStrainsSynteny1to2AdjustedHistogram[strainset.size()] += (*group)->syntenyScoreAdjusted;
            groupStrainsSynteny1to2FastHistogram[strainset.size()] += (*group)->syntenyScoreFast;
    
            groupScoreSynteny1to2SophisticatedHistogram[(int)(*group)->syntenyScoreSophisticated]++;
            groupScoreSynteny1to2OldHistogram[(int)(*group)->syntenyScoreOld]++;
            groupScoreSynteny1to2SimpleHistogram[(int)(*group)->syntenyScoreSimple]++;
            groupScoreSynteny1to2AdjustedHistogram[(int)(*group)->syntenyScoreAdjusted]++;
            groupScoreSynteny1to2FastHistogram[(int)(*group)->syntenyScoreFast]++;
        }
        
        if ((mappedGenes*100) / (mappedGenes+unmappedGenes) > maxmappedGenesSize)
            maxmappedGenesSize = (mappedGenes*100) / (mappedGenes+unmappedGenes);
        if ((unmappedGenes *100) / (mappedGenes+unmappedGenes) > maxunmappedGenesSize)
            maxunmappedGenesSize = (unmappedGenes *100) / (mappedGenes+unmappedGenes);
        if (strainset.size() > maxstrainsSize)
            maxstrainsSize = strainset.size();
        if (groupset.size() > maxgroupsetSize)
            maxgroupsetSize = groupset.size();

        if ((*group)->syntenyScoreSophisticated > maxscoreSize)
            maxscoreSize = (*group)->syntenyScoreSophisticated;
        if ((*group)->syntenyScoreOld > maxscoreSize)
            maxscoreSize = (*group)->syntenyScoreOld;
        if ((*group)->syntenyScoreSimple > maxscoreSize)
            maxscoreSize = (*group)->syntenyScoreSimple;
        if ((*group)->syntenyScoreAdjusted > maxscoreSize)
            maxscoreSize = (*group)->syntenyScoreAdjusted;
        if ((*group)->syntenyScoreFast > maxscoreSize)
            maxscoreSize = (*group)->syntenyScoreFast;
        
        if (strainset.size() == 0)
        {
            printf("%s has empty strainset!\n", (*group)->id.c_str());
        }
    }
    
    //printf("maxmappedGenesSize: %i\n", maxmappedGenesSize);
    //printf("maxunmappedGenesSize: %i\n", maxunmappedGenesSize);
    //printf("maxstrainsSize: %i\n", maxstrainsSize);
    //printf("maxgroupsetSize: %i\n", maxgroupsetSize);
    //printf("superstrainset.size(): %i\n", superstrainset.size());
    //printf("database1->groups.size(): %i\n", database1->groups.size());
    
    ofstream out( filename, ifstream::out );
    
    for (int i = 0; i < GROUP_HISTOGRAM_SIZE; i++)
        out << '\t' << i;
    
    out << "\nMapped Genes in Group";
    for (int i = 0; i <= maxmappedGenesSize; i++)
        out << '\t' << groupMappedGenesHistogram[i];
    
    out << "\nUnmapped Genes in Group";
    for (int i = 0; i <= maxunmappedGenesSize; i++)
        out << '\t' << groupUnmappedGenesHistogram[i];
    
    out << "\nStrains in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << (double)(100.0f*groupStrainsHistogram[i]) / (double)dataset1->groups.size();

    out << "\nStrains in Group 1 to 1";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << (double)(100.0f*groupStrains1to1Histogram[i]) / (double)groups1to1;

    out << "\nStrains in Group 1 to !1";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << (double)(100.0f*groupStrains1to2Histogram[i]) / (double)groups1to2;

    out << "\nStrains Synteny Sophisticated in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSyntenySophisticatedHistogram[i] / (double)groupStrainsHistogram[i];
    out << "\nStrains Synteny Old in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSyntenyOldHistogram[i]/ (double)groupStrainsHistogram[i];
    out << "\nStrains Synteny Simple in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSyntenySimpleHistogram[i]/ (double)groupStrainsHistogram[i];
    out << "\nStrains Synteny Adjusted in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSyntenyAdjustedHistogram[i]/ (double)groupStrainsHistogram[i];
    out << "\nStrains Synteny Fast in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSyntenyFastHistogram[i]/ (double)groupStrainsHistogram[i];
    
    out << "\nScore Synteny Sophisticated";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSyntenySophisticatedHistogram[i];
    out << "\nScore Synteny Old";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSyntenyOldHistogram[i];
    out << "\nScore Synteny Simple";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSyntenySimpleHistogram[i];
    out << "\nScore Synteny Adjusted";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSyntenyAdjustedHistogram[i];
    out << "\nScore Synteny Fast";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSyntenyFastHistogram[i];

    out << "\nGenes in Group 1 to 1";
    out << '\t' << genes1to1;
    
    out << "\nGenes in Group 1 to !1";
    out << '\t' << genes1to2;

    out << "\nStrains Synteny 1to1 Sophisticated in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to1SophisticatedHistogram[i] / (double)groupStrains1to1Histogram[i];
    out << "\nStrains Synteny 1to1 Old in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to1OldHistogram[i]/ (double)groupStrains1to1Histogram[i];
    out << "\nStrains Synteny 1to1 Simple in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to1SimpleHistogram[i]/ (double)groupStrains1to1Histogram[i];
    out << "\nStrains Synteny 1to1 Adjusted in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to1AdjustedHistogram[i]/ (double)groupStrains1to1Histogram[i];
    out << "\nStrains Synteny 1to1 Fast in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to1FastHistogram[i]/ (double)groupStrains1to1Histogram[i];
    
    out << "\nStrains Synteny 1to!1 Sophisticated in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to2SophisticatedHistogram[i] / (double)groupStrains1to2Histogram[i];
    out << "\nStrains Synteny 1to!1 Old in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to2OldHistogram[i]/ (double)groupStrains1to2Histogram[i];
    out << "\nStrains Synteny 1to!1 Simple in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to2SimpleHistogram[i]/ (double)groupStrains1to2Histogram[i];
    out << "\nStrains Synteny 1to!1 Adjusted in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to2AdjustedHistogram[i]/ (double)groupStrains1to2Histogram[i];
    out << "\nStrains Synteny 1to!1 Fast in Group";
    for (int i = 0; i <= maxstrainsSize; i++)
        out << '\t' << groupStrainsSynteny1to2FastHistogram[i]/ (double)groupStrains1to2Histogram[i];

    out << "\nScore Synteny 1to1 Sophisticated";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to1SophisticatedHistogram[i];
    out << "\nScore Synteny 1to1 Old";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to1OldHistogram[i];
    out << "\nScore Synteny 1to1 Simple";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to1SimpleHistogram[i];
    out << "\nScore Synteny 1to1 Adjusted";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to1AdjustedHistogram[i];
    out << "\nScore Synteny 1to1 Fast";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to1FastHistogram[i];
    
    out << "\nScore Synteny 1to!1 Sophisticated";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to2SophisticatedHistogram[i];
    out << "\nScore Synteny 1to!1 Old";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to2OldHistogram[i];
    out << "\nScore Synteny 1to!1 Simple";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to2SimpleHistogram[i];
    out << "\nScore Synteny 1to!1 Adjusted";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to2AdjustedHistogram[i];
    out << "\nScore Synteny 1to!1 Fast";
    for (int i = 0; i <= maxscoreSize; i++)
        out << '\t' << groupScoreSynteny1to2FastHistogram[i];
    //
    out << "\nGroups in Group";
    for (int i = 0; i <= maxgroupsetSize; i++)
        out << '\t' << (double)(groupToGroupsHistogram[i]*100.0f) / (double)dataset1->groups.size();
    
    out << "\nMapped gene sizes";
    for (int i = 0; i <= maxmappedSizeHistogram; i++)
        out << '\t' << mappedSizeHistogram[i];
    
    out << "\nUnmapped gene sizes";
    for (int i = 0; i <= maxunmappedSizeHistogram; i++)
        out << '\t' << unmappedSizeHistogram[i];
    
    out << "\n";
    out.close();
    printf("Done Generating %s\n", filename.c_str());
}

/*
void Histogram(Database *database1, Database *database2, string filename)
{
    int maxmappedGenesSize = 0;
    int maxunmappedGenesSize = 0;
    int maxstrainsSize = 0;
    int maxgroupsetSize = 0;
    int maxmappedSizeHistogram = 0;
    int maxunmappedSizeHistogram = 0;
    int groupMappedGenesHistogram[GROUP_HISTOGRAM_SIZE];
    int groupUnmappedGenesHistogram[GROUP_HISTOGRAM_SIZE];
    int groupStrainsHistogram[GROUP_HISTOGRAM_SIZE];
    int groupToGroupsHistogram[GROUP_HISTOGRAM_SIZE];
    int mappedSizeHistogram[GROUP_HISTOGRAM_SIZE];
    int unmappedSizeHistogram[GROUP_HISTOGRAM_SIZE];
    set<Strain*> superstrainset;
    printf("Generating %s\n", filename.c_str());
    for (int i = 0; i < GROUP_HISTOGRAM_SIZE; i++)
    {
        groupMappedGenesHistogram[i] = 0;
        groupUnmappedGenesHistogram[i] = 0;
        groupStrainsHistogram[i] = 0;
        groupToGroupsHistogram[i] = 0;
        mappedSizeHistogram[i] = 0;
        unmappedSizeHistogram[i] = 0;
    }

    for (auto group = database1->groups.begin(); group != database1->groups.end(); group++)
    {
        set<Group*> groupset;
        set<Strain*> strainset;
        groupset.clear();
        strainset.clear();
        int mappedGenes = 0;
        int unmappedGenes = 0;
        //int strains = 0;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            mappedGenes += (*gene)->antiVotes - (*gene)->votes;
            unmappedGenes += (*gene)->votes;
            
            mappedSizeHistogram[(*gene)->id.size() / 100] += (*gene)->antiVotes - (*gene)->votes;
            unmappedSizeHistogram[(*gene)->id.size() / 100]  += (*gene)->votes;
            
            if ((*gene)->antiVotes - (*gene)->votes && (*gene)->id.size() / 100 > maxmappedSizeHistogram)
                maxmappedSizeHistogram = (*gene)->id.size() / 100;
            
            if ((*gene)->votes && (*gene)->id.size() / 100 > maxunmappedSizeHistogram)
                maxunmappedSizeHistogram = (*gene)->id.size() / 100;


            if ((*gene)->strains.empty())
            {
                printf("%s has gene with no strain!\n", (*group)->id.c_str());
                printf("%s\n", (*gene)->id.c_str());
            }
            strainset.insert((*gene)->strains.begin(), (*gene)->strains.end());
            superstrainset.insert((*gene)->strains.begin(), (*gene)->strains.end());
            try
            {
                Group *g = database2->grouppool2.at((*gene)->id);
                groupset.insert(g);
            }
            catch (const out_of_range & e){}
        }
        //test
        if (mappedGenes >= GROUP_HISTOGRAM_SIZE || unmappedGenes >= GROUP_HISTOGRAM_SIZE || strainset.size() >= GROUP_HISTOGRAM_SIZE || groupset.size() >= GROUP_HISTOGRAM_SIZE)
        {
            printf("GROUP_HISTOGRAM_SIZE of %i is too small!\n", GROUP_HISTOGRAM_SIZE);
            exit(1);
        }

        groupMappedGenesHistogram[(mappedGenes*100) / (mappedGenes+unmappedGenes)]++;
        groupUnmappedGenesHistogram[(unmappedGenes *100) / (mappedGenes+unmappedGenes)]++;
        groupStrainsHistogram[strainset.size()]++;
        groupToGroupsHistogram[groupset.size()]++;
        
        if ((mappedGenes*100) / (mappedGenes+unmappedGenes) > maxmappedGenesSize)
            maxmappedGenesSize = (mappedGenes*100) / (mappedGenes+unmappedGenes);
        if ((unmappedGenes *100) / (mappedGenes+unmappedGenes) > maxunmappedGenesSize)
            maxunmappedGenesSize = (unmappedGenes *100) / (mappedGenes+unmappedGenes);
        if (strainset.size() > maxstrainsSize)
            maxstrainsSize = strainset.size();
        if (groupset.size() > maxgroupsetSize)
            maxgroupsetSize = groupset.size();
        
        if (strainset.size() == 0)
        {
            printf("%s has empty strainset!\n", (*group)->id.c_str());
        }
    }
    
    //printf("maxmappedGenesSize: %i\n", maxmappedGenesSize);
    //printf("maxunmappedGenesSize: %i\n", maxunmappedGenesSize);
    //printf("maxstrainsSize: %i\n", maxstrainsSize);
    //printf("maxgroupsetSize: %i\n", maxgroupsetSize);
    //printf("superstrainset.size(): %i\n", superstrainset.size());
    //printf("database1->groups.size(): %i\n", database1->groups.size());
    
    ofstream out( filename, ifstream::out );
    
    for (int i = 0; i < GROUP_HISTOGRAM_SIZE; i++)
        out << '\t' << i;
    
    out << "\nMapped Genes in Group";
    for (int i = 0; i < maxmappedGenesSize; i++)
        out << '\t' << groupMappedGenesHistogram[i];

    out << "\nUnmapped Genes in Group";
    for (int i = 0; i < maxunmappedGenesSize; i++)
        out << '\t' << groupUnmappedGenesHistogram[i];

    out << "\nStrains in Group";
    for (int i = 0; i < maxstrainsSize; i++)
        out << '\t' << (100*groupStrainsHistogram[i]) / database1->groups.size();
    
    out << "\nGroups in Group";
    for (int i = 0; i < maxgroupsetSize; i++)
        out << '\t' << (groupToGroupsHistogram[i]*100) / database1->groups.size();

    out << "\nMapped gene sizes";
    for (int i = 0; i < maxmappedSizeHistogram; i++)
        out << '\t' << mappedSizeHistogram[i];

    out << "\nUnmapped gene sizes";
    for (int i = 0; i < maxunmappedSizeHistogram; i++)
        out << '\t' << unmappedSizeHistogram[i];

    out << "\n";
    out.close();
    printf("Done Generating %s\n", filename.c_str());
}

void Analyze()
{
    Database *db = new Database();
    Dataset *orthoMCL = new Dataset(db);
    Dataset *proteinOrtho = new Dataset(db);
    Timer timer;
    
    timer.Start();
    Parse::fnaGroups0("/Users/Camous/Desktop/Masters/group_alns_orthomcl/", orthoMCL);
    printf("Parsing took:\t%i seconds\n\n", (int)timer.Elapsed());
    timer.Start();
    Parse::fnaGroups0("/Users/Camous/Desktop/Masters/group_alns_proteinortho/", proteinOrtho);
    printf("Parsing took:\t%i seconds\n\n", (int)timer.Elapsed());
    
    unordered_map<Group*, set<Group*>*> groupmap1, groupmap2;
    
    int mappedGenes = 0;
    int unmappedGenes = 0;
    uint64_t unsplitnucleotides = 0;
    uint64_t splitnucleotides = 0;
    int splitgenes = 0;
    int unsplitgenes = 0;

    for (auto group = orthoMCL->groups.begin(); group != orthoMCL->groups.end(); group++)
    {
        set<Group*> *groupset = new set<Group*>;
        
        int nucleotides = 0;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            nucleotides += (*gene)->length;
            try
            {
                //Group *g = proteinOrtho->grouppool.at((*gene)->strain->id + '|' + (*gene)->id);
                Group *g = proteinOrtho->grouppool.at(*gene);
                groupset->insert(g);
                
                mappedGenes++;
            }
            catch (const out_of_range & e)
            {
                unmappedGenes++;
            }
        }
        groupmap1[(*group)] = groupset;
        if (groupset->size() == 1)
        {
            unsplitnucleotides += nucleotides;
            unsplitgenes++;
        }
        else
        {
            splitnucleotides += nucleotides;
            splitgenes++;
        }
    }
    printf("Mapped Genes:\t%i\t(%.2f%%)\n", mappedGenes, mappedGenes * 100.0f / orthoMCL->genecount);
    printf("Unmapped Genes:\t%i\t(%.2f%%)\n\n", unmappedGenes, unmappedGenes * 100.0f / orthoMCL->genecount);
    mappedGenes = 0;
    unmappedGenes = 0;
    for (auto group = proteinOrtho->groups.begin(); group != proteinOrtho->groups.end(); group++)
    {
        set<Group*> *groupset = new set<Group*>;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                //Group *g = orthoMCL->grouppool.at((*gene)->strain->id + '|' + (*gene)->id);
                Group *g = orthoMCL->grouppool.at(*gene);
                groupset->insert(g);
                mappedGenes++;
            }
            catch (const out_of_range & e)
            {
                unmappedGenes++;
            }
        }
        groupmap2[(*group)] = groupset;
    }
    printf("Mapped Genes:\t%i\t(%.2f%%)\n", mappedGenes, mappedGenes * 100.0f / proteinOrtho->genecount);
    printf("Unmapped Genes:\t%i\t(%.2f%%)\n\n", unmappedGenes, unmappedGenes * 100.0f / proteinOrtho->genecount);
    
    int emptyMaps = 0;
    int nonemptyMaps = 0;
    int groupsInNonemptyMaps = 0;
    int totalMaps = 0;
    int groupsInAllMaps = 0;
    
    ofstream out( "/Users/Camous/Desktop/Masters/orthomcl_to_proteinortho_groupmap.txt", ifstream::out );
    
    for (auto map = groupmap1.begin(); map != groupmap1.end(); map++)
    {
        out << map->first->id << ':';
        
        totalMaps++;
        groupsInAllMaps += map->second->size();
        
        if (map->second->size() == 0)
            emptyMaps++;
        else
        {
            nonemptyMaps++;
            groupsInNonemptyMaps += map->second->size();
        }
        
        for (auto group = map->second->begin(); group != map->second->end(); group++)
            out << ' ' << (*group)->id;
        out << "\n";
    }
    out.close();
    printf("Orthomcl groups with no mappings:\t%i\t(%.2f%%)\n", emptyMaps, emptyMaps * 100.0f / orthoMCL->groups.size() );
    printf("Average groups that each orthomcl group maps to in proteinortho:\t%f\n", (double) groupsInAllMaps / totalMaps );
    //printf("The genes in each nonempty orthomcl group map to %f groups in proteinortho on average\n\n", (double) groupsInNonemptyMaps / nonemptyMaps );
    
    emptyMaps = 0;
    nonemptyMaps = 0;
    groupsInNonemptyMaps = 0;
    totalMaps = 0;
    groupsInAllMaps = 0;
    out.open( "/Users/Camous/Desktop/Masters/proteinortho_to_orthomcl_groupmap.txt", ifstream::out );
    
    for (auto map = groupmap2.begin(); map != groupmap2.end(); map++)
    {
        out << map->first->id << ':';
        totalMaps++;
        groupsInAllMaps += map->second->size();
        
        if (map->second->size() == 0)
            emptyMaps++;
        else
        {
            nonemptyMaps++;
            groupsInNonemptyMaps += map->second->size();
        }
        
        for (auto group = map->second->begin(); group != map->second->end(); group++)
            out << ' ' << (*group)->id;
        out << "\n";
    }
    out.close();
    printf("proteinortho groups with no mappings:\t%i\t(%.2f%%)\n", emptyMaps, emptyMaps * 100.0f / proteinOrtho->groups.size() );
    printf("Average groups that each proteinortho group maps to in orthomcl:\t%f\n\n", (double) groupsInAllMaps / totalMaps );
    //printf("The genes in each nonempty proteinortho group map to %f groups in orthomcl on average\n\n", (double) groupsInNonemptyMaps / nonemptyMaps );
    
    uint64_t unmappedNucleotides = 0;
    uint64_t mappedNucleotides = 0;
    unmappedGenes = 0;
    out.open( "/Users/Camous/Desktop/Masters/group_alns_orthomcl_unmapped.txt", ifstream::out );
    for (auto group = orthoMCL->groups.begin(); group != orthoMCL->groups.end(); group++)
    {
        bool singleFound = false;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if ((*gene)->votes == 0)
            {
                if (!singleFound)
                {
                    singleFound = true;
                    out << (*group)->id << ": ";
                }
                out << (*gene)->strain->id + '|'  +(*gene)->id << ' ';
                unmappedGenes++;
                unmappedNucleotides += (*gene)->length;
            }
            else
                mappedNucleotides += (*gene)->length;
        }
        if (singleFound)
            out << '\n';
    }
    out.close();
    
    printf("average nucleotides per mapped gene:\t%f\n\n", (double) mappedNucleotides / mappedGenes);
    //printf("unmapped genes in orthomcl: %i (%.2f%%)\n", unmappedGenes, unmappedGenes*100.0f / orthomclgenes);
    printf("average nucleotides per unmapped orthoMCL gene:\t%f\n\n", (double) unmappedNucleotides / unmappedGenes);
    
    unmappedNucleotides = 0;
    unmappedGenes = 0;
    out.open( "/Users/Camous/Desktop/Masters/group_alns_proteinortho_unmapped.txt", ifstream::out );
    for (auto group = proteinOrtho->groups.begin(); group != proteinOrtho->groups.end(); group++)
    {
        bool singleFound = false;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if ((*gene)->votes == 0)
            {
                if (!singleFound)
                {
                    singleFound = true;
                    out << (*group)->id << ": ";
                }
                out << (*gene)->strain->id + '|'  +(*gene)->id << ' ';
                unmappedGenes++;
                unmappedNucleotides += (*gene)->length;
            }
        }
        if (singleFound)
            out << '\n';
    }
    out.close();
    //printf("unmapped genes in proteinortho: %i (%.2f%%)\n", unmappedGenes, unmappedGenes*100.0f / proteinorthogenes);
    printf("average nucleotides per unmapped proteinOrtho gene:\t%f\n\n", (double) unmappedNucleotides / unmappedGenes);
    
    Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_histogram.txt");
    Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho_histogram.txt");
}
*/
void Split(Dataset *dataset1, Dataset *dataset2, string OUTPUT_FILE, bool truncate_id = false)
{
    //unordered_map<Group*, set<Group*>*> groupmap;
    
    //int mappedGenes = 0;
    //int unmappedGenes = 0;
    Timer timer;
    timer.Start();
    printf("Splitting Groups:\n");
    ofstream out1( OUTPUT_FILE + "1to1.groups", ifstream::out );
    ofstream out2( OUTPUT_FILE + "1to!1.groups", ifstream::out );
    for (auto group = dataset1->groups.begin(); group != dataset1->groups.end(); group++)
    {
        set<Group*> groupset;
        ofstream *out;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                Group *g = dataset2->grouppool.at(*gene);
                groupset.insert(g);
            }
            catch (const out_of_range & e){}
        }
        //groupmap[(*group)] = groupset;
        if (groupset.size() == 1)
            out = &out1;
        else
            out = &out2;
        
        *out << (*group)->id << ':';
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            if (truncate_id)
                *out << ' ' << (*gene)->strain->id + '|' + (*gene)->id;
            else
                *out << ' ' << (*gene)->strain->id_full + '|' + (*gene)->id_full;
        *out << "\n";
        
    }
    out1.close();
    out2.close();
    printf("Splitting took: %i seconds\n\n", (int)timer.Elapsed());
    
}

void GenerateGroupmap(Dataset *dataset1, Dataset *dataset2, string OUTPUT_FILE)
{
    //unordered_map<Group*, set<Group*>*> groupmap;
    
    Timer timer;
    timer.Start();
    printf("Mapping Groups:\n");
    ofstream out( OUTPUT_FILE, ifstream::out );
    for (auto group = dataset1->groups.begin(); group != dataset1->groups.end(); group++)
    {
        out << (*group)->id << ':';

        set<Group*> groupset;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                Group *g = dataset2->grouppool.at(*gene);

                if (!groupset.count(g))
                {
                    groupset.insert(g);
                    out << ' ' << g->id;
                }
            }
            catch (const out_of_range & e){}
        }
        out << '\n';
    }
    out.close();
    printf("Mapping took: %i seconds\n\n", (int)timer.Elapsed());
    
}
/*
void Convert(string input_directory, string output_group_file)
{
    Database *db = new Database();
    Dataset *dataset = new Dataset(db);
    Timer timer;
    timer.Start();
    Parse::fnaGroups0(input_directory, dataset);
    
    ofstream out( output_group_file, ifstream::out );
    
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        out << (*group)->id << ':';
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            out << ' ' << (*gene)->strain->id << '|' << (*gene)->id;
        }
        out << "\n";
    }
    out.close();
    printf("Conversion took:\t%i seconds\n\n", (int)timer.Elapsed());
}*/
/*
void Test()
{
    int i;
    Gene gene100percent1;
    Gene gene100percent2;
    Gene gene50percent1;
    Gene gene50percent2;
    Gene gene0percent1;
    Gene gene0percent2;
    Gene geneMixedpercent1;
    Gene geneMixedpercent2;

    Group group100, group50, group0, groupMixed;
    Group groups1[SIZE_OF_NEIGHTBOURHOOD+2];
    Group groups2[SIZE_OF_NEIGHTBOURHOOD+2];
    Gene genes1[SIZE_OF_NEIGHTBOURHOOD];
    Gene genes2[SIZE_OF_NEIGHTBOURHOOD];
    Gene genes3[SIZE_OF_NEIGHTBOURHOOD];
    Gene genes4[SIZE_OF_NEIGHTBOURHOOD];
    Gene genes5[SIZE_OF_NEIGHTBOURHOOD];
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD+2; i++)
    {
        groups1[i].id = "1-" + ((i > 9) ? to_string(i) : "0" + to_string(i));
        groups2[i].id = "2-" + ((i > 9) ? to_string(i) : "0" + to_string(i));
    }

    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        genes1[i].group = &groups1[i];
        genes2[i].group = &groups1[i];
        if (i < SIZE_OF_NEIGHTBOURHOOD/2 )
            genes3[i].group = &groups1[i];
        else
            genes3[i].group = &groups2[i];
        genes4[i].group = &groups2[i];
        if (i < SIZE_OF_NEIGHTBOURHOOD/2 )
            genes5[i].group = &groups1[i];
        else
            genes5[i].group = &groups1[i+2];
    }
    //gene100percent1.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    //gene100percent2.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    group100.genes.push_back(&gene100percent1);
    group100.genes.push_back(&gene100percent2);
    
    //gene50percent1.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    //gene50percent2.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    group50.genes.push_back(&gene50percent1);
    group50.genes.push_back(&gene50percent2);

    //gene0percent1.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    //gene0percent2.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    group0.genes.push_back(&gene0percent1);
    group0.genes.push_back(&gene0percent2);

    //geneMixedpercent1.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    //geneMixedpercent2.neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
    groupMixed.genes.push_back(&geneMixedpercent1);
    groupMixed.genes.push_back(&geneMixedpercent2);

    //for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    //{
    //    gene100percent1.neighbours[i] = &genes1[i];
    //    gene100percent2.neighbours[i] = &genes2[SIZE_OF_NEIGHTBOURHOOD-1-i];
    //    gene50percent1.neighbours[i] = &genes1[i];
    //    gene50percent2.neighbours[i] = &genes3[SIZE_OF_NEIGHTBOURHOOD-1-i];
    //    gene0percent1.neighbours[i] = &genes1[i];
    //    gene0percent2.neighbours[i] = &genes4[SIZE_OF_NEIGHTBOURHOOD-1-i];
    //    geneMixedpercent1.neighbours[i] = &genes1[i];
    //    geneMixedpercent2.neighbours[i] = &genes5[SIZE_OF_NEIGHTBOURHOOD-1-i];
    //}
    
    printf("Testing:\n");
    printf("SIZE_OF_NEIGHTBOURHOOD:\t%i\n", SIZE_OF_NEIGHTBOURHOOD);
    
    printf("Sophisticated score 100%%:\t%f\n", group100.SyntenizeSophisticated());
    printf("Old score 100%%:\t%f\n", group100.SyntenizeOld());
    printf("Simple score 100%%:\t%f\n\n", group100.SyntenizeSimple());
    printf("Adjusted score 100%%:\t%f\n\n", group100.SyntenizeAdjusted());

    printf("Sophisticated score 100%%:\t%f\n", group50.SyntenizeSophisticated());
    printf("Old score 50%%:\t%f\n", group50.SyntenizeOld());
    printf("Simple score 50%%:\t%f\n\n", group50.SyntenizeSimple());
    printf("Adjusted score 50%%:\t%f\n\n", group50.SyntenizeAdjusted());

    printf("Sophisticated score 0%%:\t%f\n", group0.SyntenizeSophisticated());
    printf("Old score 0%%:\t%f\n", group0.SyntenizeOld());
    printf("Simple score 0%%:\t%f\n\n", group0.SyntenizeSimple());
    printf("Adjusted score 0%%:\t%f\n\n", group0.SyntenizeAdjusted());

    printf("Sophisticated score Mixed%%:\t%f\n", groupMixed.SyntenizeSophisticated());
    printf("Old score Mixed%%:\t%f\n", groupMixed.SyntenizeOld());
    printf("Simple score Mixed%%:\t%f\n\n", groupMixed.SyntenizeSimple());
    printf("Adjusted score Mixed%%:\t%f\n\n", groupMixed.SyntenizeAdjusted());

    printf("Done testing\n");

    //Database *db = new Database();
    //Dataset *testset = new Dataset(db);
    
    //Gene genes1[TEST_GENES];
    //Gene genes2[TEST_GENES];

    //for (int group = 0; group < TEST_GROUPS; group++)
    //{
    //    for (int gene = 0; gene < TEST_GENES; gene++)
    //    {
    //        genes1[gene].id = "Gene" + to_string(gene);
    //
    //    }
    //}
}
*/
void GenerateGeneticRelationMatrix(Dataset *dataset, string OUTPUT_FILE, bool header = true)
//void Test(Dataset *dataset, string OUTPUT_FILE, bool header)
{
    Timer timer;
    timer.Start();
    double **presence_absence_matrix = new double*[dataset->groups.size()];
    double **genetic_relation_matrix = new double*[dataset->database->strains.size()];

    for (int i = 0; i < dataset->groups.size(); i++)
        presence_absence_matrix[i] = new double[dataset->database->strains.size()];
    for (int i = 0; i < dataset->database->strains.size(); i++)
        genetic_relation_matrix[i] = new double[dataset->database->strains.size()];

    printf("Generating genetic relation matrix\n");
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        {
            bool ingroup = false;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->strain == (*strain))
                    ingroup = true;
            //out << (ingroup ? '1' : '0');
            //printf("%i,%i ", distance(dataset->groups.begin(), group), distance(dataset->database->strains.begin(), strain));
            presence_absence_matrix[distance(dataset->groups.begin(), group)][distance(dataset->database->strains.begin(), strain)] = (ingroup ? 2 : 0);
        }
    }
    
    int count = 0;
    double diagonal_mean = 0;

    for (int i = 0; i < dataset->database->strains.size(); i++)
    {
        for (int j = 0; j < dataset->database->strains.size(); j++)
        {
            genetic_relation_matrix[i][j] = 0;
            for (int k = 0; k < dataset->groups.size(); k++)
                genetic_relation_matrix[i][j] += presence_absence_matrix[k][i]*presence_absence_matrix[k][j];
            
            if (i == j)
            {
                count++;
                diagonal_mean += genetic_relation_matrix[i][j];
            }
        }
    }
    
    diagonal_mean /= count;

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    if (header)
        //for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
        {
            //string id = (*strain)->id_full;
            //string::size_type loc;
            //if ((loc = id.rfind('-')) != string::npos)
            //    id.erase(id.begin(), id.begin()+loc+1);
            //out << ';' << "SM" + id;
            out << ';' << (*strain)->id_alt;
        }
    for (int i = 0; i <  dataset->database->strains.size(); i++)
    {
        if (header || (!header && i != 0 ) )
            out << '\n';
        
        if (header)
        {
            //string id = dataset->database->strains[i]->id_full;
            //string::size_type loc;
            //if ((loc = id.rfind('-')) != string::npos)
            //    id.erase(id.begin(), id.begin()+loc+1);
            //out << "SM" + id;
            out << dataset->database->strains[i]->id_alt;
        }
        for (int j = 0; j <  dataset->database->strains.size(); j++)
        {
            //set<Strain *> strains;
            if (header || (!header && j != 0) )
                out << ';';
            out << genetic_relation_matrix[i][j] / diagonal_mean;
            //out << genetic_relation_matrix[i][j];
        }
    }

    out.close();

    printf("Generating genetic relation matrix took:\t%i seconds\n\n", (int)timer.Elapsed());
}

void GeneratePresenceAbsenceMatrix(Dataset *dataset, string OUTPUT_FILE, bool header = true)
{
    Timer timer;
    timer.Start();
    
    printf("Generating presence absence matrix into %s:\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );
    
    if (header)
    {
        out << "Strain";
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
            out << '\t' << (*strain)->id_alt;
    }
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        if (header || (!header && group != dataset->groups.begin()) )
            out << '\n';
        if (header)
            out << (*group)->id;

        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
        {
            bool ingroup = false;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->strain == (*strain))
                    ingroup = true;
            if (header || (!header && strain != dataset->database->strains.begin()) )
                out << '\t';
            out << (ingroup ? '1' : '0');
        }
    }
    
    out << '\n';
    out.close();
    
    printf("Generating presence absence matrix took:\t%i seconds\n\n", (int)timer.Elapsed());
}

/*
void GeneratePresenceAbsenceMatrix(Dataset *dataset, string OUTPUT_FILE, bool header = true)
{
    Timer timer;
    timer.Start();

    printf("Generating presence absence matrix into %s:\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    if (header)
    {
        out << "Strain";
        for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
            out << ';' << (*group)->id;
    }
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        if (header || (!header && strain != dataset->database->strains.begin()) )
            out << '\n';
        if (header)
            out << (*strain)->id_alt;

        for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        {
            bool ingroup = false;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->strain == (*strain))
                    ingroup = true;
            if (header || (!header && group != dataset->groups.begin()) )
                out << ';';
            out << (ingroup ? '1' : '0');
        }
    }
    
    out << '\n';
    out.close();

    printf("Generating presence absence matrix took:\t%i seconds\n\n", (int)timer.Elapsed());
}
*/
void ConvertStrainNames(string INPUT_DIRECTORY, Database *db)
{
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    
    timer.Start();
    printf("Parsing strains from directory %s:\n", INPUT_DIRECTORY.c_str());
    //int genes =
    while ((dp = readdir(dirp)) != NULL)
    {
        string strainID = string(dp->d_name);
        
        if (strainID.size() < 5)
            continue;
        
        if (strainID.compare(strainID.size()-4,4,".gff"))
            continue;
        
        string id1, id2;
        strainID.resize(strainID.size()-4);
        string::size_type loc;
        if ((loc = strainID.rfind('-')) != string::npos)
        {
            id1 = strainID;
            id2 = strainID;
            id1.erase(id1.begin()+loc,id1.end());
            id2.erase(id2.begin(),id2.begin()+loc+1);
            id2 = "SM" + id2;
        }
        
        Strain *strain = db->strainpool.at(id1);
        //strain->id = id2;
        strain->id = strainID;
        db->strainpool[strain->id] = strain;
    }
}

void ExportGroups(Dataset *dataset, string OUTPUT_FILE, bool trauncate_id = false, bool renumerate = false)
{
    int groups = 1;
    Timer timer;
    timer.Start();
    printf("Exporting Groups:\n");
    ofstream out( OUTPUT_FILE, ifstream::out );
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++, groups++)
    {
        if (renumerate)
            out << "group" + to_string(groups);
        else
            out << (*group)->id;

        if (!std::isnan((*group)->algebraicConnectivity))
            out << '|' << (*group)->algebraicConnectivity;
        out << ':';
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            if (trauncate_id)
                out << ' ' << (*gene)->strain->id + '|' + (*gene)->id;
            else
                out << ' ' << (*gene)->strain->id_full + '|' + (*gene)->id_full;
        out << "\n";
    }
    out.close();
    printf("Exporting took: %i seconds\n\n", (int)timer.Elapsed());
}

void MatchGenestoGroups(vector<Gene*> genes, Dataset *dataset, string OUTPUT_FILE )
{
    int i = 0;
    Timer timer;
    timer.Start();
    printf("Matching genes to groups:\n");

    Progress p(genes.size()-1);

    ofstream out( OUTPUT_FILE, ifstream::out );

    for (auto gene = genes.begin(); gene != genes.end(); gene++, i++)
    {
       out << (*gene)->strain->id + '|' + (*gene)->id + ":";// + " (" + (*gene)->group->id + ";" + to_string((*gene)->groupscores.at((*gene)->group)) + ")";
        vector<Group *> matchedgroups = (*gene)->Match(dataset);
        printf("%s: ", ((*gene)->strain->id + '|' + (*gene)->id).c_str());
        for (int i = 0; i <= SIZE_OF_NEIGHTBOURHOOD; i++)
            printf("\t%i", (*gene)->matchHistogram[i]);
        printf("\n");
        //matchedgroups.insert(matchedgroups.begin(), (*gene)->group);
        for (auto group = matchedgroups.begin(); group != matchedgroups.end(); group++)
        {
            out << " (" + (*group)->id + " - " + to_string((*gene)->groupscores.at(*group)) + ")";
        }
        out << '\n';
        //p.Update(i);
    }
    out.close();
    printf("Matching genes took:\t%i seconds\n\n", (int)timer.Elapsed());
}

double MatchGenetoGroup(string geneID, string groupID, Dataset *dataset)
{
    Gene *gene = NULL;
    Group *group = NULL;
    double matchscore;

    for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        if ((*g)->id == groupID)
        {
            group = *g;
            break;
        }

    try
    {
        gene = dataset->database->genepool.at(geneID);
    }
    catch (const out_of_range & e) {}
    
    if (group == NULL)
    {
        printf("Group %s not found!\n", groupID.c_str());
        return -1;
    }
    
    if (gene == NULL)
    {
        printf("Gene %s not found!\n", geneID.c_str());
        return -1;
    }
    
    matchscore = gene->Match(group);
    
    printf("Match score of gene %s to group %s (%f): %f\n", gene->id_full.c_str(), group->id.c_str(), group->syntenyScoreSophisticated / (double) SIZE_OF_NEIGHTBOURHOOD, matchscore);

    return
        matchscore;
}

int AnalyzeParalogs(Dataset *dataset, string OUTPUT_FILE)
{
    int pure_paralog_groups = 0;
    int unambiguous_groups = 0;
    int ambiguous_groups = 0;
    int split_groups = 0;
    int two_paralogs = 0;
    int three_paralogs = 0;
    int four_paralogs = 0;
    int five_paralogs = 0;
    int six_paralogs = 0;
    int sevenplus_paralogs = 0;
    Timer timer;
    timer.Start();
    printf("Analyzing Paralogs:\n");
    //ofstream out_groups( OUTPUT_FILE, ifstream::out );
    ofstream out_analysis( OUTPUT_FILE + ".analysis", ifstream::out );

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        bool has_paralogs = false;
        unordered_map<Strain *, vector<Gene*>*> strains;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                auto v = strains.at((*gene)->strain);
                v->push_back(*gene);
                has_paralogs = true;
            }
            catch (const out_of_range & e)
            {
                strains[(*gene)->strain] = new vector<Gene*>;
                strains[(*gene)->strain]->push_back(*gene);
            }
        }
        
        if (has_paralogs)
        {
            out_analysis << (*group)->id + ": (";

            vector <Gene*> paralogs;
            vector <Gene*> homologs;
            
            for (auto strain = strains.begin(); strain != strains.end(); strain++)
            {
                if ((*strain).second->size() == 2)
                    two_paralogs++;
                if ((*strain).second->size() == 3)
                    three_paralogs++;
                if ((*strain).second->size() == 4)
                    four_paralogs++;
                if ((*strain).second->size() == 5)
                    five_paralogs++;
                if ((*strain).second->size() == 6)
                    six_paralogs++;
                if ((*strain).second->size() >= 7)
                    sevenplus_paralogs++;

                if ((*strain).second->size() > 1)
                    paralogs.insert(paralogs.end(), (*strain).second->begin(), (*strain).second->end());
                else
                    homologs.insert(homologs.end(), (*strain).second->begin(), (*strain).second->end());
            }
            
            if (homologs.size() == 0)
            {
                pure_paralog_groups++;
                continue;
            }
            
            Strain *prev_strain = (*paralogs.begin())->strain;
            
            bool split = false;
            int total_scoresum = 0;
            int max_scoresum = 0;
            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
            {
                double scoresum = 0.0f;
                (*g1)->paralog = true;

                for (auto g2 = homologs.begin(); g2 != homologs.end(); g2++)
                    scoresum += Score::Sophisticated(*g1, *g2);
                
                if (prev_strain != (*g1)->strain)
                {
                    out_analysis << ")\t\t(";
                    prev_strain = (*g1)->strain;
                }
                scoresum /= (double)homologs.size();
                
                for (auto g2 = paralogs.begin(); g2 != paralogs.end(); g2++)
                    if (g1 != g2 && scoresum > 0 && Score::Sophisticated(*g1, *g2) == 0)
                        split = true;
                out_analysis << ' ' + (*g1)->strain->id + '|' + (*g1)->id + " - " + to_string(scoresum);
                
                total_scoresum += scoresum;
                if (scoresum > max_scoresum)
                    max_scoresum = scoresum;
            }
            if (split)
                split_groups++;
            out_analysis << " )\n";

            if (max_scoresum > 0 && total_scoresum == max_scoresum)
                unambiguous_groups++;
            else
                ambiguous_groups++;
        }
    }
    //out_groups.close();
    out_analysis.close();
    printf("Pure paralog groups: %i\n", pure_paralog_groups);
    printf("Split groups: %i\n", split_groups);
    printf("Unambiguous paralog groups: %i\n", unambiguous_groups);
    printf("Ambiguous paralog groups: %i\n", ambiguous_groups);
    printf("two paralogs: %i\n", two_paralogs);
    printf("three paralogs: %i\n", three_paralogs);
    printf("four paralogs: %i\n", four_paralogs);
    printf("five paralogs: %i\n", five_paralogs);
    printf("six paralogs: %i\n", six_paralogs);
    printf("seven or more paralogs: %i\n", sevenplus_paralogs);
    printf("Analysis took:\t%i seconds\n\n", (int)timer.Elapsed());
    return 0;
}

int CompareContigs(int contig, Dataset *dataset)
{
    unordered_map<string, int> groupID;
    int counter = 0;
    
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        //int counter = 0;
        //groupID.clear();
        printf("%s: ", (*strain)->id.c_str());
        for (auto gene = (*strain)->contigs[contig]->genes.begin(); gene != (*strain)->contigs[contig]->genes.end(); gene++)
        {
            if ((*gene)->group == NULL)
                printf("----  ");
            else
            {
                try
                {
                    groupID.at((*gene)->group->id);
                }
                catch (const out_of_range & e)
                {
                    groupID[(*gene)->group->id] = counter;
                    counter++;
                }
                printf("%04i  ", groupID[(*gene)->group->id]);
            }
        }
        printf("\n");
    }
    return 0;
}

int RearrangeContigs(Dataset *dataset, string INPUT_DIRECTORY, bool truncateid = false)
{
    //int genes = 0;
    int strains = 0;
    int total_contigs = 0;
    int shortcontigs = 0;
    int emptycontigs = 0;
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    set<Gene*> duplicatedgenes;

    timer.Start();
    printf("Parsing updated contigs from directory %s:\n", INPUT_DIRECTORY.c_str());
    //int genes = 0;

    while ((dp = readdir(dirp)) != NULL)
    {
        string filename = string(dp->d_name);
        if (filename.size() < 5)
            continue;

        if (filename.compare(filename.size()-4, 4, ".txt"))
            continue;

        ifstream input = ifstream( INPUT_DIRECTORY + filename);

        vector<Contig *> new_contigs;
        Strain *strain = NULL;

        string line;
        while (getline(input, line))
        {
            //int counter = 0;

            if (line == "")
                continue;

            stringstream sstream(line);
            string s;

            vector<string> items;
            items.clear();
            
            while (getline(sstream, s, ' '))
            {
                if (s != "")
                    items.push_back(s);
            }
            string::size_type loc;
            string strainID, scaffoldID, scaffold;
            
            strainID = scaffoldID = items[0];
            if (strain == NULL)
            {
                if ((loc = strainID.find_first_of("_")) != string::npos)
                    strainID.erase(strainID.begin()+loc,strainID.end());
                if ((loc = strainID.rfind('-')) != string::npos)
                    strainID.erase(strainID.begin()+loc,strainID.end());

                strain = dataset->database->strainpool.at(strainID);
                strains++;
            }
            scaffold = items[1].substr(sizeof("contigs__")-1);

            stringstream sstream2(scaffold);
            string tmp, cont;

            total_contigs++;
            Contig *new_contig = new Contig;
            new_contig->id = scaffoldID;
            new_contig->strain = strain;
            new_contig->length = 0;
            //string test;
            //if ((loc = new_contig->id.find("_len")) != string::npos)
            //    new_contig->length = atoi(new_contig->id.substr(loc+5, new_contig->id.length()).c_str());
                //test = new_contig->id.substr(loc+5, new_contig->id.length());

            new_contigs.push_back(new_contig);

            if (truncateid && (loc = new_contig->id.find("scaf_")) != string::npos)
            {
                new_contig->id.erase(0, loc+5);
                if ((loc = new_contig->id.find("_len")) != string::npos)
                    new_contig->id.erase(loc, new_contig->id.length());
            }
            
            //int nucleotides = 0;
            while (getline(sstream2, tmp, ':'))
            {
                stringstream sstream3(tmp);

                while (getline(sstream3, cont, '.'))
                {
                    bool contig_orientation = cont[cont.length()-1] == 'f';
                    string contigID = cont.substr(0, cont.length()-1);

                    try
                    {
                        Contig *contig = strain->contigpool.at(contigID);
                        if (contig_orientation)
                        {
                            for (auto gene = contig->genes.begin(); gene != contig->genes.end(); gene++)
                            {
                                (*gene)->copies++;
                                if ((*gene)->copies == 1)
                                {
                                    //(*gene)->start = nucleotides + (*gene)->start;
                                    (*gene)->start = new_contig->length + (*gene)->start;
                                    (*gene)->contig = new_contig;
                                    new_contig->genes.push_back(*gene);
                                    new_contig->orientations.push_back((*gene)->orientation);
                                }
                            }
                        }
                        else
                        {
                            for (auto gene = contig->genes.rbegin(); gene != contig->genes.rend(); gene++)
                            {
                                (*gene)->copies++;
                                if ((*gene)->copies == 1)
                                {
                                    //(*gene)->start = nucleotides + (contig->length - ((*gene)->start + (*gene)->length));
                                    (*gene)->start = new_contig->length + (contig->length - ((*gene)->start + (*gene)->length));
                                    (*gene)->contig = new_contig;
                                    new_contig->genes.push_back(*gene);
                                    new_contig->orientations.push_back(!(*gene)->orientation);
                                }
                            }
                        }
                        new_contig->length += contig->length;
                        //nucleotides += contig->length;
                    }
                    catch (const out_of_range & e)
                    {
                        printf("Skipping contig %s %s\n", strain->id_full.c_str(), contigID.c_str());
                    }
                }
            }
            //new_contig->length = nucleotides;
        }
        strain->contigpool.clear();
        strain->contigs.clear();

        for (auto contig = new_contigs.rbegin(); contig != new_contigs.rend(); contig++)
        {
            if ((*contig)->genes.size() == 0)
            {
                emptycontigs++;
                continue;
            }

            strain->contigpool[(*contig)->id] = *contig;
            strain->contigs.insert(strain->contigs.begin(), *contig);
            
            Gene *previousGene = (*contig)->genes.back();
            
            for (vector<Gene*>::iterator gene = (*contig)->genes.begin(); gene !=  (*contig)->genes.end(); gene++)
            {
                previousGene->right_neighbour = *gene;
                (*gene)->left_neighbour = previousGene;
                previousGene = *gene;
            }
            
            bool isContigLargeEnough;
            int sizeOfContigNeighbourhood;
            
            if ((*contig)->genes.size() < SIZE_OF_NEIGHTBOURHOOD + 1)
            {
                shortcontigs++;
                isContigLargeEnough = false;
                sizeOfContigNeighbourhood = (*contig)->genes.size()-1;
            }
            else
            {
                isContigLargeEnough = true;
                sizeOfContigNeighbourhood = SIZE_OF_NEIGHTBOURHOOD;
            }

            for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            {
                int index;
                //int n;
                //(*gene)->copies++;
                if ((*gene)->copies > 1)
                    duplicatedgenes.insert(*gene);

                for (int n = 0; n < SIZE_OF_NEIGHTBOURHOOD; n++)
                    (*gene)->neighbours[n] = NULL;

                if ((*contig)->orientations[distance((*contig)->genes.begin(), gene)])
                {
                    for (int offset = -sizeOfContigNeighbourhood/2, n = 0; offset <= sizeOfContigNeighbourhood/2; offset++, n++)
                    {
                        if (offset == 0)
                        {
                            n--;
                            continue;
                        }
                        index = ((*contig)->genes.size() + distance((*contig)->genes.begin(), gene) + offset)%(*contig)->genes.size();
                        (*gene)->neighbours[n+(SIZE_OF_NEIGHTBOURHOOD-sizeOfContigNeighbourhood)/2] = (*contig)->genes[index];
                    }
                }
                else
                {
                    for (int offset = sizeOfContigNeighbourhood/2, n = 0; offset >= -sizeOfContigNeighbourhood/2; offset--, n++)
                    {
                        if (offset == 0)
                        {
                            n--;
                            continue;
                        }
                        index = ((*contig)->genes.size() + distance((*contig)->genes.begin(), gene) + offset)%(*contig)->genes.size();
                        (*gene)->neighbours[n+(SIZE_OF_NEIGHTBOURHOOD-sizeOfContigNeighbourhood)/2] = (*contig)->genes[index];
                    }
                }
            }
        }
        input.close();
    }
    printf("Empty contigs: %i\n", emptycontigs);
    printf("Short contigs: %i\n", shortcontigs);
    printf("Duplicated genes: %lu\n", duplicatedgenes.size());

    int duplicatesum = 0;
    for (auto gene = duplicatedgenes.begin(); gene != duplicatedgenes.end(); gene++)
        duplicatesum += (*gene)->copies - 1;

    printf("Duplicate genes sum: %i\n", duplicatesum);
    printf("Parsing %i contigs in %i strains took:\t%i seconds\n\n", total_contigs, strains, (int)timer.Elapsed());

    return 0;
}

int GenerateSyntenyMap(Dataset *dataset, string OUTPUT_FILE )
{
    Group *coregroup = NULL;
    Timer timer;
    timer.Start();
    printf("Generating Synteny Map:\n");
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        bool has_paralogs = false;
        set<Strain *> strains;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if (strains.count((*gene)->strain))
                has_paralogs = true;
            else
                strains.insert((*gene)->strain);
        }
        
        if (!has_paralogs && (*group)->genes.size() == dataset->database->strains.size())
            if (coregroup == NULL || (coregroup != NULL && (*group)->syntenyScoreSophisticated > coregroup->syntenyScoreSophisticated))
            {
                bool iscore = true;
                for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                    if ((*gene)->contig != (*gene)->strain->contigs[0])
                        iscore = false;

                if (iscore)
                    coregroup = (*group);
            }
    }

    if (coregroup == NULL)
    {
        printf("Error: No core group found!\n");
        //return(1);
    }
    else
        printf("Core group: %s - %f\n", coregroup->id.c_str(), coregroup->syntenyScoreSophisticated);

    ofstream out( OUTPUT_FILE, ifstream::out );

    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            out << (*contig)->id;
            
            int offset = 0;
            
            if (coregroup != NULL && contig == (*strain)->contigs.begin())
            {
                for (offset = 0; offset < (*contig)->genes.size(); offset++)
                    if ((*contig)->genes[offset]->group == coregroup)
                        break;
                
                if (offset == (*contig)->genes.size())
                {
                    printf("Error: Core group %s not found in strain %s\n", coregroup->id.c_str(), (*strain)->id_full.c_str());
                    return(1);
                }
            }

            for (int index = 0; index < (*contig)->genes.size(); index++)
            {
                Gene *gene = (*contig)->genes[(index + offset)%(*contig)->genes.size()];
                
                out << ";" + ((gene->group != NULL) ? to_string(gene->group->syntenyScoreAdjusted / (double)SIZE_OF_NEIGHTBOURHOOD): to_string(0.0f));
            }
            out << '\n';
        }
    }
    out.close();
    printf("Generating Synteny Map took\t%i seconds\n\n", (int)timer.Elapsed());

    return 0;
}

int GenerateSyntenyMaps(Dataset *dataset, string OUTPUT_DIRECTORY )
{
    Group *coregroup = NULL;
    Timer timer;
    timer.Start();
    printf("Generating Synteny Maps:\n");
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        bool has_paralogs = false;
        set<Strain *> strains;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if (strains.count((*gene)->strain))
                has_paralogs = true;
            else
                strains.insert((*gene)->strain);
        }
        
        if (!has_paralogs && (*group)->genes.size() == dataset->database->strains.size())
            if (coregroup == NULL || (coregroup != NULL && (*group)->syntenyScoreSophisticated > coregroup->syntenyScoreSophisticated))
            {
                bool iscore = true;
                for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                    if ((*gene)->contig != (*gene)->strain->contigs[0])
                        iscore = false;
                
                if (iscore)
                    coregroup = (*group);
            }
    }
    
    if (coregroup == NULL)
    {
        printf("Error: No core group found!\n");
        //return(1);
    }
    else
        printf("Core group: %s - %f\n", coregroup->id.c_str(), coregroup->syntenyScoreSophisticated);
    
    //ofstream out( OUTPUT_FILE, ifstream::out );
    
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            size_t loc;
            string filename = (*contig)->id;
            if ((loc = filename.find("chromosome")) != string::npos)
                filename.erase(filename.begin()+loc+4, filename.end());
            if ((loc = filename.find("plasmid")) != string::npos)
                filename.erase(filename.begin()+loc+4, filename.end());
            if ((loc = filename.find("fragment")) != string::npos)
                filename.erase(filename.begin()+loc+4, filename.end());

            if (filename.size() < 5)
                continue;
            
            ofstream out( OUTPUT_DIRECTORY + /*(*strain)->id_full + "_" +*/  filename + ".syntenymap", ifstream::out );
            //out << (*contig)->id;
            
            int offset = 0;
            
            if (coregroup != NULL && contig == (*strain)->contigs.begin())
            {
                for (offset = 0; offset < (*contig)->genes.size(); offset++)
                    if ((*contig)->genes[offset]->group == coregroup)
                        break;
                
                if (offset == (*contig)->genes.size())
                {
                    printf("Error: Core group %s not found in strain %s\n", coregroup->id.c_str(), (*strain)->id_full.c_str());
                    return(1);
                }
            }
            out << "position;synteny\n";

            for (int index = 0; index < (*contig)->genes.size(); index++)
            {
                Gene *gene = (*contig)->genes[(index + offset)%(*contig)->genes.size()];
                
                out << int(gene->start + gene->length * 0.5f) << ";" << ((gene->group != NULL) ?gene->group->syntenyScoreAdjusted / (double)SIZE_OF_NEIGHTBOURHOOD: -0.5f) << '\n';
            }
            out.close();
        }
    }
    printf("Generating Synteny Maps took\t%i seconds\n\n", (int)timer.Elapsed());

    return 0;
}

int FindCoreGroups(Dataset *dataset)
{
    vector<Group *> coregroups;

    printf("Core Groups:\n");

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        bool has_paralogs = false;
        unordered_map<Strain *, vector<Gene*>*> strains;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                auto v = strains.at((*gene)->strain);
                v->push_back(*gene);
                has_paralogs = true;
            }
            catch (const out_of_range & e)
            {
                strains[(*gene)->strain] = new vector<Gene*>;
                strains[(*gene)->strain]->push_back(*gene);
            }
        }
        
        if (!has_paralogs && (*group)->genes.size() == dataset->database->strains.size() && (*group)->syntenyScoreSophisticated >= dataset->database->strains.size() * 0.75f)
            printf("%s\n", (*group)->id.c_str());
            //coregroups.push_back(*group);
    }
    printf("End of Core Groups\n");

    return 0;
}

int SyntenyScoreSubset(Dataset *dataset, string INPUT_FILE)
{
    int counter = 0;
    double scoresumSophisticated = 0.0f;
    double scoresumOld = 0.0f;
    double scoresumSimple = 0.0f;
    double scoresumAdjusted = 0.0f;
    double scoresumFast = 0.0f;
    double syntenyScoreSophisticated;
    double syntenyScoreOld;
    double syntenyScoreSimple;
    double syntenyScoreAdjusted;
    double syntenyScoreFast;
    Timer timer;
    
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("Error: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }

    timer.Start();
    string line;
    while (getline(input, line))
    {
        counter++;
        stringstream sstream(line);
        
        if (line.empty() || counter == 1)
            continue;
        
        //if (line[0] != 'g' || line[1] != 'r' || line[2] != 'o' || line[3] != 'u' || line[4] != 'p')
        //    continue;
        
        string item;
        vector<string> items;
        while (getline(sstream, item, ';'))
            items.push_back(item);

        if (items[0] == "")
            continue;

        Group *group = NULL;

        for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        {
            if ((*g)->id == items[0])
            {
                group = *g;
                break;
            }
        }
        
        if (group == NULL)
            printf("Error: Group %s not \n", items[0].c_str());
        else
        {
            scoresumSophisticated += group->syntenyScoreSophisticated;
            scoresumOld += group->syntenyScoreOld;
            scoresumSimple += group->syntenyScoreSimple;
            scoresumAdjusted += group->syntenyScoreAdjusted;
            scoresumFast += group->syntenyScoreFast;
        }
    }
    
    syntenyScoreSophisticated = scoresumSophisticated / (double)(counter-1);
    syntenyScoreOld = scoresumOld / (double)(counter-1);
    syntenyScoreSimple = scoresumSimple / (double)(counter-1);
    syntenyScoreAdjusted = scoresumAdjusted / (double)(counter-1);
    syntenyScoreFast = scoresumFast / (double)(counter-1);

    printf("Average Synteny Score Sophisticated:\t[%.6f]\n", syntenyScoreSophisticated);
    printf("Average Synteny Score Old:\t[%.5f]\n", syntenyScoreOld);
    printf("Average Synteny Score Simple:\t[%.5f]\n", syntenyScoreSimple);
    printf("Average Synteny Score Adjusted:\t[%.5f]\n", syntenyScoreAdjusted);
    printf("Average Synteny Score Fast:\t[%.5f]\n", syntenyScoreFast);

    return 0;
}

void AnalyzeDiscordantGroups(Dataset *dataset, string groupIDs)
{
    stringstream sstream(groupIDs);
    
    string id;
    while (getline(sstream, id, ';'))
    {
        Group *group = NULL;
        
        for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        {
            if ((*g)->id == id)
            {
                group = *g;
                break;
            }
        }

        if (group == NULL)
            printf("Error: Group %s not \n", id.c_str());
        else
        {
            unordered_map<string, double> speciescores_fast;
            unordered_map<string, double> speciescores_sophisticated;
            unordered_map<string, int> speciescount;
            
            group->syntenyScoreFast = group->SyntenizeFast();
            group->syntenyScoreSophisticated = group->SyntenizeSophisticated();
            
            set<string> specieset;
            printf("\n%s (%.2f):\n", group->id.c_str(),  group->syntenyScoreSophisticated);
            //for (auto g1 = group->genes.begin(); g1 != group->genes.end(); g1++)
            ///    printf("\t%s", (*g1)->id_full.c_str());
            //printf("\n");
            set<Strain *> strainset;
            for (auto g1 = group->genes.begin(); g1 != group->genes.end(); g1++)
                strainset.insert((*g1)->strain);

            for (auto g1 = group->genes.begin(); g1 != group->genes.end(); g1++)
            {
                for (auto g2 = g1+1; g2 != group->genes.end(); g2++)
                {
                    double score_fast = Score::Fast(*g1, *g2);
                    double score_sophisticated = Score::Sophisticated(*g1, *g2);
                    (*g1)->score += score_sophisticated;
                    (*g2)->score += score_sophisticated;
                    string tag1;
                    string tag2;

                    if ((*g1)->strain->species.empty() || (*g2)->strain->species.empty())
                        continue;

                    tag1 = (*g1)->strain->species + (*g2)->strain->species;
                    tag2 = (*g2)->strain->species + (*g1)->strain->species;
                    specieset.insert((*g1)->strain->species);
                    specieset.insert((*g2)->strain->species);

                    if (speciescount.count(tag1))
                    {
                        speciescores_fast[tag1] += score_fast;
                        speciescores_sophisticated[tag1] += score_sophisticated;
                        speciescount[tag1]++;
                    }
                    else
                    {
                        speciescores_fast[tag1] = score_fast;
                        speciescores_sophisticated[tag1] = score_sophisticated;
                        speciescount[tag1] = 1;
                    }
                    if (tag1 != tag2)
                    {
                        if (speciescount.count(tag2))
                        {
                            speciescores_fast[tag2] += score_fast;
                            speciescores_sophisticated[tag2] += score_sophisticated;
                            speciescount[tag2]++;
                        }
                        else
                        {
                            speciescores_fast[tag2] = score_fast;
                            speciescores_sophisticated[tag2] = score_sophisticated;
                            speciescount[tag2] = 1;
                        }
                    }
                }
                //printf("\t%.2f", (*g1)->score / double(group->genes.size()));

                printf("\t%s\t%.2f\t%s\n", (*g1)->id_full.c_str(), (*g1)->score / double(group->genes.size()-1), (*g1)->strain->species.c_str());
            }
            /*
            for (auto strain = strainset.begin(); strain != strainset.end(); strain++)
            {
                printf("\t%s\t%s", (*strain)->id_alt.c_str(), (*strain)->species.c_str());
                for (auto gene = group->genes.begin(); gene != group->genes.end(); gene++)
                    if (*strain == (*gene)->strain)
                        printf("\t%.2f",(*gene)->score / double(group->genes.size()-1) );
                printf("\n");
            }
            */
            //printf("\n\n%s (%lu):\tFast(%f)\tSophisticated(%f)\n", group->id.c_str(), group->genes.size(), group->syntenyScoreFast, group->syntenyScoreSophisticated);
            
            for (auto s1 = specieset.begin(); s1 != specieset.end(); s1++)
            {
                string tag = *s1 + *s1;
                
                //if (speciescount.count(tag))
                    //printf("%s:\t%.2f\t%.2f\t(%d)\n", tag.c_str(), speciescores_fast[tag] / (double)speciescount[tag], speciescores_sophisticated[tag] / (double)speciescount[tag], speciescount[tag]);
                
                for (auto s2 = specieset.begin(); s2 != specieset.end(); s2++)
                {
                    if (*s1 == *s2)
                        continue;

                    tag = *s1 + *s2;

                    //if (speciescount.count(tag))
                        //printf("%s:\t%.2f\t%.2f\t(%d)\n", tag.c_str(), speciescores_fast[tag] / (double)speciescount[tag], speciescores_sophisticated[tag] / (double)speciescount[tag], speciescount[tag]);
                }
                //printf("\n");
            }
            /*
            for (auto s = speciescores.begin(); s != speciescores.end(); s++)
            {
                printf("%s: %f\n", s->first.c_str(), s->second / (double)speciescount[s->first]);
            }*/
        }
    }
}

void GetGeneGroupScores(string geneID, Dataset *dataset)
{
    Gene *g1;

    try
    {
        g1 = dataset->database->genepool.at(geneID);
    }
    catch (const out_of_range & e)
    {
        printf("Gene %s not found!\n", geneID.c_str());
        return;
    }

    printf("%s %s:\n", g1->id_full.c_str(), g1->strain->species.c_str());
    for (auto g2 = g1->group->genes.begin(); g2 != g1->group->genes.end(); g2++)
        if (g1 != *g2)
            printf("\t%s\t%s\t%.2f\n", (*g2)->id_full.c_str(), (*g2)->strain->species.c_str(), Score::Simple(g1, *g2));
    printf("\n");
}

void TestSplit(Dataset *dataset)
{
    int orphans = 0;
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        if ((*group)->CountParalogs() > 0)
        {
            (*group)->Split5();
            orphans += (*group)->orphans;
        }
    
    printf("orphans %i\n", orphans);
}

void ExportPTT(Database *db, string OUTPUT_DIRECTORY)
{
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
        (*strain)->ExportPTT(OUTPUT_DIRECTORY);
}

void Disambiguate(Dataset *dataset, string OUTPUT_FILE)
{
    int groups = 0;
    Timer timer;
    timer.Start();
    printf("Disambigating groups into %s:\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    //set<Strain *> orphanstrains;
    unordered_map<Strain *, int> orphanstrains;
    unordered_map<string, int> orphancontigs;

    int orphans = 0;
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        if ((*group)->CountParalogs() > 0)
        {
            groups++;
            vector<Group*> splitgroups = (*group)->Split5();
            orphans += (*group)->orphans;
            
            for (auto group2 = splitgroups.begin(); group2 != splitgroups.end(); group2++)
            {
                if ((*group2)->genes.size() == 1)
                {
                    //orphanstrains.insert((*group2)->genes[0]->strain);
                    if (orphancontigs.count((*group2)->genes[0]->contig->id_full))
                        orphancontigs[(*group2)->genes[0]->contig->id]++;
                    else
                        orphancontigs[(*group2)->genes[0]->contig->id] = 1;
                    //orphancontigs.insert((*group2)->genes[0]->)
                    if (orphanstrains.count((*group2)->genes[0]->strain))
                        orphanstrains[(*group2)->genes[0]->strain]++;
                    else
                        orphanstrains[(*group2)->genes[0]->strain] = 1;

                    continue;
                }
                out << (*group2)->id << ':';
                //out << (*group2)->id;
                //if ((*group2)->algebraicConnectivity != NAN)
                //    out << '|' << (*group2)->algebraicConnectivity;
                //out << ':';
                for (auto gene = (*group2)->genes.begin(); gene != (*group2)->genes.end(); gene++)
                    out << ' ' << (*gene)->strain->id_full << '|' << (*gene)->id_full;
                out << '\n';
            }
        }
        else
        {
            out << (*group)->id;
            if (!std::isnan((*group)->algebraicConnectivity))
                out << '|' << (*group)->algebraicConnectivity;
            out << ':';
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                out << ' ' << (*gene)->strain->id_full << '|' << (*gene)->id_full;
            out << '\n';
        }
    printf("orphans %i\n", orphans);
    printf("unique orphan strains %lu\n\n", orphanstrains.size());

    for (auto orphanedstrain = orphanstrains.begin(); orphanedstrain != orphanstrains.end(); orphanedstrain++)
        printf("%s\t%i\n", orphanedstrain->first->id_alt.c_str(), orphanedstrain->second);
    printf("\n\n");
    for (auto orphancontig = orphancontigs.begin(); orphancontig != orphancontigs.end(); orphancontig++)
        printf("%s\t%i\n", orphancontig->first.c_str(), orphancontig->second);
    out.close();
    printf("Disambiguating %i groups took:\t%i seconds\n\n", groups, (int)timer.Elapsed());
}

void ExportStrainsChartData(string OUTPUT_DIRECTORY, Database *db)
{
    Timer timer;
    timer.Start();
    printf("Exporting gene chart data to %s:\n", OUTPUT_DIRECTORY.c_str());
    
    ofstream out( OUTPUT_DIRECTORY + "strains.csv", ifstream::out );
    out << "ID\n";
    
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
    {
        out << (*strain)->id_full << '\n';
        (*strain)->ExportChartData(OUTPUT_DIRECTORY, db);
    }
    out.close();
    printf("Exporting gene chart data took\t%i seconds\n\n", (int)timer.Elapsed());
}

void AnalyzeDiscordantGenes(Dataset *dataset, string INPUT_FILE)
{
    Timer timer;
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("Error: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }
    
    timer.Start();
    string line;
    while (getline(input, line))
    {
        stringstream sstream(line);
        
        if (line.empty())
            continue;
        
        string item;
        vector<string> items;
        while (getline(sstream, item, ';'))
            items.push_back(item);
        
        if (items[0] == "")
            continue;
        
        Group *group = NULL;
        
        for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        {
            if ((*g)->id == items[0])
            {
                group = *g;
                break;
            }
        }
        
        if (group == NULL)
            printf("Error: Group %s not \n", items[0].c_str());
        else
        {
            /*
            scoresumSophisticated += group->syntenyScoreSophisticated;
            scoresumOld += group->syntenyScoreOld;
            scoresumSimple += group->syntenyScoreSimple;
            scoresumAdjusted += group->syntenyScoreAdjusted;
            scoresumFast += group->syntenyScoreFast;
             */
        }
    }

}

//void GenerateSyntenyTree(string OUTPUT_DIRECTORY, Dataset *dataset)
void GenerateSyntenyTree(Dataset *dataset, string OUTPUT_FILE)
{
    Progress p(dataset->database->strains.size()-1);
    Timer timer;
    timer.Start();
    printf("Generating synteny tree to %s:\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "Strain1;Strain2;Distance\n";
    
    int counter = 0;
    for (auto s1 = dataset->database->strains.begin(); s1 != dataset->database->strains.end(); s1++)
    {
        //double genes = 0;
        for (auto s2 = s1+1; s2 != dataset->database->strains.end(); s2++)
        {
            //double genes = 0;
            double distance = 0;
            
            for (auto c1 = (*s1)->contigs.begin(); c1 != (*s1)->contigs.end(); c1++)
            {
                for (auto g1 = (*c1)->genes.begin(); g1 != (*c1)->genes.end(); g1++)
                {
                    //genes++;
                    if ((*g1)->group == NULL)
                        continue;

                    //if ((*g1)->group->CountUniqueStrains() < dataset->database->strains.size())
                    //    continue;

                    for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                        if ((*g2)->strain == *s2)
                        {
                            //genes++;
                            distance += Score::Sophisticated(*g1, *g2)/SIZE_OF_NEIGHTBOURHOOD;
                            //distance += (SIZE_OF_NEIGHTBOURHOOD - Score::Sophisticated(*g1, *g2));
                        }
                }
            }
            
            //out << (*s1)->id_alt << ';' << (*s2)->id_alt << ';' << (distance / genes) << '\n';
            out << (*s1)->id_alt << ';' << (*s2)->id_alt << ';' << distance << '\n';
        }
        //out << (*s1)->id_alt << ';' << (*s1)->id_alt << ';' << 1.0f << '\n';
        double genes = 0;
        for (auto c1 = (*s1)->contigs.begin(); c1 != (*s1)->contigs.end(); c1++)
            for (auto g1 = (*c1)->genes.begin(); g1 != (*c1)->genes.end(); g1++)
                genes++;

        out << (*s1)->id_alt << ';' << (*s1)->id_alt << ';' << genes << '\n';
        p.Update(++counter);
    }
    out.close();
    printf("Generating synteny tree took\t%i seconds\n\n", (int)timer.Elapsed());
}

//void GenerateSyntenyTreeCoreGenes(Dataset *dataset, string OUTPUT_DIRECTORY)
void GenerateSyntenyTreeCoreGenes(string OUTPUT_FILE, Dataset *dataset)
{
    Progress p(dataset->database->strains.size()-1);
    Timer timer;
    timer.Start();
    //printf("Generating synteny tree to %s:\n", OUTPUT_DIRECTORY.c_str());
    printf("Generating synteny tree to %s:\n", OUTPUT_FILE.c_str());

    //ofstream out( OUTPUT_DIRECTORY + "SyntenyTree.csv", ifstream::out );
    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "Strain1;Strain2;Distance\n";

    int counter = 0;
    for (auto s1 = dataset->database->strains.begin(); s1 != dataset->database->strains.end(); s1++)
    {
        for (auto s2 = s1+1; s2 != dataset->database->strains.end(); s2++)
        {
            double genes = 0;
            double distance = 0;
            
            for (auto c1 = (*s1)->contigs.begin(); c1 != (*s1)->contigs.end(); c1++)
            {
                for (auto g1 = (*c1)->genes.begin(); g1 != (*c1)->genes.end(); g1++)
                {
                    genes++;
                    if ((*g1)->group == NULL)
                        continue;

                    if ((*g1)->group->CountUniqueStrains() < dataset->database->strains.size())
                        continue;

                    for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                        if ((*g2)->strain == *s2)
                        {
                            genes++;
                            distance += (SIZE_OF_NEIGHTBOURHOOD - Score::Sophisticated(*g1, *g2));
                        }
                }
            }
            
            out << (*s1)->id_alt << ';' << (*s2)->id_alt << ';' << (distance / genes) << '\n';
        }
        out << (*s1)->id_alt << ';' << (*s1)->id_alt << ';' << 1.0f << '\n';
        p.Update(++counter);
    }
    out.close();
    printf("Generating synteny tree took\t%i seconds\n\n", (int)timer.Elapsed());
}

void ExportIzabelData(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Exporting group data for Izabel into %s:\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    out << "ID;GC3s;Synteny\n";
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        double GC3s = 0;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            GC3s += (*gene)->CalculateGC3s();
        out << (*group)->id << ';' << GC3s / (double)(*group)->genes.size() << ';' << (*group)->syntenyScoreSophisticated << '\n';
    }
    out.close();
    printf("Exporting data for %i groups took:\t%i seconds\n\n", dataset->groups.size(), (int)timer.Elapsed());
}

void ExportStatistics(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Exporting strain statistics into %s:\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );

    out << "Strain;Contigs;bp;ubp;GC;CDS;Orthologs;Paralogs;Uniques;tRNA;rRNA;tmRNA;0-10k;10-250k;250k-500k;500k-750k;750k-1000k;1000-5000k;>5000k;Rh01;Rh02;Rh03;Rh04;Rh05;Rh06;Rh07;Rh08;RhXX\n";

    //int counter = 0;
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        //int CDS = 0;

        uint64_t orthologs = 0;
        uint64_t paralogs = 0;
        uint64_t uniques = 0;
        uint64_t contigsize10k = 0;
        uint64_t contigsize250k = 0;
        uint64_t contigsize500k = 0;
        uint64_t contigsize750k = 0;
        uint64_t contigsize1000k = 0;
        uint64_t contigsize5000k = 0;
        uint64_t contigsize5000kplus = 0;

        bool Rh01 = false;
        bool Rh02 = false;
        bool Rh03 = false;
        bool Rh04 = false;
        bool Rh05 = false;
        bool Rh06 = false;
        bool Rh07 = false;
        bool Rh08 = false;
        bool RhXX = false;

        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            if ((*contig)->length < 0)
            {
                printf("Error: %s had contig with unknown contig size!\n", (*strain)->id_full.c_str());
                exit(1);
            }
                
            if ((*contig)->length < 10000)
                contigsize10k+= (*contig)->length;
                //contigsize10k++;
            else if ((*contig)->length < 250000)
                contigsize250k+= (*contig)->length;
                //contigsize250k++;
            else if ((*contig)->length < 500000)
                contigsize500k+= (*contig)->length;
                //contigsize500k++;
            else if ((*contig)->length < 750000)
                contigsize750k+= (*contig)->length;
               //contigsize750k++;
            else if ((*contig)->length < 1000000)
                contigsize1000k+= (*contig)->length;
                //contigsize1000k++;
            else if ((*contig)->length < 5000000)
                contigsize5000k+= (*contig)->length;
                //contigsize5000k++;
            else
                contigsize5000kplus+= (*contig)->length;
                //contigsize5000kplus++;

            size_t loc;
            if ((loc = (*contig)->id_full.find("Rh01")) != string::npos)
                Rh01 = true;
            if ((loc = (*contig)->id_full.find("Rh02")) != string::npos)
                Rh02 = true;
            if ((loc = (*contig)->id_full.find("Rh03")) != string::npos)
                Rh03 = true;
            if ((loc = (*contig)->id_full.find("Rh04")) != string::npos)
                Rh04 = true;
            if ((loc = (*contig)->id_full.find("Rh05")) != string::npos)
                Rh05 = true;
            if ((loc = (*contig)->id_full.find("Rh06")) != string::npos)
                Rh06 = true;
            if ((loc = (*contig)->id_full.find("Rh07")) != string::npos)
                Rh07 = true;
            if ((loc = (*contig)->id_full.find("Rh08")) != string::npos)
                Rh08 = true;
            if ((loc = (*contig)->id_full.find("Rh1")) != string::npos)
                RhXX = true;

            for (auto g1 = (*contig)->genes.begin(); g1 != (*contig)->genes.end(); g1++)
            {
                if ((*g1)->group == NULL)
                    uniques++;
                else
                {
                    int count = 0;
                    for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                        if ((*g2)->strain == *strain)
                            count++;

                    if (count == 1)
                        orthologs++;
                    else
                        paralogs++;
                }
            }
        }
        
        if ((*strain)->CDS != orthologs + paralogs + uniques)
            printf("Error!");

        out << (*strain)->id_alt << ';';
        out << (*strain)->contigs.size() << ';';
        out << (*strain)->bp << ';';
        out << (*strain)->ubp << ';';
        out << (*strain)->GC << ';';
        out << (*strain)->CDS << ';';
        out << orthologs << ';';
        out << paralogs << ';';
        out << uniques << ';';
        out << (*strain)->tRNA << ';';
        out << (*strain)->rRNA << ';';
        out << (*strain)->tmRNA << ';';

        out << contigsize10k << ';';
        out << contigsize250k << ';';
        out << contigsize500k << ';';
        out << contigsize750k << ';';
        out << contigsize1000k << ';';
        out << contigsize5000k << ';';
        out << contigsize5000kplus << ';';

        out << (Rh01 ? "YES" : "NO") << ';';
        out << (Rh02 ? "YES" : "NO") << ';';
        out << (Rh03 ? "YES" : "NO") << ';';
        out << (Rh04 ? "YES" : "NO") << ';';
        out << (Rh05 ? "YES" : "NO") << ';';
        out << (Rh06 ? "YES" : "NO") << ';';
        out << (Rh07 ? "YES" : "NO") << ';';
        out << (Rh08 ? "YES" : "NO") << ';';
        out << (RhXX ? "YES" : "NO") << '\n';
    }

    out.close();
    printf("Exporting data for %i strains and %i groups took:\t%i seconds\n\n", dataset->database->strains.size(), dataset->groups.size(), (int)timer.Elapsed());
}

void GeneratePresenceAbsenceGAPIT(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Generating presence absence GAPIT data into %s:\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "rs;Allele;chrom;Position;Strand;assembly;center;portLSID;assayLSID;panel;Qccode";
    
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
        out << ';' << (*strain)->id_alt;
    out << '\n';
    
    int position = 100;
    vector<int> positions;
    for (int i = 0; i <= dataset->database->strains.size(); i++)
        positions.push_back(1000);

    int chromosome = 0;
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        set<Strain *> strains;
         for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
             strains.insert((*gene)->strain);

        out << (*group)->id << ';';
        out << "A/G" << ';';
        out << '1' << ';';
        out << position << ';';
        //out << ++chromosome << ';';
        //out << 10000 << ';';
        //out << strains.size() << ';';
        //out << positions[strains.size()] << ';';
        out << "NA;NA;NA;NA;NA;NA;NA";
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
            out << ';' << ((*group)->HasStrain(*strain) ? "AA" : "GG");
        out << '\n';
        positions[strains.size()] += 1000;
        position += 100;
    }

    out.close();
    printf("Generating presence absence GAPIT data took\t%i seconds\n\n", (int)timer.Elapsed());
}

void ExportGeneScores(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Exporting gene scores into %s:\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );

    out << "Gene1;Gene2;Regular Score;Mirror Score\n";

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        for (auto gene1 = (*group)->genes.begin(); gene1 != (*group)->genes.end(); gene1++)
        {
            for (auto score = (*gene1)->scores.begin();  score != (*gene1)->scores.end(); score++ )
            {
                Gene *gene2 = (*score).first;
                out << (*gene1)->id_full << ';';
                out << gene2->id_full << ';';
                out << (*score).second << ';';
                out << (*gene1)->mirrorscores[gene2] << '\n';
            }
        }
    }
    out.close();
    printf("Exporting gene scores took\t%i seconds\n\n", (int)timer.Elapsed());
}

void GenerateDiscordants(Dataset *dataset, string OUTPUT_FILE1, string OUTPUT_FILE2)
{
    Timer timer;
    timer.Start();
    printf("Generating discordant genes into %s and %s:\n", OUTPUT_FILE1.c_str(), OUTPUT_FILE2.c_str());
    
    ofstream out1( OUTPUT_FILE1, ifstream::out );
    ofstream out2( OUTPUT_FILE2, ifstream::out );

    out1 << "Group;Strain;Gene;Contig;Location;Diff;Expected Genospecies;Matched Genospecies\n";
    out2 << "Group;Discordance\n";

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        double discordants = 0;
        Group genospeciesA;
        Group genospeciesB;
        Group genospeciesC;
        Group genospeciesD;
        Group genospeciesE;
        
        genospeciesA.syntenyScoreSophisticated = 0;
        genospeciesB.syntenyScoreSophisticated = 0;
        genospeciesC.syntenyScoreSophisticated = 0;
        genospeciesD.syntenyScoreSophisticated = 0;
        genospeciesE.syntenyScoreSophisticated = 0;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if ((*gene)->strain->species == "A")
                genospeciesA.genes.push_back(*gene);
            if ((*gene)->strain->species == "B")
                genospeciesB.genes.push_back(*gene);
            if ((*gene)->strain->species == "C")
                genospeciesC.genes.push_back(*gene);
            if ((*gene)->strain->species == "D")
                genospeciesD.genes.push_back(*gene);
            if ((*gene)->strain->species == "E")
                genospeciesE.genes.push_back(*gene);
        }

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            Group *bestgroup = NULL;
            double score = 0;
            double scorediff = 0;
            bool discordant = false;

            if (genospeciesA.genes.size() > 0)
            {
                score = genospeciesA.SyntenizeAgainstGene(*gene);
                if (bestgroup == NULL || bestgroup->syntenyScoreSophisticated < score)
                {
                    if (bestgroup != NULL)
                        scorediff = score - bestgroup->syntenyScoreSophisticated;

                    bestgroup = &genospeciesA;
                    bestgroup->syntenyScoreSophisticated = score;
                }
            }
            if (genospeciesB.genes.size() > 0)
            {
                score = genospeciesB.SyntenizeAgainstGene(*gene);
                if (bestgroup == NULL || bestgroup->syntenyScoreSophisticated < score)
                {
                    if (bestgroup != NULL)
                        scorediff = score - bestgroup->syntenyScoreSophisticated;
                    
                    bestgroup = &genospeciesB;
                    bestgroup->syntenyScoreSophisticated = score;
                }
            }
            
            if (genospeciesC.genes.size() > 0)
            {
                score = genospeciesC.SyntenizeAgainstGene(*gene);
                if (bestgroup == NULL || bestgroup->syntenyScoreSophisticated < score)
                {
                    if (bestgroup != NULL)
                        scorediff = score - bestgroup->syntenyScoreSophisticated;
                    
                    bestgroup = &genospeciesC;
                    bestgroup->syntenyScoreSophisticated = score;
                }
            }
            if (genospeciesD.genes.size() > 0)
            {
                score = genospeciesD.SyntenizeAgainstGene(*gene);
                if (bestgroup == NULL || bestgroup->syntenyScoreSophisticated < score)
                {
                    if (bestgroup != NULL)
                        scorediff = score - bestgroup->syntenyScoreSophisticated;
                    
                    bestgroup = &genospeciesD;
                    bestgroup->syntenyScoreSophisticated = score;
                }
            }
            if (genospeciesE.genes.size() > 0)
            {
                score = genospeciesE.SyntenizeAgainstGene(*gene);
                if (bestgroup == NULL || bestgroup->syntenyScoreSophisticated < score)
                {
                    if (bestgroup != NULL)
                        scorediff = score - bestgroup->syntenyScoreSophisticated;
                    
                    bestgroup = &genospeciesE;
                    bestgroup->syntenyScoreSophisticated = score;
                }
            }
            
            if (bestgroup == NULL)
                continue;

            if ((*gene)->strain->species == "A" && genospeciesA.genes.size() > 1 && bestgroup != &genospeciesA)
                discordant = true;
            if ((*gene)->strain->species == "B" && genospeciesB.genes.size() > 1 && bestgroup != &genospeciesB)
                discordant = true;
            if ((*gene)->strain->species == "C" && genospeciesC.genes.size() > 1 && bestgroup != &genospeciesC)
                discordant = true;
            if ((*gene)->strain->species == "D" && genospeciesD.genes.size() > 1 && bestgroup != &genospeciesD)
                discordant = true;
            if ((*gene)->strain->species == "E" && genospeciesE.genes.size() > 1 && bestgroup != &genospeciesE)
                discordant = true;

            if (discordant && scorediff > 1)
            {
                discordants++;

                out1 << (*gene)->group->id << ';';
                out1 << (*gene)->strain->id_full << ';';
                out1 << (*gene)->id_full << ';';
                out1 << (*gene)->contig->id << ';';
                out1 << (*gene)->start + ((*gene)->length / 2) << ';';
                out1 << scorediff << ';';
                out1 << (*gene)->strain->species << ';';

                if (bestgroup == &genospeciesA)
                   out1 << "A\n";
                if (bestgroup == &genospeciesB)
                   out1 << "B\n";
                if (bestgroup == &genospeciesC)
                   out1 << "C\n";
                if (bestgroup == &genospeciesD)
                   out1 << "D\n";
                if (bestgroup == &genospeciesE)
                   out1 << "E\n";
            }
        }

        if (discordants > 0)
            out2 << (*group)->id << ';' << (discordants / (double)(*group)->genes.size()) << '\n';
    }
    out1.close();
    out2.close();
    printf("Generating discordant genes took\t%i seconds\n\n", (int)timer.Elapsed());
}

int main000(int argc, const char * argv[])
{
    Timer timer;
    
    //Convert("/Users/Camous/Desktop/Masters/group_alns_orthomcl/", "/Users/Camous/Desktop/Masters/orthomcl.groups");
    //Convert("/Users/Camous/Desktop/Masters/group_alns_proteinortho/", "/Users/Camous/Desktop/Masters/proteinortho.groups");
    Database *db = new Database();
    Dataset *orthoMCL = new Dataset(db);
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;

    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/spades_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/spades_assemblies/group_fnas/";
    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/spades_assemblies/";

    Parse::Strains(settings->INPUT_DIRECTORY, db);
    /*
    Parse::Groups(settings->GROUP_FILE+"proteinortho.groups", proteinOrtho, false);
    Parse::Groups(settings->GROUP_FILE+"orthomcl.groups", orthoMCL, false);
    Split(orthoMCL, proteinOrtho, settings->GROUP_FILE+"orthomcl");
    Split(proteinOrtho, orthoMCL, settings->GROUP_FILE+"proteinortho");
    */
    /*
    Parse::Groups(settings->GROUP_FILE+"proteinortho.groups", proteinOrtho, true);
    Parse::Groups(settings->GROUP_FILE+"orthomcl.groups", orthoMCL, false);
    proteinOrtho->ScoreSynteny(settings->GROUP_FILE + "proteinortho.scores.csv", orthoMCL);
    Histogram(proteinOrtho, orthoMCL, settings->GROUP_FILE + "proteinortho.histogram.txt");
    */
    Parse::Groups(settings->GROUP_FILE+"proteinortho.groups", proteinOrtho, false);
    Parse::Groups(settings->GROUP_FILE+"orthomcl.groups", orthoMCL, true);
    orthoMCL->ScoreSynteny(settings->GROUP_FILE + "orthomcl.scores.csv", proteinOrtho);
    Histogram(orthoMCL, proteinOrtho, settings->GROUP_FILE + "orthomcl.histogram.txt");

    exit(0);

    //Parse::Strains(settings->INPUT_DIRECTORY, db);
    //Parse::Dendro("/Users/Camous/Downloads/Dendro_discordant_genes/group62.dendro", db);
    //exit(0);
    
    /*
     Parse::ProteinOrtho(settings->GROUP_FILE, proteinOrtho, true);
     
     //Disambiguate(proteinOrtho, settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups");
     
     proteinOrtho->ScoreSynteny(settings->GROUP_FILE + "_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores", NULL);
     exit(0);
     */
    //ExportGeneRelationships(settings->INPUT_DIRECTORY, db);
    
    Dataset *proteinOrtho_disambiguated = new Dataset(db);
    Parse::Groups( settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups", proteinOrtho_disambiguated, true);
    proteinOrtho_disambiguated->ExportFNA(settings->OUTPUT_DIRECTORY);
    //proteinOrtho_disambiguated->ScoreSynteny(settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores", NULL);
    
    //Histogram(proteinOrtho, proteinOrtho_disambiguated, settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram");
    
    exit(0);
    Parse::Strains("/Users/Camous/Desktop/Masters/new_assemblies/data/", db);
    Parse::Genospecies("/Users/Camous/Desktop/Masters/Heritability/growth.csv", db);
    Parse::Sequences("/Users/Camous/Desktop/Masters/new_assemblies/data/", db);
    
    Parse::ProteinOrtho("/Users/Camous/Desktop/Masters/new_assemblies/myproject.proteinortho", proteinOrtho, true);
    proteinOrtho->ExportFNA("/Users/Camous/Desktop/Masters/new_assemblies/group_fnas/");
    
    ExportStrainsChartData("/Users/Camous/Desktop/Masters/new_assemblies/generelationship/", db);
    
    exit(0);
    Parse::Strains("/Users/Camous/Desktop/Masters/split/gff/", db);
    Parse::Genospecies("/Users/Camous/Desktop/Masters/Heritability/growth.csv", db);
    
    //ExportPTT(db, "/Users/Camous/Desktop/Masters/ptt/");
    //exit(0);
    
    //RearrangeContigs(proteinOrtho, "/Users/Camous/Desktop/Masters/scaffolds/", true);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho_disambiguated_400.groups", proteinOrtho, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
    
    proteinOrtho->grouppool2["group9085"]->GenerateSyntenyMap("/Users/Camous/Desktop/Masters/Syntenymap/");
    exit(0);
    
    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    
    //ExportGeneRelationships("/Users/Camous/Desktop/Masters/Generelationship/", db);
    //db->strainpool["3396"]->PrintContigs();
    //db->strainpool["3396"]->ExportGeneRelationship("/Users/Camous/Desktop/Masters/Generelationship/");
    //db->strainpool["3328"]->PrintContigs();
    //db->strainpool["3328"]->ExportGeneRelationship("/Users/Camous/Desktop/Masters/Generelationship/");
    //TestSplit(proteinOrtho);
    
    //Disambiguate(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups");
    
    //proteinOrtho->grouppool2["group9085"]->Split3();
    //proteinOrtho->grouppool2["group9085"]->Split4();
    
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.csv", true);
    
    //GetGeneGroupScores("3322|07082", proteinOrtho);
    //GetGeneGroupScores("3320|02815", proteinOrtho);
    
    //proteinOrtho->database->genepool["3349|06260"]->CompareNeighbours(proteinOrtho->database->genepool["3349|06825"]);
    //proteinOrtho->database->genepool["3347|03984"]->CompareNeighbours(proteinOrtho->database->genepool["3349|06825"]);
    //GetGeneGroupScores("3306|03637", proteinOrtho);
    //proteinOrtho->database->genepool["3341|02842"]->CompareNeighbours(proteinOrtho->database->genepool["3336|04312"]);
    
    //AnalyzeDiscordantGroups(proteinOrtho, "group2390;group2447;group2468;group9085;group9087;group9088;group9141");
    
    //group2390;group2447;group2468;group9085;group9087;group9088;group9141
    
    exit(0);
    /*
     Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
     Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
     
     //GenerateGroupmap(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_to_orthomcl.groupmap");
     //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
     //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
     Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
     //SyntenyScoreSubset(proteinOrtho, "/Users/Camous/Desktop/Masters/split/discordant_genes.csv");
     
     AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups");
     exit(0);
     //MatchGenetoGroup("3207|04612", "group9519", proteinOrtho);
     //MatchGenetoGroup("3287|02586", "group9373", proteinOrtho);
     //MatchGenetoGroup("3242|01736", "group9532", proteinOrtho);
     //MatchGenetoGroup("3290|07501", "group13980", proteinOrtho);
     
     //GenerateSyntenyMap(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.map" );
     
     GenerateSyntenyMaps(proteinOrtho, "/Users/Camous/Desktop/Masters/syntenymaps.proteinortho/" );
     */
    
    //
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, true);
    
    GenerateGroupmap(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl_to_proteinortho.groupmap" );
    //orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho );
    Parse::Scores("/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    
    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt" );
    
    //GenerateSyntenyMap(orthoMCL, "/Users/Camous/Desktop/Masters/split/orthomcl.map" );
    GenerateSyntenyMaps(orthoMCL, "/Users/Camous/Desktop/Masters/syntenymaps.orthomcl/" );
    
    //
    exit(0);
    
    //AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho1to1.groups");
    //Test(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups", true);
    
    //GenerateGeneticRelationMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/GRM/proteinortho.grm.csv", true);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/GRM/proteinortho.presence_absence_matrix.csv", false);
    //AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups");
    
    //CompareContigs(0, proteinOrtho);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, true);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/GRM/orthomcl.presence_absence_matrix.csv", false);
    //AnalyzeParalogs(orthoMCL, "/Users/Camous/Desktop/Masters/split/orthomcl.groups");
    
    //CompareContigs(2, orthoMCL);
    
    //proteinOrtho->database->genepool["3328|00198"]->CompareNeighbours(proteinOrtho->database->genepool["3328|00626"]);
    
    //proteinOrtho->database->genepool["3384|00807"]->CompareNeighbours(proteinOrtho->database->genepool["3384|04976"]);
    
    //3384|04976
    //219078-219839
    //698440-699201
    
    exit(0);
    /*
     Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
     Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
     Split(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl");
     Split(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho");
     Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl1to1.groups", orthoMCL, false);
     Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho1to1.groups", proteinOrtho, false);
     exit(0);
     */
    
    //Parse::Strains("/Users/Camous/Desktop/Masters/asger/rhizobium/", db);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, true);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_alt.groups");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/groups.txt", orthoMCL, true);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_alt2.groups");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, true);
    //ExportGroups(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho_alt.groups");
    /*
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl_original.groups", orthoMCL, true);
     ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_original_alt.groups");
     exit(0);
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl_original.groups", orthoMCL, true);
     Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, false);
     orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_original.scores.txt", proteinOrtho);
     Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_original.histogram.txt");
     
     exit(0);
     */
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl.presence_absence_matrix.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.txt");
    //orthoMCL->database->genepool["3206|07195"]->CompareNeighbours(orthoMCL->database->genepool["3207|06023"]);
    /*
     Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
     
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to1.groups.txt", orthoMCL, true);
     orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to1.scores.txt", NULL);
     Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl1to1.histogram.txt");
     
     
     Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho1to1.groups.txt", proteinOrtho, true);
     proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to1.scores.txt", NULL);
     Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho1to1.histogram.txt");
     */
    //exit(0);
    //Test();
    //Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl.presence_absence_matrix.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.txt");
    //orthoMCL->database->genepool["3206|07195"]->CompareNeighbours(orthoMCL->database->genepool["3207|06023"]);
    /*
     Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
     
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to1.groups.txt", orthoMCL, true);
     orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to1.scores.txt", NULL);
     Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl1to1.histogram.txt");
     
     
     Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho1to1.groups.txt", proteinOrtho, true);
     proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to1.scores.txt", NULL);
     Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho1to1.histogram.txt");
     */
    //exit(0);
    
    //Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/asger/rhizobium/orthomcl.groups.txt", orthoMCL, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/asger/rhizobium/orthomcl.groups");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/asger/rhizobium/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //ExportGroups(proteinOrtho, "/Users/Camous/Desktop/Masters/asger/rhizobium/proteinOrtho.groups");
    //exit(0);
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl.scores.txt", orthoMCL);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, false);
    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl.histogram.txt");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, false);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/proteinortho.scores.txt", proteinOrtho);
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho.histogram.txt");
    /**/
    /**/
    /*
     int genes = 0;
     int groups = 0;
     
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, true);
     Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
     MatchGenestoGroups(orthoMCL->groups[45]->genes, orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".matches_G45.txt");
     MatchGenestoGroups(orthoMCL->groups[0]->genes, orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".matches_G0.txt");
     //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, false);
     //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
     exit(0);
     */
    /*
     Dataset *orthoMCLFiltered = new Dataset(db);
     for (auto group = orthoMCL->groups.begin(); group != orthoMCL->groups.end(); group++)
     {
     if ((*group)->syntenyScoreSophisticated == SIZE_OF_NEIGHTBOURHOOD || (*group)->syntenyScoreSimple == SIZE_OF_NEIGHTBOURHOOD)
     {
     groups++;
     genes += (*group)->genes.size();
     orthoMCLFiltered->groups.push_back(*group);
     }
     }
     Histogram(orthoMCLFiltered, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + "_filtered.histogram.txt");
     */
    /*
     Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl_original.groups", orthoMCL, true);
     Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
     //Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
     orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
     Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
     */
    /**/
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
}

int mainTemp(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;
    
    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/final_assemblies/final_project.poff.groups";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/final_assemblies/data/";
    //settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/final_assemblies/GWAScharts.maf5/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/final_assemblies/GWAScharts.maf1/";

    Timer timer;
    timer.Start();
    
    printf("Attempting to parse strains:\n");
    Parse::Strains(settings->INPUT_DIRECTORY, db);
    printf("Done parsing strains.\n\n");
    
    printf("Attempting to parse genospecies:\n");
    Parse::Genospecies(settings->INPUT_DIRECTORY + "metadata.csv", db);
    printf("Done parsing genospecies.\n\n");
    
    printf("Attempting to parse sequences:\n");
    Parse::Sequences(settings->INPUT_DIRECTORY, db);
    printf("Done parsing sequences.\n\n");
    
    printf("Attempting to parse ProteinOrtho groups:\n");
    Parse::Groups(settings->GROUP_FILE, proteinOrtho);
    printf("Done parsing ProteinOrtho groups.\n\n");
    
    printf("Attempting to parse group types:\n");
    //Parse::GroupTypes("/Users/Camous/Desktop/Masters/final_assemblies/GWAS.types.maf5.csv", proteinOrtho);
    Parse::GroupTypes("/Users/Camous/Desktop/Masters/final_assemblies/GWAS.types.maf1.csv", proteinOrtho);
    printf("Done parsing group types.\n\n");
    
    printf("Attempting to export charts data:\n");
    ExportStrainsChartData(settings->OUTPUT_DIRECTORY, db);
    printf("Done exporting chart data.\n\n");

    return 0;
}

int main(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;
    /*
    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/newest_assemblies/newest_project.poff.groups";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/newest_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/newest_assemblies/charts/";
     */

    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/final_assemblies/final_project.poff.groups";
    //settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/final_assemblies/final_project.poff_disambiguated.groups";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/final_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/final_assemblies/charts/";

    /*
    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/spades_assemblies/myproject.poff.groups";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/spades_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/spades_assemblies/charts/";
    */
    Timer timer;
    timer.Start();

    printf("Attempting to parse strains:\n");
    Parse::Strains(settings->INPUT_DIRECTORY, db);
    printf("Done parsing strains.\n\n");
    
    printf("Attempting to parse genospecies:\n");
    Parse::Genospecies(settings->INPUT_DIRECTORY + "metadata.csv", db);
    printf("Done parsing genospecies.\n\n");

    printf("Attempting to parse sequences:\n");
    Parse::Sequences(settings->INPUT_DIRECTORY, db);
    printf("Done parsing sequences.\n\n");

    printf("Attempting to parse ProteinOrtho groups:\n");
    Parse::Groups(settings->GROUP_FILE, proteinOrtho);
    printf("Done parsing ProteinOrtho groups.\n\n");

    printf("Attempting to parse group types:\n");
    Parse::GroupTypes(settings->GROUP_FILE + ".types.csv", proteinOrtho);
    printf("Done parsing group types.\n\n");

    printf("Attempting to export charts data:\n");
    ExportStrainsChartData(settings->OUTPUT_DIRECTORY, db);
    printf("Done exporting chart data.\n\n");

    printf("Attempting to generate discordants:\n");
    GenerateDiscordants(proteinOrtho, settings->GROUP_FILE + ".discordantgenes.csv", settings->GROUP_FILE + ".discordantgroups.csv");
    printf("Done generating discordants.\n\n");

    //exit(0);

    //for (auto group = proteinOrtho->groups.begin(); group != proteinOrtho->groups.end(); group++)
    //    (*group)->GenerateSyntenyMap(settings->OUTPUT_DIRECTORY);
    
    //printf("Attempting to parse synteny scores:\n");
    //Parse::Scores("/Users/Camous/Desktop/Masters/new_assemblies/myproject.poff.scores.csv", proteinOrtho);
    //printf("Done parsing synteny scores.\n\n");
    printf("Attempting to score Synteny:\n");
    proteinOrtho->ScoreSynteny(settings->GROUP_FILE + ".scores.csv", NULL);
    printf("Done scoring Synteny.\n\n");

    /*
    printf("Attempting export gene scores:\n");
    ExportGeneScores(proteinOrtho, settings->GROUP_FILE + ".genescores.csv");
    printf("Done exporting gene scores.\n\n");
    */
    //exit(0);
    printf("Attempting to generate Presence Absence Matrix:\n");
    GeneratePresenceAbsenceMatrix(proteinOrtho, settings->GROUP_FILE + ".presence_absence_matrix.csv", true);
    printf("Done generating Presence Absence Matrix.\n\n");

    printf("Attempting to generate GRM:\n");
    GenerateGeneticRelationMatrix(proteinOrtho, settings->GROUP_FILE + ".grm.csv", true);
    printf("Done generating GRM.\n\n");
    
    printf("Attempting to generate Synteny Tree:\n");
    GenerateSyntenyTree(proteinOrtho, settings->GROUP_FILE + ".syntenytree.csv");
    printf("Done generating Synteny Tree.\n\n");

    printf("Attempting to generate data for GAPIT:\n");
    GeneratePresenceAbsenceGAPIT(proteinOrtho, settings->GROUP_FILE + ".gapit.csv");
    printf("Done generating data for GAPIT.\n\n");

    printf("Attempting to export statistics:\n");
    ExportStatistics(proteinOrtho, settings->GROUP_FILE + ".statistics.csv");
    printf("Done exporting statistics.\n\n");

    /*
    printf("Attempting to parse discordants:\n");
    Parse::Discordants(settings->GROUP_FILE + ".discordants.csv", proteinOrtho);
    printf("Done parsing discordants scores.\n\n");
    */
    /*
    printf("Attempting to export charts data:\n");
    ExportStrainsChartData(settings->OUTPUT_DIRECTORY, db);
    printf("Done exporting chart data.\n\n");
     */
    printf("Attempting to generate synteny maps:\n");
    GenerateSyntenyMaps(proteinOrtho, settings->OUTPUT_DIRECTORY );
    printf("Done generating synteny maps.\n\n");

    printf("Everything took \t%i seconds\n\n", (int)timer.Elapsed());

    return 0;
}

int main0()
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;

    Parse::Strains("/Users/Camous/Desktop/Masters/split/gff/", db);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho_disambiguated_400.groups", proteinOrtho, true);
    proteinOrtho->grouppool2["group9085"]->GenerateSyntenyMap("/Users/Camous/Desktop/Masters/Syntenymap/");
    return 0;
}

int main4(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;
    
    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/new_assemblies/myproject.poff";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/";
    
    Timer timer;
    timer.Start();
    
    printf("Attempting to parse strains:\n");
    Parse::Strains(settings->INPUT_DIRECTORY, db);
    printf("Done parsing strains.\n\n");
    
    printf("Attempting to parse ProteinOrtho groups:\n");
    Parse::ProteinOrtho(settings->GROUP_FILE, proteinOrtho, true);
    printf("Done parsing ProteinOrtho groups.\n\n");
    
    printf("Attempting to generate Synteny Tree:\n");
    GenerateSyntenyTree(proteinOrtho, settings->GROUP_FILE + "syntenytree.csv");
    printf("Done generating Synteny Tree.\n\n");

    printf("Everything took \t%i seconds\n\n", (int)timer.Elapsed());
    
    return 0;
}

int main3(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;

    settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/new_assemblies/myproject.poff";
    settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/data/";
    settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/charts/";
    
    Timer timer;
    
    timer.Start();

    printf("Attempting to parse strains:\n");
    Parse::Strains(settings->INPUT_DIRECTORY, db);
    printf("Done parsing strains.\n\n");

    printf("Attempting to parse sequences:\n");
    Parse::Sequences(settings->INPUT_DIRECTORY, db);
    printf("Done parsing sequences.\n\n");

    printf("Attempting to parse coverages:\n");
    Parse::Coverage(settings->INPUT_DIRECTORY, db);
    printf("Done parsing coverages.\n\n");

    printf("Attempting to parse ProteinOrtho groups:\n");
    Parse::ProteinOrtho(settings->GROUP_FILE, proteinOrtho, true);
    printf("Done parsing ProteinOrtho groups.\n\n");

    printf("Attempting to generate GRM:\n");
    GenerateGeneticRelationMatrix(proteinOrtho, settings->GROUP_FILE + ".grm", true);
    printf("Done generating GRM.\n\n");
    
    //exit(0);

    printf("Attempting to parse synteny scores:\n");
    Parse::Scores(settings->GROUP_FILE + "_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores", proteinOrtho);
    printf("Done parsing synteny scores.\n\n");

    printf("Attempting to export charts data:\n");
    ExportStrainsChartData(settings->OUTPUT_DIRECTORY, db);
    printf("Done exporting chart data.\n\n");
    
    printf("Everything took \t%i seconds\n\n", (int)timer.Elapsed());

    return 0;
}

int main2(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *proteinOrtho = new Dataset(db);
    Dataset *proteinOrtho_disambiguated = new Dataset(db);
    Settings *settings =  new Settings;
    
    if (argc <= 3)
    {
        printf("Usage example:\nsyntenizer3000 /new_assemblies/myproject.poff /new_assemblies/data/ /new_assemblies/group_fnas/\n");
        exit(1);
    }
    else
    {
        settings->GROUP_FILE = argv[1];
        settings->INPUT_DIRECTORY = argv[2];
        settings->OUTPUT_DIRECTORY = argv[3];
        //printf("argv[1] %s\n", argv[1]);
        //printf("argv[2] %s\n", argv[2]);
    }
    if (settings->INPUT_DIRECTORY[settings->INPUT_DIRECTORY.size()-1] != '/')
        settings->INPUT_DIRECTORY =+ '/';
    if (settings->OUTPUT_DIRECTORY[settings->OUTPUT_DIRECTORY.size()-1] != '/')
        settings->OUTPUT_DIRECTORY =+ '/';
    
    printf("Attempting to parse strains:\n");
    Parse::Strains(settings->INPUT_DIRECTORY, db);
    printf("Done parsing strains.\n\n");

    printf("Attempting to parse sequences:\n");
    Parse::Sequences(settings->INPUT_DIRECTORY, db);
    printf("Done parsing sequences.\n\n");
    
    printf("Attempting to parse ProteinOrtho groups:\n");
    Parse::ProteinOrtho(settings->GROUP_FILE, proteinOrtho, true);
    printf("Done parsing ProteinOrtho groups.\n\n");
    
    printf("Attempting to export groups:\n");
    ExportGroups(proteinOrtho, settings->GROUP_FILE + ".groups");
    printf("Done exporting groups.\n\n");
    
    printf("Attempting to disambiguate groups:\n");
    Disambiguate(proteinOrtho, settings->GROUP_FILE + "_disambiguated.groups");
    printf("Done disambiguating groups.\n\n");
    
    printf("Attempting to parse disambiguated groups:\n");
    Parse::Groups( settings->GROUP_FILE + "_disambiguated.groups", proteinOrtho_disambiguated, true);
    printf("Done parsing disambiguated groups.\n\n");

    printf("Attempting to export fna groups:\n");
    proteinOrtho_disambiguated->ExportFNA(settings->OUTPUT_DIRECTORY);
    printf("Done exporting fna groups.\n\n");

    return 0;
}

int main1(int argc, const char * argv[])
{
    Timer timer;

    //Convert("/Users/Camous/Desktop/Masters/group_alns_orthomcl/", "/Users/Camous/Desktop/Masters/orthomcl.groups");
    //Convert("/Users/Camous/Desktop/Masters/group_alns_proteinortho/", "/Users/Camous/Desktop/Masters/proteinortho.groups");
    Database *db = new Database();
    Dataset *orthoMCL = new Dataset(db);
    Dataset *proteinOrtho = new Dataset(db);
    Settings *settings =  new Settings;

    if (argc <= 3)
    {
        printf("Usage example:\nsyntenizer3000 /new_assemblies/myproject.poff /new_assemblies/data/ /new_assemblies/group_fnas/\n");
        //exit(1);
        //settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/split/gff/";

        settings->INPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/data/";
        settings->OUTPUT_DIRECTORY = "/Users/Camous/Desktop/Masters/new_assemblies/group_fnas/";
        settings->GROUP_FILE = "/Users/Camous/Desktop/Masters/new_assemblies/myproject.poff";
    }
    else
    {
        settings->GROUP_FILE = argv[1];
        settings->INPUT_DIRECTORY = argv[2];
        //printf("argv[1] %s\n", argv[1]);
        //printf("argv[2] %s\n", argv[2]);
    }
    if (settings->INPUT_DIRECTORY[settings->INPUT_DIRECTORY.size()-1] != '/')
        settings->INPUT_DIRECTORY =+ '/';
    if (settings->OUTPUT_DIRECTORY[settings->OUTPUT_DIRECTORY.size()-1] != '/')
        settings->OUTPUT_DIRECTORY =+ '/';
    
    //Parse::Strains(settings->INPUT_DIRECTORY, db);
    //Parse::Dendro("/Users/Camous/Downloads/Dendro_discordant_genes/group62.dendro", db);
    //exit(0);

    Parse::Strains(settings->INPUT_DIRECTORY, db);
    Parse::Sequences(settings->INPUT_DIRECTORY, db);
    /*
    Parse::ProteinOrtho(settings->GROUP_FILE, proteinOrtho, true);

    //ExportGroups(proteinOrtho, settings->GROUP_FILE + ".groups");
    //Disambiguate(proteinOrtho, settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups");
    
    proteinOrtho->ScoreSynteny(settings->GROUP_FILE + "_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores", NULL);
    exit(0);
    */
    //ExportGeneRelationships(settings->INPUT_DIRECTORY, db);

    Dataset *proteinOrtho_disambiguated = new Dataset(db);
    Parse::Groups( settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups", proteinOrtho_disambiguated, true);
    proteinOrtho_disambiguated->ExportFNA(settings->OUTPUT_DIRECTORY);
    //proteinOrtho_disambiguated->ScoreSynteny(settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores", NULL);

    //Histogram(proteinOrtho, proteinOrtho_disambiguated, settings->GROUP_FILE + "_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram");

    exit(0);
    Parse::Strains("/Users/Camous/Desktop/Masters/new_assemblies/data/", db);
    Parse::Genospecies("/Users/Camous/Desktop/Masters/Heritability/growth.csv", db);
    Parse::Sequences("/Users/Camous/Desktop/Masters/new_assemblies/data/", db);

    Parse::ProteinOrtho("/Users/Camous/Desktop/Masters/new_assemblies/myproject.proteinortho", proteinOrtho, true);
    proteinOrtho->ExportFNA("/Users/Camous/Desktop/Masters/new_assemblies/group_fnas/");

    ExportStrainsChartData("/Users/Camous/Desktop/Masters/new_assemblies/generelationship/", db);

    exit(0);
    Parse::Strains("/Users/Camous/Desktop/Masters/split/gff/", db);
    Parse::Genospecies("/Users/Camous/Desktop/Masters/Heritability/growth.csv", db);

    //ExportPTT(db, "/Users/Camous/Desktop/Masters/ptt/");
    //exit(0);

    //RearrangeContigs(proteinOrtho, "/Users/Camous/Desktop/Masters/scaffolds/", true);

    //Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho_disambiguated_400.groups", proteinOrtho, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);

    proteinOrtho->grouppool2["group9085"]->GenerateSyntenyMap("/Users/Camous/Desktop/Masters/Syntenymap/");
    exit(0);

    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");

    //ExportGeneRelationships("/Users/Camous/Desktop/Masters/Generelationship/", db);
    //db->strainpool["3396"]->PrintContigs();
    //db->strainpool["3396"]->ExportGeneRelationship("/Users/Camous/Desktop/Masters/Generelationship/");
    //db->strainpool["3328"]->PrintContigs();
    //db->strainpool["3328"]->ExportGeneRelationship("/Users/Camous/Desktop/Masters/Generelationship/");
    //TestSplit(proteinOrtho);
    
    //Disambiguate(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho_disambiguated_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".groups");

    //proteinOrtho->grouppool2["group9085"]->Split3();
    //proteinOrtho->grouppool2["group9085"]->Split4();

    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.csv", true);

    //GetGeneGroupScores("3322|07082", proteinOrtho);
    //GetGeneGroupScores("3320|02815", proteinOrtho);
    
    //proteinOrtho->database->genepool["3349|06260"]->CompareNeighbours(proteinOrtho->database->genepool["3349|06825"]);
    //proteinOrtho->database->genepool["3347|03984"]->CompareNeighbours(proteinOrtho->database->genepool["3349|06825"]);
    //GetGeneGroupScores("3306|03637", proteinOrtho);
    //proteinOrtho->database->genepool["3341|02842"]->CompareNeighbours(proteinOrtho->database->genepool["3336|04312"]);

    //AnalyzeDiscordantGroups(proteinOrtho, "group2390;group2447;group2468;group9085;group9087;group9088;group9141");

    //group2390;group2447;group2468;group9085;group9087;group9088;group9141

    exit(0);
    /*
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);

    //GenerateGroupmap(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_to_orthomcl.groupmap");
    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    //SyntenyScoreSubset(proteinOrtho, "/Users/Camous/Desktop/Masters/split/discordant_genes.csv");

    AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups");
    exit(0);
    //MatchGenetoGroup("3207|04612", "group9519", proteinOrtho);
    //MatchGenetoGroup("3287|02586", "group9373", proteinOrtho);
    //MatchGenetoGroup("3242|01736", "group9532", proteinOrtho);
    //MatchGenetoGroup("3290|07501", "group13980", proteinOrtho);

    //GenerateSyntenyMap(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.map" );
    
    GenerateSyntenyMaps(proteinOrtho, "/Users/Camous/Desktop/Masters/syntenymaps.proteinortho/" );
     */

    //
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, true);

    GenerateGroupmap(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl_to_proteinortho.groupmap" );
    //orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho );
    Parse::Scores("/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);

    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl_rearranged_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt" );

    //GenerateSyntenyMap(orthoMCL, "/Users/Camous/Desktop/Masters/split/orthomcl.map" );
    GenerateSyntenyMaps(orthoMCL, "/Users/Camous/Desktop/Masters/syntenymaps.orthomcl/" );

    //
    exit(0);
    
    //AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho1to1.groups");
    //Test(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups", true);

    //GenerateGeneticRelationMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/GRM/proteinortho.grm.csv", true);

    //Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/GRM/proteinortho.presence_absence_matrix.csv", false);
    //AnalyzeParalogs(proteinOrtho, "/Users/Camous/Desktop/Masters/split/proteinortho.groups");

    //CompareContigs(0, proteinOrtho);
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, true);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/GRM/orthomcl.presence_absence_matrix.csv", false);
    //AnalyzeParalogs(orthoMCL, "/Users/Camous/Desktop/Masters/split/orthomcl.groups");

    //CompareContigs(2, orthoMCL);

    //proteinOrtho->database->genepool["3328|00198"]->CompareNeighbours(proteinOrtho->database->genepool["3328|00626"]);
    
    //proteinOrtho->database->genepool["3384|00807"]->CompareNeighbours(proteinOrtho->database->genepool["3384|04976"]);
    
    //3384|04976
    //219078-219839
    //698440-699201

    exit(0);
    /*
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
    Split(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/split/orthomcl");
    Split(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho");
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl1to1.groups", orthoMCL, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho1to1.groups", proteinOrtho, false);
    exit(0);
     */

    //Parse::Strains("/Users/Camous/Desktop/Masters/asger/rhizobium/", db);

    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, true);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_alt.groups");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/groups.txt", orthoMCL, true);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_alt2.groups");

    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, true);
    //ExportGroups(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho_alt.groups");
    /*
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl_original.groups", orthoMCL, true);
    ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_original_alt.groups");
    exit(0);
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl_original.groups", orthoMCL, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, false);
    orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_original.scores.txt", proteinOrtho);
    Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_original.histogram.txt");

    exit(0);
     */
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl.presence_absence_matrix.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.txt");
    //orthoMCL->database->genepool["3206|07195"]->CompareNeighbours(orthoMCL->database->genepool["3207|06023"]);
    /*
     Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
     
     Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to1.groups.txt", orthoMCL, true);
     orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to1.scores.txt", NULL);
     Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl1to1.histogram.txt");
     
     
     Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho1to1.groups.txt", proteinOrtho, true);
     proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to1.scores.txt", NULL);
     Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho1to1.histogram.txt");
     */
    //exit(0);
    //Test();
    //Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl.presence_absence_matrix.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //GeneratePresenceAbsenceMatrix(proteinOrtho, "/Users/Camous/Desktop/Masters/proteinortho.presence_absence_matrix.txt");
    //orthoMCL->database->genepool["3206|07195"]->CompareNeighbours(orthoMCL->database->genepool["3207|06023"]);
    /*
    Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
    
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to1.groups.txt", orthoMCL, true);
    orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to1.scores.txt", NULL);
    Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl1to1.histogram.txt");
    

    Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho1to1.groups.txt", proteinOrtho, true);
    proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to1.scores.txt", NULL);
    Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho1to1.histogram.txt");
    */
    //exit(0);

    //Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);

    //Parse::Groups("/Users/Camous/Desktop/Masters/asger/rhizobium/orthomcl.groups.txt", orthoMCL, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //ExportGroups(orthoMCL, "/Users/Camous/Desktop/Masters/asger/rhizobium/orthomcl.groups");

    //Parse::Groups("/Users/Camous/Desktop/Masters/asger/rhizobium/proteinortho.groups.txt", proteinOrtho, true);
    //ConvertStrainNames("/Users/Camous/Desktop/Masters/gff/", db);
    //ExportGroups(proteinOrtho, "/Users/Camous/Desktop/Masters/asger/rhizobium/proteinOrtho.groups");
    //exit(0);
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl.scores.txt", orthoMCL);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, false);
    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl.histogram.txt");

    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL, false);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", proteinOrtho, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/proteinortho.scores.txt", proteinOrtho);
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho.histogram.txt");
    /**/
    /**/
    /*
    int genes = 0;
    int groups = 0;

    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, true);
    Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    MatchGenestoGroups(orthoMCL->groups[45]->genes, orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".matches_G45.txt");
    MatchGenestoGroups(orthoMCL->groups[0]->genes, orthoMCL, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".matches_G0.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, false);
    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    exit(0);
    */
    /*
    Dataset *orthoMCLFiltered = new Dataset(db);
    for (auto group = orthoMCL->groups.begin(); group != orthoMCL->groups.end(); group++)
    {
        if ((*group)->syntenyScoreSophisticated == SIZE_OF_NEIGHTBOURHOOD || (*group)->syntenyScoreSimple == SIZE_OF_NEIGHTBOURHOOD)
        {
            groups++;
            genes += (*group)->genes.size();
            orthoMCLFiltered->groups.push_back(*group);
        }
    }
    Histogram(orthoMCLFiltered, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + "_filtered.histogram.txt");
    */
    /*
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl_original.groups", orthoMCL, true);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, false);
    //Parse::Scores("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_new_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    */
    /**/
    Parse::Groups("/Users/Camous/Desktop/Masters/split/orthomcl.groups", orthoMCL, false);
    Parse::Groups("/Users/Camous/Desktop/Masters/split/proteinortho.groups", proteinOrtho, true);
    //Parse::Scores("/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/split/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    /**/
     /*
    Dataset *proteinOrthoFiltered = new Dataset(db);
    for (auto group = proteinOrtho->groups.begin(); group != proteinOrtho->groups.end(); group++)
    {
        if ((*group)->syntenyScoreSophisticated == SIZE_OF_NEIGHTBOURHOOD || (*group)->syntenyScoreSimple == SIZE_OF_NEIGHTBOURHOOD)
        {
            groups++;
            genes += (*group)->genes.size();
            proteinOrthoFiltered->groups.push_back(*group);
        
        }
    }
    Histogram(proteinOrthoFiltered, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + "_filtered.histogram.txt");
    */
    /*
    printf("groups: %i\n", groups);
    printf("genes: %i\n", genes);
    printf("genes per group: %f\n", (double)genes / (double)groups);
     
    exit(0);
    */
    //Parse::Strains("/Users/Camous/Desktop/Masters/gff/", db);
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, true);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, false);
    //orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", proteinOrtho);
    //Histogram(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, false);
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups", proteinOrtho, true);
    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
//
    /**/
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", orthoMCL, false);
    //Parse::Groups("/Users/Camous/Desktop/Masters/asger/rhizobium/proteinortho.groups_disambiguated", proteinOrtho, true);
    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho_disambiguated2_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".scores.txt", orthoMCL);
    //Histogram(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho_disambiguated2_" + to_string(SIZE_OF_NEIGHTBOURHOOD) + ".histogram.txt");
    /**/
//
    //Split(orthoMCL, proteinOrtho, "/Users/Camous/Desktop/Masters/orthomcl");
    //Split(proteinOrtho, orthoMCL, "/Users/Camous/Desktop/Masters/proteinortho");
    /*
    Dataset *orthoMCL1to1 = new Dataset(db);
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to1.groups.txt", orthoMCL1to1);
    orthoMCL1to1->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to1.scores.txt");
    Dataset *orthoMCL1to2 = new Dataset(db);
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl1to!1.groups.txt", orthoMCL1to2);
    orthoMCL1to2->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl1to!1.scores.txt");

    Dataset *proteinOrtho1to1 = new Dataset(db);
    Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho1to1.groups.txt", proteinOrtho1to1);
    proteinOrtho1to1->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to1.scores.txt");
    Dataset *proteinOrtho1to2 = new Dataset(db);
    Parse::Groups("/Users/Camous/Desktop/Masters/proteinOrtho1to!1.groups.txt", proteinOrtho1to2);
    proteinOrtho1to2->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho1to!1.scores.txt");
    */
    /*
    Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups", dataset);
    
    ofstream out( "/Users/Camous/Desktop/Masters/orthomcl.scores", ifstream::out );
    double scoresum = 0;
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        scoresum += (*group)->Syntenize();
        out << (*group)->id.c_str() << ": Genes: " << (*group)->genes.size() << " Synteny Score: " << (*group)->syntenyScore << "\n";
    }
    printf("Average Synteny Score:\t%.2f\n", scoresum / dataset->groups.size());
    out << "Average Synteny Score: " << scoresum / dataset->groups.size();
    out.close();
     */
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl.groups.txt", orthoMCL);
    //orthoMCL->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl.scores.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho.groups.txt", dataset);
    //dataset->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho.scores");
    //Parse::Groups("/Users/Camous/Desktop/Masters/orthomcl_unmapped.groups.txt", dataset);
    //dataset->ScoreSynteny("/Users/Camous/Desktop/Masters/orthomcl_unmapped.scores.txt");
    //Parse::Groups("/Users/Camous/Desktop/Masters/proteinortho_unmapped.groups.txt", proteinOrtho);
    //proteinOrtho->ScoreSynteny("/Users/Camous/Desktop/Masters/proteinortho_unmapped.scores.txt");
    /*
    Database *database1 = new Database();
    Database *database2 = new Database();
    Timer timer;
    timer.Start();
    Parse::fnaGroups2("/Users/Camous/Desktop/Masters/group_alns_orthomcl/", database1);
    printf("Parsing took: %i seconds\n\n", (int)timer.Elapsed());
    timer.Start();
    Parse::fnaGroups2("/Users/Camous/Desktop/Masters/group_alns_proteinortho/", database2);
    printf("Parsing took: %i seconds\n\n", (int)timer.Elapsed());
    exit(0);
    timer.Start();
    int orthomclgenes = Parse::fnaGroups("/Users/Camous/Desktop/Masters/group_alns_orthomcl/", database1);
    printf("Parsing took: %i seconds\n\n", (int)timer.Elapsed());
    timer.Start();
    int proteinorthogenes = Parse::fnaGroups("/Users/Camous/Desktop/Masters/group_alns_proteinortho/", database2);
    printf("Parsing took: %i seconds\n\n", (int)timer.Elapsed());

    unordered_map<Group*, set<Group*>*> groupmap1, groupmap2;

    int mappedGenes = 0;
    int unmappedGenes = 0;
    for (auto group = database1->groups.begin(); group != database1->groups.end(); group++)
    {
        set<Group*> *groupset = new set<Group*>;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                Group *g = database2->grouppool2.at((*gene)->id);
                groupset->insert(g);

                Gene *gene2 = database2->genepool[(*gene)->id];
                while ((*gene)->votes > 0 && gene2->votes > 0 )
                {
                    (*gene)->votes--;
                    gene2->votes--;
                    mappedGenes++;
                }
                unmappedGenes += (*gene)->votes + gene2->votes;

            }
            catch (const out_of_range & e)
            {
                unmappedGenes += (*gene)->votes;
            }
        }
        groupmap1[(*group)] = groupset;
    }
    printf("Mapped Genes: %i\n", mappedGenes);
    //printf("Unmapped Genes: %i\n\n", unmappedGenes);

    for (auto group = database2->groups.begin(); group != database2->groups.end(); group++)
    {
        set<Group*> *groupset = new set<Group*>;
        
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                Group *g = database1->grouppool2.at((*gene)->id);
                groupset->insert(g);
            }
            catch (const out_of_range & e) {}
        }
        groupmap2[(*group)] = groupset;
    }

    int emptyMaps = 0;
    int nonemptyMaps = 0;
    int groupsInNonemptyMaps = 0;
    int totalMaps = 0;
    int groupsInAllMaps = 0;
    
    ofstream out( "/Users/Camous/Desktop/Masters/orthomcl_to_proteinortho_groupmap.txt", ifstream::out );

    for (auto map = groupmap1.begin(); map != groupmap1.end(); map++)
    {
        out << map->first->id << ':';
        
        totalMaps++;
        groupsInAllMaps += map->second->size();

        if (map->second->size() == 0)
            emptyMaps++;
        else
        {
            nonemptyMaps++;
            groupsInNonemptyMaps += map->second->size();
        }
        
        for (auto group = map->second->begin(); group != map->second->end(); group++)
            out << ' ' << (*group)->id;
        out << "\n";
    }
    out.close();
    printf("Orthomcl groups with no mappings: %i (%.2f%%)\n", emptyMaps, emptyMaps * 100.0f / database1->groups.size() );
    printf("The genes in each orthomcl group map to %f groups in proteinortho on average\n", (double) groupsInAllMaps / totalMaps );
    printf("The genes in each nonempty orthomcl group map to %f groups in proteinortho on average\n\n", (double) groupsInNonemptyMaps / nonemptyMaps );

    emptyMaps = 0;
    nonemptyMaps = 0;
    groupsInNonemptyMaps = 0;
    totalMaps = 0;
    groupsInAllMaps = 0;
    out.open( "/Users/Camous/Desktop/Masters/proteinortho_to_orthomcl_groupmap.txt", ifstream::out );
    
    for (auto map = groupmap2.begin(); map != groupmap2.end(); map++)
    {
        out << map->first->id << ':';
        totalMaps++;
        groupsInAllMaps += map->second->size();

        if (map->second->size() == 0)
            emptyMaps++;
        else
        {
            nonemptyMaps++;
            groupsInNonemptyMaps += map->second->size();
        }
        
        for (auto group = map->second->begin(); group != map->second->end(); group++)
            out << ' ' << (*group)->id;
        out << "\n";
    }
    out.close();
    printf("proteinortho groups with no mappings: %i (%.2f%%)\n", emptyMaps, emptyMaps * 100.0f / database2->groups.size() );
    printf("The genes in each proteinortho group map to %f groups in orthomcl on average\n", (double) groupsInAllMaps / totalMaps );
    printf("The genes in each nonempty proteinortho group map to %f groups in orthomcl on average\n\n", (double) groupsInNonemptyMaps / nonemptyMaps );

    int64_t unmappedNucleotides = 0;
    int64_t mappedNucleotides = 0;
    unmappedGenes = 0;
    out.open( "/Users/Camous/Desktop/Masters/group_alns_orthomcl_unmapped.txt", ifstream::out );
    for (auto group = database1->groups.begin(); group != database1->groups.end(); group++)
    {
        bool singleFound = false;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if ((*gene)->votes > 0)
            {
                if (!singleFound)
                {
                    singleFound = true;
                    out << '>' << (*group)->id << '\n';
                }
                out << (*gene)->votes << " x " <<(*gene)->id << '\n';
                unmappedGenes += (*gene)->votes;
                unmappedNucleotides += (*gene)->votes * (*gene)->id.size();
            }
            
            mappedNucleotides += ((*gene)->antiVotes - (*gene)->votes) * (*gene)->id.size();
        }
        if (singleFound)
            out << '\n';
    }
    out.close();

    printf("average nucleotides per mapped gene: %f\n\n", (double) mappedNucleotides / mappedGenes);
    printf("unmapped genes in orthomcl: %i (%.2f%%)\n", unmappedGenes, unmappedGenes*100.0f / orthomclgenes);
    printf("average nucleotides per unmapped gene: %f\n\n", (double) unmappedNucleotides / unmappedGenes);

    unmappedNucleotides = 0;
    unmappedGenes = 0;
    out.open( "/Users/Camous/Desktop/Masters/group_alns_proteinortho_unmapped.txt", ifstream::out );
    for (auto group = database2->groups.begin(); group != database2->groups.end(); group++)
    {
        bool singleFound = false;
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if ((*gene)->votes > 0)
            {
                if (!singleFound)
                {
                    singleFound = true;
                    out << '>' << (*group)->id << '\n';
                }
                out << (*gene)->votes << " x " <<(*gene)->id << '\n';
                unmappedGenes += (*gene)->votes;
                unmappedNucleotides += (*gene)->votes * (*gene)->id.size();
            }
        }
        if (singleFound)
            out << '\n';
    }
    out.close();
    printf("unmapped genes in proteinortho: %i (%.2f%%)\n", unmappedGenes, unmappedGenes*100.0f / proteinorthogenes);

    printf("average nucleotides per unmapped gene: %f\n\n", (double) unmappedNucleotides / unmappedGenes);

    Histogram(database1, database2, "/Users/Camous/Desktop/Masters/orthomcl_histogram.txt");
    Histogram(database2, database1, "/Users/Camous/Desktop/Masters/proteinortho_histogram.txt");
    */
    //
    exit(0);
    {
    Settings *settings =  new Settings;

    if (argc <= 1)
        settings->INPUT_DIRECTORY = "/Users/Camous/Dev/Syntenizer2000/data/ygob/";
    else
        settings->INPUT_DIRECTORY = argv[1];

    unordered_map<int64_t, double> *similaritypool = new unordered_map<int64_t, double>;
    unordered_map<string, Gene *> *genepool = new unordered_map<string, Gene *>;
    unordered_map<Gene *, Group *> *grouppool = new unordered_map<Gene *, Group *>;
    vector<Group*> *groups = new vector<Group*>;
    vector<Strain*> *strains = new vector<Strain*>;

    ifstream input(settings->INPUT_DIRECTORY+SIMILAR_SEQUENCE_FILE, ifstream::in);
    string line;
    int i;
    string str[8];
    //Timer timer;

    //vector<Group> groups, newGroups;
    strains->clear();
    //groups.clear();
    timer.Start();
    int64_t counter = 0;
    while (getline(input, line))
    {
        stringstream sstream(line);
        for (i = 0; i < 8; i++)
            getline(sstream, str[i], '\t');

        if (str[0] != str[1])
        {
            double score = (1.0 / THRESHOLD) * pow(stod(str[4]), stod(str[5]));
            if (score > 1.0)
                continue;

            Gene *gene[2];

            //printf("tuple: \"%s\" \"%s\" score: %f\n",  str[0].c_str(), str[1].c_str(), score);
            for (int j = 0; j < 2; j++)
            {
                /*
                size_t loc = str[j].find("|");
                string strainID = str[j].substr(0, loc);
                string geneID = str[j].substr(loc+1);
                Strain *strain = GetStrain(strainID, &strains);
                
                if (strain == NULL)
                {
                    printf("Strain \"%s\" not found!\n", strainID.c_str());
                    exit(1);
                }
                 */
                
                //printf("Got strain: %s\n", strain->id.c_str());
                //gene[j] = strain->GetGene(geneID);
                try
                {
                    gene[j] = genepool->at(str[j]);
                    /*
                     gene[j] = strain->genemap.at(geneID);
                    if (gene[j]->id.empty() || gene[j]->id == "")
                    {
                        printf("Gene \"%s\" in strain \"%s\" not found!\n", geneID.c_str(), strainID.c_str());
                        exit(1);
                    }
                     */
                }
                catch (...)
                {
                    Strain *strain = new Strain;
                    size_t loc = str[j].find("|");
                    strain->id = str[j].substr(0, loc);
                    strain->Parse(settings->INPUT_DIRECTORY, genepool);
                    //string geneID = str[j].substr(loc+1);
                    
                    //strain->id =

                    //ParseStrain(strain, INPUT_DIRECTORY);
                    strains->push_back(strain);
                    
                    try
                    {
                        //gene[j] = (*genepool)[str[j]];
                        gene[j] = genepool->at(str[j]);
                    }
                    catch(...)
                    {
                        printf("Gene \"%s\" in strain \"%s\" not found!\n", str[j].substr(loc+1).c_str(), strain->id.c_str());
                        exit(1);
                    }
                }
                /*
                if (gene[j] == NULL)
                {
                    printf("Gene \"%s\" not found!\n", geneID.c_str());
                    exit(1);
                }*/
                //gene->group = group;
                //group->genes.push_back(gene);
            }
            /**/
            int64_t adr_1 = (int64_t)gene[0];
            int64_t adr_2 = (int64_t)gene[1];
            //pair<int64_t,int64_t> pair = adr_1 < adr_2 ? make_pair(adr_1, adr_2) : make_pair(adr_2, adr_1);
            int64_t pair = adr_1 < adr_2 ? (adr_1 ^ adr_2) : (adr_2 ^ adr_1);
            (*similaritypool)[pair] = 1.0-score;
            /**/
            Relation *r1 = new Relation(gene[1], score);
            Relation *r2 = new Relation(gene[0], score);
            gene[0]->relations.push_back(r1);
            gene[1]->relations.push_back(r2);
        }
        counter++;
        //if (timer.Elapsed() > RUN_LENGTH)
        //    break;
    }
    input.close();
    printf("Parsed %lli similarities in %i secs.\n", counter, (int)timer.Elapsed());

    int64_t genes = 0;
    counter = 0;
    for (auto it = genepool->begin(); it != genepool->end(); it++)
    {
        genes++;
        Gene *gene = it->second;
        counter += gene->relations.size();
    }
    printf("total relations: %lli, total genes: %lli, avg raw group size %lli.\n", counter, genes, genes > 0 ? counter/genes: 0);
    /*
    timer.Start();
    printf("Populating Groups\n");
    counter = PopulateGroups(grouppool, groups, genepool);
    printf("Populated %lli groups in %i secs.\n", counter, (int)timer.Elapsed());
    
    timer.Start();
    printf("Populating Neighbours\n");
    PopulateNeighbours(strains);
    printf("Populated neighbours in %i secs.\n", (int)timer.Elapsed());
        */
    /*
    double avg = 0;
    for (auto pool = grouppool.begin(); pool !=  grouppool.end(); pool++)
    {
        Group *group = pool->second;
        avg += 1.0f / group->genes.size();
    }
    */
    printf("average number of genes per group: %f\n", (double)genepool->size()/(double)counter);
    printf("size of genepool: %lu\n", genepool->size());
    printf("size of grouppool: %lu\n", grouppool->size());

    {
        ofstream output( string(settings->INPUT_DIRECTORY) + "raw.groups", ifstream::out );
        int64_t numGenesFiltered = 0;
        int64_t numGenesUnfiltered = 0;
        int64_t numGroups = 0;
        
        for (auto group = groups->begin(); group != groups->end(); group++)
        {
            int members = 0;
            stringstream filtered;
            //double highestScore = 0.0f;
            
            filtered << (*group)->id << ':';
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            {
                filtered << " " << (*gene)->strain->id << '|' << (*gene)->id;
                members++;
                numGenesUnfiltered++;
            }
            filtered << '\n';
            
            if (members > 1)
            {
                numGenesFiltered += members;
                numGroups++;
                output << filtered.str();
            }
            /*
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->score > highestScore)
                    highestScore = (*gene)->score;
            
            for (gene = group->genes.begin(); gene != group->genes.end(); gene++)
            {
                if (!(*gene)->excluded && ((*gene)->score / highestScore) * 100.0f >= treshold && (*gene)->score > 0.0f)
                {
                    filtered << " " << (*gene)->strain << '|' << (*gene)->id;
                    members++;
                }
            }
            filtered << '\n';
            
            //if (highestScore > 0.0f)
            //{
            if (members > 1)
            {
                numGenes += members;
                numGroups++;
                output << filtered.str();
            }
            //else
            //    counter++;
            //}
            //else counter++;
             */
        }
        //printf("Threshold: %i\n", treshold);
        printf("Total number of filtered groups: %lld\n", numGroups);
        printf("Total number of genes in filtered groups: %lld\n\n", numGenesFiltered);
        printf("Total number of genes in unfiltered groups: %lld\n\n", numGenesUnfiltered);

        //printf("counter: %i\n", counter);
        output.close();
    }
    /*
    //Group *g = groups.front();
    double score = 0;
    for (auto g = groups->begin(); g != groups->end(); g++)
    {
        score = Syntenize3000(*g, similaritypool);
        if (score != 0)
        {
            printf("%s: synteny score %f\n", (*g)->id.c_str(), score);
            return 0;
        }
    }
     */
    return 1;
    printf("Syntenizing\n");
    timer.Start();
    //Syntenize(&groups);
    printf("Syntenizing took %i seconds\n", (int)timer.Elapsed());
    //
    long double sumd = 0.0f;
    int64_t sumi = 0;
    genes = 0;
    ofstream output( string(settings->INPUT_DIRECTORY) + "scores_" + settings->GROUP_FILE, ifstream::out );
    
    vector<int> range;
    
    range.resize(SIZE_OF_NEIGHTBOURHOOD);
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        range[i] = 0;
    
    int s;
    vector<Group*>::iterator group;
    for (group = groups->begin(); group != groups->end(); group++)
    {
        vector<Gene*>::iterator gene;
        
        output << (*group)->id << ':';
        
        for (s = 0, gene = (*group)->genes.begin(); gene != (*group)->genes.end(); s++, gene++)
        {
            if ((*group)->genes.size() > 1)
            {
                sumd += (*gene)->score/((long double)(*group)->genes.size()-1.0f);
                sumi += (*gene)->score;
                
                int r = (*gene)->score/((long double)(*group)->genes.size()-1.0f);
                
                if (r == SIZE_OF_NEIGHTBOURHOOD)
                    r--;
                
                range[r]++;
            }
            else
                range[0]++;
            
            genes++;
            output << " " << (*gene)->score;
        }
        output << "\n";
    }
    output.close();
    printf("Weighed average: %Lf\n", sumd / (long double)genes);
    printf("Raw average: %lld\n", sumi / genes);
    printf("number of groups: %lu\n", groups->size());
    printf("number of genes: %lld\n", genes);
    printf("range:");

    return 0;
    }
}
