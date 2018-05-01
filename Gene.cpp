//
//  Gene.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Gene.hpp"
//#include "nw.hpp"

vector<Group *> Gene::Match(Dataset *dataset)
{
    double highestscoresum = 0.0f;
    vector<Group*> bestgroups;
    
    for (int i = 0; i <= SIZE_OF_NEIGHTBOURHOOD; i++)
        matchHistogram[i] = 0;

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        double scoresum = 0.0f;
        int count = 0;
        
        if (this->group == *group)
            continue;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            if (*gene == this)
                continue;

            count++;
            int score = Score::Fast(this, *gene);
            matchHistogram[score]++;
            scoresum += score;
        }
        scoresum = scoresum / (double)count;
        groupscores[*group] = scoresum;
        if (scoresum == highestscoresum)
            bestgroups.push_back(*group);
        else if (scoresum > highestscoresum)
        {
            bestgroups.clear();
            highestscoresum = scoresum;
            bestgroups.push_back(*group);
        }
    }
    return(bestgroups);
}

double Gene::Match(Group *group)
{
    double matchscore = 0.0f;
    int count = 0;
    
    for (auto gene = group->genes.begin(); gene != group->genes.end(); gene++)
    {
        if (*gene == this)
            continue;

        count++;
        matchscore += Score::Sophisticated(this, *gene);
    }
    matchscore /= (double)count;

    return(matchscore/(double) SIZE_OF_NEIGHTBOURHOOD);
}

void Gene::CompareNeighbours(Gene *gene)
{
    double score1 = Score::Sophisticated(this, gene);
    double score2 = Score::Simple(this, gene);
    double score3 = Score::Naive(this, gene);
    double score4 = Score::Old(this, gene);
    double score5 = Score::Adjusted(this, gene);
    double score6 = Score::Fast(this, gene);
    printf("Neighbourhood Comparison [%s] %s (%lu) - [%s] %s (%lu)\n", (orientation ? "+" : "-"), (strain->id + '|' + id).c_str(), contig->genes.size(), (gene->orientation ? "+" : "-"), (gene->strain->id + '|' + gene->id).c_str(), gene->contig->genes.size());
    printf("Scores: Sophisticated(%.2f)\tSimple(%.2f)\tNaive(%.2f)\tOld(%.2f)\tAdjusted(%.2f)\tFast(%.2f)\n", score1, score2, score3, score4, score5, score6);

    unordered_map<string, int> groupID;
    int counter = 0;
    for (int j = 0; j < 2; j++)
    {
        Gene *g;
        if (j == 0)
            g = this;
        else
            g = gene;
        for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        {
            if (g->neighbours[i] == NULL || g->neighbours[i]->group == NULL)
                printf("---  ");
            else
            {
                try
                {
                    if (g->neighbours[i]->group != NULL)
                        groupID.at(g->neighbours[i]->group->id);
                }
                catch (const out_of_range & e)
                {
                    groupID[g->neighbours[i]->group->id] = counter;
                    counter++;
                }
                printf("%03i  ", groupID[g->neighbours[i]->group->id]);
            }
        }
        printf("\n");
    }
}

// TAA TAG TGA TGG ATG
double Gene::CalculateGC3s()
{
    double TA = 0;
    double GC = 0;
    for (int codon = 0; codon < sequence.length(); codon+=3)
    {
        if (sequence[codon] == 'T' && (sequence[codon+1] == 'G' || sequence[codon+1] == 'A') && (sequence[codon+2] == 'G' || sequence[codon+2] == 'A'))
            continue;

        if (sequence[codon] == 'A' && sequence[codon+1] == 'T' && sequence[codon+2] == 'G')
            continue;

        if (sequence[codon+2] == 'G' || sequence[codon+2] == 'C')
            GC++;
        if (sequence[codon+2] == 'T' || sequence[codon+2] == 'A')
            TA++;
    }

    return GC / (GC + TA);
}
