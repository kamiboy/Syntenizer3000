//
//  Group.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Group.hpp"

//using namespace std;

bool Group::HasGene(Gene* g)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (*gene == g)
            return true;
    
    return false;
}

void Group::RemoveGene(Gene* g)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (*gene == g)
            genes.erase(gene);
}

void Group::InsertGenes(vector<Gene*> newGenes, unordered_map<Gene *, Group *> *grouppool)
{
    for (auto gene = newGenes.begin(); gene != newGenes.end(); gene++)
        if (!this->HasGene(*gene))
        {
            //printf("moving gene %s into %s.\n", (*gene)->id.c_str(), this->id.c_str());
            genes.push_back(*gene);
            (*gene)->group = this;
            (*grouppool)[*gene] = this;
        }
}

int Group::SharedGenes(vector<Gene*> newGenes)
{
    int shared = 0;
    vector<Gene*>::iterator gene;
    
    for (gene = newGenes.begin(); gene != newGenes.end(); gene++)
        if (!this->HasGene(*gene))
            shared++;
    
    return shared;
}

double Group::SyntenizeSophisticated()
{
    int count = 0;
    double scoresum = 0.0f;
    
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
            if (*g1 != *g2)
            {
                count++;
                double score = Score::Sophisticated(*g1, *g2);

                (*g1)->score += score;
                (*g2)->score += score;
                scoresum += score;
                /*
                double score1 = nw_align((*g1)->neighbours, (*g2)->neighbours, 1, false);
                double score2 = nw_align((*g1)->neighbours, (*g2)->neighbours, 1, true);
                //printf("score1: %f\n", score1);
                //printf("score2: %f\n", score2);
                count++;
                scoresum += (score1 > score2) ? score1 : score2;
                 */
            }

    syntenyScoreSophisticated = scoresum / (double)count;
    return syntenyScoreSophisticated;
}
/*
double Group::SyntenizeNaive()
{
    int count = 0;
    int scoresum = 0;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
            if (*g1 != *g2)
            {
                count++;
                for (int i = 0; i < (*g1)->neighbours.size(); i++)
                    for (int j = 0; j < (*g2)->neighbours.size(); j++)
                        if ((*g1)->neighbours[i]->group == (*g2)->neighbours[j]->group)
                            scoresum++;
            }
    syntenyScoreNaive = (double)scoresum / (double)count;
    
    return syntenyScoreNaive;

}
*/
/*
 double Group::SyntenizeNaive()
 {
 int count = 0;
 int scoresum = 0;
 for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
 for (auto g2 = g1+1; g2 != genes.end(); g2++)
 if (*g1 != *g2)
 {
 count++;
 for (int i = 0; i < (*g1)->neighbours.size(); i++)
 for (int j = 0; j < (*g2)->neighbours.size(); j++)
 if ((*g1)->neighbours[i]->group == (*g2)->neighbours[j]->group)
 scoresum++;
 }
 syntenyScoreNaive = (double)scoresum / (double)count;
 
 return syntenyScoreNaive;
 
 }
 */
double Group::SyntenizeOld()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Old(*g1, *g2);
            }
        }
    }
    
    syntenyScoreOld = (double)scoresum / (double)count;
    
    return syntenyScoreOld;
    
}

double Group::SyntenizeSimple()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Simple(*g1, *g2);
            }
        }
    }
    
    syntenyScoreSimple = (double)scoresum / (double)count;
    
    return syntenyScoreSimple;
}

double Group::SyntenizeAdjusted()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Adjusted(*g1, *g2);
            }
        }
    }
    
    syntenyScoreAdjusted = (double)scoresum / (double)count;
    
    return syntenyScoreAdjusted;
    
}

double Group::SyntenizeTest()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Test(*g1, *g2);
            }
        }
    }
    
    syntenyScoreTest = (double)scoresum / (double)count;
    
    return syntenyScoreTest;
    
}

double Group::SyntenizeFast()
{
    int count = 0;
    double scoresum = 0.0f;
    
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
            if (*g1 != *g2)
            {
                count++;
                double score = Score::Fast(*g1, *g2);
                //(*g1)->score += score;
                //(*g2)->score += score;
                scoresum += score;
            }
    
    syntenyScoreFast = scoresum / (double)count;
    return syntenyScoreFast;
}

double Group::SyntenizeAgainstGene(Gene *g1)
{
    double count = 0;
    double scoresum = 0;

    for (auto g2 = genes.begin(); g2 != genes.end(); g2++)
        if (g1 != *g2)
        {
            count++;
            scoresum += Score::Sophisticated(g1, *g2);
            //scoresum += Score::Simple(g1, *g2);
        }

    return(scoresum / count);
}

int Group::CountUniqueStrains()
{
    set<Strain *> uniquestrains;
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        uniquestrains.insert((*gene)->strain);

    return uniquestrains.size();
}

int Group::CountOrthologs()
{
    int orthologs = 0;
    vector<Strain *> strains;
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        strains.push_back((*gene)->strain);

    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (count(strains.begin(), strains.end(), (*gene)->strain) == 1)
            orthologs++;
    }

    return orthologs;
}

int Group::CountParalogs()
{
    int paralogs = 0;
    vector<Strain *> strains;
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        strains.push_back((*gene)->strain);
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (count(strains.begin(), strains.end(), (*gene)->strain) > 1)
            paralogs++;
    }
    
    return paralogs;
    /*
    int paralogs = 0;
    set<Strain *> strains;
    set<Strain *> paralogstrains;
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);
    }
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs++;
    
    return paralogs;
     */
}

void Group::Prune()
{
    int pure_paralog_groups = 0;
    int unambiguous_groups = 0;
    int ambiguous_groups = 0;
    int split_groups = 0;

    bool has_paralogs = false;
    unordered_map<Strain *, vector<Gene*>*> strains;

    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
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
        vector <Gene*> paralogs;
        vector <Gene*> homologs;
        
        for (auto strain = strains.begin(); strain != strains.end(); strain++)
        {
            if ((*strain).second->size() > 1)
                paralogs.insert(paralogs.end(), (*strain).second->begin(), (*strain).second->end());
            else
                homologs.insert(homologs.end(), (*strain).second->begin(), (*strain).second->end());
        }
        
        if (homologs.size() == 0)
        {
            pure_paralog_groups++;
            return;
        }
        
        Strain *prev_strain = (*paralogs.begin())->strain;
        
        bool split = false;
        int total_scoresum = 0;
        int max_scoresum = 0;
        for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
        {
            double scoresum = 0.0f;
            
            for (auto g2 = homologs.begin(); g2 != homologs.end(); g2++)
                scoresum += Score::Sophisticated(*g1, *g2);
            
            if (prev_strain != (*g1)->strain)
            {
                prev_strain = (*g1)->strain;
            }
            scoresum /= (double)homologs.size();
            
            for (auto g2 = paralogs.begin(); g2 != paralogs.end(); g2++)
                if (g1 != g2 && scoresum > 0 && Score::Sophisticated(*g1, *g2) == 0)
                    split = true;
            
            total_scoresum += scoresum;
            if (scoresum > max_scoresum)
                max_scoresum = scoresum;
        }
        if (split)
            split_groups++;
        
        if (max_scoresum > 0 && total_scoresum == max_scoresum)
            unambiguous_groups++;
        else
            ambiguous_groups++;
    }
}

bool Group::HasStrain(Strain *strain)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if ((*gene)->strain == strain)
            return true;

    return false;
}

vector<Group *> Group::Split2()
{
    vector<Group *> splitgroups;
    vector<Gene *> candidates;
    
    Gene *candidate = NULL;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            double score = Score::Sophisticated(*g1, *g2);
            (*g1)->score += score;
            (*g2)->score += score;
            
            if (candidate == NULL || candidate->score < (*g1)->score)
                candidate = *g1;

            if (candidate == NULL || candidate->score < (*g2)->score)
                candidate = *g2;
        }
    }
    
    candidate->group = NULL;
    candidates.push_back(candidate);
    genes.erase(find(genes.begin(), genes.end(), candidate));

    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        Gene *bestcandidate = NULL;
        for (auto g2 = candidates.begin(); g2 != candidates.end(); g2++)
        {
            double score;
            if ((*g2)->group == NULL)
                score = Score::Sophisticated(*g1, *g2);
            else
                score = (*g2)->group->SyntenizeAgainstGene(*g1);

            if (bestcandidate == NULL || bestcandidate->score < score)
            {
                bestcandidate = *g2;
                bestcandidate->score = score;
            }
        }

        if (bestcandidate->score < 3)
        {
            (*g1)->group = NULL;
            candidates.push_back(*g1);
            //genes.erase(g1);
            continue;
        }

        if (bestcandidate->group == NULL)
        {
            Group *newgroup = new Group;
            newgroup->id = id + "_split" + to_string(splitgroups.size());
            bestcandidate->group = newgroup;
            bestcandidate->group->genes.push_back(bestcandidate);
            splitgroups.push_back(newgroup);
        }
        (*g1)->group = bestcandidate->group;
        bestcandidate->group->genes.push_back(*g1);
    }

    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
    {
        printf("%s\t(%lu)\t[%.2f]\t{%i}\n", (*group)->id.c_str(), (*group)->genes.size(), (*group)->SyntenizeSophisticated(), (*group)->CountParalogs());
        /*
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            printf(" %s", (*gene)->id_full.c_str());
            //printf(" %s", (*gene)->strain->species.c_str());

        printf("\n\n");
         */
        printf("\n");
    }
    
    return splitgroups;
}

vector<Group *> Group::Split3()
{
    vector<Group *> splitgroups;
    vector<Gene *> candidates;
    vector<Gene *> paralogs;
    vector<Gene *> orthologs;
    set<Strain *> strains;
    set<Strain *> paralogstrains;

    printf("%s [%lu] (%i):", id.c_str(), genes.size(), CountParalogs());

    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);
        
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs.push_back(*gene);
        else
            orthologs.push_back(*gene);

    if (paralogs.size() < 2)
        printf("!");

    printf(" (%lu)", paralogs.size() );
    Gene *candidate = NULL;
    for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != paralogs.end(); g2++)
        {
            double score = Score::Sophisticated(*g1, *g2);
            (*g1)->score += score;
            (*g2)->score += score;
            
            if (candidate == NULL || candidate->score < (*g1)->score)
                candidate = *g1;
            
            if (candidate == NULL || candidate->score < (*g2)->score)
                candidate = *g2;
        }
    }

    candidate->group = new Group;
    candidate->group->id = id + "_split" + to_string(splitgroups.size());
    candidate->group->genes.push_back(candidate);
    splitgroups.push_back(candidate->group);
    candidates.push_back(candidate);
    paralogs.erase(find(paralogs.begin(), paralogs.end(), candidate));

    for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
    {
        Gene *bestcandidate = NULL;
        for (auto g2 = candidates.begin(); g2 != candidates.end(); g2++)
        {
            double score = (*g2)->group->SyntenizeAgainstGene(*g1);
            
            if (!(*g2)->group->HasStrain((*g1)->strain) && (bestcandidate == NULL || bestcandidate->score < score ))
            {
                bestcandidate = *g2;
                bestcandidate->score = score;
            }
        }

        if (bestcandidate == NULL || bestcandidate->score == 0)
        {
            (*g1)->group = new Group;
            (*g1)->group->id = id + "_split" + to_string(splitgroups.size());
            (*g1)->group->genes.push_back(*g1);
            splitgroups.push_back((*g1)->group);
            candidates.push_back(*g1);

            continue;
        }

        (*g1)->group = bestcandidate->group;
        bestcandidate->group->genes.push_back(*g1);
    }

    orphans = 0;
    printf(" {%lu}", splitgroups.size());
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
    {
        //printf("%s (%lu) [%.2f] {%i}:\n", (*group)->id.c_str(), (*group)->genes.size(), (*group)->SyntenizeSophisticated(), (*group)->CountParalogs());
        if ((*group)->genes.size() > 1)
            //printf("%s\t%lu\t%.2f\t%i\n", (*group)->id.c_str(), (*group)->genes.size(), (*group)->SyntenizeSophisticated(), (*group)->CountParalogs());
            printf(" (%lu %.2f)", (*group)->genes.size(), (*group)->SyntenizeSophisticated());
        else
        {
            //printf("orphan\n");
            orphans++;
        }
        /*
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            printf(" %s", (*gene)->id_full.c_str());
        //printf(" %s", (*gene)->strain->species.c_str());
        
        printf("\n\n");
        */
    }
    printf(" - %i\n", orphans);

    return splitgroups;
}

vector<Group *> Group::Split4()
{
    vector<Group *> splitgroups;
    vector<Gene *> candidates;
    vector<Gene *> paralogs;
    vector<Gene *> orthologs;
    set<Strain *> strains;
    set<Strain *> paralogstrains;

    printf("%s [%lu] (%i):", id.c_str(), genes.size(), CountParalogs());

    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);

    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs.push_back(*gene);
        else
            orthologs.push_back(*gene);

    printf(" (%lu)", paralogs.size() );

    do
    {
        Gene *candidate = NULL;
        for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
        {
            for (auto g2 = g1+1; g2 != paralogs.end(); g2++)
            {
                double score = Score::Sophisticated(*g1, *g2);
                //double score = Score::Simple(*g1, *g2);
                (*g1)->score += score;
                (*g2)->score += score;
                
                if (candidate == NULL || candidate->score < (*g1)->score)
                    candidate = *g1;
                
                if (candidate == NULL || candidate->score < (*g2)->score)
                    candidate = *g2;
            }
        }
        
        if (paralogs.size() == 1)
            candidate = paralogs[0];
        
        candidate->group = new Group;
        candidate->group->id = id + "_split" + to_string(splitgroups.size());
        candidate->group->genes.push_back(candidate);
        splitgroups.push_back(candidate->group);
        candidates.push_back(candidate);
        paralogs.erase(find(paralogs.begin(), paralogs.end(), candidate));
        
        bool restart;
        do
        {
            restart = false;
            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
            {
                if (!candidate->group->HasStrain((*g1)->strain) && candidate->group->SyntenizeAgainstGene(*g1) > 0 )
                {
                    candidate->group->genes.push_back(*g1);
                    paralogs.erase(g1);
                    restart = true;
                    break;
                }
            }
        }
        while(restart);
    }
    while (paralogs.size() > 0);

    Group *largestgroup = NULL;
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        if (largestgroup == NULL || (*group)->genes.size() > largestgroup->genes.size())
            largestgroup = *group;
        else if ((*group)->genes.size() == largestgroup->genes.size())
            largestgroup = (*group)->SyntenizeSophisticated() > largestgroup->SyntenizeSophisticated() ? *group : largestgroup;

    largestgroup->genes.insert(largestgroup->genes.end(), orthologs.begin(), orthologs.begin());

    orphans = 0;
    printf(" {%lu}", splitgroups.size());
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
    {
        if ((*group)->genes.size() > 1)
            printf(" (%lu %.2f)", (*group)->genes.size(), (*group)->SyntenizeSophisticated());
            //printf(" (%lu %.2f)", (*group)->genes.size(), (*group)->SyntenizeSimple());
        else
        {
            orphans++;
        }
    }
    printf(" - %i\n", orphans);
    
    return splitgroups;
}

vector<Group *> Group::Split5()
{
    vector<Group *> splitgroups;
    vector<Gene *> candidates;
    vector<Gene *> paralogs;
    //vector<Gene *> orthologs;
    set<Strain *> strains;
    set<Strain *> paralogstrains;
    Group *orthologs = new Group;

    //printf("%s [%lu] (%i):", id.c_str(), genes.size(), CountParalogs());
    //printf("%s:", id.c_str());
    //size_t loc;
    //if ((loc = id.find(' ')) != string::npos)
    //    id.erase(id.begin()+loc,id.end());

    printf("%s\t%lu\t%.2f\t%i", id.c_str(), genes.size(), SyntenizeSophisticated(), CountParalogs());

    orthologs->id = id;
    splitgroups.push_back(orthologs);

    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);
    
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs.push_back(*gene);
        else
            orthologs->genes.push_back(*gene);

    //printf(" (%lu)", paralogs.size() );

    Gene *bestcandidate;
    do
    {
        bestcandidate = NULL;
        for (auto gene = paralogs.begin(); gene != paralogs.end(); gene++)
        {
            if (orthologs->HasStrain((*gene)->strain))
                continue;

            (*gene)->score = orthologs->SyntenizeAgainstGene(*gene);
            if (bestcandidate == NULL || (*gene)->score > bestcandidate->score)
                bestcandidate = *gene;
        }

        if (bestcandidate != NULL)
        {
            orthologs->genes.push_back(bestcandidate);
            paralogs.erase(find(paralogs.begin(), paralogs.end(), bestcandidate));
        }
    }
    while(bestcandidate != NULL);

    do
    {
        Gene *candidate = NULL;
        for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
        {
            for (auto g2 = g1+1; g2 != paralogs.end(); g2++)
            {
                double score = Score::Sophisticated(*g1, *g2);
                //double score = Score::Simple(*g1, *g2);
                (*g1)->score += score;
                (*g2)->score += score;
                
                if (candidate == NULL || candidate->score < (*g1)->score)
                    candidate = *g1;

                if (candidate == NULL || candidate->score < (*g2)->score)
                    candidate = *g2;
            }
        }

        if (paralogs.size() == 1)
            candidate = paralogs[0];

        candidate->group = new Group;
        candidate->group->id = id + "-" + to_string(splitgroups.size()+1);
        candidate->group->genes.push_back(candidate);
        splitgroups.push_back(candidate->group);
        candidates.push_back(candidate);
        paralogs.erase(find(paralogs.begin(), paralogs.end(), candidate));

        bool restart;
        do
        {
            restart = false;
            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
            {
                if (!candidate->group->HasStrain((*g1)->strain) && candidate->group->SyntenizeAgainstGene(*g1) > 0 )
                {
                    candidate->group->genes.push_back(*g1);
                    paralogs.erase(g1);
                    restart = true;
                    break;
                }
            }
        }
        while(restart);
    }
    while (paralogs.size() > 0);

    vector<Gene *> orphans;

    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
    {
        if ((*group)->genes.size() == 1)
            orphans.push_back((*group)->genes[0]);
    }

    //splitgroups.erase(remove_if(splitgroups.begin(), splitgroups.end(), [](Group *group) { return group->genes.size() == 1; }), splitgroups.end());

    int unorphaned = 0;
    for (auto gene = orphans.begin(); gene != orphans.end(); gene++)
    {
        Group *bestgroup = NULL;

        for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        {
            if ((*group)->HasStrain((*gene)->strain) || (*gene)->group == *group)
                continue;
            if (bestgroup == NULL || (*group)->genes.size() > bestgroup->genes.size())
                bestgroup = *group;
        }
        
        if (bestgroup != NULL)
        {
            bestgroup->genes.push_back(*gene);
            unorphaned++;

            splitgroups.erase(find(splitgroups.begin(), splitgroups.end(), (*gene)->group));
        }
    }
    
    this->orphans = 0;
    printf("\t%lu\t", splitgroups.size());
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        if ((*group)->genes.size() > 1)
            //int nop;
            printf("{%lu %.2f} ", (*group)->genes.size(), (*group)->SyntenizeSophisticated());
        else
        {
            this->orphans++;
            //printf("\t%s\t%lu", (*group)->genes[0]->strain->id_alt.c_str(), (*group)->genes[0]->contig->genes.size());
            printf("\t%s\t%.2f", (*group)->genes[0]->strain->id_alt.c_str(), orthologs->SyntenizeAgainstGene((*group)->genes[0]));
        }

    this->orphans = orphans.size() - unorphaned;
    //printf(" - %d\n", this->orphans);
    printf("\n");
    return splitgroups;
}

vector<Group *> Group::Split()
{
    vector<Group *> splitgroups;
    int pure_paralog_groups = 0;
    int unambiguous_groups = 0;
    int ambiguous_groups = 0;
    int split_groups = 0;
    
    bool has_paralogs = false;
    unordered_map<Strain *, vector<Gene*>*> strains;
    
    for (auto g1 = this->genes.begin(); g1 != this->genes.end(); g1++)
    {
        //vector<Gene*> pruned;
        Gene *bestparalog = NULL;
        if ((*g1)->paralog)
        {
            (*g1)->score = 0;
            for (auto g2 = this->genes.begin(); g2 != this->genes.end(); g2++)
                if ((*g1) != (*g2))
                    (*g1)->score += Score::Sophisticated(*g1, *g2);
            (*g1)->score /= double(this->genes.size()-1);
            if ((*g1)->score == 0)
                this->genes.erase(g1);
            else
            {
                if (bestparalog == NULL || (*g1)->score > bestparalog->score)
                    bestparalog = *g1;
            }
        }
        /*
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
         */
    }
    /*
     if (has_paralogs)
     {
     vector <Gene*> paralogs;
     vector <Gene*> homologs;
     
     for (auto strain = strains.begin(); strain != strains.end(); strain++)
     {
     if ((*strain).second->size() > 1)
     paralogs.insert(paralogs.end(), (*strain).second->begin(), (*strain).second->end());
     else
     homologs.insert(homologs.end(), (*strain).second->begin(), (*strain).second->end());
     }
     
     if (homologs.size() == 0)
     {
     pure_paralog_groups++;
     return splitgroups;
     }
     
     Strain *prev_strain = (*paralogs.begin())->strain;
     
     bool split = false;
     int total_scoresum = 0;
     int max_scoresum = 0;
     for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
     {
     double scoresum = 0.0f;
     
     for (auto g2 = homologs.begin(); g2 != homologs.end(); g2++)
     scoresum += Score::Sophisticated(*g1, *g2);
     
     if (prev_strain != (*g1)->strain)
     {
     prev_strain = (*g1)->strain;
     }
     scoresum /= (double)homologs.size();
     
     for (auto g2 = paralogs.begin(); g2 != paralogs.end(); g2++)
     if (g1 != g2 && scoresum > 0 && Score::Sophisticated(*g1, *g2) == 0)
     split = true;
     
     total_scoresum += scoresum;
     if (scoresum > max_scoresum)
     max_scoresum = scoresum;
     }
     if (split)
     split_groups++;
     
     if (max_scoresum > 0 && total_scoresum == max_scoresum)
     unambiguous_groups++;
     else
     ambiguous_groups++;
     }
     */
    return splitgroups;
}

void Group::GenerateSyntenyMap(string OUTPUT_DIRECTORY )
{
    ofstream out( OUTPUT_DIRECTORY + id + ".syntenymap.csv", ifstream::out );

    printf("Generating Synteny map for %s\n", id.c_str());
    SyntenizeSophisticated();
    //SyntenizeFast();

    set<Group *> groups;
    double bestscore = -1;
    Gene *candidate = NULL;
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
    {
        for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
            if ((*gene)->neighbours[i] != NULL )
                groups.insert((*gene)->neighbours[i]->group);

        if (candidate == NULL || (*gene)->score > bestscore)
        {
            bestscore = (*gene)->score;
            candidate = *gene;
        }
    }

    printf("Best candidate %s: %.2f\n", candidate->id_full.c_str(), candidate->score );

    unordered_map<Group *, int> groupcolours;
    int colour = 0;
    for (auto group = groups.begin(); group != groups.end(); group++)
        groupcolours[*group] = ++colour;

    vector<Gene *> remaining;
    vector<Gene *> sorted;
    sorted.push_back(candidate);
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if ((*gene) != candidate)
            remaining.push_back(*gene);

    do
    {
        bestscore = -1;
        Gene *bestmatch = NULL;
        for (auto gene = remaining.begin(); gene != remaining.end(); gene++)
        {
            double score = Score::Sophisticated(candidate, *gene);
            //double score = Score::Fast(candidate, *gene);
            if (bestmatch == NULL || score > bestscore)
            {
                bestscore = score;
                bestmatch = *gene;
            }
        }
        printf("Best match %s: %.2f\n", bestmatch->id_full.c_str(), bestscore );

        sorted.push_back(bestmatch);
        remaining.erase(find(remaining.begin(), remaining.end(), bestmatch));
        candidate = bestmatch;
    }
    while (remaining.size() > 0);

    for (auto gene = sorted.begin(); gene != sorted.end(); gene++)
    {
        out << (*gene)->strain->id_alt;
        if (*gene != sorted.back())
            out << ';';
        else
            //out << ";colour\n";
            out << "\n";
    }

    for (int i = 0; i < sorted.size(); i++)
        //out << "0" << (i < sorted.size()-1 ? ";" : ";" + to_string(groupcolours.size()+1) + "\n");
        out << "0" << (i < sorted.size()-1 ? ";" : "\n");
    for (int i = 0; i < sorted.size(); i++)
        //out << SIZE_OF_NEIGHTBOURHOOD+2 << (i < sorted.size()-1 ? ";" : ";" + to_string(groupcolours.size()+2) + "\n");
        out << SIZE_OF_NEIGHTBOURHOOD+2 << (i < sorted.size()-1 ? ";" : "\n");
    for (int i = 0; i < sorted.size(); i++)
        //out << (SIZE_OF_NEIGHTBOURHOOD/2)+1 << (i < sorted.size()-1 ? ";" : ";" + to_string(groupcolours.size()+2) + "\n");
        out << (SIZE_OF_NEIGHTBOURHOOD/2)+1 << (i < sorted.size()-1 ? ";" : "\n");

    for (auto line = groupcolours.begin(); line != groupcolours.end(); line++)
    {
        Group *group = line->first;
        int colour = line->second;

        for (int s = 0; s < sorted.size(); s++)
        {
            Gene * g1 = sorted[s];
            for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
                if (g1->neighbours[i] != NULL && g1->neighbours[i]->group == group)
                {
                    if (i+1 <= (SIZE_OF_NEIGHTBOURHOOD/2))
                        out << i+1;
                    else
                        out << i+2;
                    break;
                }

            if (s < sorted.size()-1)
                out << ';';
            else
                out << '\n';
                //out << ';' << colour << '\n';
        }

        /*
        for (int s = 0; s < sorted.size()-1; s++)
        {
            Gene * g1 = sorted[s];
            Gene * g2 = sorted[s+1];
            
            //if ((*g1)->group != NULL && (*g2)->group != NULL)
            bool connects = false;

            for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
                for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
                    if (g1->neighbours[i] != NULL && g2->neighbours[j] != NULL && g1->neighbours[i]->group == group && g2->neighbours[j]->group == group)
                    {
                        out << i+1 << ';' << j+1;
                        connects = true;
                    }
            if (!connects)
                out << ';';
            
            if (s < sorted.size()-2)
                out << ';';
            else
                out << ';' << colour << '\n';
        }
        */
    }
    out.close();
}

void Group::GenerateSyntenyHistogram(string OUTPUT_DIRECTORY )
{
 /*
    Timer timer;
    timer.Start();
    printf("Generating group synteny histogram to %s:\n", OUTPUT_DIRECTORY.c_str());
    
    ofstream out( OUTPUT_DIRECTORY + id + ".syntenyhistogram.csv", ifstream::out );
    out << "ID\n";
    
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
    {
        out << (*strain)->id_full << '\n';
        (*strain)->ExportGeneRelationship(db, OUTPUT_DIRECTORY);
    }
    out.close();
    printf("Generating Synteny Maps took\t%i seconds\n\n", (int)timer.Elapsed());
 */
}

Group::Group(void)
{
    //id = "";
    //genes.clear();
    //scoresSum.clear();
    coorthologs = -1;
    //discordant = false;
    algebraicConnectivity = NAN;
    syntenyScoreSophisticated = NAN;
    syntenyScoreOld = NAN;
    syntenyScoreSimple = NAN;
    syntenyScoreAdjusted = NAN;
    syntenyScoreFast = NAN;
    syntenyScoreTest = NAN;
}
/*
Group::~Group(void)
{
    printf("Group %s is being deconstructed.\n", id.c_str());
}
*/
