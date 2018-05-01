//
//  Strain.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Strain.hpp"

Strain::Strain()
{
    id.clear();
    contigs.clear();
    species.clear();
    //genes.clear();
    //genemap.clear();
    bp = 0;
    ubp = 0;
    CDS = 0;
    tRNA = 0;
    tmRNA = 0;
    rRNA = 0;
    GC = NAN;
}

void Strain::ExportPTT(string OUTPUT_DIRECTORY)
{
    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        ofstream out( OUTPUT_DIRECTORY + id_full + "_" + (*contig)->id + ".ptt", ifstream::out );
        
        out << id_full << "_" << (*contig)->id << ", complete sequence - " << '\n' << (*contig)->genes.size() << " proteins\n" ;
        out << "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
        {
            out << (*gene)->start << ".." << (*gene)->start + (*gene)->length << '\t';
            out << ((*gene)->orientation ? '+' : '-') << '\t';
            out << int((*gene)->length / 3) << '\t';
            out << (*gene)->id_full << '\t';
            out << '-' << '\t';
            out << '-' << '\t';
            out << '-' << '\t';
            out << '-' << '\t';
            out << (*gene)->product << '\n';
        }
        out.close();
    }
}

void Strain::PrintContigs()
{
    /*
    double score1 = Score::Sophisticated(this, gene);
    double score2 = Score::Simple(this, gene);
    double score3 = Score::Naive(this, gene);
    double score4 = Score::Old(this, gene);
    double score5 = Score::Adjusted(this, gene);
    double score6 = Score::Fast(this, gene);
    printf("Neighbourhood Comparison [%s] %s (%lu) - [%s] %s (%lu)\n", (orientation ? "+" : "-"), (strain->id + '|' + id).c_str(), contig->genes.size(), (gene->orientation ? "+" : "-"), (gene->strain->id + '|' + gene->id).c_str(), gene->contig->genes.size());
    printf("Scores: Sophisticated(%.2f)\tSimple(%.2f)\tNaive(%.2f)\tOld(%.2f)\tAdjusted(%.2f)\tFast(%.2f)\n", score1, score2, score3, score4, score5, score6);
    
    */
    
    unordered_map<string, int> groupID;
    int counter = 0;
    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        printf("%s [%i]:\t\t\t", (*contig)->id.c_str(), (*contig)->genes.size());
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
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
}

void Strain::ExportChartData(string OUTPUT_DIRECTORY, Database *db)
{
    ofstream out0( OUTPUT_DIRECTORY + id_full + "_contigs.csv", ifstream::out );
    out0 << "ID;Bound;Colour\n";

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        string colour = "white";
        for (auto cmap = db->contigcolours.begin(); cmap != db->contigcolours.end(); cmap++ )
        {
            if ((*contig)->id_full.find(cmap->first) != string::npos)
            {
                colour = cmap->second;
                break;
            }
        }
        /*
        if ((*contig)->id_full.find("chromosome") != string::npos)
            colour = "palegreen";
        else if ((*contig)->id_full.find("plasmid") != string::npos)
            colour = "skyblue";
        else if ((*contig)->id_full.find("fragment") != string::npos)
            colour = "red";
        else
            colour = "white";
        */

        out0 << (*contig)->id << ";0;" << colour << '\n';
        out0 << (*contig)->id << ';' << (*contig)->length << ';' << colour << '\n';
    }
    out0.close();

    ofstream out1( OUTPUT_DIRECTORY + id_full + "_placements.csv", ifstream::out );
    out1 << "ID;Location;Synteny;Abundance;Coverage;GC3s;Type\n";

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            //out1 << (*contig)->id << ";"<< (*gene)->start + ((*gene)->length / 2) << ';' << ((*gene)->group != NULL ? (*gene)->group->SyntenizeAgainstGene((*gene))/(double)SIZE_OF_NEIGHTBOURHOOD : 0.0f) << ';' << ((*gene)->group != NULL ? (double)(*gene)->group->CountUniqueStrains() / (double)db->strains.size() : -0.01f / (double)db->strains.size()) << ';' << (*gene)->coverage << ';' << (*gene)->GC3s << ';' << (*gene)->type << '\n';
            out1 << (*contig)->id << ";"<< (*gene)->start + ((*gene)->length / 2) << ';' << 0.0f << ';' << ((*gene)->group != NULL ? (double)(*gene)->group->CountUniqueStrains() / (double)db->strains.size() : -0.01f / (double)db->strains.size()) << ';' << (*gene)->coverage << ';' << (*gene)->GC3s << ';' << (*gene)->type << '\n';
    }
    out1.close();

    unordered_map<Group *, set<Gene *>> links;

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            if ((*gene)->group != NULL)
                links[(*gene)->group].insert(*gene);

    ofstream out2( OUTPUT_DIRECTORY + id_full + "_links.csv", ifstream::out );
    out2 << "ID1;Location1;ID2;Location2\n";

    for (auto link = links.begin(); link != links.end(); link++)
        if (link->second.size() > 1)
            for (auto g1 = link->second.begin(); g1 != link->second.end(); g1++)
                for (auto g2 = g1; g2 != link->second.end(); g2++)
                    if (*g1 != *g2)
                        out2 << (*g1)->contig->id << ";" << (*g1)->start + ((*g1)->length / 2)  << ";" << (*g2)->contig->id << ";" << (*g2)->start + ((*g2)->length / 2) << '\n';

    out2.close();
    /*
    ofstream out3( OUTPUT_DIRECTORY + id_full + "_colours.csv", ifstream::out );
    out3 << "ID;Value\n";
    
    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        out3 << (*contig)->id << ";";

        if ((*contig)->id_full.find("chromosome") != string::npos)
            out3 << "palegreen\n";
        else if ((*contig)->id_full.find("plasmid") != string::npos)
            out3 << "skyblue\n";
        else if ((*contig)->id_full.find("fragment") != string::npos)
            out3 << "red\n";
        else
            out3 << "white\n";
    }
    out3.close();
     */
}

/*
void Strain::Parse(string INPUT_DIRECTORY, unordered_map<string, Gene *> *genepool)
{
    string type;
    
    //FILE test = FILE(INPUT_DIRECTORY + strainID + "." + type).exis
    
    ifstream input;
    
    if ((input = ifstream( INPUT_DIRECTORY + this->id + "." +  TYPE_GFF)).good())
    {
        type = TYPE_GFF;
    }
    else if ((input = ifstream( INPUT_DIRECTORY + this->id + "." +  TYPE_TAB)).good())
    {
        type = TYPE_TAB;
    }
    else
    {
        printf("Valid input file type not found; valid file types: \"TAB\" or \"GFF\"\n");
        printf("input directory was: %s\n", INPUT_DIRECTORY.c_str());
        printf("strainID was: %s\n", this->id.c_str());
        exit(1);
    }
    
    string line;
    while (getline(input, line))
    {
        int counter = 0;
        
        if (type == TYPE_GFF)
        {
            if (line[0] == '#')
                continue;
            
            stringstream sstream(line);
            string s;
            
            vector<string> items;
            items.clear();
            
            while (getline(sstream, s, '\t'))
            {
                items.push_back(s);
                counter++;
                if (counter >= 9)
                    break;
            }
            
            if (items[2] == "CDS")
            {
                Gene *gene = new Gene;
                
                size_t loc1 = items[8].find_first_of("ID=");
                size_t loc2 = items[8].find_first_of(';', loc1);
                
                gene->id = items[8].substr(loc1+3, loc2-3);
                gene->start = stoi(items[3]);
                gene->end = stoi(items[4]);
                gene->orientation = (items[6] == "-");
                //gene->strain = strainID;
                gene->strain = this;
                gene->left_neighbour = NULL;
                gene->right_neighbour = NULL;
                gene->group = NULL;
                gene->votes = 0;
                gene->antiVotes = 0;
                gene->relations.clear();
                gene->scoresMatrix.clear();
                gene->excluded = false;
                //gene->neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
                //strain->genemap[gene->id] = gene;
                this->genes.push_back(gene);
                (*genepool)[ this->id + "|" + gene->id] = gene;
            }
        }
        else if (type == TYPE_TAB)
        {
            if (line.empty())
                continue;
            
            stringstream sstream(line);
            string s;
            
            vector<string> items;
            items.clear();
            
            while (getline(sstream, s, '\t'))
            {
                items.push_back(s);
                counter++;
                if (counter >= 4)
                    break;
            }
            
            Gene *gene = new Gene;
            
            //strain->genes.resize(strain->genes.size()+1);
            //gene = &strain->genes.back();
            
            gene->id = items[0];
            gene->orientation = (items[1] == "1");
            gene->start = stoi(items[2]);
            gene->end = stoi(items[3]);
            //gene->strain = strainID;
            gene->strain = this;
            gene->left_neighbour = NULL;
            gene->right_neighbour = NULL;
            gene->group = NULL;
            gene->votes = 0;
            gene->antiVotes = 0;
            gene->scoresMatrix.clear();
            gene->excluded = false;
            //gene->neighbours.resize(SIZE_OF_NEIGHTBOURHOOD);
            //strain->genemap[gene->id] = gene;
            this->genes.push_back(gene);
            (*genepool)[this->id + "|" + gene->id] = gene;
        }
    }
    
    Gene *previousGene = this->genes.back();
    
    for (vector<Gene*>::iterator gene = this->genes.begin(); gene !=  this->genes.end(); gene++)
    {
        previousGene->right_neighbour = *gene.base();
        (*gene)->left_neighbour = previousGene;
        previousGene = *gene.base();
    }
    input.close();
}
*/
