//
//  Dataset.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 11/09/2017.
//
//

#include "Dataset.hpp"

Dataset::Dataset(Database *db)
{
    database = db;
    genecount = 0;
    syntenyScoreSophisticated = NAN;
    syntenyScoreOld = NAN;
    syntenyScoreAdjusted = NAN;
    syntenyScoreTest = NAN;
    syntenyScoreFast = NAN;
}

void Dataset::ScoreSynteny(string OUTPUT_FILE, Dataset *dataset2)
{
    ofstream out( OUTPUT_FILE, ifstream::out );

    int counter = 0;
    int groups1to1 = 0;
    int groups1to2 = 0;
    double scoresumSophisticated = 0.0f;
    double scoresumOld = 0.0f;
    double scoresumSimple = 0.0f;
    double scoresumAdjusted = 0.0f;
    //double scoresumTest = 0.0f;
    double scoresumFast = 0.0f;
    double scoresum1to1Sophisticated = 0.0f;
    double scoresum1to2Sophisticated = 0.0f;
    double scoresum1to1Old = 0.0f;
    double scoresum1to2Old = 0.0f;
    double scoresum1to1Simple = 0.0f;
    double scoresum1to2Simple = 0.0f;
    double scoresum1to1Adjusted = 0.0f;
    double scoresum1to2Adjusted = 0.0f;
    //double scoresum1to1Test = 0.0f;
    //double scoresum1to2Test = 0.0f;
    double scoresum1to1Fast = 0.0f;
    double scoresum1to2Fast = 0.0f;
    int genes1to1 = 0;
    int genes1to2 = 0;
    uintmax_t nucleotides1to1 = 0;
    uintmax_t nucleotides1to2 = 0;

    Progress p(groups.size()-1);
    printf("Scoring Syntenies for:\t%s\n", OUTPUT_FILE.c_str());
    //out << "Group\tGenes\tSynteny Score Sophisticated\tSynteny Score Old\tSynteny Score Simple\tSynteny Score Adjusted\tSynteny Score Fast\tSynteny Score Test\n";
    //out << "Group|Connectivity\tGenes\tSynteny Score Sophisticated\tSynteny Score Old\tSynteny Score Simple\tSynteny Score Adjusted\tSynteny Score Fast\n";
    out << "Group\tConnectivity\tGenes\tStrains\tOrthologs\tParalogs\tChromosome Genes\tPlasmid Genes\tFragment Genes\tProtein Products\tSynteny Score Old\tSynteny Score Adjusted\tSynteny Score Fast\tSynteny Score Simple\tSynteny Score Sophisticated\n";
    for (auto group = groups.begin(); group != groups.end(); group++)
    {
        set<Group*> groupset;
        if ((*group)->genes.size() > 1)
        {
            scoresumSophisticated += (*group)->SyntenizeSophisticated();
            scoresumOld += (*group)->SyntenizeOld();
            scoresumSimple += (*group)->SyntenizeSimple();
            scoresumAdjusted += (*group)->SyntenizeAdjusted();
            //scoresumTest += (*group)->SyntenizeTest();
            scoresumFast += (*group)->SyntenizeFast();
            counter++;
            //scoresumSophisticated += (*group)->SyntenizeSophisticated() / (double)(*group)->genes.size();
            //scoresumNaive += (*group)->SyntenizeNaive() / (double)(*group)->genes.size();
            //scoresumSimple += (*group)->SyntenizeSimple() / (double)(*group)->genes.size();

            /*
            out << (*group)->id.c_str();
            if (!std::isnan((*group)->algebraicConnectivity))
                out << '|' << (*group)->algebraicConnectivity;
            out << "\t" << (*group)->genes.size() << "\t" << (*group)->syntenyScoreSophisticated << "\t" << (*group)->syntenyScoreOld << "\t" << (*group)->syntenyScoreSimple << "\t" << (*group)->syntenyScoreAdjusted << "\t" << (*group)->syntenyScoreFast << "\n";
            */
            
            set<string> types;
            int chromosome = 0;
            int plasmid = 0;
            int fragment = 0;

            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            {
                if ((*gene)->contig->id_full.find("chromosome") != string::npos)
                    chromosome++;
                if ((*gene)->contig->id_full.find("plasmid") != string::npos)
                    plasmid++;
                if ((*gene)->contig->id_full.find("fragment") != string::npos)
                    fragment++;

                types.insert((*gene)->product);
            }

            out << (*group)->id.c_str() << "\t";
            out << (*group)->algebraicConnectivity << "\t";
            out << (*group)->genes.size() << "\t";
            out << (*group)->CountUniqueStrains() << "\t";
            out << (*group)->CountOrthologs() << "\t";
            out << (*group)->CountParalogs() << "\t";

            out << chromosome << "\t";
            out << plasmid << "\t";
            out << fragment << "\t";

            int counter = 1;
            for (auto type = types.begin(); type != types.end(); type++, counter++)
            {
                out << (*type);
                if (counter < types.size())
                    out << ';';
                else
                    out << '\t';
            }

            out << (*group)->syntenyScoreOld << "\t";
            out << (*group)->syntenyScoreAdjusted << "\t";
            out << (*group)->syntenyScoreFast << "\t";
            out << (*group)->syntenyScoreSimple << "\t";
            out << (*group)->syntenyScoreSophisticated << "\n";

            p.Update(distance(groups.begin(), group));
            if (dataset2 != NULL)
            {
                int nucleotides = 0;
                for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                {
                    nucleotides += (*gene)->length;
                    try
                    {
                        Group *g = dataset2->grouppool.at(*gene);
                        groupset.insert(g);
                    }
                    catch (const out_of_range & e){}
                }
                if (groupset.size() == 0)
                    continue;
                if (groupset.size() == 1)
                {
                    groups1to1++;
                    scoresum1to1Sophisticated += (*group)->syntenyScoreSophisticated;
                    scoresum1to1Old += (*group)->syntenyScoreOld;
                    scoresum1to1Simple += (*group)->syntenyScoreSimple;
                    scoresum1to1Adjusted += (*group)->syntenyScoreAdjusted;
                    //scoresum1to1Test += (*group)->syntenyScoreTest;
                    scoresum1to1Fast += (*group)->syntenyScoreFast;
                    genes1to1 += (*group)->genes.size();
                    nucleotides1to1 += nucleotides;
                }
                else
                {
                    groups1to2++;
                    scoresum1to2Sophisticated += (*group)->syntenyScoreSophisticated;
                    scoresum1to2Old += (*group)->syntenyScoreOld;
                    scoresum1to2Simple += (*group)->syntenyScoreSimple;
                    scoresum1to2Adjusted += (*group)->syntenyScoreAdjusted;
                    //scoresum1to2Test += (*group)->syntenyScoreTest;
                    scoresum1to2Fast += (*group)->syntenyScoreFast;
                    genes1to2 += (*group)->genes.size();
                    nucleotides1to2 += nucleotides;
                }
            }
        }
    }
//Score:Genes
    //6:2
    //7:5
    //10:4
    //3:8
    
    //((6 + 7 + 10 +3)/(2 + 5 + 4 +8))-(((6/2)+(7/5)+(10/4)+(3/8))/4)
    syntenyScoreSophisticated = scoresumSophisticated / (double)counter;
    syntenyScoreOld = scoresumOld / (double)counter;
    syntenyScoreSimple = scoresumSimple / (double)counter;
    syntenyScoreAdjusted = scoresumAdjusted / (double)counter;
    //syntenyScoreTest = scoresumTest / (double)counter;
    syntenyScoreFast = scoresumFast / (double)counter;
    //syntenyScoreSophisticated = scoresumSophisticated / (double)groups.size();
    //syntenyScoreNaive = scoresumNaive / (double)groups.size();
    //syntenyScoreSimple = scoresumSimple / (double)groups.size();
    printf("\nAverage Synteny Score Sophisticated:\t[%.6f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreSophisticated, scoresum1to1Sophisticated / (double) groups1to1, scoresum1to2Sophisticated / (double) groups1to2);
    printf("\nAverage Synteny Score Old:\t[%.5f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreOld, scoresum1to1Old / (double) groups1to1, scoresum1to2Old / (double) groups1to2);
    printf("\nAverage Synteny Score Simple:\t[%.5f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreSimple, scoresum1to1Simple / (double) groups1to1, scoresum1to2Simple / (double) groups1to2);
    printf("\nAverage Synteny Score Adjusted:\t[%.5f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreAdjusted, scoresum1to1Adjusted / (double) groups1to1, scoresum1to2Adjusted / (double) groups1to2);
    printf("\nAverage Synteny Score Fast:\t[%.5f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreFast, scoresum1to1Fast / (double) groups1to1, scoresum1to2Fast / (double) groups1to2);
    //printf("\nAverage Synteny Score Test:\t[%.5f]\t(1:1)[%.6f]\t(1:!1)[%.6f]", syntenyScoreTest, scoresum1to1Test / (double) groups1to1, scoresum1to2Test / (double) groups1to2);
    /*
    out << "\n\nAverage Synteny Score Sophisticated\t" << syntenyScoreSophisticated;
    out << "\nAverage 1:1 Synteny Score Sophisticated\t" << (double)scoresum1to1Sophisticated / (double) groups1to1;
    out << "\nAverage 1:!1 Synteny Score Sophisticated\t" << (double)scoresum1to2Sophisticated / (double) groups1to2;
    
    out << "\n\nAverage Synteny Score Old\t" << syntenyScoreOld;
    out << "\nAverage 1:1 Synteny Score Old\t" << (double)scoresum1to1Old / (double) groups1to1;
    out << "\nAverage 1:!1 Synteny Score Old\t" << (double)scoresum1to2Old / (double) groups1to2;
    
    out << "\n\nAverage Synteny Score Simple\t" << syntenyScoreSimple;
    out << "\nAverage 1:1 Synteny Score Simple\t" << (double)scoresum1to1Simple / (double) groups1to1;
    out << "\nAverage 1:!1 Synteny Score Simple\t" << (double)scoresum1to2Simple / (double) groups1to2;

    out << "\n\nAverage Synteny Score Adjusted\t" << syntenyScoreAdjusted;
    out << "\nAverage 1:1 Synteny Score Adjusted\t" << (double)scoresum1to1Adjusted / (double) groups1to1;
    out << "\nAverage 1:!1 Synteny Score Adjusted\t" << (double)scoresum1to2Adjusted / (double) groups1to2;

    out << "\n\nAverage Synteny Score Fast\t" << syntenyScoreFast;
    out << "\nAverage 1:1 Synteny Score Fast\t" << (double)scoresum1to1Fast / (double) groups1to1;
    out << "\nAverage 1:!1 Synteny Score Fast\t" << (double)scoresum1to2Fast / (double) groups1to2;

    //out << "\n\nAverage Synteny Score Test\t" << syntenyScoreTest;
    //out << "\nAverage 1:1 Synteny Score Test\t" << (double)scoresum1to1Test / (double) groups1to1;
    //out << "\nAverage 1:!1 Synteny Score Test\t" << (double)scoresum1to2Test / (double) groups1to2;

    out << "\n\nAverage 1:1 nucleotides per gene\t" << (double)nucleotides1to1 / (double) genes1to1;
    out << "\nAverage 1:!1 nucleotides per gene\t" << (double)nucleotides1to2 / (double) genes1to2;

    out << "\n\nAverage 1:1 genes per group\t" << (double)genes1to1 / (double) groups1to1;
    out << "\nAverage 1:!1 genes per group\t" << (double)genes1to2 / (double) groups1to2;

    out << "\n";
    */
    out.close();
}

void Dataset::ExportFNA(string OUTPUT_DIRECTORY)
{
    Timer timer;
    timer.Start();
    printf("Exporting fna groups to directory %s:\n", OUTPUT_DIRECTORY.c_str());
    for (auto group = groups.begin(); group != groups.end(); group++)
    {
        ofstream out( OUTPUT_DIRECTORY + (*group)->id + ".fna", ifstream::out );
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            out << '>' << (*group)->id << '|' << (*gene)->strain->id_full << '|' << (*gene)->id_full << '\n';
            out << (*gene)->sequence << '\n';
        }
        out.close();
    }
    printf("Exporting %lu fna groups took:\t%i seconds\n\n", groups.size(), (int)timer.Elapsed());
}

void Dataset::PrintStats()
{
    //database->strains
    //unordered_map<int, Strain *>
    for (int strainsize = database->strains.size(); strainsize > 1; strainsize-- )
    {
        vector<Strain *> strains;
        for (auto group = groups.begin() ; group < groups.end(); group++ )
        {
            set<Strain *> strainset;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                strainset.insert((*gene)->strain);
            
            if (strainset.size() != strainsize)
                continue;

            vector<Strain *> missing;
        }
    }
}
