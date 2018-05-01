//
//  Score.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 19/09/2017.
//
//

#include "Score.hpp"
/*
double Score::Synteny(Gene *gene1, Gene *gene2, Type type )
{
    double score = NAN;

    switch (type)
    {
        case FAST:
        {
            score = Score::Fast(gene1, gene2);
        }
        break;
        
        case SIMPLE:
        {
            score = Score::Simple(gene1, gene2);
        }
        break;
        
        case NAIVE:
        {
            score = Score::Naive(gene1, gene2);
        }
        break;
            
        case ADJUSTED:
        {
            score = Score::Adjusted(gene1, gene2);
        }
        break;
        
        case OLD:
        {
            score = Score::Old(gene1, gene2);
        }
        break;

        case SOPHISTICATED:
        {
            //gene1->group->
            double s1 = Score::Sophisticated(gene1, gene2);
            double s2 = Score::SophisticatedMirrored(gene1, gene2) * 0.75;

            score = (s1 > s2) ? s1 : s2;
        }
        break;
    }

    return score;
}
*/
int Score::Fast(Gene *gene1, Gene *gene2)
{
    int score = 0;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2)-1; i >= 0; i--)
        if (gene1->neighbours[i] == NULL || gene2->neighbours[i] == NULL || gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL || gene1->neighbours[i]->group != gene2->neighbours[i]->group)
            break;
        else
            score++;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2); i < SIZE_OF_NEIGHTBOURHOOD; i++)
        if (gene1->neighbours[i] == NULL || gene2->neighbours[i] == NULL || gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL || gene1->neighbours[i]->group != gene2->neighbours[i]->group)
            break;
        else
            score++;
    
    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    double     mirrorscore = 0.0f;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    try
    {
        score = gene1->scores.at(gene2);
        mirrorscore = gene1->mirrorscores[gene2];
        
        return (score > mirrorscore*0.75) ? score : mirrorscore*0.75;
    }
    catch (exception e){}

    try
    {
        score = gene2->scores.at(gene1);
        mirrorscore = gene2->mirrorscores[gene1];
        
        return (score > mirrorscore*0.75) ? score : mirrorscore*0.75;
    }
    catch (exception e){}

    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[ j-1 ] != NULL && gene2->neighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[ j-1 ] != NULL && gene2->neighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
    }
    
    mirrorscore = SophisticatedMirrored(gene1, gene2);

    gene1->scores[gene2] = score;
    gene1->mirrorscores[gene2] = mirrorscore;
    //gene2->scores[gene1] = score;
    //gene2->mirrorscores[gene1] = mirrorscore;

    /*
    uint64_t g1 = (uint64_t) gene1;
    uint64_t g2 = (uint64_t) gene2;
    //uint
    //uint64_t key = (gene1 > gene2) ? gene1 | : gene2;

    Score::scorepool[(g1 > g2) ? g1: g2] = score;
    Score::mirrorscorepool[(g1 > g2) ? g1: g2] = Score::SophisticatedMirrored(gene1, gene2);

    //scorepool[make_pair((gene1 > gene2) ? gene1: gene2 , (gene1 > gene2) ? gene2: gene1)] = score;
    //mirrorscorepool[make_pair((gene1 > gene2) ? gene1: gene2 , (gene1 > gene2) ? gene2: gene1)] = Score::SophisticatedMirrored(gene1, gene2);

    //double mirroredscore = SophisticatedMirrored(gene1, gene2) * 0.75;
    //score = (score > mirroredscore) ? score : mirroredscore;
     */
    //return score;
   return (score > mirrorscore*0.75) ? score : mirrorscore*0.75;
}

double Score::SophisticatedMirrored(Gene *gene1, Gene *gene2)
{
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    Gene* mirroredNeighbours[SIZE_OF_NEIGHTBOURHOOD];
    
    for (int g = 0, m = SIZE_OF_NEIGHTBOURHOOD-1; g < SIZE_OF_NEIGHTBOURHOOD; g++, m-- )
        mirroredNeighbours[m] = gene2->neighbours[g];
    
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[ j-1 ] != NULL && mirroredNeighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && mirroredNeighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == mirroredNeighbours[ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[ j-1 ] != NULL && mirroredNeighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && mirroredNeighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == mirroredNeighbours[ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
    }
    return score;
}

double Score::Naive(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
                score++;
    
    return score;
}

double Score::Simple(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = false;
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
            {
                found = true;
                break;
            }
        }
        if (found)
            score++;
    }
    
    return score;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
            {
                partialmatch = 1.0f / (double)(1 + abs(i - j));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        
        double distance = (i < (SIZE_OF_NEIGHTBOURHOOD / 2)) ? (SIZE_OF_NEIGHTBOURHOOD / 2) - i : i+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
    }
    
    return matches;
}

double Score::Old(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        if (gene1->neighbours[i] != NULL &&  gene2->neighbours[i] != NULL && gene1->neighbours[i]->group != NULL &&  gene2->neighbours[i]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[i]->group)
            matches++;
        else
        {
            double maxpartialmatch = 0.0;
            double partialmatch;
            
            for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            {
                if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
                {
                    partialmatch = 1.0f / double(abs(i-j)+1.0f);
                    
                    if (partialmatch > maxpartialmatch)
                        maxpartialmatch = partialmatch;
                }
            }
            
            matches += maxpartialmatch;
        }
    }
    return matches;
}

double Score::Test(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = rand() % 2;
        //for (int j = 0; j < gene2->neighbours.size(); j++)
        // {
        // if (true)
        // break;
        // }
        if (found)
            score++;
    }
    
    return score;
}
/*
static void StoreScore(Gene *gene1, Gene *gene2, double score)
{
    Score::scorepool[make_pair((gene1 > gene2) ? gene1: gene2 , (gene1 > gene2) ? gene2: gene1)] = score;
}

static double FetchScore(Gene *gene1, Gene *gene2)
{
    double score = NAN;
    
    try
    {
        score = Score::scorepool.at(make_pair((gene1 > gene2) ? gene1: gene2 , (gene1 > gene2) ? gene2: gene1));
    }
    catch (const out_of_range & e){}

    return score;
}
*/
//--------------------------------------------------------------------------------------------------------------------------------------------
/*
int Score::Fast(Gene *gene1, Gene *gene2)
{
    int score1 = Score::Fast(gene1, gene2, false);
    int score2 = Score::Fast(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

int Score::Fast(Gene *gene1, Gene *gene2, bool mirror)
{
    int score = 0;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2)-1; i >= 0; i--)
        if (gene1->neighbours[0][i] == NULL || gene2->neighbours[mirror][i] == NULL || gene1->neighbours[0][i]->group == NULL || gene2->neighbours[mirror][i]->group == NULL || gene1->neighbours[0][i]->group != gene2->neighbours[mirror][i]->group)
            break;
        else
            score++;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2); i < SIZE_OF_NEIGHTBOURHOOD; i++)
        if (gene1->neighbours[0][i] == NULL || gene2->neighbours[mirror][i] == NULL || gene1->neighbours[0][i]->group == NULL || gene2->neighbours[mirror][i]->group == NULL || gene1->neighbours[0][i]->group != gene2->neighbours[mirror][i]->group)
            break;
        else
            score++;
    
    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    double score1 = Score::Sophisticated(gene1, gene2, false);
    double score2 = Score::Sophisticated(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2, bool mirror)
{
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[0][ j-1 ] != NULL && gene2->neighbours[mirror][ i-1 ] != NULL && gene1->neighbours[0][ j-1 ]->group != NULL && gene2->neighbours[mirror][ i-1 ]->group != NULL && gene1->neighbours[0][ j-1 ]->group == gene2->neighbours[mirror][ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[0][ j-1 ] != NULL && gene2->neighbours[mirror][ i-1 ] != NULL && gene1->neighbours[0][ j-1 ]->group != NULL && gene2->neighbours[mirror][ i-1 ]->group != NULL && gene1->neighbours[0][ j-1 ]->group == gene2->neighbours[mirror][ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
    }
    return score;
}

double Score::Naive(Gene *gene1, Gene *gene2)
{
    double score1 = Score::Naive(gene1, gene2, false);
    double score2 = Score::Naive(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

double Score::Naive(Gene *gene1, Gene *gene2, bool mirror)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            if (gene1->neighbours[0][i] != NULL && gene2->neighbours[mirror][j] != NULL && gene1->neighbours[0][i]->group != NULL && gene2->neighbours[mirror][j]->group != NULL && gene1->neighbours[0][i]->group == gene2->neighbours[mirror][j]->group)
                score++;
    
    return score;
}

double Score::Simple(Gene *gene1, Gene *gene2)
{
    double score1 = Score::Simple(gene1, gene2, false);
    double score2 = Score::Simple(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

double Score::Simple(Gene *gene1, Gene *gene2, bool mirror)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = false;
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[0][i] != NULL && gene2->neighbours[mirror][j] != NULL && gene1->neighbours[0][i]->group != NULL && gene2->neighbours[mirror][j]->group != NULL && gene1->neighbours[0][i]->group == gene2->neighbours[mirror][j]->group)
            {
                found = true;
                break;
            }
        }
        if (found)
            score++;
    }
    
    return score;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    double score1 = Score::Adjusted(gene1, gene2, false);
    double score2 = Score::Adjusted(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

double Score::Adjusted(Gene *gene1, Gene *gene2, bool mirror)
{
    int i, j;
    double matches = 0.0f;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[0][i] != NULL && gene2->neighbours[mirror][j] != NULL && gene1->neighbours[0][i]->group != NULL && gene2->neighbours[mirror][j]->group != NULL && gene1->neighbours[0][i]->group == gene2->neighbours[mirror][j]->group)
            {
                partialmatch = 1.0f / (double)(1 + abs(i - j));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        
        double distance = (i < (SIZE_OF_NEIGHTBOURHOOD / 2)) ? (SIZE_OF_NEIGHTBOURHOOD / 2) - i : i+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
    }
    
    return matches;
}

double Score::Old(Gene *gene1, Gene *gene2)
{
    double score1 = Score::Old(gene1, gene2, false);
    double score2 = Score::Old(gene1, gene2, true);
    return score1 > score2 ? score1 : score2;
}

double Score::Old(Gene *gene1, Gene *gene2, bool mirror)
{
    int i, j;
    double matches = 0.0f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        if (gene1->neighbours[0][i] != NULL &&  gene2->neighbours[mirror][i] != NULL && gene1->neighbours[0][i]->group != NULL &&  gene2->neighbours[mirror][i]->group != NULL && gene1->neighbours[0][i]->group == gene2->neighbours[mirror][i]->group)
            matches++;
        else
        {
            double maxpartialmatch = 0.0;
            double partialmatch;
            
            for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            {
                if (gene1->neighbours[0][i] != NULL && gene2->neighbours[mirror][j] != NULL && gene1->neighbours[0][i]->group != NULL && gene2->neighbours[mirror][j]->group != NULL && gene1->neighbours[0][i]->group == gene2->neighbours[mirror][j]->group)
                {
                    partialmatch = 1.0f / double(abs(i-j)+1.0f);
                    
                    if (partialmatch > maxpartialmatch)
                        maxpartialmatch = partialmatch;
                }
            }
            
            matches += maxpartialmatch;
        }
    }
    return matches;
}

double Score::Test(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = rand() % 2;
        //for (int j = 0; j < gene2->neighbours.size(); j++)
        // {
        // if (true)
        // break;
        // }
        if (found)
            score++;
    }
    
    return score;
}
*/
//--------------------------------------------------------------------------------------------------------------------------------------------
/*
int Score::Fast(Gene *gene1, Gene *gene2)
{
    int score = 0;

    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2)-1; i >= 0; i--)
        if (gene1->neighbours[gene1->orientation][i] == NULL || gene2->neighbours[gene2->orientation][i] == NULL || gene1->neighbours[gene1->orientation][i]->group == NULL || gene2->neighbours[gene2->orientation][i]->group == NULL || gene1->neighbours[gene1->orientation][i]->group != gene2->neighbours[gene2->orientation][i]->group)
            break;
        else
            score++;

    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2); i < SIZE_OF_NEIGHTBOURHOOD; i++)
        if (gene1->neighbours[gene1->orientation][i] == NULL || gene2->neighbours[gene2->orientation][i] == NULL || gene1->neighbours[gene1->orientation][i]->group == NULL || gene2->neighbours[gene2->orientation][i]->group == NULL || gene1->neighbours[gene1->orientation][i]->group != gene2->neighbours[gene1->orientation][i]->group)
            break;
        else
            score++;

    //for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    //    if (gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL)
    //        continue;
    //    else if (gene1->neighbours[i]->group == gene2->neighbours[i]->group)
    //        score++;
    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[gene1->orientation][ j-1 ] != NULL && gene2->neighbours[gene2->orientation][ i-1 ] != NULL && gene1->neighbours[gene1->orientation][ j-1 ]->group != NULL && gene2->neighbours[gene2->orientation][ i-1 ]->group != NULL && gene1->neighbours[gene1->orientation][ j-1 ]->group == gene2->neighbours[gene2->orientation][ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[gene1->orientation][ j-1 ] != NULL && gene2->neighbours[gene2->orientation][ i-1 ] != NULL && gene1->neighbours[gene1->orientation][ j-1 ]->group != NULL && gene2->neighbours[gene2->orientation][ i-1 ]->group != NULL && gene1->neighbours[gene1->orientation][ j-1 ]->group == gene2->neighbours[gene2->orientation][ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
        //k++ ;
    }
    return score;
}

double Score::Naive(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            if (gene1->neighbours[gene1->orientation][i] != NULL && gene2->neighbours[gene2->orientation][j] != NULL && gene1->neighbours[gene1->orientation][i]->group != NULL && gene2->neighbours[gene2->orientation][j]->group != NULL && gene1->neighbours[gene1->orientation][i]->group == gene2->neighbours[gene2->orientation][j]->group)
                score++;
    
    return score;
}

double Score::Simple(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = false;
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[gene1->orientation][i] != NULL && gene2->neighbours[gene2->orientation][j] != NULL && gene1->neighbours[gene1->orientation][i]->group != NULL && gene2->neighbours[gene2->orientation][j]->group != NULL && gene1->neighbours[gene1->orientation][i]->group == gene2->neighbours[gene2->orientation][j]->group)
            {
                found = true;
                break;
            }
        }
        if (found)
            score++;
    }
    
    return score;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[gene1->orientation][i] != NULL && gene2->neighbours[gene2->orientation][j] != NULL && gene1->neighbours[gene1->orientation][i]->group != NULL && gene2->neighbours[gene2->orientation][j]->group != NULL && gene1->neighbours[gene1->orientation][i]->group == gene2->neighbours[gene2->orientation][j]->group)
            {
                partialmatch = 1.0f / (double)(1 + abs(i - j));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        
        double distance = (i < (SIZE_OF_NEIGHTBOURHOOD / 2)) ? (SIZE_OF_NEIGHTBOURHOOD / 2) - i : i+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
    }
    
    return matches;
}

double Score::Old(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        if (gene1->neighbours[gene1->orientation][i] != NULL &&  gene2->neighbours[gene2->orientation][i] != NULL && gene1->neighbours[gene1->orientation][i]->group != NULL &&  gene2->neighbours[gene2->orientation][i]->group != NULL && gene1->neighbours[gene1->orientation][i]->group == gene2->neighbours[gene2->orientation][i]->group)
            matches++;
        else
        {
            double maxpartialmatch = 0.0;
            double partialmatch;
            
            for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            {
                if (gene1->neighbours[gene1->orientation][i] != NULL && gene2->neighbours[gene2->orientation][j] != NULL && gene1->neighbours[gene1->orientation][i]->group != NULL && gene2->neighbours[gene2->orientation][j]->group != NULL && gene1->neighbours[gene1->orientation][i]->group == gene2->neighbours[gene2->orientation][j]->group)
                {
                    partialmatch = 1.0f / double(abs(i-j)+1.0f);
                    
                    if (partialmatch > maxpartialmatch)
                        maxpartialmatch = partialmatch;
                }
            }
            
            matches += maxpartialmatch;
        }
    }
    return matches;
}

double Score::Test(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = rand() % 2;
        //for (int j = 0; j < gene2->neighbours.size(); j++)
        // {
        // if (true)
        // break;
        // }
        if (found)
            score++;
    }
    
    return score;
}
*/
/*
int Score::Fast(Gene *gene1, Gene *gene2)
{
    bool flipped = gene1->orientation != gene2->orientation;
    int score = 0;

    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
            if (gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL)
                continue;
            else if (!flipped && gene1->neighbours[i]->group == gene2->neighbours[i]->group)
                score++;
            else if (flipped && gene1->neighbours[i]->group == gene2->neighbours[SIZE_OF_NEIGHTBOURHOOD-1 - i]->group)
                score++;

    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    //int        k = 0;
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0, jj;
    double     score = 0.0f;
    //string     seq_1_al, seq_2_al;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            if (gene1->orientation == gene2->orientation)
                jj = j-1;
            else
                jj = (SIZE_OF_NEIGHTBOURHOOD-1) - (j - 1);

            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[ jj ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ jj ]->group == gene2->neighbours[ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        if (gene1->orientation == gene2->orientation)
            jj = j-1;
        else
            jj = (SIZE_OF_NEIGHTBOURHOOD-1) - (j - 1);
        
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[ jj ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ jj ]->group == gene2->neighbours[ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
        //k++ ;
    }
    return score;
}
*/
//--------------------------------------------------------------------------------------------------------------------------------------------
/*
double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    double score1 = Sophisticated(gene1->neighbours, gene2->neighbours, false);
    double score2 = Sophisticated(gene1->neighbours, gene2->neighbours, true);

    return(score1 > score2 ? score1 : score2);
}

double Score::Sophisticated(vector<Gene*>&seq_1, vector<Gene*>&seq_2, bool mirror)
{
    int        k = 0;
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    //string     seq_1_al, seq_2_al;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    if (mirror)
        reverse( seq_1.begin(), seq_1.end() );
    
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (seq_1[ j-1 ]->group != NULL && seq_2[ i-1 ]->group != NULL && seq_1[ j-1 ]->group == seq_2[ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                //seq_1_al += "\t-\t\t" ;
                //seq_2_al += seq_2[ i-1 ]->group->id + "\t" ;
                //score -= d;
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                //seq_1_al += seq_1[ j-1 ]->group->id + "\t"  ;
                //seq_2_al += seq_2[ i-1 ]->group->id + "\t"  ;
                //score += (seq_1[ j-1 ]->group == seq_2[ i-1 ]->group) ? 2 : -1;
                //score += (seq_1[ j-1 ]->group == seq_2[ i-1 ]->group) ? 1 : -0.25f;
                score += (seq_1[ j-1 ]->group != NULL && seq_2[ i-1 ]->group != NULL && seq_1[ j-1 ]->group == seq_2[ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                //seq_1_al += seq_1[ j-1 ]->group->id + "\t"  ;
                //seq_2_al += "\t-\t\t" ;
                //score -= d;
                score -= 0.5f;
                j-- ;
        }
        k++ ;
    }
    
    //cout << seq_1_al << "\n";
    //cout << seq_2_al << "\n";
    
    //return  score/2.0f;
    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2, bool mirror)
{
    double fU, fD, fL ;
    char ptr;
    int i = 0, j = 0;
    double score = 0;
    double F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    //vector<Gene*> seq_1_al, seq_2_al;
    vector<string> seq_1_al, seq_2_al;

    seq_1_al.clear();
    seq_2_al.clear();

    F[ 0 ][ 0 ] =  0;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  j * GAP_PENALTY ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  i * GAP_PENALTY ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] + GAP_PENALTY;
            if (mirror)
                fD = F[ i-1 ][ j-1 ] + ((gene1->neighbours[SIZE_OF_NEIGHTBOURHOOD - j]->group == gene2->neighbours[ i-1 ]->group) ? MATCH_SCORE :GAP_PENALTY);
            else
                fD = F[ i-1 ][ j-1 ] + ((gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group) ? MATCH_SCORE : GAP_PENALTY);
            
            fL = F[ i ][ j-1 ] + GAP_PENALTY;
            F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;
            
            traceback[ i ][ j ] =  ptr ;
        }
    }
    i-- ; j-- ;
    int counter = 0;
    while( i > 0 || j > 0 )
    {
        if (i < 0 || j < 0 )
            printf("!\n");
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                //if (mirror)
                //    seq_1_al.insert(seq_1_al.begin(), NULL);
                //else
                //    seq_1_al.push_back(NULL);
                
                seq_1_al.push_back("groupNULL");
                seq_2_al.push_back(gene2->neighbours[ i-1 ]->group->id) ;
                i-- ;
                score += GAP_SCORE;
                counter++;
                break ;
                
            case '\\':
                if (mirror)
                {
                    //seq_1_al.insert(seq_1_al.begin(), gene1->neighbours[ i-1 ]);
                    seq_1_al.push_back(gene1->neighbours[ SIZE_OF_NEIGHTBOURHOOD - j]->group->id);
                    score += (gene1->neighbours[ SIZE_OF_NEIGHTBOURHOOD-j ]->group == gene2->neighbours[ i-1 ]->group) ? MATCH_SCORE : GAP_SCORE;
                }
                else
                {
                    seq_1_al.push_back(gene1->neighbours[ j-1 ]->group->id);
                    score += (gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group) ? MATCH_SCORE : GAP_SCORE;
                }
     
                seq_2_al.push_back(gene2->neighbours[ i-1 ]->group->id) ;
                i-- ;  j-- ;
                counter++;
                break ;
                
            case '-' :
                if (mirror)
                    //seq_1_al.insert(seq_1_al.begin(), gene1->neighbours[ i-1 ]) ;
                    seq_1_al.push_back(gene1->neighbours[ SIZE_OF_NEIGHTBOURHOOD - j]->group->id) ;
                else
                    seq_1_al.push_back(gene1->neighbours[ j-1 ]->group->id) ;
                seq_2_al.push_back("groupNULL");
                j-- ;
                score += GAP_SCORE;
                counter++;
                break;
        }
    }
    if (!mirror)
        reverse( seq_1_al.begin(), seq_1_al.end() );
    reverse( seq_2_al.begin(), seq_2_al.end() );
     if (counter != seq_1_al.size() || counter != seq_2_al.size())
         printf("Warning: counter (%d) != seq_1_al.size(%lu) || seq_2_al.size(%lu)\n", counter, seq_1_al.size(), seq_2_al.size());
     else
         printf("Counter:\t%d\n", counter);
    
     double testscore = 0;
     for (int i = 0; i < seq_1_al.size(); i++)
     {
         if (seq_1_al[i] == "groupNULL" || seq_2_al[i] == "groupNULL")
             testscore += GAP_SCORE;
         else
             testscore += (seq_1_al[i] == seq_2_al[i]) ? MATCH_SCORE : GAP_SCORE;
     }
     printf("\nOriginal:\n");
     for (j = 0; j < gene1->neighbours.size(); j++)
         printf("%s\t", gene1->neighbours[j]->group->id.c_str());
     printf("\n");
     for (i = 0; i < gene2->neighbours.size(); i++)
         printf("%s\t", gene2->neighbours[i]->group->id.c_str());
     
     printf("\nModified:\n");
     for (int i = 0; i < seq_1_al.size(); i++)
        printf("%s\t", seq_1_al[i].c_str());
     printf("\n");
     for (int i = 0; i < seq_2_al.size(); i++)
             printf("%s\t", seq_2_al[i].c_str());
     printf("\n");
     
     if (testscore != score)
         printf("Warning: score (%f) != testscore (%f)\n", score, testscore);
    else
        printf("Score:\t%f\n\n", score);
    return score;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    int i, j, jj;
    double matches = 0.0f;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;

    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;

    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
    
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->orientation == gene2->orientation)
                jj = j;
            else
                jj = (SIZE_OF_NEIGHTBOURHOOD - 1) - j;

            if (gene1->neighbours[i]->group != NULL && gene2->neighbours[jj]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[jj]->group)
            {
                partialmatch = 1.0f / (double)(1 + abs(i - j));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        double distance = (i < (SIZE_OF_NEIGHTBOURHOOD / 2)) ? i - (SIZE_OF_NEIGHTBOURHOOD / 2) : i+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
    }

    return matches;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    double score1 = Adjusted(gene1, gene2, false);
    double score2 = Adjusted(gene1, gene2, true);

    return(score1 > score2 ? score1 : score2);
}

double Score::Adjusted(Gene *gene1, Gene *gene2, bool mirrored)
{
    int i, j;
    double matches = 0.0f;
    
    //if (gene1 == gene2)
    //    return 0.0f;

    //double weight = 1.0;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;

    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;

    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
            {
                
                int iloc = (i / 2)+1;
                int jloc = (j / 2)+1;
                
                if (i%2)
                    iloc *= -1;
                if (j%2)
                    jloc *= -1;
                
                if (mirrored)
                    jloc *= -1;
                
                partialmatch = 1.0f / (double)(abs(iloc-jloc)+1.0f);
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        //0..99
        //0..49
        //50..99
        double distance = (j < 50) ? j - (SIZE_OF_NEIGHTBOURHOOD / 2) : j+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
        //matches += maxpartialmatch * weight;
        
        //if (i%2)
        //    weight *= 0.9;
    }
    
    return matches;
}

double Score::Old(Gene *gene1, Gene *gene2)
{
    int i, j, ii, jj;
    double matches = 0.0f;
    
    //if (g1 == g2)
    //    return 0.0f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        if (gene1->orientation == gene2->orientation)
            ii = i;
        else
            ii = (SIZE_OF_NEIGHTBOURHOOD - 1) - i;

        if (gene1->neighbours[i]->group != NULL &&  gene2->neighbours[ii]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[ii]->group)
            matches++;
        else
        {
            double maxpartialmatch = 0.0;
            double partialmatch;
            
            for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            {
                if (gene1->orientation == gene2->orientation)
                    jj = i;
                else
                    jj = (SIZE_OF_NEIGHTBOURHOOD - 1) - j;

                if (gene1->neighbours[i]->group != NULL && gene2->neighbours[jj]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[jj]->group)
                {
                    partialmatch = (i/2 == j/2) ? 1.0 : 1.0f / double(abs((i/2)-(j/2))+1.0f);
                    
                    if (partialmatch > maxpartialmatch)
                        maxpartialmatch = partialmatch;
                }
            }
            
            matches += maxpartialmatch;
        }
    }
    return matches;
}
*/
