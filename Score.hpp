//
//  Score.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 19/09/2017.
//
//

#ifndef Score_hpp
#define Score_hpp

//#define GAP_PENALTY     -1.0f
//#define MATCH_SCORE     1.0f
//#define GAP_SCORE       -0.25f

#include <stdio.h>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include "Gene.hpp"
//enum Type{FAST, SIMPLE, NAIVE, SOPHISTICATED, ADJUSTED, OLD};

class Score
{
public:
    //static void StoreScore(Gene *, Gene *, double);
    //static double FetchScore(Gene *, Gene *);
    //static double Synteny(Gene*,Gene*, Type);
    static int Fast(Gene *, Gene *);
    static double Simple(Gene *, Gene *);
    static double Naive(Gene *, Gene *);
    static double Sophisticated(Gene *, Gene *);
    static double SophisticatedMirrored(Gene *, Gene *);
    static double Adjusted(Gene *, Gene *);
    static double Old(Gene *, Gene *);
    static double Test(Gene *, Gene *);
    //static map<u_int64_t, double> scorepool;
    //static map<u_int64_t, double> mirrorscorepool;
    //static unordered_map<pair<Gene*,Gene*>, double> scorepool;
    //static unordered_map<pair<Gene*,Gene*>, double> mirrorscorepool;
private:
    /*
    static double Old(Gene *, Gene *, bool );
    static double Simple(Gene *, Gene *, bool );
    static int Fast(Gene *, Gene *, bool);
    static double Adjusted(Gene *, Gene *, bool );
    static double Sophisticated(Gene *, Gene *, bool );
    static double Naive(Gene *, Gene *, bool );*/
};

#endif /* Score_hpp */
