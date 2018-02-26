#ifndef _MT_MOVIE_H_
#define _MT_MOVIE_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>

#include "mt_propagule_orig.h"


LIBEA_ANALYSIS_TOOL(movie) {
    
    int update_max = 500;
    typename EA::individual_type best_founder = *ea.begin();
    
//    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
//        
//        // not preserving location
//        typename EA::individual_type tmp(*i->traits().founder());
//        
//        for (int j=0; j<=update_max; ++j) {
//            tmp.update();
//        }
//        recalculate_fitness(tmp, ea);
//        
//        
//        
//        
//        // Calc fitness for each subpop
//        //                eval_two_stripes(tmp);
//        
//        // copy the stripe fit to the accumulator and also the subpop
//        double sf =get<STRIPE_FIT>(tmp);
//        
//        if (sf > max_fit) {
//            max_fit = sf;
//            best = *i;
//        }
//    }
//    
    
    datafile df("movie.dat");
    df.write(get<SPATIAL_X>(ea));
    df.write(get<SPATIAL_Y>(ea));
    df.endl();
    
    //typename EA::individual_type best_founder(*best.traits().founder());
    //typename EA::individual_type best_founder = *best.traits().founder();
    
    for (int j=0; j<=update_max; ++j) {
        best_founder.update();
        df.write(j);
        df.write(get<MULTICELL_REP_TIME>(best_founder,0));
        df.write(get<GROUP_RESOURCE_UNITS>(best_founder,0));
        // grab info based on location...
        for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
            for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                
                //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
                typename EA::individual_type::environment_type::location_type* l = &best_founder.env().location(x,y);
                
                
                if (l->occupied()) {
                    std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
                    
                    if(lt == "not") {
                        df.write("1");
                    }
                    if (lt == "nand") {
                        df.write("2");
                    }
                    if (lt == "or") {
                        df.write("3");
                    }
                    if (lt == "ornot") {
                        df.write("4");
                    }
                    if (lt == "and") {
                        df.write("5");
                    }
                    if (lt == "andnot") {
                        df.write("6");
                    }
                    if (lt == "nor") {
                        df.write("7");
                    }
                    if (lt == "xor") {
                        df.write("8");
                    }
                    if (lt == "equals") {
                        df.write("9");
                    }
                    if (lt == "") {
                        df.write("0");
                    }
                    
                    
                } else {
                    df.write("-1");
                }
                
            }
        }
        df.endl();
        
    }
    
    df.endl();
    
}







#endif

