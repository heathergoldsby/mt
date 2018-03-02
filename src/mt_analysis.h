
#ifndef _MT_ANALYSIS_H_
#define _MT_ANALYSIS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
//#include <ea/analysis/tool.h>
#include <ea/analysis.h>
//#include <ea/functional.h>
#include <ea/digital_evolution/instruction_set.h>
//#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/digital_evolution/environment.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>


namespace ealib {
    namespace analysis {
        LIBEA_ANALYSIS_TOOL(task_profile) {
            
            // find the best...
            typename EA::individual_type best_founder = *ea.begin();
            int update_max = 2000;
            int best_update = 2000;
            int ind_count = 0;
            int best_ind = 0;

            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                typename EA::individual_ptr_type test_ea = ea.make_individual(*i->traits().founder());
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*test_ea,0) < get<GROUP_REP_THRESHOLD>(*test_ea)) &&
                       (cur_update < update_max)){
                    test_ea->update();
                    ++cur_update;
                }
                
                if (cur_update < best_update)  {
                    int germ_count = 0;
                    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                            
                            typename EA::individual_type::environment_type::location_type l = test_ea->env().location(x,y);
                            if (l.occupied()) {
                                if (get<GERM_STATUS>(*l.inhabitant(), true)) {
                                    germ_count++;
                                }
                            }
                        }
                        
                    }

                    if (germ_count) {
                        best_update = cur_update;
                        best_ind = ind_count;
                    }
                }
                ind_count++;
            }
            
            int cur_ind_count = 0;
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                if (cur_ind_count == best_ind) {
                    best_founder = *i;
                    break;
                }
                cur_ind_count++;
            }
            
            datafile df("task_profile.dat");
            df.add_field("x_position")
            .add_field("y_position")
            .add_field("gs")
            .add_field("workload")
            .add_field("task_profile");
            
            
            typename EA::individual_ptr_type control_ea = ea.make_individual(*best_founder.traits().founder());
            control_ea->resources().reset();
            
            // replay! till the group amasses the right amount of resources
            // or exceeds its window...
            int cur_update = 0;
            
            // and run till the group amasses the right amount of resources
            while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                   (cur_update < update_max)){
                control_ea->update();
                ++cur_update;
            }
            accumulator_set<double, stats<tag::mean, tag::variance> > workload;

            for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                    df.write(x);
                    df.write(y);
                    
                    typename EA::individual_type::environment_type::location_type l = control_ea->env().location(x,y);
                    if (l.occupied()) {
                        df.write(get<GERM_STATUS>(*l.inhabitant(), true))
                        .write(get<WORKLOAD>(*l.inhabitant(), 0))
                        .write(get<TASK_PROFILE>(*l.inhabitant(),""));
                        workload(get<WORKLOAD>(*l.inhabitant(), 0));
                    } else {
                        df.write("2")
                        .write("0")
                        .write("-");
                    }
                    df.endl();
                    
                }
                
            }
            df.endl();
            df.write(cur_update);
            df.write(control_ea->population().size());
            df.write(best_ind);
            df.write(variance(workload));

            df.endl();
            
        }
        
        LIBEA_ANALYSIS_TOOL(movie_gs) {
            
            // find the best...
            typename EA::individual_type best_founder = *ea.begin();
            int update_max = 2000;
            int best_update = 2000;
            int ind_count = 0;
            int best_ind = 0;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                typename EA::individual_ptr_type test_ea = ea.make_individual(*i->traits().founder());
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*test_ea,0) < get<GROUP_REP_THRESHOLD>(*test_ea)) &&
                       (cur_update < update_max)){
                    test_ea->update();
                    ++cur_update;
                }
                
                if (cur_update < best_update)  {
                    int germ_count = 0;
                    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                            
                            typename EA::individual_type::environment_type::location_type l = test_ea->env().location(x,y);
                            if (l.occupied()) {
                                if (get<GERM_STATUS>(*l.inhabitant(), true)) {
                                    germ_count++;
                                }
                            }
                        }
                        
                    }
                    
                    if (germ_count) {
                        best_update = cur_update;
                        best_ind = ind_count;
                    }
                }
                ind_count++;
            }
            
            int cur_ind_count = 0;
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                if (cur_ind_count == best_ind) {
                    best_founder = *i;
                    break;
                }
                cur_ind_count++;
            }
            
            
            
            datafile df("movie.dat");

            df.write(get<SPATIAL_X>(ea));
            df.write(get<SPATIAL_Y>(ea));
            df.endl();

            int cur_update = 0;
            typename EA::individual_ptr_type control_ea = ea.make_individual(*best_founder.traits().founder());
            
            while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                   (cur_update < update_max)){
                control_ea->update();
                cur_update++;
                
                df.write(cur_update);
//                df.write(get<GROUP_RESOURCE_UNITS>(*control_ea,0) );
//                df.write(get<MULTICELL_REP_TIME>(best_founder,0));
//                df.write(get<GROUP_RESOURCE_UNITS>(best_founder,0));
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                        
                        //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
                        typename EA::individual_type::environment_type::location_type l = control_ea->env().location(x,y);
                        if (l.occupied()) {
                            df.write(get<GERM_STATUS>(*l.inhabitant(), true))
                            .write(get<WORKLOAD>(*l.inhabitant(),0));
                            
                        } else {
                            df.write("2")
                            .write("0");
                        }
                        
                    }
                }
                df.endl();
                
            }
            
            df.endl();
            
        }


    }
}

#endif
