
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
#include <ea/math/information.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace ealib::math;


namespace ealib {
    namespace analysis {
        
        
        LIBEA_ANALYSIS_TOOL(temporal_poly) {
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > not_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > nand_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > and_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > ornot_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > or_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > andnot_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > nor_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > xor_age;
            accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > equals_age;

            
            typename EA::individual_type best_founder = *ea.begin();
            int update_max = 2000;
            int best_update = 2000;
            int ind_count = 0;
            int best_ind = 0;
            std::vector<std::string> tps;
            
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
            
        
            
            typename EA::individual_ptr_type control_ea = ea.make_individual(*best_founder.traits().founder());
            control_ea->resources().reset();
            
            int cur_update = 0; 
            // and run till the group amasses the right amount of resources
            while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                    
                    
                    // see what the org did...
                    // reset its counts to 0
                    // add on it's age to the appropriate spot
                    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                            
                            typename EA::individual_type::environment_type::location_type l = control_ea->env().location(x,y);
                            if (l.occupied()) {
                                // org age?
                                int age = cur_update - get<IND_BIRTH_UPDATE>(*l.inhabitant(),0);
                                while(get<TASK_NOT>(*l.inhabitant(),0.0) > 0) {
                                    not_age(age);
                                    get<TASK_NOT>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_NAND>(*l.inhabitant(),0.0) > 0) {
                                    nand_age(age);
                                    get<TASK_NAND>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_AND>(*l.inhabitant(),0.0) > 0) {
                                    and_age(age);
                                    get<TASK_AND>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_ORNOT>(*l.inhabitant(),0.0) > 0) {
                                    ornot_age(age);
                                    get<TASK_ORNOT>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_OR>(*l.inhabitant(),0.0) > 0) {
                                    or_age(age);
                                    get<TASK_OR>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_ANDNOT>(*l.inhabitant(),0.0) > 0) {
                                    andnot_age(age);
                                    get<TASK_ANDNOT>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_NOR>(*l.inhabitant(),0.0) > 0) {
                                    nor_age(age);
                                    get<TASK_NOR>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_XOR>(*l.inhabitant(),0.0) > 0) {
                                    xor_age(age);
                                    get<TASK_XOR>(*l.inhabitant(),0.0)--;
                                }
                                while(get<TASK_EQUALS>(*l.inhabitant(),0.0) > 0) {
                                    equals_age(age);
                                    get<TASK_EQUALS>(*l.inhabitant(),0.0)--;
                                }
                                
                            }
                        }
                        
                    }
                }
            
            datafile _df("temporal_poly.dat");
            _df.add_field("mean_not_age")
            .add_field("var_not_age")
            .add_field("mean_nand_age")
            .add_field("var_nand_age")
            .add_field("mean_and_age")
            .add_field("var_and_age")
            .add_field("mean_ornot_age")
            .add_field("var_ornot_age")
            .add_field("mean_or_age")
            .add_field("var_or_age")
            .add_field("mean_andnot_age")
            .add_field("var_andnot_age")
            .add_field("mean_nor_age")
            .add_field("var_nor_age")
            .add_field("mean_xor_age")
            .add_field("var_xor_age")
            .add_field("mean_equals_age")
            .add_field("var_equals_age");
            ;
            
            if (count(not_age)) {
                _df.write(mean(not_age))
                .write(variance(not_age));
            } else {
                _df.write(0)
                .write(0);
            }
            if (count(nand_age)) {
                _df.write(mean(nand_age))
                .write(variance(nand_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(and_age)) {
                _df.write(mean(and_age))
                .write(variance(and_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(ornot_age)) {
                _df.write(mean(ornot_age))
                .write(variance(ornot_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(or_age)) {
                _df.write(mean(or_age))
                .write(variance(or_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(andnot_age)) {
                _df.write(mean(andnot_age))
                .write(variance(andnot_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(nor_age)) {
                _df.write(mean(nor_age))
                .write(variance(nor_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(xor_age)) {
                _df.write(mean(xor_age))
                .write(variance(xor_age));
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (count(equals_age)) {
                _df.write(mean(equals_age))
                .write(variance(equals_age));
            } else {
                _df.write(0)
                .write(0);
            }
            _df.endl();
            
        }
        
        LIBEA_ANALYSIS_TOOL(task_profile) {
            
            // find the best...
            typename EA::individual_type best_founder = *ea.begin();
            int update_max = 2000;
            int best_update = 2000;
            int ind_count = 0;
            int best_ind = 0;
            std::vector<std::string> tps;

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
                        tps.push_back(get<TASK_PROFILE>(*l.inhabitant(),""));
                    } else {
                        df.write("2")
                        .write("0")
                        .write("-");
                    }
                    df.endl();
                    
                }
                
            }
            df.endl();
            
            float shannon_sum = 0;
            for (int m = 0; m < tps.size(); m++) {
                for (int n = m+1; n < tps.size(); n++){
                    shannon_sum += mutual_information(tps[m], tps[n]);
                }
            }
            
            df.write(cur_update);
            df.write(control_ea->population().size());
            df.write(best_ind);
            df.write(variance(workload));
            df.write(shannon_sum);
            df.write(shannon_sum / tps.size());
            

            df.endl();
            
        }
        
        LIBEA_ANALYSIS_TOOL(task_profile2) {
            
            // find the best...
            typename EA::individual_type best_founder = *ea.begin();
            int update_max = 2000;
            int best_update = 2000;
            int ind_count = 0;
            int best_ind = 0;
            std::vector<std::string> tps;
            
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
            while ((cur_update < update_max)){
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
                        tps.push_back(get<TASK_PROFILE>(*l.inhabitant(),""));
                    } else {
                        df.write("2")
                        .write("0")
                        .write("-");
                    }
                    df.endl();
                    
                }
                
            }
            df.endl();
            int max_size = 200;
            for (int m = 0; m < tps.size(); m++) {
                tps[m].resize(max_size);
            }

            
            float shannon_sum = 0;
            for (int m = 0; m < tps.size(); m++) {
                for (int n = m+1; n < tps.size(); n++){
                    shannon_sum += mutual_information(tps[m], tps[n]);
                }
            }
            
            df.write(cur_update);
            df.write(control_ea->population().size());
            df.write(best_ind);
            df.write(variance(workload));
            df.write(shannon_sum);
            df.write(shannon_sum / tps.size());
            
            
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
        LIBEA_ANALYSIS_TOOL(dom_mutational_analysis) {
            
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
            
            
        
            int cur_update = 0;
            typename EA::individual_ptr_type control_ea = ea.make_individual(*best_founder.traits().founder());
            
            
            datafile df("mutational_analysis");
            df.add_field("control_fit")
            .add_field("control_size")
            .add_field("num_ko_inviable")
            .add_field("num_ko_uni")
            .add_field("num_ko_multi")
            .add_field("fit_uni")
            .add_field("fit_multi")
            .add_field("size_multi")
            .add_field("num_neg_mut")
            .add_field("num_neutral_mut")
            .add_field("num_pos_mut")
            ;
            
            
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                df.write(control_ea->population().size());
                
                float control_fit = cur_update;
                
            
                //
                float num_inviable = 0;
                float num_uni = 0;
                float num_multi = 0;
                float fit_uni = 0;
                float fit_multi = 0;
                float size_multi = 0;
                float num_pos = 0;
                float num_neg = 0;
                float num_neutral = 0;
                
                
                
                // ok we need to iterate through size...
                // fixed size 100 genome...
                for (int z =0; z < 100; z++) {
                    for (int q = 0; q < control_ea->isa().size(); q++) {
                        typename EA::individual_ptr_type knockout_loc = ea.make_individual(*best_founder.traits().founder());
                        put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                        
                        
                       
                        knockout_loc->population()[0]->genome()[z] = q;
                        
                        int cur_update = 0;
                        int update_max = 2000;
                        // and run till the group amasses the right amount of resources
                        while ((get<GROUP_RESOURCE_UNITS>(*knockout_loc,0) < get<GROUP_REP_THRESHOLD>(*knockout_loc)) &&
                               (cur_update < update_max)){
                            knockout_loc->update();
                            ++cur_update;
                        }
                        // assess:
                        
                        if (cur_update == update_max){
                            num_inviable++;
                        } else {
                            
                            if (knockout_loc->population().size() < 1.5) {
                                num_uni ++;
                                fit_uni += cur_update;
                            } else {
                                num_multi++;
                                fit_multi += cur_update;
                                size_multi += (knockout_loc->population().size());
                            }
                            
                            if (cur_update == control_fit) {
                                num_neutral ++;
                            }
                            if (cur_update < control_fit) {
                                num_pos ++;
                            }
                            if (cur_update > control_fit) {
                                num_neg ++;
                            }
                        }
                        
                        
                    }
                }
                df.write(num_inviable)
                .write(num_uni)
                .write(num_multi);
                if (num_uni > 0) {
                    df.write(fit_uni / num_uni);
                } else {
                    df.write(0.0);
                }
                
                if (num_multi > 0) {
                    df.write(fit_multi / num_multi)
                    .write(size_multi / num_multi);
                } else {
                    df.write(0.0)
                    .write(0.0);
                }
                df.write (num_neg)
                .write(num_neutral)
                .write(num_pos);
                
                
                
        }
        


    }
}

#endif
