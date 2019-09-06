
#ifndef _EALIFE_LOD_KNOCKOUTS_FITNESS_H_
#define _EALIFE_LOD_KNOCKOUTS_FITNESS_H_

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
        
        
        LIBEA_ANALYSIS_TOOL(lod_fitness) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;

            int num_rep = 100;

            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.end(); --i;
            
            for (int nr = 0; nr < num_rep; nr++) {
                
            
                // should define checkpoint + analysis input
                //ea is the thing loaded from the checkpoint; EA is its type
                EA metapop; // a new EA
                //typedef metadata md_type;
            
                typename EA::md_type md(ea.md());
                // override md settings (pop size, geo, etc)
                
                
                metapop.initialize(md);
                put<METAPOPULATION_SIZE>(32, metapop);
                put<RUN_UPDATES>(10000, metapop);
                put<RNG_SEED>(nr, metapop);

                if (nr != 0) {
                    metapop.reset_rng(nr);
                }

                typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
                typename EA::population_type init_mc;
                init_mc.insert(init_mc.end(),control_mc);
                
                std::swap(metapop.population(), init_mc);
                
                add_event<mt_gls_propagule>(metapop);
                
                int max_size = 32;
                int max_update = 50000;
                int cur_update = 0;
                
                while ((metapop.size() < max_size) &&
                       (cur_update < max_update)){
                    metapop.update();
                    ++cur_update;
                }
                
                // get workload
                float total_workload = 0;
                typedef typename EA::subpopulation_type::population_type subpop_type;

                for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                    for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;
                        total_workload += get<WORKLOAD>(org, 0.0);
                    }
                }
            
            
                df.write("final")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
            }
            

            
            
            // **i is the EA, AS OF THE TIME THAT IT DIED!
            
            // To replay, need to create new eas for each knockout exper.
            // setup the population (really, an ea):

            typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
            int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);

            // ok we need to iterate through size...
            // fixed size 100 genome...
            int uni_count = 0;
            int num_uni_viable = 0;
            int num_uni_inviable = 0;
            int num_uni = 0;


            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {
                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());

                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    
                    knockout_loc->population()[0]->genome()[z] = q;
                    typename EA::individual_ptr_type knockout_loc2 (knockout_loc);
                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }
                    


                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }
                        
                        num_uni_viable++;
                        

                        for (int nr = 0; nr < num_rep; nr++) {
                            // should define checkpoint + analysis input
                            //ea is the thing loaded from the checkpoint; EA is its type
                            EA metapop; // a new EA
                            //typedef metadata md_type;
                            typename EA::md_type md(ea.md());
                            // override md settings (pop size, geo, etc)
                            
                            
                            metapop.initialize(md);
                            put<METAPOPULATION_SIZE>(32, metapop);
                            put<RUN_UPDATES>(50000, metapop);
                            put<RNG_SEED>(nr, metapop);
                            if (nr != 0) {
                                metapop.reset_rng(nr);
                            }
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
                            while ((metapop.size() < max_size) &&
                                   (cur_update < max_update)){
                                metapop.update();
                                ++cur_update;
                            }
                            
                            // get workload
                            float total_workload = 0;
                            typedef typename EA::subpopulation_type::population_type subpop_type;
                            
                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                    typename EA::subpopulation_type::individual_type& org=**m;
                                    total_workload += get<WORKLOAD>(org, 0.0);
                                }
                            }
                            
                            df.write("final")
                            .write("uni")
                            .write(uni_count)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;

                    }
                    
                }
                
            }
            df2.write("final")
            .write(num_uni)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();

            
        }
        
        LIBEA_ANALYSIS_TOOL(lod_fitness_no_mutations) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("update")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;
            
            int num_rep = 100;
            
            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.end(); --i;
            
            for (int nr = 0; nr < num_rep; nr++) {
                
                
                // should define checkpoint + analysis input
                //ea is the thing loaded from the checkpoint; EA is its type
                EA metapop; // a new EA
                //typedef metadata md_type;
                
                typename EA::md_type md(ea.md());
                // override md settings (pop size, geo, etc)
                
                
                metapop.initialize(md);
                put<METAPOPULATION_SIZE>(32, metapop);
                put<RUN_UPDATES>(10000, metapop);
                put<RNG_SEED>(nr, metapop);
                put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                put<MUTATION_PER_SITE_P>(0, metapop);
                put<GERM_MUTATION_PER_SITE_P>(0, metapop);

                
                if (nr != 0) {
                    metapop.reset_rng(nr);
                }
                
                typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
                int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);

                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
                put<TASK_MUTATION_PER_SITE_P>(0, *control_mc);
                put<MUTATION_PER_SITE_P>(0, *control_mc);

                
                typename EA::population_type init_mc;
                init_mc.insert(init_mc.end(),control_mc);
                
                std::swap(metapop.population(), init_mc);
                
                add_event<mt_gls_propagule>(metapop);
                
                int max_size = 32;
                int max_update = 50000;
                int cur_update = 0;
                
                while ((metapop.size() < max_size) &&
                       (cur_update < max_update)){
                    metapop.update();
                    ++cur_update;
                }
                
                // get workload
                float total_workload = 0;
                typedef typename EA::subpopulation_type::population_type subpop_type;
                
                for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                    for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;
                        total_workload += get<WORKLOAD>(org, 0.0);
                    }
                }
                
                
                df.write("final")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
            }
            
            
            
            
            // **i is the EA, AS OF THE TIME THAT IT DIED!
            
            // To replay, need to create new eas for each knockout exper.
            // setup the population (really, an ea):
            
            typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
            int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);

            // ok we need to iterate through size...
            // fixed size 100 genome...
            int uni_count = 0;
            int num_uni_viable = 0;
            int num_uni_inviable = 0;
            int num_uni = 0;
            
            
            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {
                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                    
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc);


                    knockout_loc->population()[0]->genome()[z] = q;
                    typename EA::individual_ptr_type knockout_loc2 (knockout_loc);
                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc2);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc2);
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);

                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }
                    
                    
                    
                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }
                        
                        num_uni_viable++;
                        
                        
                        for (int nr = 0; nr < num_rep; nr++) {
                            // should define checkpoint + analysis input
                            //ea is the thing loaded from the checkpoint; EA is its type
                            EA metapop; // a new EA
                            //typedef metadata md_type;
                            typename EA::md_type md(ea.md());
                            // override md settings (pop size, geo, etc)
                            
                            
                            metapop.initialize(md);
                            put<METAPOPULATION_SIZE>(32, metapop);
                            put<RUN_UPDATES>(50000, metapop);
                            put<RNG_SEED>(nr, metapop);
                            put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                            put<MUTATION_PER_SITE_P>(0, metapop);
                            put<GERM_MUTATION_PER_SITE_P>(0, metapop);

                            if (nr != 0) {
                                metapop.reset_rng(nr);
                            }
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
                            while ((metapop.size() < max_size) &&
                                   (cur_update < max_update)){
                                metapop.update();
                                ++cur_update;
                            }
                            
                            // get workload
                            float total_workload = 0;
                            typedef typename EA::subpopulation_type::population_type subpop_type;
                            
                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                    typename EA::subpopulation_type::individual_type& org=**m;
                                    total_workload += get<WORKLOAD>(org, 0.0);
                                }
                            }
                            
                            df.write("final")
                            .write("uni")
                            .write(uni_count)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;
                        
                    }
                    
                }
                
            }
            df2.write("final")
            .write(num_uni)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();
            
            
        }

        
        LIBEA_ANALYSIS_TOOL(lod_fitness_at_trans) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;
            
            int num_rep = 100;
            
            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.begin(); i++;
            
            // find the first to transition
            for( ; i!=lod.end(); i++) {
                if (i->size() > 2) {
                    break;
                }
            }
            
            
            
            for (int nr = 0; nr < num_rep; nr++) {
                
                
                // should define checkpoint + analysis input
                //ea is the thing loaded from the checkpoint; EA is its type
                EA metapop; // a new EA
                //typedef metadata md_type;
                
                typename EA::md_type md(ea.md());
                // override md settings (pop size, geo, etc)
                
                
                metapop.initialize(md);
                put<METAPOPULATION_SIZE>(32, metapop);
                put<RUN_UPDATES>(10000, metapop);
                put<RNG_SEED>(nr, metapop);
                
                if (nr != 0) {
                    metapop.reset_rng(nr);
                }
                
                typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
                
                typename EA::population_type init_mc;
                init_mc.insert(init_mc.end(),control_mc);
                
                std::swap(metapop.population(), init_mc);
                
                add_event<mt_gls_propagule>(metapop);
                
                int max_size = 32;
                int max_update = 50000;
                int cur_update = 0;
                
                while ((metapop.size() < max_size) &&
                       (cur_update < max_update)){
                    metapop.update();
                    ++cur_update;
                }
                
                // get workload
                float total_workload = 0;
                typedef typename EA::subpopulation_type::population_type subpop_type;
                
                for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                    for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;
                        total_workload += get<WORKLOAD>(org, 0.0);
                    }
                }
                
                
                df.write("trans")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
            }
            
            
            
            
            // **i is the EA, AS OF THE TIME THAT IT DIED!
            
            // To replay, need to create new eas for each knockout exper.
            // setup the population (really, an ea):
            
            typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
            int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);

            // ok we need to iterate through size...
            // fixed size 100 genome...
            int uni_count = 0;
            int num_uni_viable = 0;
            int num_uni_inviable = 0;
            int num_uni = 0;
            
            
            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {
                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                    
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    
                    knockout_loc->population()[0]->genome()[z] = q;
                    typename EA::individual_ptr_type knockout_loc2 (knockout_loc);
                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }

                    
                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }

                        
                        num_uni_viable++;
                        
                        
                        for (int nr = 0; nr < num_rep; nr++) {
                            // should define checkpoint + analysis input
                            //ea is the thing loaded from the checkpoint; EA is its type
                            EA metapop; // a new EA
                            //typedef metadata md_type;
                            typename EA::md_type md(ea.md());
                            // override md settings (pop size, geo, etc)
                            
                            
                            metapop.initialize(md);
                            put<METAPOPULATION_SIZE>(32, metapop);
                            put<RUN_UPDATES>(10000, metapop);
                            put<RNG_SEED>(nr, metapop);
                            if (nr != 0) {
                                metapop.reset_rng(nr);
                            }
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
                            while ((metapop.size() < max_size) &&
                                   (cur_update < max_update)){
                                metapop.update();
                                ++cur_update;
                            }
                            
                            // get workload
                            float total_workload = 0;
                            typedef typename EA::subpopulation_type::population_type subpop_type;
                            
                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                    typename EA::subpopulation_type::individual_type& org=**m;
                                    total_workload += get<WORKLOAD>(org, 0.0);
                                }
                            }
                            
                            df.write("trans")
                            .write("uni")
                            .write(uni_count)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;
                        
                    }
                    
                }
                
            }
            df2.write("final")
            .write(uni_count)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();
            
            
        }
        

        
        
    }
}


#endif
