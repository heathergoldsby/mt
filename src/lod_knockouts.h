
#ifndef _EALIFE_LOD_KNOCKOUTS_H_
#define _EALIFE_LOD_KNOCKOUTS_H_

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
        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         
         using namespace ealib;
         //    using namespace ealib::analysis;
         
         line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
         typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
         
         for( ; i!=lod.end(); ++i) {
         typename EA::individual_ptr_type control_ea = ea.make_individual();
         typename EA::individual_type::individual_ptr_type o = i->make_individual(i->founder().repr());
         o->hw().initialize();
         control_ea->append(o);
         }
         
         */
        LIBEA_ANALYSIS_TOOL(lod_knockouts) {
            //        template <typename EA>
            //        struct lod_knockouts : public unary_function<EA> {
            //        struct lod_knockouts : public ealib::analysis::unary_function<EA> {
            
            //            static const char* name() { return "lod_knockouts"; }
            //
            //            virtual void operator()(EA& ea) {
            //                using namespace ealib;
            //                using namespace ealib::analysis;
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i; ++i;
            
            datafile df("lod_knockouts.dat");
            df.add_field("lod_depth")
            .add_field("control_fit")
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
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
//            for( ; i!=lod.end(); ++i) {
            for( ; i!=lod.end(); i++) {
                if ((lod_depth % 10) != 0) {
                    lod_depth ++;
                    continue;
                }
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                
                typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());

                
                
                
                //typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual();
                //knockout<instructions::rx_msg,instructions::nop_x>(*knockout_rx_ea);
                
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 2000;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                df.write(control_ea->population().size());
                
                float control_fit = cur_update;

                // get the genome...
//                typename EA::individual_type::individual_ptr_type j = *i->traits().founder()->population().begin();
//                typename EA::subpopulation_type::genome_type r(j->genome().begin(), j->genome().begin()+j->hw().original_size());
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
                        typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                        put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);

                        
                        //typename EA::subpopulation_type::genome_type r2(knockout_loc->population()[0]->genome());
                        //r);
                        knockout_loc->population()[0]->genome()[z] = q;
                        //r2[z] = q;
                    
                        //(*(knockout_loc->population().begin()))->genome()[z] = q;

                    
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
            

                
                
                
//                // setup location knockout founder
//                typename EA::individual_type::individual_ptr_type ko_l = i->copy_individual(**((*i).traits().founder()->population().begin()));
//                ko_l->hw().initialize();                
//                knockout_location_ea->insert(knockout_location_ea->end(), ko_l);
//
                

                
//                cur_update = 0;
//                while ((get<GROUP_RESOURCE_UNITS>(*knockout_rx_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_rx_ea)) &&
//                       (cur_update < update_max)){
//                    knockout_rx_ea->update();
//                    ++cur_update;
//                }
//                df.write(cur_update);
                
//                cur_update = 0;
//                while ((get<GROUP_RESOURCE_UNITS>(*knockout_location_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_location_ea)) &&
//                       (cur_update < update_max)){
//                    knockout_location_ea->update();
//                    ++cur_update;
//                }
//                df.write(cur_update);
                
                
                df.endl();
                
                ++lod_depth;
                
            }
        }
        
        LIBEA_ANALYSIS_TOOL(lod_knockouts_capabilities) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i; ++i;
            
            datafile df("lod_capability_knockout.dat");
            df.add_field("lod_depth")
            .add_field("control_fit")
            .add_field("control_size")
            .add_field("rx_ko_fit")
            .add_field("rx_ko_size")
            .add_field("neighbor_ko_fit")
            .add_field("neighbor_ko_size")
            .add_field("gs_sense_ko_fit")
            .add_field("gs_sense_ko_size")
            .add_field("res_sense_ko_fit")
            .add_field("res_sense_ko_size")
            .add_field("soma_ko_fit")
            .add_field("soma_ko_size")
            ;
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); i++) {
                if ((lod_depth % 10) != 0) {
                    lod_depth ++;
                    continue;
                }
                
                df.write(lod_depth);
                
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *control_ea);
             
                
                typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual(*i->traits().founder());
                knockout<instructions::rx_msg,instructions::nop_x>(*knockout_rx_ea);
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_rx_ea);

                
                typename EA::individual_ptr_type knockout_neighbor_ea = ea.make_individual(*i->traits().founder());
                knockout<instructions::is_neighbor,instructions::nop_x>(*knockout_neighbor_ea);
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_neighbor_ea);

                
                typename EA::individual_ptr_type knockout_sense_gs_ea = ea.make_individual(*i->traits().founder());
                knockout<if_germ,instructions::nop_x>(*knockout_sense_gs_ea);
                knockout<if_soma,instructions::nop_x>(*knockout_sense_gs_ea);
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_sense_gs_ea);
                
                typename EA::individual_ptr_type knockout_sense_res_ea = ea.make_individual(*i->traits().founder());
                knockout<if_res_more_than_thresh,instructions::nop_x>(*knockout_sense_res_ea);
                knockout<if_res_less_than_thresh,instructions::nop_x>(*knockout_sense_res_ea);
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_sense_res_ea);
                
                
                typename EA::individual_ptr_type knockout_soma_ea = ea.make_individual(*i->traits().founder());
                knockout<become_soma,instructions::nop_x>(*knockout_soma_ea);
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_soma_ea);
                
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 2000;

                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                df.write(control_ea->population().size());
                
                
                cur_update = 0;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_rx_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_rx_ea)) &&
                       (cur_update < update_max)){
                    knockout_rx_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(knockout_rx_ea->population().size());
                
                
                cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_neighbor_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_neighbor_ea)) &&
                       (cur_update < update_max)){
                    knockout_neighbor_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(knockout_neighbor_ea->population().size());

                cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_sense_gs_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_sense_gs_ea)) &&
                       (cur_update < update_max)){
                    knockout_sense_gs_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(knockout_sense_gs_ea->population().size());
                cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_sense_res_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_sense_res_ea)) &&
                       (cur_update < update_max)){
                    knockout_sense_res_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(knockout_sense_res_ea->population().size());
                
                cur_update = 0;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_soma_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_soma_ea)) &&
                       (cur_update < update_max)){
                    knockout_soma_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(knockout_soma_ea->population().size());
                

                
                
                df.endl();
                
                ++lod_depth;
            }

        }
        
                

        
        
        LIBEA_ANALYSIS_TOOL(lod_transition) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_transition.dat");
            df.add_field("lod_depth")
            .add_field("fit")
            .add_field("size")
            .add_field("not")
            .add_field("nand")
            .add_field("and")
            .add_field("ornot")
            .add_field("or")
            .add_field("andnot")
            .add_field("nor")
            .add_field("xor ")
            .add_field("equals")
            ;
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 2000;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(control_ea->population().size());
                
                df.write(get<TASK_NOT>(*control_ea, 0.0))
                .write(get<TASK_NAND>(*control_ea, 0.0))
                .write(get<TASK_AND>(*control_ea, 0.0))
                .write(get<TASK_ORNOT>(*control_ea, 0.0))
                .write(get<TASK_OR>(*control_ea, 0.0))
                .write(get<TASK_ANDNOT>(*control_ea, 0.0))
                .write(get<TASK_NOR>(*control_ea, 0.0))
                .write(get<TASK_XOR>(*control_ea, 0.0))
                .write(get<TASK_EQUALS>(*control_ea, 0.0))
                .endl();
                
                ++lod_depth;
            }
            
        }

            

        
        LIBEA_ANALYSIS_TOOL(lod_report_gs) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_report_gs.dat");
            df.add_field("lod_depth")
            .add_field("fit")
            .add_field("size")
            .add_field("num_germ")
            .add_field("num_soma")
            .add_field("germ_workload")
            .add_field("germ_workload_var")
            .add_field("soma_workload")
            .add_field("soma_workload_var")
            ;
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 2000;
                
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                
                df.write(cur_update);
                df.write(control_ea->population().size());
                
                int germ_count = 0;
                int soma_count = 0;
                accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc;
                accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc;
                
                for(typename EA::subpopulation_type::population_type::iterator j=control_ea->population().begin(); j!=control_ea->population().end(); ++j) {
                    
                    typename EA::subpopulation_type::individual_type& org=**j;
                    if (get<GERM_STATUS>(org, true)) {
                        ++germ_count;
                        germ_workload_acc(get<WORKLOAD>(org, 0.0));
                    } else {
                        soma_workload_acc(get<WORKLOAD>(org, 0.0));
                        ++soma_count;
                    }
                }
                
                
                if (germ_count) {
                    df.write(germ_count)
                    .write(soma_count)
                    .write(mean(germ_workload_acc))
                    .write(variance(germ_workload_acc))
                    ;
                } else {
                    df.write(0)
                    .write(0)
                    .write(0)
                    .write(0);
                }
                if (soma_count) {
                df.write(mean(soma_workload_acc))
                .write(variance(soma_workload_acc));
                } else {
                    df.write(0)
                    .write(0);
                }
                
                df.endl();
                
                ++lod_depth;
            }
            
        }
        
        
        LIBEA_ANALYSIS_TOOL(lod_gls_circle_square_plot) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_gls_circle_square_plot.dat");
            df.add_field("lod_depth");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
                control_ea->resources().reset();

                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 2000;
                

                // and run till the group amasses the right amount of resources
                //while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                while ((get<DIVIDE_REMOTE>(*control_ea,0) == 0) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                //df.write(cur_update);
                
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
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
                
                ++lod_depth;
            }
        }

        LIBEA_ANALYSIS_TOOL(lod_archive_reversion) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            int lod_length = lod.size();
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            
            
            int lod_depth = 0;
            int revert = false;
            int flag_one = -1;
            int flag_two = lod_length - 1;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                int bupdate = get<IND_BIRTH_UPDATE>(*i);
                // if birth update > 1000000 look for reversion.
                ++lod_depth;

                if (bupdate > 1000000) {
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
                    put<RNG_SEED>(get<RNG_SEED>(*i->traits().founder()), *control_ea);

                    // replay! till the group amasses the right amount of resources
                    // or exceeds its window...
                    int cur_update = 0;
                    int update_max = 2000;
                
                
                    // and run till the group amasses the right amount of resources
                    //while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                    while ((get<DIVIDE_REMOTE>(*control_ea,0) == 0) &&
                       (cur_update < update_max)){
                        control_ea->update();
                        ++cur_update;
                    }
                
                    
                
                    if ((control_ea->size() == 1) && (revert==false)) {
                        revert = true;
                        typename EA::population_type output;
                        std::string fname = "archive_revert.xml";
                        archive::load_if(fname, output, ea);
                    
                        int archive_mark = get<ARCHIVE_MARK>(ea,0);
                        // copy the population:
                        typename EA::individual_ptr_type arch_ind = ea.make_individual(*i->traits().founder());

                        get<ARCHIVE_MARK>(*arch_ind,0) = archive_mark;

                        for(int k=0; k < get<METAPOPULATION_SIZE>(ea); k++) {
                            output.insert(output.end(), ea.copy_individual(*arch_ind));
                        }
                    
                        // save the output archive:
                        archive::save(fname, output, ea);
                        flag_one = floor((lod_length - lod_depth)/2) + lod_depth;
                    
                    }
                
                    if ((lod_depth == flag_one) || (lod_depth == flag_two)) {
                        typename EA::population_type output;
                        std::string fname = "archive_revert_";

                        if (lod_depth == flag_one) {
                            fname += "1.xml";
                        } else {
                            fname += "2.xml";
                        }
                   
                    
                        int archive_mark = get<ARCHIVE_MARK>(ea,0);
                        // copy the population:
                        get<ARCHIVE_MARK>(*i,0) = archive_mark;
                    
                        for(int k=0; k < get<METAPOPULATION_SIZE>(ea); k++) {
                            output.insert(output.end(), ea.copy_individual(*i));
                        }
                    
                        // save the output archive:
                        archive::save(fname, output, ea);
                    }
                }
                
            }
        }

        
        

    }
}


#endif
