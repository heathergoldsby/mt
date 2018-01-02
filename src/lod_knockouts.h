
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
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_knockouts.dat");
            df.add_field("lod_depth")
            .add_field("control_fit")
            .add_field("control_size")
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
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual();
                control_ea->initialize(ea.md());
                control_ea->reset_rng(ea.rng().seed());
                //control_ea->rng().reset(get<RNG_SEED>(*i));
                
                typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual();
                knockout_rx_ea->initialize(ea.md());
                knockout_rx_ea->rng().reset(get<RNG_SEED>(ea));
                knockout<instructions::rx_msg,instructions::nop_x>(*knockout_rx_ea);
                
//                typename EA::individual_ptr_type knockout_location_ea = ea.make_individual();
//                knockout_location_ea->initialize(ea.md());
//                knockout_location_ea->rng().reset(get<RNG_SEED>(ea));
//                knockout<instructions::get_xy,instructions::nop_x>(*knockout_location_ea);
//                
                // setup the founder
                //typename EA::individual_type::individual_ptr_type o=i->copy_individual(**((*i).traits().founder()->population().begin()));
                typename EA::individual_type::individual_ptr_type o=i->copy_individual(**(i->population().begin()));

                
                o->hw().initialize();
                control_ea->insert(control_ea->end(), o);
                
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
//                // setup send knockout founder
////                typename EA::individual_type::individual_ptr_type ko_s =i->copy_individual(**((*i).traits().founder()->population().begin()));
//                typename EA::individual_type::individual_ptr_type ko_s=i->copy_individual(**(i->population().begin()));
//
//                
//                ko_s->hw().initialize();
//                knockout_rx_ea->insert(knockout_rx_ea->end(), ko_s);
//
//                // get the genome...
                typename EA::individual_type::individual_ptr_type j = *(i->population().begin());
                typename EA::subpopulation_type::genome_type r(j->genome().begin(), j->genome().begin()+j->hw().original_size());
//
                
                float num_uni = 0;
                float num_multi = 0;
                float fit_uni = 0;
                float fit_multi = 0;
                float size_multi = 0;
                float num_pos = 0;
                float num_neg = 0;
                float num_neutral = 0;
                
//                // for each location, knockout!
//                r[0] =  knockout_rx_ea->isa()["nop_x"];
                
                
                // ok we need to iterate through size...
                for (int z =0; z < r.size(); z++) {
                    for (int q = 0; q < control_ea->isa().size(); q++) {
                        typename EA::individual_ptr_type knockout_loc = ea.make_individual();
                        knockout_loc->initialize(ea.md());
                        knockout_rx_ea->rng().reset(get<RNG_SEED>(ea));
                        typename EA::subpopulation_type::genome_type r2(r);
                        r2[z] = q;
                    
                        //typename EA::individual_type::individual_ptr_type ko_l = i->copy_individual(r2);
                        typename EA::individual_type::individual_ptr_type ko_l = knockout_loc->make_individual(r2);
                        inherits_from(*j, *ko_l, *i);


                        ko_l->hw().initialize();
                        knockout_loc->insert(knockout_loc->end(), ko_l);
                    
                    
                        int cur_update = 0;
                        int update_max = 2000;
                        // and run till the group amasses the right amount of resources
                        while ((get<GROUP_RESOURCE_UNITS>(*knockout_loc,0) < get<GROUP_REP_THRESHOLD>(*knockout_loc)) &&
                               (cur_update < update_max)){
                            knockout_loc->update();
                            ++cur_update;
                        }
                        // assess:
                        
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
                df.write(num_uni)
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
            //            }
            
            //        };
        }
    }
}
#endif
