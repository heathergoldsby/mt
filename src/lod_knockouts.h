
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
            .add_field("no_knockouts")
            .add_field("rx_knockedout")
            .add_field("location_knockedout");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual();
                control_ea->rng().reset(get<RNG_SEED>(*i));
                
                typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual();
                knockout_rx_ea->rng().reset(get<RNG_SEED>(*i));
                knockout<instructions::rx_msg,instructions::nop_x>(*knockout_rx_ea);
                
                typename EA::individual_ptr_type knockout_location_ea = ea.make_individual();
                knockout_location_ea->rng().reset(get<RNG_SEED>(*i));
                knockout<instructions::get_xy,instructions::nop_x>(*knockout_location_ea);
                
                // setup the founder
                typename EA::individual_type::individual_ptr_type o=i->copy_individual(**((*i).traits().founder()->population().begin()));
                
                o->hw().initialize();
                control_ea->insert(control_ea->end(), o);
                
                // setup send knockout founder
                typename EA::individual_type::individual_ptr_type ko_s =i->copy_individual(**((*i).traits().founder()->population().begin()));

                ko_s->hw().initialize();
                knockout_rx_ea->insert(knockout_rx_ea->end(), ko_s);

                
                // setup location knockout founder
                typename EA::individual_type::individual_ptr_type ko_l = i->copy_individual(**((*i).traits().founder()->population().begin()));
                ko_l->hw().initialize();                
                knockout_location_ea->insert(knockout_location_ea->end(), ko_l);

                
                // replay! till the group amasses the right amount of resources
                // or exceeds its window...
                int cur_update = 0;
                int update_max = 1000;
                // and run till the group amasses the right amount of resources
                while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                       (cur_update < update_max)){
                    control_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                
                cur_update = 0;
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_rx_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_rx_ea)) &&
                       (cur_update < update_max)){
                    knockout_rx_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                
                cur_update = 0;
                while ((get<GROUP_RESOURCE_UNITS>(*knockout_location_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_location_ea)) &&
                       (cur_update < update_max)){
                    knockout_location_ea->update();
                    ++cur_update;
                }
                df.write(cur_update);
                
                
                df.endl();
                
                ++lod_depth;
            }
            //            }
            
            //        };
        }
    }
}
#endif
