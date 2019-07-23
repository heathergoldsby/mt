
#ifndef _EALIFE_LOD_CONTROL_ANALYSIS_H_
#define _EALIFE_LOD_CONTROL_ANALYSIS_H_

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

LIBEA_MD_DECL(ARCHIVE_OUTPUT_SIZE, "ea.mt.archive_output_size", int);


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
        LIBEA_ANALYSIS_TOOL(lod_control_analysis) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i; ++i;
            typename line_of_descent<EA>::iterator iend=lod.end(); --iend;
            
            datafile df("lod_control.dat");
            datafile df2("lod_control_genomes.dat");

            df.add_field("lod_depth")
            .add_field("flag_0")
            .add_field("flag_1")
            .add_field("germ")
            .add_field("size")
            ;
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            //            for( ; i!=lod.end(); ++i) {
            
            for( ; i!=lod.end(); i++) {
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                float flag_0 = 0;
                float flag_1 = 0;
                float num_germ = 0;
                int count = 0;
                for(typename EA::subpopulation_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                    if (get<FLAG>(**j,-1) == 0){
                        flag_0++;
                    }
                    if (get<FLAG>(**j,-1) == 1){
                        flag_1++;
                    }
                    if (get<GERM_STATUS>(**j, -1) == true) {
                        num_germ++;
                    }
                    df2.write(lod_depth)
                    .write(count)
                    .write(get<FLAG>(**j,-1))
                    .write(get<GERM_STATUS>(**j, -1));
                    
                    for(typename EA::subpopulation_type::genome_type::iterator k=(*j)->genome().begin(); k!=(*j)->genome().end(); ++k) {
                        df2.write(*k)
                        .write(" ");
                       
                    }
                    df2.endl();
                    
                    count++;
                }

                df.write(lod_depth)
                .write(flag_0)
                .write(flag_1)
                .write(num_germ)
                .write(i->population().size());
                df.endl();
                
                ++lod_depth;
                
            }
        }
        
        
        
    }
}


#endif
