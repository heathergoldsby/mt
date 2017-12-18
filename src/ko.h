



LIBEA_ANALYSIS_TOOL(ko) {
    
//    int update_max = 1000;
//    //typename EA::individual_type best_founder = *ea.begin();
//    
//
//    datafile df("ko.dat");
////    df.write(get<SPATIAL_X>(ea));
////    df.write(get<SPATIAL_Y>(ea));
////    df.endl();
////    
//    //typename EA::individual_type best_founder(*best.traits().founder());
//    //typename EA::individual_type best_founder = *best.traits().founder();
//    
//    typename EA::individual_ptr_type best_founder;
//    int best_time = 1000;
//    
//    // go through population... rerun, find best
//    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
//        
//        // setup the population (really, an ea):
//        typename EA::individual_ptr_type p = ea.make_individual();
//        p->initialize(ea.md());
//        p->reset_rng(i->rng().seed());
//        
//        for(typename EA::subpopulation_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
//            
//            typename EA::subpopulation_type::individual_ptr_type q = i->copy_individual(**j);
//            p->insert(p->end(), q);
//        }
//        
//        
//        int cur_update = 0;
//        typename EA::individual_ptr_type test_mc = p;
//        
//        while ((get<GROUP_RESOURCE_UNITS>(*test_mc,0) < get<GROUP_REP_THRESHOLD>(*test_mc)) &&
//               (cur_update < update_max)){
//            test_mc->update();
//            ++cur_update;
//        }
//        
//        if (cur_update < best_time) {
//            best_founder = test_mc;
//            best_time = cur_update;
//        }
//
//    }
////    
////    typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(), (*j)->genome().begin()+(*j)->hw().original_size());
////    
//    
//    /* given the best... knockouts!  replace each position with nop-x*/
//    /*std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
//    
//    // Must use representation size of 100.
//    assert(repr.size() == 100);
//    
//    repr[0] =  ea.isa()["h_alloc"];*/
//    // record mc size and time to replicate
//    
//    typename EA::subpopulation_type::individual_type::genome_type r(best_founder->traits()->founder()->population().begin().genome());
////
////    
////    for(std::size_t i=0; i<repr.size(); ++i) {
////    }
//
    
    
    
}


