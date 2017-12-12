
LIBEA_ANALYSIS_TOOL(ko) {
    
    int update_max = 1000;
    //typename EA::individual_type best_founder = *ea.begin();
    

    datafile df("ko.dat");
//    df.write(get<SPATIAL_X>(ea));
//    df.write(get<SPATIAL_Y>(ea));
//    df.endl();
//    
    //typename EA::individual_type best_founder(*best.traits().founder());
    //typename EA::individual_type best_founder = *best.traits().founder();
    
    typename EA::individual_ptr_type best_founder;
    int best_time = 1000;
    
    // go through population... rerun, find best
    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
        
        int cur_update = 0;
        typename EA::individual_ptr_type test_mc = copy_individual(*(i->traits().founder()));
        
        while ((get<GROUP_RESOURCE_UNITS>(*test_mc,0) < get<GROUP_REP_THRESHOLD>(*test_mc)) &&
               (cur_update < update_max)){
            test_mc->update();
            ++cur_update;
        }
        
        if (cur_update < best_time) {
            best_founder = test_mc;
            best_time = cur_update;
        }


    }
    /*
    for (int j=0; j<=update_max; ++j) {
        best_founder.update();
        df.write(j);
        df.write(get<MULTICELL_REP_TIME>(best_founder,0));
        df.write(get<GROUP_RESOURCE_UNITS>(best_founder,0));
        // grab info based on location...
        for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
            for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                
               
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
     */
    
}


