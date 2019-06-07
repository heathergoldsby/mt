//
//  ps_simple.h
//  mt
//
//  Created by Heather Goldsby on 5/23/19.
//  Copyright © 2019 Michigan State University. All rights reserved.
//

#ifndef ps_simple_h
#define ps_simple_h

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/mutation.h>



using namespace ealib;

LIBEA_MD_DECL(START_PROPAGULE_SIZE, "ea.mt.start_propagule_size", int);
LIBEA_MD_DECL(MEMBER_START_PROPAGULE, "ea.mt.member_start_propagule", int);
LIBEA_MD_DECL(FLAG, "ea.mt.flag", int);
LIBEA_MD_DECL(FLAG_LOCK, "ea.mt.flag_lock", int);


/*! Execute the next instruction if the cell was part of the propagule
 */
DIGEVO_INSTRUCTION_DECL(if_member_start_propagule) {
    if(!get<MEMBER_START_PROPAGULE>(*p,0)) {
        hw.advanceHead(Hardware::IP);
    }
}


/*! Execute the next instruction if the cell was not part of the propagule
 */
DIGEVO_INSTRUCTION_DECL(if_not_member_start_propagule) {
    if(get<MEMBER_START_PROPAGULE>(*p,0)) {
        hw.advanceHead(Hardware::IP);
    }
}

DIGEVO_INSTRUCTION_DECL(flag_0) {
    if (get<FLAG_LOCK>(*p,0) == 0) {
        get<FLAG>(*p,0) = 0;
    }
}

DIGEVO_INSTRUCTION_DECL(flag_1) {
    if (get<FLAG_LOCK>(*p,0) == 0) {
        get<FLAG>(*p,0) = 1;
    }
}

DIGEVO_INSTRUCTION_DECL(unlock_flag) {
    get<FLAG_LOCK>(*p,0) = 0;
}

DIGEVO_INSTRUCTION_DECL(lock_flag) {
    get<FLAG_LOCK>(*p,0) = 1;
}

// make sure resources are moved to multi. check gls for example

//! Performs multicell replication using germ lines. One cells is selected, mutated, and then used to create the appropriate number of cells. Thus, the starting multicell offspring is clonal.
template <typename MEA>
struct mt_ps_propagule : end_of_update_event<MEA> {
    //! Constructor.
    mt_ps_propagule(MEA& mea) : end_of_update_event<MEA>(mea), _df("ps_simple.dat") {
        _df.add_field("update")
        .add_field("mean_rep_time")
        .add_field("mean_res")
        .add_field("mean_multicell_size")
        .add_field("mean_prop_size")
        .add_field("replication_count")
        .add_field("mean_generation");
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~mt_ps_propagule() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(MEA& mea) {
        
        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
        accumulator_set<double, stats<tag::mean> > gen;
        
        
        // Replicate!
        int ru = 1;
        if ((mea.current_update() % ru) == 0) {
            
            // See if any subpops have exceeded the resource threshold
            typename MEA::population_type offspring;
            for(typename MEA::iterator i=mea.begin(); i!=mea.end(); ++i) {
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                gen(get<IND_GENERATION>(*i));
                // figure out which individuals from the parent comprise the propagule:
                typedef typename MEA::subpopulation_type::population_type propagule_type;
                
                // track multicells (even those that don't replicate)
                if ((mea.current_update() % 100) == 0) {
                    
                    multicell_rep.push_back(get<MULTICELL_REP_TIME>(*i,0));
                    multicell_res.push_back(get<GROUP_RESOURCE_UNITS>(*i,0));
                }
                
                
                
                
                if (get<DIVIDE_REMOTE>(*i,0)){
                    
                    
                    // get a new subpopulation:
                    typename MEA::individual_ptr_type p = mea.make_individual();
                    p->initialize(mea.md());
                    p->reset_rng(mea.rng().seed());
                    
                    int num_moved = 0;
                    int num_alive = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            if (get<GERM_STATUS>(**j, true)) {

                            typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                            (*j)->genome().begin()+(*j)->hw().original_size());
                            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                            
                            inherits_from(**j, *q, *p);
                            mutate(*q,m,*p);
                            put<MEMBER_START_PROPAGULE>(1,*q);
                            
                            p->insert(p->end(), q);
                            
                            ++num_moved;
                            }
                            ++num_alive;
                            
                        }
                    }
                    
                    if (num_moved == 0) { continue; }
                    multicell_prop_size.push_back(num_moved);
                    multicell_size.push_back(num_alive);


                    put<START_PROPAGULE_SIZE>(num_moved, *p);
                    

                    // track last replication state
                    //int rep_size = i->population().size();
                    
                    
                    offspring.insert(offspring.end(),p);
                    
                    
                    // replication
                    ++num_rep;
                    
                    // reset parent multicell
                    i->resources().reset();
                    put<GROUP_RESOURCE_UNITS>(0,*i);
                    put<MULTICELL_REP_TIME>(0,*i);
                    put<DIVIDE_REMOTE>(0,*i);
                    
                    // i == parent individual;
                    typename MEA::population_type parent_pop, offspring_pop;
                    parent_pop.push_back(*i.base());
                    offspring_pop.push_back(p);
                    inherits(parent_pop, offspring_pop, mea);
                    
                    
                    
                }
            }
            
            
            // select surviving parent groups
            if (offspring.size() > 0) {
                int n = get<METAPOPULATION_SIZE>(mea) - offspring.size();
                
                typename MEA::population_type survivors;
                select_n<selection::random< > >(mea.population(), survivors, n, mea);
                
                // add the offspring to the list of survivors:
                survivors.insert(survivors.end(), offspring.begin(), offspring.end());
                
                // and swap 'em in for the current population:
                std::swap(mea.population(), survivors);
            }
        }
        
        if ((mea.current_update() % 100) == 0) {
            
            if (multicell_rep.size() > 0) {
                _df.write(mea.current_update())
                .write(std::accumulate(multicell_rep.begin(), multicell_rep.end(), 0.0)/multicell_rep.size());
                if (multicell_res.size() > 0) {
                    _df.write(std::accumulate(multicell_res.begin(), multicell_res.end(), 0.0) /multicell_res.size())
                    .write(std::accumulate(multicell_size.begin(), multicell_size.end(), 0.0) /multicell_size.size());
                } else {
                    _df.write(0.0)
                    .write(0.0);
                }
                if (multicell_prop_size.size() > 0) {
                    _df.write(std::accumulate(multicell_prop_size.begin(), multicell_prop_size.end(), 0.0)/multicell_prop_size.size());
                } else {
                    _df.write(0.0);
                }

                _df.write(num_rep)
                .write(mean(gen))
                .endl();
                num_rep = 0;
                multicell_rep.clear();
                multicell_res.clear();
                multicell_size.clear();
                multicell_prop_size.clear();

            } else {
                _df.write(mea.current_update())
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(num_rep)
                .write(mean(gen))
                .endl();
            }
        }
        
    }
    
    datafile _df;
    std::deque<double> multicell_rep;
    std::deque<double> multicell_res;
    std::deque<double> multicell_size;
    std::deque<double> multicell_prop_size;
    
    int num_rep;
    
};

template <typename EA>
struct size_based_resources : end_of_update_event<EA> {
    //! Constructor.
    size_based_resources(EA& ea) : end_of_update_event<EA>(ea), _df("size_based_res.dat") {

    }
    
    
    //! Destructor.
    virtual ~size_based_resources() {
    }
    
    //! Give resources to populations
    virtual void operator()(EA& ea) {
        
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            float start_size = get<START_PROPAGULE_SIZE>(*i, 1);
            float current_size = i->population().size();
            float amt = current_size/start_size;
            get<GROUP_RESOURCE_UNITS>(*i, 0) += (current_size/start_size);
        }
    }
    datafile _df;

};
/*
template <typename EA>
struct germ_soma_based_resources : end_of_update_event<EA> {
    //! Constructor.
    germ_soma_based_resources(EA& ea) : end_of_update_event<EA>(ea), _df("size_based_res.dat") {
        
    }
    
    
    //! Destructor.
    virtual ~germ_soma_based_resources() {
    }
    
    //! Give resources to populations
    virtual void operator()(EA& ea) {
        
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            float num_germ = i->population().size() - get<CURR_SOMA_SIZE>(*i,0.0);
            float reward = 0;
            if (num_germ) {
                reward =  i->population().size() / (num_germ);
            }
            get<GROUP_RESOURCE_UNITS>(*i, 0) += reward;
        }
    }
    datafile _df;
    
};*/


template <typename EA>
struct flag_based_resources : end_of_update_event<EA> {
    //! Constructor.
    flag_based_resources(EA& ea) : end_of_update_event<EA>(ea), _df("size_based_res.dat") {
        
    }
    
    
    //! Destructor.
    virtual ~flag_based_resources() {
    }
    
    //! Give resources to populations
    virtual void operator()(EA& ea) {
        
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            float reward = 1;
            float num_0 = 0;
            float num_1 = 0;
            for(typename EA::subpopulation_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                if (get<FLAG>(**j,0) == 0){
                    num_0++;
                } else {
                    num_1++;
                }
                
            }
            if (num_0) {
                reward += ((num_1 + 1)/num_0);
            }
            get<GROUP_RESOURCE_UNITS>(*i, 0) += reward;
        }
    }
    datafile _df;
    
};



#endif /* ps_simple_h */
