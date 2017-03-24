//
//  ts.h
//  ealife
//
//  Created by Heather Goldsby on 8/23/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_TS_H_
#define _EALIFE_TS_H_



#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>

#include "selfrep_not_ancestor.h"
#include "repro_not_ancestor.h"
#include "resource_consumption.h"
#include "configurable_mutation.h"


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ea;
using namespace boost::accumulators;

LIBEA_MD_DECL(TASK_SWITCHING_COST, "ea.ts.task_switching_cost", int);
LIBEA_MD_DECL(LAST_TASK, "ea.ts.last_task", std::string);
LIBEA_MD_DECL(NUM_SWITCHES, "ea.ts.num_switches", int);

LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.ts.germ_mutation_per_site_p", double);




/*! If an organism changes tasks, then it incurs a task-switching cost.
 */

template <typename EA>
struct task_switching_cost : task_performed_event<EA> {
    
    task_switching_cost(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_switching_cost() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r, // amount of resource consumed
                            EA& ea) {
        
        if (exists<LAST_TASK>(ind) && 
            (task->name() != get<LAST_TASK>(ind, ""))) {
            
            ind.hw().add_cost(get<TASK_SWITCHING_COST>(ea)); 
            get<NUM_SWITCHES>(ind, 0) += 1; 
        }
        put<LAST_TASK>(task->name(), ind); 
        
    }
};


/*! Prints information about the mean number of task-switches
 */


template <typename EA>
struct task_switch_tracking : end_of_update_event<EA> {
    task_switch_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("ts.dat") { 
        _df.add_field("update")
        .add_field("sub_pop_size")
        .add_field("pop_size")
        .add_field("mean ts");
        
    }
    
    //! Destructor.
    virtual ~task_switch_tracking() {
    }
    
    //! Track how many task-switches are being performed!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            double ts = 0;
            double org = 0;
            
            int sub_pop_size = 0;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                ++sub_pop_size;
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j){
                    
                    typename EA::individual_type::individual_type& ind=**j;
                    if (ind.alive()) {
                        ts += get<NUM_SWITCHES>(ind, 0);
                        ++org;
                    }
                }
            }
            ts /= org;
            _df.write(ea.current_update())
            .write(sub_pop_size)
            .write(org)
            .write(ts)
            .endl();
        }
        
    }
    datafile _df;    
    
};



//! Performs group replication using germ lines.
template <typename EA>
struct ts_replication : end_of_update_event<EA> {
    //! Constructor.
    ts_replication(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ts_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < 2) continue; 
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) && 
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                // why can't I just grab the founder? Does this make a copy?
                // grab a copy of the founder!
                
                typename EA::individual_type::individual_type prop = (*i).founder();
                prop.repr().resize((*i).founder().hw().original_size());

                //typename EA::individual_type::individual_type j = **(i->population().begin());
                //typename EA::individual_type::individual_type prop = j;
                //prop.repr().resize(j.hw().original_size());
                prop.hw().initialize();
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea)); 
                mutate(prop,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::individual_ptr_type o=p->make_individual(prop.repr());
                p->append(o);
                offspring.push_back(p);
                
                // reset resource units
                i->env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*i.base());
                offspring_pop.push_back(p);
                inherits(parent_pop, offspring_pop, ea);
            }
        }
        
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<META_POPULATION_SIZE>(ea) - offspring.size(); 
            
            typename EA::population_type survivors;
            select_n<selection::random>(ea.population(), survivors, n, ea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(ea.population(), survivors);
        }
       
    }

    
    
};
 



/*! An organism rotates to face its parent....
 */
template <typename EA>
struct ts_birth_event : birth_event<EA> {
    
    //! Constructor.
    ts_birth_event(EA& ea) : birth_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~ts_birth_event() {
    }
    
    /*! Called for every inheritance event. We are using the orientation of the first parent...
     */
    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
                            typename EA::individual_type& parent, // individual parent
                            EA& ea) {
        ea.env().face_org(parent, offspring);
        //get<GERM_STATUS>(offspring, true) = get<GERM_STATUS>(ind(parents.begin(),ea), true);
        
    }
};



#endif



