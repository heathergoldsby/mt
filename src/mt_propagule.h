

#ifndef _EALIFE_MTPROPAGULE_H_
#define _EALIFE_MTPROPAGULE_H_

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


// RES_UPDATE, MULTICELL_REP_TIME, DIVIDE_REMOTE

LIBEA_MD_DECL(DIVIDE_REMOTE, "ea.mt.divide_remote", int); // 0 = no divide; 1 divide
LIBEA_MD_DECL(DIVIDE_ALT, "ea.mt.divide_alt", int); // 0 = remote; 1 local
LIBEA_MD_DECL(MULTICELL_REP_TIME, "ea.mt.mcreptime", int);

/* Divide remote only works if there are enough resources... */
DIGEVO_INSTRUCTION_DECL(h_divide_remote) {
    if(hw.age() >= (0.8 * hw.original_size())) {
        typename Hardware::genome_type& r=hw.repr();
        
        // Check to see if the offspring would be a good length.
        int divide_pos = hw.getHeadLocation(Hardware::RH);
        int extra_lines = r.size() - hw.getHeadLocation(Hardware::WH);
        
        int child_size = r.size() - divide_pos - extra_lines;
        int parent_size = r.size() - child_size - extra_lines;
        double ratio = 2.0;
        
        if ((child_size < (hw.original_size()/ratio)) ||
            (child_size > (hw.original_size()*ratio)) ||
            (parent_size < (hw.original_size()/ratio)) ||
            (parent_size > (hw.original_size()*ratio))){
            // fail!
            return;
        }
    
        typename Hardware::genome_type::iterator f=r.begin(),l=r.begin();
        std::advance(f, hw.getHeadLocation(Hardware::RH));
        std::advance(l, hw.getHeadLocation(Hardware::WH));
//        typename Hardware::genome_type offr(f, l);
        
        r.resize(parent_size);
//        replicate(p, offr, ea);
        hw.replicated_soft_reset();
        
        if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<GROUP_REP_THRESHOLD>(ea, 0.0)) {
            // raise flag
            put<DIVIDE_REMOTE>(1, ea);
        }
    }
}



/* Divide remote only works if there are enough resources... */
DIGEVO_INSTRUCTION_DECL(h_alt_divide) {
    if(hw.age() >= (0.8 * hw.original_size())) {
        typename Hardware::genome_type& r=hw.repr();
        
        // Check to see if the offspring would be a good length.
        int divide_pos = hw.getHeadLocation(Hardware::RH);
        int extra_lines = r.size() - hw.getHeadLocation(Hardware::WH);
        
        int child_size = r.size() - divide_pos - extra_lines;
        int parent_size = r.size() - child_size - extra_lines;
        double ratio = 2.0;
        
        if ((child_size < (hw.original_size()/ratio)) ||
            (child_size > (hw.original_size()*ratio)) ||
            (parent_size < (hw.original_size()/ratio)) ||
            (parent_size > (hw.original_size()*ratio))){
            // fail!
            return;
        }
        
        
        typename Hardware::genome_type::iterator f=r.begin(),l=r.begin();
        std::advance(f, hw.getHeadLocation(Hardware::RH));
        std::advance(l, hw.getHeadLocation(Hardware::WH));
        typename Hardware::genome_type offr(f, l);
        
        
        r.resize(parent_size);
        hw.replicated_soft_reset();
        
        // remote = 0; local = 1.
        if(get<DIVIDE_ALT>(ea, 0) == 0) {
            if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<GROUP_REP_THRESHOLD>(ea, 0.0)) {
                // set rest to zero                // raise flag
                put<DIVIDE_REMOTE>(1, ea);
            }
            put<DIVIDE_ALT>(1,ea);
        } else {
            replicate(p, offr, ea);
            put<DIVIDE_ALT>(0,ea);
        }
        
    }

    
    
}



// make sure resources are moved to multi. check gls for example

//! Performs multicell replication using germ lines. One cells is selected, mutated, and then used to create the appropriate number of cells. Thus, the starting multicell offspring is clonal.
template <typename MEA>
struct mt_propagule : end_of_update_event<MEA> {
    //! Constructor.
    mt_propagule(MEA& mea) : end_of_update_event<MEA>(mea), _df("mt.dat") {
        _df.add_field("update")
        .add_field("mean_rep_time")
        .add_field("mean_res")
        .add_field("mean_multicell_size")
        .add_field("replication_count");
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~mt_propagule() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(MEA& mea) {


        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));

        
        // Replicate!
        
            // See if any subpops have exceeded the resource threshold
            typename MEA::population_type offspring;
            for(typename MEA::iterator i=mea.begin(); i!=mea.end(); ++i) {
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename MEA::subpopulation_type::population_type propagule_type;
//                propagule_type propagule;
                
                // track multicells (even those that don't replicate)
                if ((mea.current_update() % 100) == 0) {
                    
                    int alive_count = 0;
                    
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            alive_count++;
                        }
                        
                    }
                    
                    multicell_rep.push_back(get<MULTICELL_REP_TIME>(*i,0));
                    multicell_res.push_back(get<GROUP_RESOURCE_UNITS>(*i,0));
                    multicell_size.push_back(alive_count);
                }
                
            
            
                
                if (get<DIVIDE_REMOTE>(*i,0)){
                    
                    
                    // get a new subpopulation:
                    typename MEA::individual_ptr_type p = mea.make_individual();
                    p->initialize(mea.md());
                    p->reset_rng(mea.rng().seed());
                    
                    int num_moved = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                            (*j)->genome().begin()+(*j)->hw().original_size());
                            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                            
                            inherits_from(**j, *q, *p);
                            mutate(*q,m,*p);

                            
                            p->insert(p->end(), q);
                            
                            ++num_moved;
                            
                        }
                        if (num_moved > 0) {
                            break;
                        }
                    }
                    
                    offspring.insert(offspring.end(),p);
                    
                    
                    // replication
                    ++num_rep;
                    
                    // reset parent multicell
                    i->resources().reset();
                    int res_amt = get<GROUP_RESOURCE_UNITS>(*i) - get<GROUP_REP_THRESHOLD>(*i, 0.0);
                    put<GROUP_RESOURCE_UNITS>(res_amt,*i);
                    put<MULTICELL_REP_TIME>(0,*i);
                    put<DIVIDE_REMOTE>(0,*i);
                    
                    // i == parent individual;
                    typename MEA::population_type parent_pop, offspring_pop;
                    parent_pop.push_back(*i.base());
                    offspring_pop.push_back(p);
                    inherits(parent_pop, offspring_pop, mea);
                    
                    
                    
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
                .write(std::accumulate(multicell_rep.begin(), multicell_rep.end(), 0.0)/multicell_rep.size())
                .write(std::accumulate(multicell_res.begin(), multicell_res.end(), 0.0)/multicell_res.size())
                .write(std::accumulate(multicell_size.begin(), multicell_size.end(), 0.0)/multicell_size.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
                multicell_rep.clear();
                multicell_res.clear();
                multicell_size.clear();
            } else {
                _df.write(mea.current_update())
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(num_rep)
                .endl();
            }
        }
        
    }
    
    datafile _df;
    std::deque<double> multicell_rep;
    std::deque<double> multicell_res;
    std::deque<double> multicell_size;
    
    int num_rep;
    
    
};


#endif
