

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

#include "gls.h"

//#include "stripes.h"

using namespace ealib;


// RES_UPDATE, MULTICELL_REP_TIME, DIVIDE_REMOTE

LIBEA_MD_DECL(DIVIDE_REMOTE, "ea.mt.divide_remote", int); // 0 = no divide; 1 divide
LIBEA_MD_DECL(DIVIDE_ALT, "ea.mt.divide_alt", int); // 0 = remote; 1 local
LIBEA_MD_DECL(MULTICELL_REP_TIME, "ea.mt.mcreptime", int);
LIBEA_MD_DECL(IND_REP_THRESHOLD, "ea.mt.ind_rep_threshold", int); // 0 = no divide; 1 divide
LIBEA_MD_DECL(COST_START_UPDATE, "ea.mt.cost_start_update", int);
LIBEA_MD_DECL(COST_RAMP, "ea.mt.cost_ramp", int);


//! Execute the next instruction if group resources exceed threshold.
DIGEVO_INSTRUCTION_DECL(if_res_more_than_thresh) {
    if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) < get<GROUP_REP_THRESHOLD>(ea, 0.0)) {
        hw.advanceHead(Hardware::IP);
    }
}

//! Execute the next instruction if group resources is less than the threshold.
DIGEVO_INSTRUCTION_DECL(if_res_less_than_thresh) {
    if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<GROUP_REP_THRESHOLD>(ea, 0.0)) {
        hw.advanceHead(Hardware::IP);
    }
}


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
            // set rest to zero
            int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<GROUP_REP_THRESHOLD>(ea, 0.0);
            put<GROUP_RESOURCE_UNITS>(res_amt,ea);
            // raise flag
            put<DIVIDE_REMOTE>(1, ea);
        }
    }
}


/* changed to using ramped costs... */
/*
DIGEVO_INSTRUCTION_DECL(h_divide_local) {
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
        
        // defaults to no cost... ramps up the cost one step per update till max.
        int local_cost = 0;
        int birth_update = get<IND_BIRTH_UPDATE>(ea);
        int start_update = get<COST_START_UPDATE>(ea);
        if (get<IND_REP_THRESHOLD>(ea, 0.0) > 0) {
            int cu =ea.current_update() + birth_update;
            local_cost = floor((cu - start_update)/get<COST_RAMP>(ea,1));
            if (local_cost < 0) { local_cost = 0; }
            if (local_cost > get<IND_REP_THRESHOLD>(ea, 0.0)) { local_cost = get<IND_REP_THRESHOLD>(ea, 0.0); }
        }
        
        typename Hardware::genome_type::iterator f=r.begin(),l=r.begin();
        std::advance(f, hw.getHeadLocation(Hardware::RH));
        std::advance(l, hw.getHeadLocation(Hardware::WH));
        typename Hardware::genome_type offr(f, l);
        
        r.resize(parent_size);
        hw.replicated_soft_reset();
        
        if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) >= local_cost) {
            // raise flag
            int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - local_cost;
            put<GROUP_RESOURCE_UNITS>(res_amt,ea);
            replicate(p, offr, ea);
            
        }
        
        
    }
}*/



DIGEVO_INSTRUCTION_DECL(h_divide_local) {
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
        
        if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<IND_REP_THRESHOLD>(ea, 0.0)) {
            // raise flag
            int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<IND_REP_THRESHOLD>(ea, 0.0);
            put<GROUP_RESOURCE_UNITS>(res_amt,ea);
            replicate(p, offr, ea);
            
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
                // set rest to zero
                int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<GROUP_REP_THRESHOLD>(ea, 0.0);
                put<GROUP_RESOURCE_UNITS>(res_amt,ea);
                // raise flag
                put<DIVIDE_REMOTE>(1, ea);
            }
            put<DIVIDE_ALT>(1,ea);
        } else {
            if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<IND_REP_THRESHOLD>(ea, 0.0)) {
                // raise flag
                int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<IND_REP_THRESHOLD>(ea, 0.0);
                put<GROUP_RESOURCE_UNITS>(res_amt,ea);
                replicate(p, offr, ea);
                
            }
            
            put<DIVIDE_ALT>(0,ea);
        }
        
    }
    
    
    
}



/* Divide remote only works if there are enough resources... */
DIGEVO_INSTRUCTION_DECL(h_divide_multicell) {
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
        
        if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<GROUP_REP_THRESHOLD>(ea, 0.0)) {
            // set rest to zero
            int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<GROUP_REP_THRESHOLD>(ea, 0.0);
            put<GROUP_RESOURCE_UNITS>(res_amt,ea);
            put<DIVIDE_REMOTE>(1, ea);
        } else {
            
            if (get<GROUP_RESOURCE_UNITS>(ea, 0.0) > get<IND_REP_THRESHOLD>(ea, 0.0)) {
                // raise flag
                int res_amt = get<GROUP_RESOURCE_UNITS>(ea) - get<IND_REP_THRESHOLD>(ea, 0.0);
                put<GROUP_RESOURCE_UNITS>(res_amt,ea);
                replicate(p, offr, ea);
                
            }
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
        int ru = 1;
        if ((mea.current_update() % ru) == 0) {
            
            // See if any subpops have exceeded the resource threshold
            typename MEA::population_type offspring;
            for(typename MEA::iterator i=mea.begin(); i!=mea.end(); ++i) {
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename MEA::subpopulation_type::population_type propagule_type;
                
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


// make sure resources are moved to multi. check gls for example

//! Performs multicell replication using germ lines. One cells is selected, mutated, and then used to create the appropriate number of cells. Thus, the starting multicell offspring is clonal.
template <typename MEA>
struct mt_gls_propagule : end_of_update_event<MEA> {
    //! Constructor.
    mt_gls_propagule(MEA& mea) : end_of_update_event<MEA>(mea), _df("mt_gls.dat") {
        _df.add_field("update")
        .add_field("mean_rep_time")
        .add_field("mean_res")
        .add_field("mean_multicell_size")
        .add_field("mean_germ_num")
        .add_field("mean_pop_num")
        .add_field("mean_germ_percent")
        .add_field("mean_germ_workload")
        .add_field("mean_germ_workload_var")
        .add_field("mean_soma_workload")
        .add_field("mean_soma_workload_var")
        .add_field("replication_count")
        .add_field("mean_uni_rep_time")
        .add_field("mean_uni_workload")
        .add_field("mean_mc_rep_time")
        .add_field("mean_mc_workload");
        
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~mt_gls_propagule() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(MEA& mea) {
        
        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));

        
        // Replicate!
        int ru = 1;
        if ((mea.current_update() % ru) == 0) {
            
            // See if any subpops have exceeded the resource threshold
            typename MEA::population_type offspring;
            for(typename MEA::iterator i=mea.begin(); i!=mea.end(); ++i) {
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename MEA::subpopulation_type::population_type propagule_type;
                
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
                    typename MEA::subpopulation_type::individual_type germ;
                    int germ_present = false;
                    
                    // If so, setup a new replicate pop.
                    // Find a germ...
                    std::random_shuffle(i->population().begin(), i->population().end(), mea.rng());
                    
                    int germ_count = 0;
                    int pop_count = 0;
                    accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc;
                    accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc;

                    
                    
                    
                    // get a new subpopulation:
                    typename MEA::individual_ptr_type p = mea.make_individual();
                    p->initialize(mea.md());
                    p->reset_rng(mea.rng().seed());
                    
                    int total_workload = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        typename MEA::subpopulation_type::individual_type& org=**j;
                        if (get<GERM_STATUS>(org, true)) {
                            ++germ_count;
                            germ_workload_acc(get<WORKLOAD>(org, 0.0));
                            if (!germ_present){
                                typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(), (*j)->genome().begin()+(*j)->hw().original_size());
                                typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                                
                                inherits_from(**j, *q, *p);
                                
                                mutate(*q,m,*p);
                                
                                p->insert(p->end(), q);
                                
                                germ_present = true;
                            }
                        } else {
                            soma_workload_acc(get<WORKLOAD>(org, 0.0));
                        }
                        ++pop_count;
                        total_workload += get<WORKLOAD>(org, 0.0);
                        
                    }
                    
                    if (!germ_present) continue;
                    
                    if (pop_count == 1) { // track as uni
                        uni_rep_time_acc.push_back(get<MULTICELL_REP_TIME>(*i));
                        uni_workload_acc.push_back(total_workload);
                    } else {
                        mc_rep_time_acc.push_back(get<MULTICELL_REP_TIME>(*i));
                        mc_workload_acc.push_back(total_workload);
                    }
                    
                    pop_num.push_back(pop_count);
                    germ_num.push_back(germ_count);
                    germ_percent.push_back(germ_count/((double) i->population().size())*100.0);
                    germ_workload.push_back(mean(germ_workload_acc));
                    germ_workload_var.push_back(variance(germ_workload_acc));
                    
                    if (germ_count != pop_count) {
                        soma_workload.push_back(mean(soma_workload_acc));
                        soma_workload_var.push_back(variance(soma_workload_acc));
                    } else {
                        soma_workload.push_back(0);
                        soma_workload_var.push_back(0);
                    }


                    
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
            _df.write(mea.current_update());
            
            if (multicell_rep.size() > 0) {
                _df.write(std::accumulate(multicell_rep.begin(), multicell_rep.end(), 0.0)/multicell_rep.size())
                .write(std::accumulate(multicell_res.begin(), multicell_res.end(), 0.0)/multicell_res.size())
                .write(std::accumulate(multicell_size.begin(), multicell_size.end(), 0.0)/multicell_size.size());
            } else {
                _df.write(0.0)
                .write(0.0)
                .write(0.0);
            }
            
            
            if (germ_num.size() > 0) {
                _df.write(std::accumulate(germ_num.begin(), germ_num.end(), 0.0)/germ_num.size())
                .write(std::accumulate(pop_num.begin(), pop_num.end(), 0.0)/pop_num.size())
                .write(std::accumulate(germ_percent.begin(), germ_percent.end(), 0.0)/germ_percent.size())
                .write(std::accumulate(germ_workload.begin(), germ_workload.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(germ_workload_var.begin(), germ_workload_var.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(soma_workload.begin(), soma_workload.end(), 0.0)/soma_workload.size())
                .write(std::accumulate(soma_workload_var.begin(), soma_workload_var.end(), 0.0)/soma_workload.size());
            } else {
                _df.write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .write(0);
                
            }
            
            
            _df.write(num_rep);
            
            if (uni_rep_time_acc.size() > 0) {
                _df.write(std::accumulate(uni_rep_time_acc.begin(), uni_rep_time_acc.end(), 0.0)/uni_rep_time_acc.size())
                .write(std::accumulate(uni_workload_acc.begin(), uni_workload_acc.end(), 0.0)/uni_workload_acc.size());
            } else {
                _df.write(0)
                .write(0);
            }
            
            if (mc_rep_time_acc.size() > 0) {
                _df.write(std::accumulate(mc_rep_time_acc.begin(), mc_rep_time_acc.end(), 0.0)/mc_rep_time_acc.size())
                .write(std::accumulate(mc_workload_acc.begin(), mc_workload_acc.end(), 0.0)/mc_workload_acc.size());
            } else {
                _df.write(0)
                .write(0);
            }
            /*         .add_field("mean_uni_rep_time")
             .add_field("mean_mc_rep_time")
             .add_field("mean_uni_workload")
             .add_field("mean_mc_workload");*/
            _df.endl();
            num_rep = 0;
            multicell_rep.clear();
            multicell_res.clear();
            multicell_size.clear();
            germ_num.clear();
            germ_percent.clear();
            pop_num.clear();
            germ_workload.clear();
            germ_workload_var.clear();
            soma_workload.clear();
            soma_workload_var.clear();
            uni_workload_acc.clear();
            mc_workload_acc.clear();
            uni_rep_time_acc.clear();
            mc_rep_time_acc.clear();
            
        }
        
    }
    
    datafile _df;
    std::deque<double> multicell_rep;
    std::deque<double> multicell_res;
    std::deque<double> multicell_size;
    std::deque<double> germ_num;
    std::deque<double> germ_percent;
    std::deque<double> pop_num;
    std::deque<double> germ_workload;
    std::deque<double> germ_workload_var;
    std::deque<double> soma_workload;
    std::deque<double> soma_workload_var;
    std::deque<double> uni_workload_acc;
    std::deque<double> mc_workload_acc;
    std::deque<double> uni_rep_time_acc;
    std::deque<double> mc_rep_time_acc;

   
    int num_rep;
    
    
};




template <typename MEA>
struct dol_tracking : end_of_update_event<MEA> {
    dol_tracking(MEA& ea) : end_of_update_event<MEA>(ea), _df("dol.dat") {
        _df.add_field("update")
        .add_field("mean_shannon_sum")
        .add_field("mean_shannon_norm")
        .add_field("mean_active_pop")
        .add_field("mean_pop_count");
        
    }
    
    //! Destructor.
    virtual ~dol_tracking() {
    }
    
    //! Track how many task-switches are being performed!
    virtual void operator()(MEA& ea) {
        
        if ((ea.current_update() % 100) == 0) {
            _df.write(ea.current_update());
            double shannon_sum_all = 0;
            double shannon_norm_all = 0;
            double active_pop_all = 0;
            double pop_count_all = 0;
            double num_multis = 0;
            
            for(typename MEA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                num_multis++;
                std::vector< std::vector<double> > pij;
                std::vector<double> pj (9);
                double pop_count = 0;
                double active_pop = 0;
                
                // cycle through orgs and create matrix for shannon mutual information.
                for(typename MEA::subpopulation_type::iterator j=i->begin(); j!=i->end(); ++j) {
                    typename MEA::subpopulation_type::individual_type& org=*j;
                    ++pop_count;
                    std::vector<double> porg (9);
                    porg[0] = get<TASK_NOT>(org,0.0);
                    porg[1] = get<TASK_NAND>(org,0.0);
                    porg[2] = get<TASK_AND>(org,0.0);
                    porg[3] = get<TASK_ORNOT>(org,0.0);
                    porg[4] = get<TASK_OR>(org,0.0);
                    porg[5] = get<TASK_ANDNOT>(org,0.0);
                    porg[6] = get<TASK_NOR>(org,0.0);
                    porg[7] = get<TASK_XOR>(org,0.0);
                    porg[8] = get<TASK_EQUALS>(org,0.0);
                    
                    double total_num_tasks = std::accumulate(porg.begin(), porg.end(), 0);
                    
                    // Normalize the tasks and add to matrix
                    if(total_num_tasks > 0) {
                        for (unsigned int k=0; k<porg.size(); ++k) {
                            porg[k] /= total_num_tasks;
                        }
                        ++active_pop;
                        pij.push_back(porg);
                    }
                }
                
                double shannon_sum = 0.0;
                double shannon_norm = 0.0;
                if (active_pop > 1) {
                    // figure out pj
                    for (unsigned int k=0; k<pj.size(); ++k) {
                        for (int m=0; m<active_pop; ++m) {
                            pj[k] += pij[m][k];
                        }
                        pj[k] /= active_pop;
                    }
                    
                    // compute shannon mutual information based on matrix...
                    double shannon_change = 0.0;
                    double t_pij = 0.0;
                    double t_pi = 1.0/active_pop;
                    double t_pj = 0;
                    double pij_sum = 0.0;
                    // calculate shannon mutual information
                    for (unsigned int i=0; i<active_pop; i++) {
                        for (int j=0; j<pj.size(); j++) {
                            t_pij = pij[i][j]/active_pop;
                            t_pj = pj[j];
                            pij_sum += t_pij;
                            if (t_pi && t_pj && t_pij) {
                                shannon_change= (t_pij * log(t_pij / (t_pi * t_pj)));
                                shannon_sum += shannon_change;
                            }
                        }
                    }
                }
                if ((shannon_sum > 0 ) &&(active_pop > 0)) {
                    shannon_norm = shannon_sum / log((double)active_pop);
                } 
                
                shannon_sum_all += shannon_sum;
                shannon_norm_all += shannon_norm;
                active_pop_all += active_pop;
                pop_count_all += pop_count;
                
                
            }
            _df.write(shannon_sum_all / num_multis)
            .write(shannon_norm_all / num_multis)
            .write(active_pop_all / num_multis)
            .write(pop_count_all / num_multis);
            _df.endl();
        }
    }
    datafile _df;
    
};


#endif
