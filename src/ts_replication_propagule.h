//
//  ts_replication_propagule.h
//  ps
//
//  Created by Heather Goldsby on 1/25/17.
//  Copyright © 2017 Michigan State University. All rights reserved.
//

#ifndef ts_replication_propagule_h
#define ts_replication_propagule_h





#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/configurable_mutation.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>
#include <ea/metadata.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution.h>

using namespace ealib;
using namespace boost::accumulators;


template <typename EA>
struct ts_replication_propagule_hetero : end_of_update_event<EA> {
    //! Constructor.
    ts_replication_propagule_hetero(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ts_replication_propagule_hetero() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& mea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator k=mea.begin(); k!=mea.end(); ++k) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (k->population().size() < 2) continue;
            
            if (exists<GROUP_RESOURCE_UNITS>(*k) &&
                (get<GROUP_RESOURCE_UNITS>(*k) > get<GROUP_REP_THRESHOLD>(*k))){
                
                // the number of parents selected is the propagule size or 1, if
                // the propagule's composition is clonal.
                std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
                assert(prop_size > 0);
                if (prop_size > k->size()) {
                    prop_size = k->size();
                }
                
                
                // get a new subpopulation:
                typename EA::individual_ptr_type p = mea.make_individual();
                p->initialize(mea.md());
                p->reset_rng(mea.rng().seed());
                
                
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename EA::subpopulation_type::population_type propagule_type;
                propagule_type propagule;
                mea.rng().sample_without_replacement(k->population().begin(),
                                                     k->population().end(),
                                                     std::back_inserter(propagule),
                                                     prop_size);
                
                
                int pop_size = get<POPULATION_SIZE>(mea);
                
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
                
                std::vector<int> used_pos;
                std::vector<int> used_pos_with_avail_neighbors;
                
                
                int max_x =get<SPATIAL_X>(mea);
                
                // now add a new individual built from each of the propagules to the
                // subpopulation:
                for(typename propagule_type::iterator i=propagule.begin(); i!=propagule.end(); ++i) {
                    // grab the original part of the propagule's genome; note that it could have been
                    // changed (implicit-like mutations):
                    typename EA::subpopulation_type::genome_type r((*i)->genome().begin(),
                                                                    (*i)->genome().begin()+(*i)->hw().original_size());
                    typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                    
                    inherits_from(**i, *q, *p);
                    
                    mutate(*q,m,*p);
                    
                    
                    if (used_pos.size() == 0) {
                        int pos = mea.rng().uniform_integer(0,pop_size);
                        used_pos.push_back(pos);
                        used_pos_with_avail_neighbors.push_back(pos);
                        p->insert_at(p->end(),q, p->env().location(pos).position());
                    } else {
                        
                        //
                        //            int max_x =get<SPATIAL_X>(mea);
                        //
                        
                        typename EA::subpopulation_type::individual_ptr_type o(q);
                        
                        bool not_placed = true;
                        int place = -1;
                        
                        while (not_placed) {
                            int n = mea.rng().uniform_integer(0,used_pos_with_avail_neighbors.size());
                            int neighbor = used_pos_with_avail_neighbors[n];
                            int dir_try = 0;
                            int dir = mea.rng().uniform_integer(0,3);
                            
                            
                            
                            while (dir_try <= 4) {
                                
                                dir = (dir + 1);
                                if (dir > 3) { dir = 0; }
                                dir_try++;
                                
                                // N
                                if (dir == 0) {
                                    if ((neighbor - max_x) < 0) { continue; }
                                    place = neighbor - max_x;
                                } else if (dir == 1) { // E
                                    if ((neighbor % max_x) == (max_x -1)) { continue; }
                                    place = neighbor + 1;
                                } else if (dir == 2) { // S
                                    if ((neighbor + max_x) >= get<POPULATION_SIZE>(mea)) { continue; }
                                    place = neighbor + max_x;
                                } else if (dir == 3) { // W
                                    if ((neighbor % max_x) == 0) { continue; }
                                    place = neighbor - 1;
                                }
                                
                                // already used
                                if (std::find(used_pos.begin(), used_pos.end(), place) != used_pos.end()) {
                                    place = -1;
                                    continue;
                                } else {
                                    break;
                                }
                                
                            }
                            
                            // not placed, then try again
                            if (place == -1) {
                                used_pos_with_avail_neighbors.erase(used_pos_with_avail_neighbors.begin()+n);
                                continue;
                            }
                            
                            
                            
                            p->insert_at(p->end(),q, p->env().location(place).position());
                            used_pos.push_back(place);
                            used_pos_with_avail_neighbors.push_back(place);
                            not_placed = false;
                        }
                        
                        //}
                    }
                }
                
                
                // reset resource units
                k->resources().reset();
                put<GROUP_RESOURCE_UNITS>(0,*k);
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*k.base());
                offspring.insert(offspring.end(),p);
                inherits(parent_pop, offspring_pop, mea);

            }

        }
        
        get<NUM_GROUP_REPLICATIONS>(mea,0) += offspring.size();
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<METAPOPULATION_SIZE>(mea) - offspring.size();
            
            typename EA::population_type survivors;
            select_n<selection::random< > >(mea.population(), survivors, n, mea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(mea.population(), survivors);
        }
        
    }
    
    
    
};

template <typename EA>
struct ts_replication_propagule : end_of_update_event<EA> {
    //! Constructor.
    ts_replication_propagule(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ts_replication_propagule() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& mea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator k=mea.begin(); k!=mea.end(); ++k) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (k->population().size() < 2) continue;
            
            if (exists<GROUP_RESOURCE_UNITS>(*k) &&
                (get<GROUP_RESOURCE_UNITS>(*k) > get<GROUP_REP_THRESHOLD>(*k))){
                
                // the number of parents selected is the propagule size or 1, if
                // the propagule's composition is clonal.
                std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
                assert(prop_size > 0);
                if (prop_size > k->size()) {
                    prop_size = k->size();
                }
                
                
                // get a new subpopulation:
                typename EA::individual_ptr_type p = mea.make_individual();
                p->initialize(mea.md());
                p->reset_rng(mea.rng().seed());
                
                
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename EA::subpopulation_type::population_type propagule_type;
                propagule_type propagule;
                mea.rng().sample_without_replacement(k->population().begin(),
                                                     k->population().end(),
                                                     std::back_inserter(propagule),
                                                     prop_size);
                
                
                int pop_size = get<POPULATION_SIZE>(mea);
                
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
                
                std::vector<int> used_pos;
                std::vector<int> used_pos_with_avail_neighbors;
                
                
                int max_x =get<SPATIAL_X>(mea);
                
                // now add a new individual built from each of the propagules to the
                // subpopulation:
                for(typename propagule_type::iterator i=propagule.begin(); i!=propagule.end(); ++i) {
                    // grab the original part of the propagule's genome; note that it could have been
                    // changed (implicit-like mutations):
                    typename EA::subpopulation_type::genome_type r((*i)->genome().begin(),
                                                                   (*i)->genome().begin()+(*i)->hw().original_size());
                    typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                    
                    inherits_from(**i, *q, *p);
                    
                    mutate(*q,m,*p);
                    
                    
                        int pos = mea.rng().uniform_integer(0,pop_size);
                        used_pos.push_back(pos);
                        used_pos_with_avail_neighbors.push_back(pos);
                        p->insert_at(p->end(),q, p->env().location(pos).position());

                        for (int k=1; k<get<NUM_PROPAGULE_GERM>(mea); ++k) {

                            typename EA::subpopulation_type::individual_ptr_type o(q);
                        
                            bool not_placed = true;
                            int place = -1;
                        
                            while (not_placed) {
                                int n = mea.rng().uniform_integer(0,used_pos_with_avail_neighbors.size());
                                int neighbor = used_pos_with_avail_neighbors[n];
                                int dir_try = 0;
                                int dir = mea.rng().uniform_integer(0,3);
                            
                            
                            
                                while (dir_try <= 4) {
                                
                                    dir = (dir + 1);
                                    if (dir > 3) { dir = 0; }
                                    dir_try++;
                                
                                    // N
                                    if (dir == 0) {
                                        if ((neighbor - max_x) < 0) { continue; }
                                        place = neighbor - max_x;
                                    } else if (dir == 1) { // E
                                        if ((neighbor % max_x) == (max_x -1)) { continue; }
                                        place = neighbor + 1;
                                    } else if (dir == 2) { // S
                                        if ((neighbor + max_x) >= get<POPULATION_SIZE>(mea)) { continue; }
                                        place = neighbor + max_x;
                                    } else if (dir == 3) { // W
                                        if ((neighbor % max_x) == 0) { continue; }
                                        place = neighbor - 1;
                                    }
                                
                                    // already used
                                    if (std::find(used_pos.begin(), used_pos.end(), place) != used_pos.end()) {
                                        place = -1;
                                        continue;
                                    } else {
                                        break;
                                    }
                                
                                }
                            
                                // not placed, then try again
                                if (place == -1) {
                                    used_pos_with_avail_neighbors.erase(used_pos_with_avail_neighbors.begin()+n);
                                    continue;
                                }
                            
                            
                            
                                p->insert_at(p->end(),q, p->env().location(place).position());
                                used_pos.push_back(place);
                            used_pos_with_avail_neighbors.push_back(place);
                            not_placed = false;
                            } // end while
                        
                        } // end for
                }
                
                
                // reset resource units
                k->resources().reset();
                put<GROUP_RESOURCE_UNITS>(0,*k);
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*k.base());
                offspring.insert(offspring.end(),p);
                inherits(parent_pop, offspring_pop, mea);
                
            }
            
        }
        
        get<NUM_GROUP_REPLICATIONS>(mea,0) += offspring.size();
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<METAPOPULATION_SIZE>(mea) - offspring.size();
            
            typename EA::population_type survivors;
            select_n<selection::random< > >(mea.population(), survivors, n, mea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(mea.population(), survivors);
        }
        
    }
    
    
    
};


#endif
