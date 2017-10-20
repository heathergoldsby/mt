/* gls.cpp
 *
 * This file is part of EALife.
 *
 * Copyright 2012 David B. Knoester, Heather J. Goldsby.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "gls.h"
//#include "subpopulation_founder.h"
#include <ea/digital_evolution/ancestors/multi_birth_selfrep_not_ancestor.h>
#include <ea/subpopulation_founder.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/digital_evolution.h>
#include <ea/cmdline_interface.h>
#include <ea/line_of_descent.h>
#include <ea/generational_models/periodic_competition.h>
#include <ea/generational_models/moran_process.h>
#include <ea/datafiles/fitness.h>
#include <ea/line_of_descent.h>

//! Configuration object for an EA.
struct lifecycle : public default_lifecycle {
    /*! Called after EA initialization.
     
     This is a good place to handle programmatic setup tasks.  E.g., adding
     instructions to a digital evolution ISA, loading external data files,
     and the like.
     */
    template <typename EA>
    void after_initialization(EA& ea) {
        using namespace ealib::instructions;
        append_isa<nop_a>(0,ea);
        append_isa<nop_b>(0,ea);
        append_isa<nop_c>(0,ea);
        append_isa<nop_x>(ea);
        append_isa<mov_head>(ea);
        append_isa<if_label>(ea);
        append_isa<h_search>(ea);
        append_isa<nand>(ea);
        append_isa<push>(ea);
        append_isa<pop>(ea);
        append_isa<swap>(ea);
        append_isa<inc>(ea);
        append_isa<dec>(ea);
        append_isa<tx_msg_check_task>(ea);
        append_isa<tx_msg>(ea);
        append_isa<rx_msg>(ea);
        append_isa<bc_msg>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea);
        append_isa<h_alloc>(ea);
        append_isa<h_copy>(ea);
        append_isa<h_divide>(ea);
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<get_xy>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<become_soma>(ea);
        append_isa<if_germ>(ea);
        append_isa<if_soma>(ea);
        
        
        add_event<task_resource_consumption>(ea);
        //add_event<task_switching_cost>(this, ea);
        //add_event<ts_birth_event>(this,ea);
        add_event<task_mutagenesis>(ea);
        add_event<gs_inherit_event>(ea);

        typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
        typedef typename EA::resource_ptr_type resource_ptr_type;
        
        // Add tasks
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<0> >("and", ea);
        task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<0> >("ornot", ea);
        task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<0> >("or", ea);
        task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<0> >("andnot", ea);
        task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<0> >("nor", ea);
        task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<0> >("xor", ea);
        task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<0> >("equals", ea);
        
        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
        double init_amt = get<RES_INITIAL_AMOUNT>(ea, 0);
        double inflow = get<RES_INFLOW_AMOUNT>(ea,0);
        double outflow = get<RES_OUTFLOW_FRACTION>(ea,0);
        double frac = get<RES_FRACTION_CONSUMED>(ea,0);
//        
//        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
//        resource_ptr_type resA = make_resource("resA", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resB = make_resource("resB", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resC = make_resource("resC", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resD = make_resource("resD", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resE = make_resource("resE", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resF = make_resource("resF", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resG = make_resource("resG", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resH = make_resource("resH", 100.0, 1.0, 0.01, 0.05, ea);
//        resource_ptr_type resI = make_resource("resI", 100.0, 1.0, 0.01, 0.05, ea);
        
        
        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
        resource_ptr_type resA = make_resource("resA", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resB = make_resource("resB", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resC = make_resource("resC", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resD = make_resource("resD", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resE = make_resource("resE", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resF = make_resource("resF", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resG = make_resource("resG", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resH = make_resource("resH", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resI = make_resource("resI", init_amt, inflow, outflow, frac, ea);
        
        put<TASK_MUTATION_MULT>(get<NOT_MUTATION_MULT>(ea), *task_not);
        put<TASK_MUTATION_MULT>(get<NAND_MUTATION_MULT>(ea), *task_nand);
        put<TASK_MUTATION_MULT>(get<AND_MUTATION_MULT>(ea), *task_and);
        put<TASK_MUTATION_MULT>(get<ORNOT_MUTATION_MULT>(ea), *task_ornot);
        put<TASK_MUTATION_MULT>(get<OR_MUTATION_MULT>(ea), *task_or);
        put<TASK_MUTATION_MULT>(get<ANDNOT_MUTATION_MULT>(ea), *task_andnot);
        put<TASK_MUTATION_MULT>(get<NOR_MUTATION_MULT>(ea), *task_nor);
        put<TASK_MUTATION_MULT>(get<XOR_MUTATION_MULT>(ea), *task_xor);
        put<TASK_MUTATION_MULT>(get<EQUALS_MUTATION_MULT>(ea), *task_equals);
        
        task_not->consumes(resA);
        task_nand->consumes(resB);
        task_and->consumes(resC);
        task_ornot->consumes(resD);
        task_or->consumes(resE);
        task_andnot->consumes(resF);
        task_nor->consumes(resG);
        task_xor->consumes(resH);
        task_equals->consumes(resI);
        

        
    }
    
};

template <typename T>
struct subpop_trait : subpopulation_founder_trait<T>, fitness_trait<T> {
    typedef subpopulation_founder_trait<T> parent1_type;
    typedef fitness_trait<T> parent2_type;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp("subpopulation_founder_trait", boost::serialization::base_object<parent1_type>(*this));
        ar & boost::serialization::make_nvp("fitness_trait", boost::serialization::base_object<parent2_type>(*this));
    }
};

/*! Artificial life simulation definition.
 */
typedef digital_evolution
< lifecycle
, recombination::asexual
, round_robin
, multibirth_selfrep_not_ancestor
, empty_neighbor
, dont_stop
, generate_single_ancestor
> sea_type;

typedef metapopulation
< sea_type
, quiet_nan
, mutation::operators::no_mutation
, recombination::no_recombination
, generational_models::isolated_subpopulations
, ancestors::default_subpopulation
, dont_stop
, fill_metapopulation
, default_lifecycle
, subpop_trait
> mea_type;



/*!
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    
    
    virtual void gather_options() {
        add_option<SPATIAL_X>(this);
        add_option<SPATIAL_Y>(this);
        add_option<METAPOPULATION_SIZE>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        add_option<SCHEDULER_RESOURCE_SLICE>(this);
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<MUTATION_UNIFORM_INT_MIN>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        // ts specific options
        add_option<GROUP_REP_THRESHOLD>(this);
        add_option<TASK_SWITCHING_COST>(this);
        add_option<LAST_TASK>(this);
        add_option<NUM_SWITCHES>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<START_UPDATE>(this);

        
        // gls specific options
        add_option<TASK_MUTATION_PER_SITE_P>(this);
//        add_option<GERM_MUTATION_PER_SITE_P>(this);
//        add_option<GROUP_REP_THRESHOLD>(this);
        
        // gls mutation multipliers
        add_option<NOT_MUTATION_MULT>(this);
        add_option<NAND_MUTATION_MULT>(this);
        add_option<AND_MUTATION_MULT>(this);
        add_option<ORNOT_MUTATION_MULT>(this);
        add_option<OR_MUTATION_MULT>(this);
        add_option<ANDNOT_MUTATION_MULT>(this);
        add_option<NOR_MUTATION_MULT>(this);
        add_option<XOR_MUTATION_MULT>(this);
        add_option<EQUALS_MUTATION_MULT>(this);
        
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<gls_replication>(ea);

//        add_event<ts_replication>(this,ea);
        add_event<task_performed_tracking>(ea);
        add_event<task_switch_tracking>(ea);
        add_event<subpopulation_founder_event>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
