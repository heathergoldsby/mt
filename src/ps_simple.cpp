#include <ea/digital_evolution.h>
#include <ea/cmdline_interface.h>
#include <ea/subpopulation_founder.h>
#include <ea/line_of_descent.h>
#include <ea/analysis/archive.h>

#include "gls.h"
#include "mt_propagule_orig.h"
#include "ps_simple.h"
#include "multi_birth_selfrep_not_remote_ancestor.h"





using namespace ealib;

//#include "movie_evo_plane.h"


//! Configuration object for an EA.
struct lifecycle : public default_lifecycle {
    
    //! Called as the final step of EA construction (must not depend on configuration parameters)
    template <typename EA>
    void after_initialization(EA& ea) {
        if(ea.isa().size()) {
            return;
        }
        
        using namespace instructions;
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
        append_isa<tx_msg>(ea);
        append_isa<rx_msg>(ea);
        append_isa<bc_msg>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea);
        append_isa<h_alloc>(ea);
        append_isa<h_copy>(ea);
        append_isa<h_divide_soft_parent_reset>(ea);
        append_isa<input>(ea);
        append_isa<output>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<is_neighbor>(ea);
        append_isa<h_divide_remote>(ea);
        append_isa<h_alt_divide>(ea);
        append_isa<become_soma>(ea);
        append_isa<if_germ>(ea);
        append_isa<if_soma>(ea);
        append_isa<get_germ_size>(ea);
        append_isa<if_member_start_propagule>(ea);
        append_isa<if_not_member_start_propagule>(ea);
        
        add_event<gs_inherit_event>(ea);

        
        
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




typedef digital_evolution
< lifecycle
, recombination::asexual
, round_robin
, multibirth_selfrep_not_remote_unfixed_ancestor
, empty_facing_neighbor
, dont_stop
, generate_single_ancestor
> sea_type;


typedef metapopulation
< sea_type
, quiet_nan
, mutation::operators::no_mutation
, quiet_nan
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
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<SCHEDULER_RESOURCE_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<MUTATION_UNIFORM_INT_MIN>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this);
        
        add_option<ANALYSIS_INPUT>(this);
        
        
        // ts specific options
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<GROUP_REP_THRESHOLD>(this);
        
        
    }
    
    virtual void gather_tools() {
        
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<mt_ps_propagule>(ea);
        add_event<germ_soma_based_resources>(ea);
        
//        add_event<size_based_resources>(ea);
        
    }
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
