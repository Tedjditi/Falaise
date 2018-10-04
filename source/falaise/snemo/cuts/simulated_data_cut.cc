// falaise/snemo/cuts/simulated_data_cut.cc

// Ourselves:
#include <falaise/snemo/cuts/simulated_data_cut.h>

/// Standard library:
#include <sstream>
#include <stdexcept>

// Third party:
// - Bayeux/datatools:
#include <bayeux/datatools/properties.h>
#include <bayeux/datatools/things.h>

// This project :
#include <falaise/snemo/datamodels/data_model.h>

namespace snemo {

namespace cut {

// Registration instantiation macro :
CUT_REGISTRATION_IMPLEMENT(simulated_data_cut, "snemo::cut::simulated_data_cut")

void simulated_data_cut::_set_defaults() {
  _SD_label_ = "";
  _mode_ = MODE_UNDEFINED;
  _flag_name_ = "";
  _hit_category_ = "";
  _hit_category_range_min_ = -1;
  _hit_category_range_max_ = -1;
  _hit_property_logic_ = "";
  _gamma_epsilonZ_ = 0;
  _gamma_epsilonY_ = 0;
  _first_ = -1;
  _last_ = -1;
  _wall_ = 0;
  return;
}

void simulated_data_cut::set_SD_label(const std::string& SD_label_) {
  _SD_label_ = SD_label_;
  return;
}

const std::string& simulated_data_cut::get_SD_label() const { return _SD_label_; }

uint32_t simulated_data_cut::get_mode() const { return _mode_; }

bool simulated_data_cut::is_mode_flag() const { return _mode_ & MODE_FLAG; }

bool simulated_data_cut::is_mode_range_hit_category() const {
  return _mode_ & MODE_RANGE_HIT_CATEGORY;
}

bool simulated_data_cut::is_mode_has_hit_category() const { return _mode_ & MODE_HAS_HIT_CATEGORY; }

bool simulated_data_cut::is_mode_has_hit_property() const { return _mode_ & MODE_HAS_HIT_PROPERTY; }

bool simulated_data_cut::is_mode_gamma_has_gone_straight() const { return _mode_ & MODE_GAMMA_HAS_GONE_STRAIGHT; }

  bool simulated_data_cut::is_mode_has_hit_on_wall() const { return _mode_ & MODE_HAS_HIT_ON_WALL; }
  
void simulated_data_cut::set_flag_name(const std::string& flag_name_) {
  _flag_name_ = flag_name_;
  return;
}

const std::string& simulated_data_cut::get_flag_name() const { return _flag_name_; }

simulated_data_cut::simulated_data_cut(datatools::logger::priority logger_priority_)
    : cuts::i_cut(logger_priority_) {
  _set_defaults();
  return;
}

simulated_data_cut::~simulated_data_cut() {
  if (is_initialized()) this->simulated_data_cut::reset();
  return;
}

void simulated_data_cut::reset() {
  _set_defaults();
  this->i_cut::_reset();
  this->i_cut::_set_initialized(false);
  return;
}

void simulated_data_cut::initialize(const datatools::properties& configuration_,
                                    datatools::service_manager& /* service_manager_ */,
                                    cuts::cut_handle_dict_type& /* cut_dict_ */) {
  DT_THROW_IF(is_initialized(), std::logic_error,
              "Cut '" << get_name() << "' is already initialized ! ");

  this->i_cut::_common_initialize(configuration_);

  if (_SD_label_.empty()) {
    if (configuration_.has_key("SD_label")) {
      set_SD_label(configuration_.fetch_string("SD_label"));
    } else {
      set_SD_label(snemo::datamodel::data_info::default_simulated_data_label());
    }
  }

  if (_mode_ == MODE_UNDEFINED) {
    if (configuration_.has_flag("mode.flag")) {
      _mode_ |= MODE_FLAG;
    }

    if (configuration_.has_flag("mode.has_hit_category")) {
      _mode_ |= MODE_HAS_HIT_CATEGORY;
    }

    if (configuration_.has_flag("mode.range_hit_category")) {
      _mode_ |= MODE_RANGE_HIT_CATEGORY;
    }

    if (configuration_.has_flag("mode.has_hit_property")) {
      _mode_ |= MODE_HAS_HIT_PROPERTY;
    }

    if (configuration_.has_flag("mode.gamma_has_gone_straight")) {
      _mode_ |= MODE_GAMMA_HAS_GONE_STRAIGHT;
    }

    if (configuration_.has_flag("mode.has_hit_on_wall")) {
      _mode_ |= MODE_HAS_HIT_ON_WALL;
    }
    
    DT_THROW_IF(_mode_ == MODE_UNDEFINED, std::logic_error,
                "Missing at least a 'mode.XXX' property ! ");

    // mode FLAG:
    if (is_mode_flag()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using FLAG mode...");
      DT_THROW_IF(!configuration_.has_key("flag.name"), std::logic_error,
                  "Missing 'flag.name' property !");
      set_flag_name(configuration_.fetch_string("flag.name"));
    }  // end if is_mode_flag

    // mode HAS_HIT_CATEGORY:
    if (is_mode_has_hit_category()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using HAS_HIT_CATEGORY mode...");
      DT_THROW_IF(!configuration_.has_key("has_hit_category.category"), std::logic_error,
                  "Missing 'has_hit_category.category' property !");
      _hit_category_ = configuration_.fetch_string("has_hit_category.category");
    }  // end if is_mode_has_hit_category

    // mode RANGE_HIT_CATEGORY:
    if (is_mode_range_hit_category()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using RANGE_HIT_CATEGORY mode...");
      DT_THROW_IF(!configuration_.has_key("range_hit_category.category"), std::logic_error,
                  "Missing 'range_hit_category.category' property !");
      _hit_category_ = configuration_.fetch_string("range_hit_category.category");

      int count = 0;
      if (configuration_.has_key("range_hit_category.min")) {
        const int nmin = configuration_.fetch_integer("range_hit_category.min");
        DT_THROW_IF(nmin < 0, std::logic_error, "Invalid min number of hits(" << nmin << ") !");
        _hit_category_range_min_ = nmin;
        count++;
      }
      if (configuration_.has_key("range_hit_category.max")) {
        const int nmax = configuration_.fetch_integer("range_hit_category.max");
        DT_THROW_IF(nmax < 0, std::logic_error, "Invalid max number of hits(" << nmax << ") !");
        _hit_category_range_max_ = nmax;
        count++;
      }
      DT_THROW_IF(count == 0, std::logic_error,
                  "Missing 'range_hit_category.min' or 'range_hit_category.max' property !");
      if (count == 2 && _hit_category_range_min_ >= 0 && _hit_category_range_max_ >= 0) {
        DT_THROW_IF(_hit_category_range_min_ > _hit_category_range_max_, std::logic_error,
                    "Invalid 'range_hit_category.min' > 'range_hit_category.max' values !");
      }
    }  // end if is_mode_range_hit_category

    // MODE GAMMA_HAS_GONE_STRAIGHT
    if (is_mode_gamma_has_gone_straight()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using GAMMA_HAS_GONE_STRAIGHT mode...");
      DT_THROW_IF(!configuration_.has_key("gamma_has_gone_straight_category.category"), std::logic_error, "Missing 'gamma_has_gone_straight_category.category' property !");
      _hit_category_=configuration_.fetch_string("gamma_has_gone_straight_category.category");

      if (configuration_.has_key("gamma_has_gone_straight_category.epsilonZ")) {
	const double epsilonZ = configuration_.fetch_real("gamma_has_gone_straight_category.epsilonZ");
	DT_THROW_IF(epsilonZ < 0, std::logic_error, "Missing 'gamma_has_gone_straight_category.epsilonZ' property !");
	_gamma_epsilonZ_ = epsilonZ;
      }

      if (configuration_.has_key("gamma_has_gone_straight_category.epsilonY")) {
	const double epsilonY = configuration_.fetch_real("gamma_has_gone_straight_category.epsilonY");
	DT_THROW_IF(epsilonY < 0, std::logic_error, "Missing 'gamma_has_gone_straight_category.epsilonY' property !");
	_gamma_epsilonY_ = epsilonY;
      }
      
    }

    // MODE_HAS_HIT_ON_WALL
    if (is_mode_has_hit_on_wall()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using HAS_HIT_ON_WALL mode ...");
      
      if (configuration_.has_key("has_hit_on_wall_category.wall")) {
	const double wall = configuration_.fetch_real("has_hit_on_wall_category.wall");
	DT_THROW_IF(wall == 0, std::logic_error, "Missing 'has_hit_on_wall_category.wall' property !");
	_wall_ = wall;
      }

    }
    
    if (is_mode_has_hit_property()) {
      DT_LOG_DEBUG(get_logging_priority(), "Using HAS_HIT_PROPERTY mode...");
      DT_THROW_IF(!configuration_.has_key("has_hit_property.category"), std::logic_error,
                  "Missing 'has_hit_property.category' property !");
      _hit_category_ = configuration_.fetch_string("has_hit_property.category");

      DT_THROW_IF(!configuration_.has_key("has_hit_property.keys"), std::logic_error,
                  "Missing 'has_hit_property.keys' property !");
      std::vector<std::string> keys;
      configuration_.fetch("has_hit_property.keys", keys);

      // Look for the logical relation between properties
      if (keys.size() > 1) {
        DT_THROW_IF(!configuration_.has_key("has_hit_property.logic"), std::logic_error,
                    "Missing 'has_hit_property.logic' property !");
        _hit_property_logic_ = configuration_.fetch_string("has_hit_property.logic");
        DT_THROW_IF(_hit_property_logic_ != "and" && _hit_property_logic_ != "or", std::logic_error,
                    "Logic operation '" << _hit_property_logic_ << "' is unkown or not supported!");
      }

      for (std::vector<std::string>::const_iterator ikey = keys.begin(); ikey != keys.end();
           ++ikey) {
        const std::string check_key = "has_hit_property." + *ikey + ".values";
        DT_THROW_IF(!configuration_.has_key(check_key), std::logic_error,
                    "Missing '" << check_key << "' property !");
        // 2012-11-25 XG: It is maybe better to store things in an
        // 'unordered' way (thus store property_keys in a
        // independent list) to let the user define the best way to
        // look for property keys. For example, it is more efficient
        // to look for a given 'creator_process' and then
        // 'creator_category' since there can be a lot of hits
        // fulfilling 'creator_category' condition and only few for
        // 'creator_process'.
        configuration_.fetch(check_key, _hit_property_values_[*ikey]);
      }
    }  // end if is_mode_hit_has_property
  }

  this->i_cut::_set_initialized(true);
  return;
}

  bool simulated_data_cut::is_inside_volume(const geomtools::vector_3d& Position) {

    if ( Position.getX() < 425 && Position.getX() > -425 ) {
      if ( Position.getY() < 2463 && Position.getY() > -2463) {
	if (Position.getZ() < 1682 && Position.getZ() > -1682) {
	  return true;
	}
      }
    }
    return false;
  }

bool simulated_data_cut::is_colinear(mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_min,mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_max) {
  const geomtools::vector_3d& P_1 = it_t_min->get().get_position_stop() - it_t_min->get().get_position_start();
  const geomtools::vector_3d& P_2 = it_t_max->get().get_position_start() - it_t_min->get().get_position_start();


  if (P_1.x() == 0) return false;
  
  const double coeff = P_2.x() / P_1.x() ;

  double Y;
  double Z;
  
  Y = coeff * P_1.getY();
  Z = coeff * P_1.getZ();

  
  if ( abs( P_2.getY() - Y ) > _gamma_epsilonY_ ) return false;
  if ( abs( P_2.getZ() - Z ) > _gamma_epsilonZ_ ) return false;

  return true; 
  
}
  
  bool simulated_data_cut::is_colinear(const mctools::simulated_data::hit_handle_collection_type& step_handle_track,const mctools::simulated_data::hit_handle_collection_type& step_handle_calo, int track_id) {

  mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_track = step_handle_track.begin();

  mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_calo = step_handle_calo.begin();

  mctools::simulated_data::hit_handle_collection_type::const_iterator it_step = step_handle_track.begin();

  mctools::simulated_data::hit_handle_collection_type::const_iterator it_step_calo = step_handle_calo.begin();
  
  bool check_is_colinear = false;
  double t_track = 10000;
  double t_calo = 10000;
  double t_check = 10000;
  //unsigned int count_ = 1;
  
  DT_LOG_DEBUG(get_logging_priority(), "Entering Search of First step in Track VOLUME...");

  //std::cout << "Number of Step = " << step_handle.size() << std::endl << std::endl;
  
  while (it_step != step_handle_track.end()){

    DT_LOG_DEBUG(get_logging_priority(), "Getting time of step...");

    //std::cout << "STEP NUMBER " << count_ << std::endl << std::endl;
    //count_++;
      
    t_check = it_step->get().get_time_start();

    DT_LOG_DEBUG(get_logging_priority(), "Checkign Track_id of step...");

    if (it_step->get().get_track_id() == track_id) {

      DT_LOG_DEBUG(get_logging_priority(), "Entering Check is inside Volume...");
      if (simulated_data_cut::is_inside_volume(it_step->get().get_position_start())) {
	
	if (t_check < t_track) {
	  t_track = t_check;
	  it_t_track = it_step;
	}

	if (t_check == 0 ) break;
	
      }
      
      DT_LOG_DEBUG(get_logging_priority(), "Exiting Check is inside Volume...");
      
    }

    it_step++;
    //if (count_ > step_handle.size() ) break;

  }

  DT_LOG_DEBUG(get_logging_priority(), "End of search of First Step...");
  DT_LOG_DEBUG(get_logging_priority(), "Entering Search of First step in CALO...");

  
  while (it_step_calo != step_handle_calo.end()){

    DT_LOG_DEBUG(get_logging_priority(), "CALO: Getting time of step...");

    //std::cout << "STEP NUMBER " << count_ << std::endl << std::endl;
    //count_++;
      
    t_check = it_step_calo->get().get_time_start();

    	
	if (t_check < t_calo) {
	  t_calo = t_check;
	  it_t_calo = it_step_calo;
	}
	
    
      
    it_step_calo++;
      //if (count_ > step_handle.size() ) break;

  }

  DT_LOG_DEBUG(get_logging_priority(), "End of search of CALO First Step...");
  
  if ( t_track == 10000 || t_calo == 10000) {
    DT_LOG_DEBUG(get_logging_priority(), "NO CANDIDATE FOUND...");
    return check_is_colinear;
  }

  DT_LOG_DEBUG(get_logging_priority(), "Found First and last step...");

   // DT_LOG_DEBUG(get_logging_priority(), "ID of First and Last step selected...");

   // std::cout << "FIRST = HIT n° " << it_t_track->get().get_hit_id() << "  t_min = " << it_t_track->get().get_time_start() << std::endl << "LAST = HIT n° " << it_t_calo->get().get_hit_id() << "  t_max = " << it_t_calo->get().get_time_start() << std::endl << std::endl;
  
  DT_LOG_DEBUG(get_logging_priority(), "Entering Check of colinearity...");
  check_is_colinear = is_colinear(it_t_track,it_t_calo);

  DT_LOG_DEBUG(get_logging_priority(), "Exiting Check of colinearity...");

  //std::cout << "IS_COLINEAR ? " << check_is_colinear << std::endl;
  
  return check_is_colinear;
}

  bool simulated_data_cut::is_on_wall( mctools::simulated_data::hit_handle_collection_type::const_iterator it_step, double wall) {

    DT_LOG_DEBUG(get_logging_priority(),"Entering IS_ON_WALL...");
    
    const geomtools::vector_3d& Hit_calo = it_step->get().get_position_start();

    if (Hit_calo.x() == 0) return false;

    DT_LOG_DEBUG(get_logging_priority(),"Calcul of X_check and Wall_check...");
    
    double X_check = round( Hit_calo.x() / abs(Hit_calo.x()) );
    double Wall_check = round ( wall / abs(wall) );

    //std::cout <<" X_check = " << X_check << "  Wall_check = " << Wall_check << "   Wall value = " << wall << std::endl << std::endl;
    
    
    if ( X_check == Wall_check ) return true;

    return false;


  }
  
  bool simulated_data_cut::has_hit_on_target_wall(const mctools::simulated_data::hit_handle_collection_type& step_handle, double wall) {

    mctools::simulated_data::hit_handle_collection_type::const_iterator it_step = step_handle.begin();
    mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_calo = step_handle.begin();
    
    DT_LOG_DEBUG(get_logging_priority(), " Search of first CALO HIT...");

    double t_check = 10000;
    double t_calo = 10000;
    
    while (it_step != step_handle.end()){

      DT_LOG_DEBUG(get_logging_priority(), "CALO: Getting time of step...");
      
      t_check = it_step->get().get_time_start();
      
      if (t_check < t_calo) {
	t_calo = t_check;
	it_t_calo = it_step;
      }
	
      it_step++;
      //if (count_ > step_handle.size() ) break;

    }

    DT_LOG_DEBUG(get_logging_priority(), "Check is_on_wall...");

    return is_on_wall(it_t_calo, wall);
    
  }
  
int simulated_data_cut::_accept() {
  int cut_returned = cuts::SELECTION_INAPPLICABLE;

  // Get event record
  const datatools::things& ER = get_user_data<datatools::things>();

  if (!ER.has(_SD_label_)) {
    DT_LOG_DEBUG(get_logging_priority(), "Event record has no '" << _SD_label_ << "' bank !");
    return cut_returned;
  }

  // Get simulated data bank
  const mctools::simulated_data& SD = ER.get<mctools::simulated_data>(_SD_label_);

  // Check if the simulated data has a property flag with a specific name :
  bool check_flag = true;
  if (is_mode_flag()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running FLAG mode...");
    const bool check = SD.get_properties().has_flag(_flag_name_);
    if (!check) {
      check_flag = false;
    }
  }

  // Check if the simulated data has some specific category of hits :
  bool check_has_hit_category = true;
  if (is_mode_has_hit_category()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running HAS_HIT_CATEGORY mode...");
    const bool check = SD.has_step_hits(_hit_category_);
    if (!check) {
      check_has_hit_category = false;
    }
  }

  // Check if the simulated data has some specific category of hits :
  bool check_range_hit_category = true;
  if (is_mode_range_hit_category()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running RANGE_HIT_CATEGORY mode...");
    if (!SD.has_step_hits(_hit_category_)) {
      DT_LOG_DEBUG(get_logging_priority(),
                   "Simulated data has no step hit from '" << _hit_category_ << "' category");
      return cuts::SELECTION_INAPPLICABLE;
    }
    const size_t nhits = SD.get_number_of_step_hits(_hit_category_);
    bool check = true;
    if (_hit_category_range_min_ >= 0) {
      if (nhits < (size_t)_hit_category_range_min_) {
        check = false;
      }
    }
    if (_hit_category_range_max_ >= 0) {
      if (nhits > (size_t)_hit_category_range_max_) {
        check = false;
      }
    }
    if (!check) {
      check_range_hit_category = false;
    }
  }

  // Check if the simulated data has gamma with no deviation
  bool check_gamma_has_gone_straight = true;
  if (is_mode_gamma_has_gone_straight()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running GAMMA_HAS_GONE_STRAIGHT mode...");
    const mctools::simulated_data::primary_event_type& primary_event_ = SD.get_primary_event();

    const bool check_calo_hit = SD.has_step_hits("calo");
    
    if (!check_calo_hit) {
      DT_LOG_DEBUG(get_logging_priority(), "NO CALO HIT...");
      check_gamma_has_gone_straight = false;
    }else {
    
      DT_LOG_DEBUG(get_logging_priority(), "Getting all step...");
      // std::vector< std::string > the_categories;
      // SD.get_step_hits_categories(the_categories);
      // for(auto &categ : the_categories)
      //   std::cout<<categ<<std::endl;
      //the_categories.clear();

      
      const mctools::simulated_data::hit_handle_collection_type& _step_hits_handle_track =
	SD.get_step_hits(_hit_category_);
      
      const mctools::simulated_data::hit_handle_collection_type& _step_hits_handle_calorimeter = SD.get_step_hits("calo");
      
    
      DT_LOG_DEBUG(get_logging_priority(), "Entering Loop on ALL particle simulated..."); 
      
      for (unsigned int it_pp = 0 ;it_pp < primary_event_.get_number_of_particles();it_pp++){
	
	DT_LOG_DEBUG(get_logging_priority(),"Checking Primary event type ");

	if (!primary_event_.get_particle(it_pp).is_gamma()) continue;

	DT_LOG_DEBUG(get_logging_priority(),"Entering Is Colinear");

	check_gamma_has_gone_straight = simulated_data_cut::is_colinear(_step_hits_handle_track,_step_hits_handle_calorimeter,it_pp+1);

	DT_LOG_DEBUG(get_logging_priority(),"Exiting Selection or Looping");

	if (!check_gamma_has_gone_straight) break;

      }
    }
    
    DT_LOG_DEBUG(get_logging_priority(), "Exiting GAMMA_HAS_GONE_STRAIGHT mode...");
  }

  // check if there is a hit on target wall
  bool check_has_hit_on_wall = true;
  if (is_mode_has_hit_on_wall()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running HAS_HIT_ON_WALL mode...");
    const bool check_calo = SD.has_step_hits("calo");

    //std::cout<<"NUMBER OF HIT = " << SD.get_step_hits("__visu.tracks").size() << std::endl << " Wall value = " << _wall_ << std::endl;

    // std::vector< std::string > the_categories;
    // SD.get_step_hits_categories(the_categories);
    // for(auto &categ : the_categories)
    //    std::cout<<categ<<std::endl;
    // the_categories.clear();
    // std::cout << std::endl;
    
    if (!check_calo) {
      DT_LOG_DEBUG(get_logging_priority(), "NO CALO HIT...");
      check_has_hit_on_wall = false;
    } else {

      DT_LOG_DEBUG(get_logging_priority(), "CALO HIT FOUND...");
      
      const mctools::simulated_data::hit_handle_collection_type& _step_hits_handle_calo = SD.get_step_hits("calo");
      DT_LOG_DEBUG(get_logging_priority(), "EnteringG HAS_HIT_ON_WALL check...");
      
      check_has_hit_on_wall = simulated_data_cut::has_hit_on_target_wall(_step_hits_handle_calo,_wall_);

      DT_LOG_DEBUG(get_logging_priority(), "Exiting HAS_HIT_ON_WALL check...");
      
    }

    
    DT_LOG_DEBUG(get_logging_priority(), "Exiting HAS_HIT_ON_WALL mode...");
  }
 
  
  // Check if the simulated step hits have some property value :
  bool check_has_hit_property = true;
  if (is_mode_has_hit_property()) {
    DT_LOG_DEBUG(get_logging_priority(), "Running HAS_HIT_PROPERTY mode...");
    if (!SD.has_step_hits(_hit_category_)) {
      DT_LOG_DEBUG(get_logging_priority(),
                   "Simulated data has no step hit from '" << _hit_category_ << "' category");
      return cuts::SELECTION_INAPPLICABLE;
    }
    const mctools::simulated_data::hit_handle_collection_type& the_step_hits =
        SD.get_step_hits(_hit_category_);

    // Iterators
    mctools::simulated_data::hit_handle_collection_type::const_iterator istart =
        the_step_hits.begin();
    mctools::simulated_data::hit_handle_collection_type::const_iterator istop = the_step_hits.end();
    property_values_dict_type::const_iterator iprop = _hit_property_values_.begin();

    while (iprop != _hit_property_values_.end()) {
      const std::string& a_key = iprop->first;
      const std::vector<std::string>& the_values = iprop->second;
      DT_LOG_TRACE(get_logging_priority(), "Checking property '" << a_key);

      geomtools::base_hit::has_string_property_predicate str_pred(a_key, the_values);
      // wrapper predicates :
      datatools::mother_to_daughter_predicate<geomtools::base_hit, mctools::base_step_hit> pred(
          str_pred);
      datatools::handle_predicate<mctools::base_step_hit> pred_via_handle(pred);

      // Update iterator position
      mctools::simulated_data::hit_handle_collection_type::const_iterator ifound =
          std::find_if(istart, istop, pred_via_handle);
      if (ifound == the_step_hits.end()) {
        check_has_hit_property = false;
        if (_hit_property_logic_ == "and") break;
        // Go to the next property (OR mode)
        istart = the_step_hits.begin();
        istop = the_step_hits.end();
        iprop++;
      } else if (ifound == istop) {
        // This hit does not fulfill requirement : restart to first
        // property but at different iterator position
        istart = ++ifound;
        istop = the_step_hits.end();
        iprop = _hit_property_values_.begin();
      } else {
        // Found one hit : check next property for this hit
        check_has_hit_property = true;
        if (_hit_property_logic_ == "or") break;
        istart = ifound;
        istop = ++ifound;
        iprop++;
      }
    }
  }

  cut_returned = cuts::SELECTION_REJECTED;
  if (check_flag && check_has_hit_category && check_range_hit_category && check_has_hit_property && check_gamma_has_gone_straight && check_has_hit_on_wall) {
    cut_returned = cuts::SELECTION_ACCEPTED;
  }

  return cut_returned;
}

}  // end of namespace cut

}  // end of namespace snemo

DOCD_CLASS_IMPLEMENT_LOAD_BEGIN(snemo::cut::simulated_data_cut, ocd_) {
  ocd_.set_class_name("snemo::cut::simulated_data_cut");
  ocd_.set_class_description(
      "Cut based on criteria applied to the simulated data bank stored in the event record");
  ocd_.set_class_library("falaise");
  // ocd_.set_class_documentation("");

  cuts::i_cut::common_ocd(ocd_);

  {
    // Description of the 'SD_label' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("SD_label")
        .set_terse_description("The name of the Simulated Data bank")
        .set_traits(datatools::TYPE_STRING)
        .set_default_value_string(snemo::datamodel::data_info::default_simulated_data_label())
        .add_example(
            "Set the default value::                          \n"
            "                                                 \n"
            "  SD_label : string = \"SD\"                     \n"
            "                                                 \n");
  }

  {
    // Description of the 'mode.flag' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("mode.flag")
        .set_terse_description("Mode with a special request flag")
        .set_traits(datatools::TYPE_BOOLEAN)
        .add_example(
            "Activate the requested flag mode::               \n"
            "                                                 \n"
            "  mode.flag : boolean = true                     \n"
            "                                                 \n");
  }

  {
    // Description of the 'mode.has_hit_category' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("mode.has_hit_category")
        .set_terse_description("Mode with a special requested hit category")
        .set_traits(datatools::TYPE_BOOLEAN)
        .add_example(
            "Activate the requested hit category mode::       \n"
            "                                                 \n"
            "  mode.has_hit_category : boolean = true         \n"
            "                                                 \n");
  }

  {
    // Description of the 'mode.range_hit_category' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("mode.range_hit_category")
        .set_terse_description("Mode with a special requested hit category to be ranged")
        .set_traits(datatools::TYPE_BOOLEAN)
        .add_example(
            "Activate the requested ranged hit category mode::\n"
            "                                                 \n"
            "  mode.range_hit_category : boolean = true       \n"
            "                                                 \n");
  }

  {
    // Description of the 'mode.has_hit_property' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("mode.has_hit_property")
        .set_terse_description("Mode with a special requested hit property")
        .set_traits(datatools::TYPE_BOOLEAN)
        .add_example(
            "Activate the requested hit property mode::       \n"
            "                                                 \n"
            "  mode.has_hit_property : boolean = true         \n"
            "                                                 \n");
  }

  {
    // Description of the 'flag.name' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("flag.name")
        .set_terse_description("Name of the requested flag")
        .set_triggered_by_flag("mode.flag")
        .set_traits(datatools::TYPE_STRING)
        .add_example(
            "Set a specific requested flag name::                      \n"
            "                                                          \n"
            "  flag.name : string = \"high_energy\"                    \n"
            "                                                          \n");
  }

  {
    // Description of the 'has_hit_category.category' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("has_hit_category.category")
        .set_terse_description("Name of the requested hit category")
        .set_triggered_by_flag("mode.has_hit_category")
        .set_traits(datatools::TYPE_STRING)
        .add_example(
            "Set a specific hit category::                             \n"
            "                                                          \n"
            "  has_hit_category.category : string = \"calo\"           \n"
            "                                                          \n");
  }

  {
    // Description of the 'range_hit_category.category' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("range_hit_category.category")
        .set_terse_description("Name of the requested ranged hit category")
        .set_triggered_by_flag("mode.range_hit_category")
        .set_traits(datatools::TYPE_STRING)
        .add_example(
            "Set a specific hit category to be ranged::                \n"
            "                                                          \n"
            "  range_hit_category.category : string = \"gg\"           \n"
            "                                                          \n");
  }

  {
    // Description of the 'range_hit_category.min' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("range_hit_category.min")
        .set_terse_description("Minimum number of hits of the requested ranged hit category")
        .set_triggered_by_flag("mode.range_hit_category")
        .set_traits(datatools::TYPE_INTEGER)
        .add_example(
            "Set a specific minimum number of hits::                   \n"
            "                                                          \n"
            "  range_hit_category.min : integer = 3                    \n"
            "                                                          \n");
  }

  {
    // Description of the 'range_hit_category.max' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("range_hit_category.max")
        .set_terse_description("Maximum number of hits of the requested ranged hit category")
        .set_triggered_by_flag("mode.range_hit_category")
        .set_traits(datatools::TYPE_INTEGER)
        .add_example(
            "Set a specific maximum number of hits::                   \n"
            "                                                          \n"
            "  range_hit_category.max : integer = 20                   \n"
            "                                                          \n");
  }

  {
    // Description of the 'has_hit_property.category' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("has_hit_property.category")
        .set_terse_description("Name of the hit category to search for a specific property")
        .set_triggered_by_flag("mode.has_hit_property")
        .set_traits(datatools::TYPE_STRING)
        .add_example(
            "Set a specific hit category to be investigated::          \n"
            "                                                          \n"
            "  has_hit_property.category : string = \"gg\"             \n"
            "                                                          \n");
  }

  {
    // Description of the 'has_hit_property.keys' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("has_hit_property.keys")
        .set_terse_description("Names of property keys to be tested in a specific hit category")
        .set_triggered_by_flag("mode.has_hit_property")
        .set_traits(datatools::TYPE_STRING, datatools::configuration_property_description::ARRAY)
        .add_example(
            "Set the lists of keys to be tested::                      \n"
            "                                                          \n"
            "  has_hit_property.keys : string[2] = \\                  \n"
            "     \"noisy\" \\                                         \n"
            "     \"bad\"                                              \n"
            "                                                          \n");
  }

  {
    // Description of the 'has_hit_property.${has_hit_property.keys}.values' configuration property
    // :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("has_hit_property.${has_hit_property.keys}.values")
        .set_traits(datatools::TYPE_STRING, datatools::configuration_property_description::ARRAY)
        .add_example(
            "Set the lists of values to be associated to tested keys::   \n"
            "                                                            \n"
            "  has_hit_property.keys : string[2] = \\                    \n"
            "     \"creator_process\" \\                                 \n"
            "     \"g4_volume\"                                          \n"
            "  has_hit_property.creator_process.values : string[2] = \\  \n"
            "     \"brems\" \\                                           \n"
            "     \"compton\"                                            \n"
            "  has_hit_property.g4_volume.values : string[2] = \\        \n"
            "     \"drift_cell_core.log\" \\                             \n"
            "     \"tracker_chamber.log\"                                \n"
            "                                                            \n");
  }

  // Additional configuration hints :
  ocd_.set_configuration_hints(
      "Here is a full configuration example in the     \n"
      "``datatools::properties`` ASCII format::        \n"
      "                                                \n"
      "   SD_label : string = \"SD\"                   \n"
      "   mode.flag : boolean = false                  \n"
      "   # flag.name : string = \"high_energy\"       \n"
      "   mode.has_hit_category : boolean = false      \n"
      "   # has_hit_category.category : string = \"gg\"\n"
      "   mode.range_hit_category : boolean = true     \n"
      "   range_hit_category.category : string = \"gg\"\n"
      "   range_hit_category.min : integer = 5         \n"
      "   range_hit_category.max : integer = 20        \n"
      "   mode.has_hit_property : boolean = true       \n"
      "   has_hit_property.category : string = \"gg\"  \n"
      "   has_hit_property.logic : string = \"or\"     \n"
      "   has_hit_property.keys : string[2] = \\                   \n"
      "     \"creator_process\" \\                                 \n"
      "     \"g4_volume\"                                          \n"
      "   has_hit_property.creator_process.values : string[1] = \\ \n"
      "         \"brems\"                                          \n"
      "   has_hit_property.g4_volume.values : string[2] = \\       \n"
      "         \"tracking_chamber.log\"  \"drift_cell.log\"       \n"
      "                                                            \n");

  ocd_.set_validation_support(true);
  ocd_.lock();
  return;
}
DOCD_CLASS_IMPLEMENT_LOAD_END()  // Closing macro for implementation

// Registration macro for class 'snemo::cut::simulated_data_cut' :
DOCD_CLASS_SYSTEM_REGISTRATION(snemo::cut::simulated_data_cut, "snemo::cut::simulated_data_cut")

/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** End: --
*/
