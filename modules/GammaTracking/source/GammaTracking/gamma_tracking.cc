// Ourselves:
#include <GammaTracking/gamma_tracking.h>
// This project
#include <GammaTracking/tof_computing.h>

// Standard library:
#include <algorithm>
#include <iostream>

// Third party:
// - GSL:
#include <gsl/gsl_cdf.h>
// - Boost:
#include <boost/next_prior.hpp>
// #include <boost/fusion/iterator/next.hpp>
// - Bayeux/datatools:
#include <datatools/properties.h>

namespace gt {

  gamma_tracking::gamma_tracking() {
    _set_defaults();
    _initialized_ = false;
    return;
  }

  gamma_tracking::gamma_tracking(const gamma_tracking &gt_) {
    _initialized_ = gt_._initialized_;
    _initialized_combi_ = gt_._initialized_combi_;
    _logging_priority_ = gt_._logging_priority_;
    _absolute_ = gt_._absolute_;
    _extern_ = gt_._extern_;
    _max_ = gt_._max_;
    _min_prob_ = gt_._min_prob_;
    _starts_ = gt_._starts_;
    _pair_calo_ = gt_._pair_calo_;
    _single_calo_ = gt_._single_calo_;
    _combination_ = gt_._combination_;
    _all_combi_ = gt_._all_combi_;
    _min_chi2_ = gt_._min_chi2_;

    for (std::map<const list_type *, double>::const_iterator mit = gt_._chi2_.begin();
	 mit != gt_._chi2_.end(); ++mit) {
      solution_type::iterator it = std::find(_pair_calo_.begin(), _pair_calo_.end(), *(mit->first));
      if (it != _pair_calo_.end()) _chi2_[&(*it)] = mit->second;
    }

    for (std::map<const list_type *, double>::const_iterator mit = gt_._proba_.begin();
	 mit != gt_._proba_.end(); ++mit) {
      solution_type::iterator it = std::find(_pair_calo_.begin(), _pair_calo_.end(), *(mit->first));
      if (it != _pair_calo_.end()) _proba_[&(*it)] = mit->second;
    }
    return;
  }

  gamma_tracking::~gamma_tracking() {
    if (is_initialized()) reset();
    return;
  }

  void gamma_tracking::set_logging_priority(datatools::logger::priority priority_) {
    _logging_priority_ = priority_;
    return;
  }

  datatools::logger::priority gamma_tracking::get_logging_priority() const {
    return _logging_priority_;
  }

  bool gamma_tracking::is_initialized() const { return _initialized_; }

  void gamma_tracking::set_initialized(bool initialized_) {
    _initialized_ = initialized_;
    return;
  }
  void gamma_tracking::set_initialized_combi(bool initialized_) {
    _initialized_combi_ = initialized_;
    return;
  }
  void gamma_tracking::initialize(const datatools::properties &config_) {
    DT_THROW_IF(is_initialized(), std::logic_error, "Already initialized !");

    datatools::logger::priority p = datatools::logger::extract_logging_configuration(config_, datatools::logger::PRIO_UNDEFINED, true);

    if (p != datatools::logger::PRIO_UNDEFINED) {
      set_logging_priority(p);
    }

    if (config_.has_flag("use_absolute")) {
      set_absolute(true);
    }

    if (config_.has_flag("use_extern")) {
      set_extern(true);
    }
    if (config_.has_flag("print_all")) {
      set_print(true);
    }
    if (config_.has_flag("print_size")) {
      set_print_size(true);
    }

    // if (config_.has_key("maximal_gamma_size")) {
    //   _max_ = config_.fetch_integer("maximal_gamma_size");
    // }

    if (config_.has_key("minimal_probability")) {
      _min_prob_ = config_.fetch_real("minimal_probability");
    }
    set_initialized_combi(false);
    set_initialized(true);
    return;
  }

  void gamma_tracking::_set_defaults() {
    _logging_priority_ = datatools::logger::PRIO_WARNING;
    _max_ = 0;
    _min_prob_ = 1e-5;
    _min_chi2_.insert(std::make_pair(1, gsl_cdf_chisq_Qinv(_min_prob_, 1)));
    _absolute_ = false;
    _extern_ = false;
    _print_ = false;
    _print_size_ = false;
    return;
  }

  bool gamma_tracking::has_tracks() { return _pair_calo_.size() > 0; }

  void gamma_tracking::add(int number_) {
    list_type tamp;
    tamp.push_back(number_);

    if (is_inside_serie(tamp)) return;
    _pair_calo_.push_back(tamp);
    _proba_[&(_pair_calo_.back())] = 1.0;
    _chi2_[&(_pair_calo_.back())] = 0.0;
    ++_max_;
  }

  void gamma_tracking::add(int number_,solution_type& _selected_gamma_) {
    list_type tamp;
    tamp.push_back(number_);

    if (is_inside_serie(tamp)) return;
    _selected_gamma_.push_back(tamp);
    _proba_[&(_selected_gamma_.back())] = 1.0;
    _chi2_[&(_selected_gamma_.back())] = 0.0;
    ++_max_;
  }
  void gamma_tracking::add_probability(int number1_, int number2_, double proba_) {
   
    if (number1_ == number2_) return;

    DT_LOG_TRACE(get_logging_priority(),
		 "Probability(" << number1_ << "," << number2_ << ") = " << proba_);
    if (proba_ < _min_prob_) {
      DT_LOG_TRACE(get_logging_priority(),
		   "Probability below threshold (" << proba_ << "<" << _min_prob_ << ")");
      return;
    }
    

    list_type tamp;
    tamp.push_back(number1_);
    tamp.push_back(number2_);
    if (is_inside_serie(tamp)) return;
    _pair_calo_.push_back(tamp);

    const double chi2 = gsl_cdf_chisq_Qinv(proba_, 1);
    _chi2_[&(_pair_calo_.back())] = chi2;
    _proba_[&(_pair_calo_.back())] = proba_;
  }

  void gamma_tracking::add_chi2(int number1_, int number2_, double chi2_) {
    if (number1_ == number2_ || chi2_ > get_chi_limit(1)) {
      DT_LOG_TRACE(get_logging_priority(),
		   "X² value below minimal value (" << chi2_ << "<" << get_chi_limit(1));
      return;
    }
    const double proba = gsl_cdf_chisq_Q(chi2_, 1);
    add_probability(number1_, number2_, proba);
  }

  void gamma_tracking::add_start(int number_) {
    bool its_inside = false;
    for (solution_type::iterator rit = _pair_calo_.begin(); rit != _pair_calo_.end() && !its_inside; ++rit) {
      its_inside = is_inside((*rit), number_);
    }

    if (!its_inside) {
      DT_LOG_WARNING(
		     get_logging_priority(),
		     "Calorimeter " << number_ << " is not in the collection. You should add it before.");
      return;
    }
    _starts_.push_back(number_);
    return;
  }

  void gamma_tracking::dump(std::ostream &out_) const {
    for (solution_type::const_iterator it = _pair_calo_.begin(); it != _pair_calo_.end(); ++it) {
      if (boost::next(it) == _pair_calo_.end()) {
	out_ << "`- ";
      } else {
	out_ << "|- ";
      }
      out_ << "for list ";
      for (list_type::const_iterator iit = it->begin(); iit != it->end(); ++iit) {
	out_ << *iit;
	if (boost::next(iit) != it->end()) {
	  out_ << "->";
	}
      }
      out_ << " - probability is " << _proba_.at(&(*it)) << std::endl;
    }
    return;
  }

  bool gamma_tracking::is_inside_serie(const list_type &list_) const {
    if (std::find(_pair_calo_.begin(), _pair_calo_.end(), list_) != _pair_calo_.end()) {
      return true;
    }
    return false;
  }
  bool gamma_tracking::is_inside(const list_type &check_, int value_) const {
    if (std::find(check_.begin(), check_.end(), value_) != check_.end()) {
      return true;
    }
    return false;
  }
  bool gamma_tracking::is_inside(const list_type &check_, const list_type &values_) const {
    for (list_type::const_iterator it = values_.begin(); it != values_.end(); ++it) {
      if (std::find(check_.begin(), check_.end(), (*it)) != check_.end()) {
	return true;
      }
    }
    return false;
  }

  void gamma_tracking::extract(list_type &source_, const list_type &values_) {
    for (list_type::const_iterator it = values_.begin(); it != values_.end(); ++it) {
      list_type::iterator to_extract = std::find(source_.begin(), source_.end(), *it);
      if (to_extract != source_.end()) source_.erase(to_extract);
    }
    return;
  }

  void gamma_tracking::put_inside(const list_type &from_, list_type &to_) {
    for (list_type::const_iterator it = from_.begin(); it != from_.end(); ++it) to_.push_back(*it);
    to_.sort();
    to_.unique();
    return;
  }

  void gamma_tracking::set_absolute(bool absolute_) {
    _absolute_ = absolute_;
    return;
  }
  void gamma_tracking::set_print(bool print_) {
    _print_ = print_;
    return;
  }
  void gamma_tracking::set_print_size(bool print_) {
    _print_size_ = print_;
    return;
  }
  bool gamma_tracking::is_to_print() { return _print_; }
  bool gamma_tracking::is_to_print_size() { return _print_size_;}
  bool gamma_tracking::is_absolute() { return _absolute_; }

  void gamma_tracking::set_extern(bool extern_) {
    _extern_ = extern_;
    return;
  }

  bool gamma_tracking::is_extern() { return _extern_; }

  void gamma_tracking::set_probability_min(double min_prob_) {
    _min_prob_ = min_prob_;
    for (std::map<int, double>::iterator it = _min_chi2_.begin(); it != _min_chi2_.end(); ++it) {
      it->second = gsl_cdf_chisq_Qinv(_min_prob_, it->first);
    }
  }

  const gamma_tracking::solution_type &gamma_tracking::get_all() const { return _all_combi_; }

  double gamma_tracking::get_probability(const list_type &scin_ids_) const {
    double probability = datatools::invalid_real();
    solution_type::const_iterator it = std::find(_pair_calo_.begin(), _pair_calo_.end(), scin_ids_);
    if (it != _pair_calo_.end()) {
      probability = _proba_.at(&(*it));
    }
    return probability;
  }
  double gamma_tracking::get_probability(int scin_id1_, int scin_id2_) const {
    list_type l1;
    l1.push_back(scin_id1_);
    l1.push_back(scin_id2_);
    return get_probability(l1);
  }

  double gamma_tracking::get_chi2(const list_type &scin_ids_) const {
    double chi2 = datatools::invalid_real();
    solution_type::const_iterator it = std::find(_pair_calo_.begin(), _pair_calo_.end(), scin_ids_);
    if (it != _pair_calo_.end()) {
      chi2 = _chi2_.at(&(*it));
    }
    return chi2;
  }
  double gamma_tracking::get_chi2(int scin_id1_, int scin_id2_) const {
    list_type l1;
    l1.push_back(scin_id1_);
    l1.push_back(scin_id2_);
    return get_chi2(l1);
  }
  double gamma_tracking::get_chi_limit(unsigned int freedom_) {
    if (!_min_chi2_.count(freedom_)) {
      _min_chi2_[freedom_] = gsl_cdf_chisq_Qinv(_min_prob_, freedom_);
    }
    return _min_chi2_[freedom_];
  }

  const event &gamma_tracking::get_event() const { return _event_; }

  event &gamma_tracking::grab_event() { return _event_; }

  gamma_tracking::list_type gamma_tracking::fuse_calo_list(const list_type list1, const list_type list2){

    list_type fused_list;

    for (auto &it: list1) {

      if (it == list2.front()) continue;
      fused_list.push_back(it);
    }
    
    for (auto &it: list2)
      fused_list.push_back(it);

    
    return fused_list;
  }

  gamma_tracking::list_type gamma_tracking::convert_int_to_list(int value){
    list_type list;
    list.push_back(value);
    return list;
  }
    
  void gamma_tracking::prepare_process() {

    
    const event::calorimeter_collection_type &the_gamma_calos = _event_.get_calorimeters();

    for (event::calorimeter_collection_type::const_iterator icalo = the_gamma_calos.begin();
	 icalo != the_gamma_calos.end(); ++icalo) {

      _single_calo_.push_back(icalo->first);
      
      for (event::calorimeter_collection_type::const_iterator jcalo = the_gamma_calos.begin();
	   jcalo != the_gamma_calos.end(); ++jcalo) {
	
	if (icalo == jcalo) continue;
	
	event::calorimeter_collection_type::const_iterator it1 = icalo;
	event::calorimeter_collection_type::const_iterator it2 = jcalo;
	  
	const double tof_chi2 = tof_computing::get_chi2(it1->second, it2->second);
	const double tof_prob = tof_computing::get_internal_probability(tof_chi2);
	
	DT_LOG_DEBUG(get_logging_priority(), "X²(" << it1->first << "->" << it2->first
		     << ") = " << tof_chi2 << ", P = " << tof_prob);
	add_probability(it1->first, it2->first, tof_prob);
      }
    }
    return;
  }
  
  void gamma_tracking::initialize_combi(){
	  
    list_type fused_list;
    bool check_is_already;
    
      for(auto &Pair1: _pair_calo_) {

	for(auto &Pair2: _pair_calo_) {

	  if (Pair1.back() == Pair2.front()){

	    check_is_already = false;

	    for (auto &check : Pair1){

	      if (check == Pair2.back()) check_is_already = true;
	    }
	    
	    if (check_is_already) continue;

	    
	    const int ndf = Pair1.size() -1  + Pair2.size() -1;
	    const double chi2_combi = _chi2_[&(Pair1)] + _chi2_[&(Pair2)];
	    
	    DT_LOG_TRACE(get_logging_priority(), "X²[" << &(Pair1) << "] = " << _chi2_[&(Pair1)]);
	    DT_LOG_TRACE(get_logging_priority(), "X²[" << &(Pair2) << "] = " << _chi2_[&(Pair2)]);
	    DT_LOG_TRACE(get_logging_priority(), "X² = " << chi2_combi);
	    
	    if (chi2_combi < get_chi_limit(ndf)) {
	      
	      const double prob_combi = tof_computing::get_internal_probability(chi2_combi,ndf);
	      fused_list = fuse_calo_list(Pair1,Pair2);
	      
	      _combination_.push_back(fused_list);
	      _chi2_[&(_combination_.back())] = chi2_combi;
	      _proba_[&(_combination_.back())] = prob_combi;

	    }
	  }
	}	  
      }
    
    
    gamma_tracking::set_initialized_combi(true);
  }
  void gamma_tracking::fill_combi(){

	  
    list_type fused_list;
    bool check_is_already;
    
    DT_LOG_TRACE(get_logging_priority(), "Combination Size " << _combination_.size());

    if (_initialized_combi_ == false) {
      DT_LOG_TRACE(get_logging_priority(), "Combination not initialized");
      
      return;
      }

      
    for (auto &Combi: _combination_) {

      for (auto &Pair: _pair_calo_) {

	if (Combi.back() == Pair.front()){

	  check_is_already = false;
	
	  for (auto &check : Combi){

	    if (check == Pair.back()) check_is_already = true;
	  }

	  if (check_is_already) continue;
	  
	  const int ndf = Combi.size() -1 + Pair.size() - 1;
	  const double chi2_combi = _chi2_[&(Combi)] + _chi2_[&(Pair)];
	  
	  DT_LOG_TRACE(get_logging_priority(), "X²[" << &(Combi) << "] = " << _chi2_[&(Combi)]);
	  DT_LOG_TRACE(get_logging_priority(), "X²[" << &(Pair) << "] = " << _chi2_[&(Pair)]);
	  DT_LOG_TRACE(get_logging_priority(), "X² = " << chi2_combi);

	  if (chi2_combi < get_chi_limit(ndf)){
	    const double prob_combi = tof_computing::get_internal_probability(chi2_combi,ndf);
	    fused_list = fuse_calo_list(Combi,Pair);

	    _combination_.push_back(fused_list);
	    _chi2_[&(_combination_.back())] = chi2_combi;
	    _proba_[&(_combination_.back())] = prob_combi;
	  }
	}
      }

    }
  }

  void gamma_tracking::sort_list(solution_type &list_) {

    if (list_.size() <= 1) return;

    bool _combi_has_changed = true;
    while (_combi_has_changed) {
      _combi_has_changed = false;
      solution_type::iterator it1 = list_.begin();
      solution_type::iterator it2 = list_.begin();
      it2++;
      
      while (it2 != list_.end() && !_combi_has_changed) {

	if ( it1->size() < it2->size()) {

	  _combi_has_changed = true;
	  list_.splice(it1, list_, it2);
	  
	}
	else if (it1->size() == it2->size() && _proba_[&(*it1)] < _proba_[&(*it2)] ){

	  _combi_has_changed = true;
	  list_.splice(it1, list_, it2);
	}
	else {
	  it1++;
	  it2++;
	}
      }
    }
    return;

 
  }
  void gamma_tracking::sort_all_combi(){

    gamma_tracking::sort_list(_pair_calo_);
    gamma_tracking::sort_list(_combination_);
    
    return;
  }
  
  void gamma_tracking::gamma_selection(solution_type &selected_gamma_){

    bool has_all_empty=false;
    if (_combination_.empty() && _pair_calo_.empty()) has_all_empty=true;
    
    const double UProb = _min_prob_;
    
    selected_gamma_.clear();

    if (!has_all_empty){
      solution_type::iterator it_list_pair = _pair_calo_.begin();
    
      while (it_list_pair != _pair_calo_.end()){
      
	if(_proba_[&(*it_list_pair)] < UProb)
	  it_list_pair = _pair_calo_.erase(it_list_pair);
	else it_list_pair++;
      
      }

      solution_type::iterator it_list_combi = _combination_.begin();
    
      while (it_list_combi != _combination_.end()){
	
	if(_proba_[&(*it_list_combi)] < UProb)
	  it_list_combi = _combination_.erase(it_list_combi);
	else it_list_combi++;
      
      }    
      
      while (!_combination_.empty() ) {

	selected_gamma_.push_back(_combination_.front());
     
	for (auto &number : selected_gamma_.back()) {

	  if (_combination_.empty() && _pair_calo_.empty()) break;

	  solution_type::iterator it_list = _combination_.begin();

	  while (it_list != _combination_.end()){

	    if(std::find(it_list->begin(), it_list->end(), number) != it_list->end()){
	      it_list = _combination_.erase(it_list);
	    
	    }
	    else it_list++;

	  }

	  it_list = _pair_calo_.begin();
	  while (it_list != _pair_calo_.end()){

	    if(std::find(it_list->begin(), it_list->end(), number) != it_list->end()){
	      it_list = _pair_calo_.erase(it_list);
	    
	    }
	    else it_list++;

	
	  }
	}
      }
    
      while (!_pair_calo_.empty() ) {

	selected_gamma_.push_back(_pair_calo_.front());
     
	for (auto &number : selected_gamma_.back()) {
	  
	  if (_pair_calo_.empty()) break;
	
	  solution_type::iterator it_list = _pair_calo_.begin();
	  
	  while (it_list != _pair_calo_.end()){

	    if(std::find(it_list->begin(), it_list->end(), number) != it_list->end()){
	      it_list = _pair_calo_.erase(it_list);
	    
	    }
	    else it_list++;

	  }

	}
      }  
    }
    
    for (auto &it_selec : selected_gamma_) {
      for (auto &Num : it_selec) {

	list_type::iterator it_single = std::find(_single_calo_.begin(), _single_calo_.end(),Num);
	
	if(it_single != _single_calo_.end()) _single_calo_.erase(it_single); 

	if(_single_calo_.empty()) break;
      }
      if(_single_calo_.empty()) break;
    }
    
    if(!_single_calo_.empty()) {
      for( auto &Single : _single_calo_ ) {
	selected_gamma_.push_back(gamma_tracking::convert_int_to_list(Single));
      }
    }
  }
  
  
  void gamma_tracking::fuse_all_combi(){

    for (auto &it: _combination_)
      _all_combi_.push_back(it);

    for (auto &it: _pair_calo_)
      _all_combi_.push_back(it);

    for (auto &it: _single_calo_)
      _all_combi_.push_back(gamma_tracking::convert_int_to_list(it));
    
    return;

  }
  void gamma_tracking::print_all_proba(){

    for (auto &it : _combination_) {
      std::cout <<"C"<< _proba_[&(it)] << std::endl;
      std::cout <<"D"<< it.size() << std::endl;
    }
    for (auto &it : _pair_calo_){
      std::cout << "P" << _proba_[&(it)] << std::endl;
      std::cout <<"D"<< it.size() << std::endl;
    }
    
  }
  void gamma_tracking::print_size_combi(){
    int R3=0,R4=0,R5=0,R6=0,R7=0,R8=0,R9=0,R10=0;

    for (auto &it: _combination_){

      if(it.size() == 3) R3++;
      if(it.size() == 4) R4++;
      if(it.size() == 5) R5++;
      if(it.size() == 6) R6++;
      if(it.size() == 7) R7++;
      if(it.size() == 8) R8++;
      if(it.size() == 9) R9++;
      if(it.size() == 10) R10++;
      
    }
    std::cout<<_pair_calo_.size()<<" "<<R3<<" "<<R4<<" "<<R5<<" "<<R6<<" "<<R7<<" "<<R8<<" "<<R9<<" "<<R10<<std::endl;
  }

  void gamma_tracking::process(solution_type &selected_gamma_) {

    gamma_tracking::initialize_combi();
    gamma_tracking::fill_combi();
    gamma_tracking::sort_all_combi();
    if (is_to_print()) gamma_tracking::print_all_proba();
    if (is_to_print_size()) gamma_tracking::print_size_combi();
    gamma_tracking::fuse_all_combi(); //temporary
    
    
    gamma_tracking::gamma_selection(selected_gamma_);    
      
  }

  void gamma_tracking::reset() {
    _single_calo_.clear();
    _pair_calo_.clear();
    _combination_.clear();
    _all_combi_.clear();
    _chi2_.clear();
    _proba_.clear();
    _starts_.clear();
    _event_.reset();
    return;
  }

}  // namespace gt
