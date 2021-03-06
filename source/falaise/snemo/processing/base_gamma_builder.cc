/// \file falaise/snemo/processing/base_gamma_builder.cc

// Ourselves:
#include "falaise/snemo/processing/base_gamma_builder.h"

// Third party:
// - GSL:
#include <gsl/gsl_cdf.h>
// - Bayeux/datatools:
#include <bayeux/datatools/properties.h>
// - Bayeux/geomtools:
#include <bayeux/geomtools/manager.h>

// This project:
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/xcalo_locator.h>

namespace snemo {

namespace processing {

datatools::logger::priority base_gamma_builder::get_logging_priority() const {
  return _logging_priority;
}

void base_gamma_builder::set_logging_priority(datatools::logger::priority priority_) {
  _logging_priority = priority_;
  return;
}

const std::string& base_gamma_builder::get_id() const { return _id_; }

const snemo::geometry::calo_locator& base_gamma_builder::get_calo_locator() const {
  DT_THROW_IF(!is_initialized(), std::logic_error,
              "Driver '" << get_id() << "' is not initialized !");
  return _locator_plugin_->get_calo_locator();
}

const snemo::geometry::xcalo_locator& base_gamma_builder::get_xcalo_locator() const {
  DT_THROW_IF(!is_initialized(), std::logic_error,
              "Driver '" << get_id() << "' is not initialized !");
  return _locator_plugin_->get_xcalo_locator();
}

const snemo::geometry::gveto_locator& base_gamma_builder::get_gveto_locator() const {
  DT_THROW_IF(!is_initialized(), std::logic_error,
              "Driver '" << get_id() << "' is not initialized !");
  return _locator_plugin_->get_gveto_locator();
}

bool base_gamma_builder::is_initialized() const { return _initialized_; }

void base_gamma_builder::_set_initialized(bool i_) {
  _initialized_ = i_;
  return;
}

void base_gamma_builder::_initialize(const datatools::properties& setup_) {
  DT_THROW_IF(is_initialized(), std::logic_error,
              "Driver '" << get_id() << "' is already initialized !");

  DT_THROW_IF(!has_geometry_manager(), std::logic_error, "Missing geometry manager !");
  DT_THROW_IF(!get_geometry_manager().is_initialized(), std::logic_error,
              "Geometry manager is not initialized !");

  // Extract the setup of the base gamma builder :
  datatools::properties bgb_setup;
  setup_.export_and_rename_starting_with(bgb_setup, "BGB.", "");

  // Logging priority
  datatools::logger::priority lp = datatools::logger::extract_logging_configuration(bgb_setup);
  DT_THROW_IF(lp == datatools::logger::PRIO_UNDEFINED, std::logic_error,
              "Invalid logging priority level for base gamma builder !");
  set_logging_priority(lp);

  // Get geometry locator plugin
  const geomtools::manager& geo_mgr = get_geometry_manager();
  std::string locator_plugin_name;
  if (bgb_setup.has_key("locator_plugin_name")) {
    locator_plugin_name = bgb_setup.fetch_string("locator_plugin_name");
  } else {
    // If no locator plugin name is set, then search for the first one
    const geomtools::manager::plugins_dict_type& plugins = geo_mgr.get_plugins();
    for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
         ip != plugins.end(); ip++) {
      const std::string& plugin_name = ip->first;
      if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(plugin_name)) {
        DT_LOG_DEBUG(get_logging_priority(), "Find locator plugin with name = " << plugin_name);
        locator_plugin_name = plugin_name;
        break;
      }
    }
  }
  // Access to a given plugin by name and type :
  DT_THROW_IF(!geo_mgr.has_plugin(locator_plugin_name) ||
                  !geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name),
              std::logic_error, "Found no locator plugin named '" << locator_plugin_name << "'");
  _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);

  // Select calorimeter hits based on associated tags
  if (bgb_setup.has_key("select_calorimeter_hits")) {
    _select_calorimeter_hits_ = bgb_setup.fetch_boolean("select_calorimeter_hits");
  }

  if (_select_calorimeter_hits_) {
    DT_THROW_IF(!bgb_setup.has_key("select_calorimeter_hits.tags"), std::logic_error,
                "Missing 'select_calorimeter_hits.tags' field!");
    bgb_setup.fetch("select_calorimeter_hits.tags", _select_calorimeter_hits_tags_);
  }

  // Extrapolation on the source foil given charged particle
  if (bgb_setup.has_key("add_foil_vertex_extrapolation")) {
    _add_foil_vertex_extrapolation_ = bgb_setup.fetch_boolean("add_foil_vertex_extrapolation");
  }

  if (_add_foil_vertex_extrapolation_) {
    if (bgb_setup.has_key("add_foil_vertex_extrapolation.minimal_probability")) {
      _add_foil_vertex_minimal_probability_ =
          bgb_setup.fetch_real("add_foil_vertex_extrapolation.minimal_probability");
      if (!bgb_setup.has_explicit_unit("add_foil_vertex_extrapolation.minimal_probability")) {
        _add_foil_vertex_minimal_probability_ *= CLHEP::perCent;
      }
    }
  }

  // Search for gamma from e+/e- annihilation
  if (bgb_setup.has_key("add_gamma_from_annihilation")) {
    _add_gamma_from_annihilation_ = bgb_setup.fetch_boolean("add_gamma_from_annihilation");
  }

  if (_add_gamma_from_annihilation_) {
    if (bgb_setup.has_key("add_gamma_from_annihilation.minimal_probability")) {
      _add_gamma_from_annihilation_minimal_probability_ =
          bgb_setup.fetch_real("add_gamma_from_annihilation.minimal_probability");
      if (!bgb_setup.has_explicit_unit("add_gamma_from_annihilation.minimal_probability")) {
        _add_gamma_from_annihilation_minimal_probability_ *= CLHEP::perCent;
      }
    }
  }
  return;
}

void base_gamma_builder::_clear_working_arrays() {
  _ignored_hits_.clear();
  _used_hits_.clear();
  return;
}

void base_gamma_builder::_reset() {
  _set_initialized(false);
  this->base_gamma_builder::_set_defaults();
  this->base_gamma_builder::_clear_working_arrays();
  return;
}

void base_gamma_builder::set_geometry_manager(const geomtools::manager& gmgr_) {
  DT_THROW_IF(is_initialized(), std::logic_error, "Already initialized/locked !");
  _geometry_manager_ = &gmgr_;
  return;
}

const geomtools::manager& base_gamma_builder::get_geometry_manager() const {
  DT_THROW_IF(!has_geometry_manager(), std::logic_error, "No geometry manager is setup !");
  return *_geometry_manager_;
}

bool base_gamma_builder::has_geometry_manager() const { return _geometry_manager_ != 0; }

void base_gamma_builder::_set_defaults() {
  _logging_priority = datatools::logger::PRIO_WARNING;
  _geometry_manager_ = 0;
  _locator_plugin_ = 0;
  _add_foil_vertex_extrapolation_ = true;
  _add_foil_vertex_minimal_probability_ = 1.0 * CLHEP::perCent;
  _add_gamma_from_annihilation_ = false;
  _add_gamma_from_annihilation_minimal_probability_ = 1.0 * CLHEP::perCent;
  _select_calorimeter_hits_ = false;
  _select_calorimeter_hits_tags_.clear();
  return;
}

// Constructor
base_gamma_builder::base_gamma_builder(const std::string& id_) {
  _id_ = id_;
  _set_initialized(false);
  _set_defaults();
  return;
}

base_gamma_builder::~base_gamma_builder() {
  if (is_initialized()) {
    _reset();
  }
  return;
}

void base_gamma_builder::tree_dump(std::ostream& out_, const std::string& title_,
                                   const std::string& indent_, bool inherit_) const {
  std::string indent;
  if (!indent_.empty()) {
    indent = indent_;
  }
  if (!title_.empty()) {
    out_ << indent << title_ << std::endl;
  }

  out_ << indent << datatools::i_tree_dumpable::tag << "Logging          : '"
       << datatools::logger::get_priority_label(_logging_priority) << "'" << std::endl;
  out_ << indent << datatools::i_tree_dumpable::tag << "Initialized      : " << is_initialized()
       << std::endl;
  out_ << indent << datatools::i_tree_dumpable::tag << "Geometry manager : " << _geometry_manager_
       << std::endl;
  if (has_geometry_manager()) {
    out_ << indent << datatools::i_tree_dumpable::tag << "Geometry setup label   : '"
         << _geometry_manager_->get_setup_label() << "'" << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Geometry setup version : '"
         << _geometry_manager_->get_setup_version() << "'" << std::endl;
  }

  out_ << indent << datatools::i_tree_dumpable::tag
       << "Foil vertex extrapolation : " << _add_foil_vertex_extrapolation_ << std::endl;
  if (_add_foil_vertex_extrapolation_) {
    out_ << indent << datatools::i_tree_dumpable::skip_tag << datatools::i_tree_dumpable::last_tag
         << "Minimal TOF probability : " << _add_foil_vertex_minimal_probability_ / CLHEP::perCent
         << "%" << std::endl;
  }
  out_ << indent << datatools::i_tree_dumpable::tag
       << "Search for gamma from e+/e- annihilation : " << _add_gamma_from_annihilation_
       << std::endl;
  if (_add_gamma_from_annihilation_) {
    out_ << indent << datatools::i_tree_dumpable::skip_tag << datatools::i_tree_dumpable::last_tag
         << "Minimal TOF probability : "
         << _add_gamma_from_annihilation_minimal_probability_ / CLHEP::perCent << "%" << std::endl;
  }
  out_ << indent << datatools::i_tree_dumpable::tag
       << "Selection of calorimeter hits : " << _select_calorimeter_hits_ << std::endl;
  if (_select_calorimeter_hits_) {
    for (size_t i = 0; i < _select_calorimeter_hits_tags_.size(); i++) {
      out_ << indent << datatools::i_tree_dumpable::skip_tag;
      if (i + 1 == _select_calorimeter_hits_tags_.size()) {
        out_ << datatools::i_tree_dumpable::last_tag;
      } else {
        out_ << datatools::i_tree_dumpable::tag;
      }
      out_ << "tag[" << i << "] = " << _select_calorimeter_hits_tags_[i] << std::endl;
    }
  }
  out_ << indent << datatools::i_tree_dumpable::inherit_tag(inherit_) << "End." << std::endl;
  return;
}

int base_gamma_builder::process(const base_gamma_builder::hit_collection_type& calo_hits_,
                                snemo::datamodel::particle_track_data& ptd_) {
  int status = 0;
  DT_THROW_IF(!is_initialized(), std::logic_error,
              "Gamma builder '" << _id_ << "' is not initialized !");

  status = _prepare_process(calo_hits_, ptd_);
  if (status != 0) {
    DT_LOG_ERROR(get_logging_priority(), "Pre-processing '" << get_id() << "' has failed !");
    return status;
  }

  status = _process_algo(_used_hits_, ptd_);
  if (status != 0) {
    DT_LOG_ERROR(get_logging_priority(),
                 "Processing of '" << get_id() << "' algorithm has failed !");
    return status;
  }

  status = _post_process(calo_hits_, ptd_);
  if (status != 0) {
    DT_LOG_ERROR(get_logging_priority(),
                 "Post-processing of '" << get_id() << "' algorithm has failed !");
    return status;
  }
  return status;
}

int base_gamma_builder::_prepare_process(const base_gamma_builder::hit_collection_type& calo_hits_,
                                         snemo::datamodel::particle_track_data& /*ptd_*/) {
  this->base_gamma_builder::_clear_working_arrays();
  _used_hits_.reserve(calo_hits_.size());
  _ignored_hits_.reserve(calo_hits_.size());
  for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator ihit =
           calo_hits_.begin();
       ihit != calo_hits_.end(); ++ihit) {
    const snemo::datamodel::calibrated_calorimeter_hit& a_calo_hit = ihit->get();
    bool use_hit = false;
    if (_select_calorimeter_hits_) {
      const datatools::properties& the_auxiliaries = a_calo_hit.get_auxiliaries();
      for (size_t i = 0; i < _select_calorimeter_hits_tags_.size(); i++) {
        if (the_auxiliaries.has_flag(_select_calorimeter_hits_tags_[i])) {
          use_hit = true;
          break;
        }
      }
    } else {
      use_hit = true;
    }
    if (use_hit) {
      _used_hits_.push_back(*ihit);
    } else {
      _ignored_hits_.push_back(*ihit);
    }
  }
  DT_LOG_DEBUG(get_logging_priority(),
               "Number of calorimeter hits used: " << _used_hits_.size() << " ("
                                                   << _ignored_hits_.size() << " skipped)");
  return 0;
}

int base_gamma_builder::_post_process(const base_gamma_builder::hit_collection_type& /*calo_hits_*/,
                                      snemo::datamodel::particle_track_data& ptd_) {
  // Add the ignored hits to the list of non associated calorimeters
  // ptd_.reset_non_associated_calorimeters();
  snemo::datamodel::calibrated_calorimeter_hit::collection_type& calos =
      ptd_.grab_non_associated_calorimeters();
  calos.assign(_ignored_hits_.begin(), _ignored_hits_.end());

  // Given charged particle then process gammas
  snemo::datamodel::particle_track_data::particle_collection_type gamma_particles;
  ptd_.fetch_particles(gamma_particles, snemo::datamodel::particle_track::NEUTRAL);
  if (gamma_particles.empty()) {
    DT_LOG_DEBUG(get_logging_priority(), "No gamma particles have been found !");
    return 0;
  }

  snemo::datamodel::particle_track_data::particle_collection_type charged_particles;
  ptd_.fetch_particles(charged_particles, snemo::datamodel::particle_track::NEGATIVE |
                                              snemo::datamodel::particle_track::POSITIVE |
                                              snemo::datamodel::particle_track::UNDEFINED);
  if (charged_particles.empty()) {
    DT_LOG_DEBUG(get_logging_priority(), "No charged particles have been found !");
    return 0;
  }

  for (snemo::datamodel::particle_track_data::particle_collection_type::const_iterator iparticle =
           charged_particles.begin();
       iparticle != charged_particles.end(); ++iparticle) {
    const snemo::datamodel::particle_track& a_particle = iparticle->get();

    // No calorimeter associated, no TOF computation
    if (!a_particle.has_associated_calorimeter_hits()) continue;

    for (snemo::datamodel::particle_track_data::particle_collection_type::iterator igamma =
             gamma_particles.begin();
         igamma != gamma_particles.end(); ++igamma) {
      snemo::datamodel::particle_track& a_gamma = igamma->grab();

      // snemo::datamodel::particle_track::vertex_collection_type & vertices =
      // a_gamma.grab_vertices();
      snemo::datamodel::particle_track::vertex_collection_type the_vertices_2;
      a_gamma.fetch_vertices(the_vertices_2,
                             snemo::datamodel::particle_track::VERTEX_ON_MAIN_CALORIMETER |
                                 snemo::datamodel::particle_track::VERTEX_ON_X_CALORIMETER |
                                 snemo::datamodel::particle_track::VERTEX_ON_GAMMA_VETO);
      if (the_vertices_2.empty()) {
        DT_LOG_DEBUG(get_logging_priority(), "Gamma track has no vertices associated !");
        continue;
      }

      for (snemo::datamodel::particle_track::vertex_collection_type::iterator ivtx =
               the_vertices_2.begin();
           ivtx != the_vertices_2.end(); ++ivtx) {
        const geomtools::blur_spot& a_spot = ivtx->get();
        // Get associated calorimeter:
        if (!a_gamma.has_associated_calorimeter_hits()) {
          DT_LOG_DEBUG(get_logging_priority(),
                       "Gamma track is not associated to any calorimeter block !");
          continue;
        }
        const snemo::datamodel::calibrated_calorimeter_hit::collection_type& hits =
            a_gamma.get_associated_calorimeter_hits();
        geomtools::base_hit::has_geom_id_predicate hit_pred(a_spot.get_geom_id());
        datatools::mother_to_daughter_predicate<geomtools::base_hit,
                                                snemo::datamodel::calibrated_calorimeter_hit>
            pred_M2D(hit_pred);
        datatools::handle_predicate<snemo::datamodel::calibrated_calorimeter_hit> pred_via_handle(
            pred_M2D);
        snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator found =
            std::find_if(hits.begin(), hits.end(), pred_via_handle);
        DT_THROW_IF(
            found == hits.end(), std::logic_error,
            "Calibrated calorimeter hit with id " << a_spot.get_geom_id() << " can not be found");
        const snemo::datamodel::calibrated_calorimeter_hit& a_calo_hit_2 = found->get();
        const double gamma_time = a_calo_hit_2.get_time();
        const double gamma_sigma_time = a_calo_hit_2.get_sigma_time();
        // const double gamma_energy       = a_calo_hit_2.get_energy();
        // const double gamma_sigma_energy = a_calo_hit_2.get_sigma_energy();

        // Perform extrapolation to the source foil given charged particles
        // and favor it wrt the search of gamma from e+/e- annihilation
        if (_add_foil_vertex_extrapolation_) {
          // Get track length
          const snemo::datamodel::base_trajectory_pattern& a_track_pattern =
              a_particle.get_trajectory().get_pattern();
          const geomtools::i_shape_1d& a_shape = a_track_pattern.get_shape();
          const double particle_track_length = a_shape.get_length();
          DT_LOG_DEBUG(get_logging_priority(),
                       "Track length = " << particle_track_length / CLHEP::mm << " mm");

          snemo::datamodel::particle_track::vertex_collection_type vertices;
          a_particle.fetch_vertices(vertices,
                                    snemo::datamodel::particle_track::VERTEX_ON_SOURCE_FOIL);
          if (vertices.empty()) continue;
          const geomtools::vector_3d& a_foil_vertex = vertices.front().get().get_position();

          // Compute theoritical time for the gamma in case it comes from the foil vertex
          const double gamma_track_length = (a_foil_vertex - a_spot.get_position()).mag();
          const double gamma_time_th = gamma_track_length / CLHEP::c_light;

          // Assume particle are electron/positron
          const snemo::datamodel::calibrated_calorimeter_hit::collection_type& the_calorimeters =
              a_particle.get_associated_calorimeter_hits();
          // Only take care of the first associated calorimeter
          const snemo::datamodel::calibrated_calorimeter_hit& a_calo_hit =
              the_calorimeters.front().get();
          const double particle_time = a_calo_hit.get_time();
          const double particle_sigma_time = a_calo_hit.get_sigma_time();
          const double particle_energy = a_calo_hit.get_energy();
          const double particle_sigma_energy = a_calo_hit.get_sigma_energy();
          const double particle_mass = CLHEP::electron_mass_c2;
          const double beta = std::sqrt(particle_energy * (particle_energy + 2. * particle_mass)) /
                              (particle_energy + particle_mass);
          const double particle_time_th = particle_track_length / CLHEP::c_light / beta;
          const double sigma_particle_time_th =
              particle_time_th * std::pow(particle_mass, 2) /
              (particle_energy * (particle_energy + particle_mass) *
               (particle_energy + 2 * particle_mass)) *
              particle_sigma_energy;

          const double dt_int = particle_time - gamma_time - (particle_time_th - gamma_time_th);
          const double sigma = std::pow(particle_sigma_time, 2) + std::pow(gamma_sigma_time, 2) +
                               std::pow(sigma_particle_time_th, 2);
          const double chi2_int = std::pow(dt_int, 2) / sigma;
          const double int_prob = gsl_cdf_chisq_Q(chi2_int, 1) * 100. * CLHEP::perCent;
          const double int_prob_limit = _add_foil_vertex_minimal_probability_;
          if (int_prob > int_prob_limit) {
            DT_LOG_DEBUG(get_logging_priority(),
                         "Adding foil vertex (with internal probability = " << int_prob << ")");

            snemo::datamodel::particle_track::handle_spot hBSv(new geomtools::blur_spot);
            // a_gamma.grab_vertices().push_back(hBSv);
            a_gamma.grab_vertices().insert(a_gamma.grab_vertices().begin(), hBSv);
            geomtools::blur_spot& spot_v = hBSv.grab();
            spot_v.set_hit_id(0);
            spot_v.grab_auxiliaries().store(
                snemo::datamodel::particle_track::vertex_type_key(),
                snemo::datamodel::particle_track::vertex_on_source_foil_label());
            spot_v.set_blur_dimension(geomtools::blur_spot::dimension_three);
            spot_v.set_position(a_foil_vertex);
            break;
          }
        }
        // Perform the search for gamma from e+/e- annihilation
        if (_add_gamma_from_annihilation_) {
          // Get calorimeter hit position associated to charged particle
          const snemo::datamodel::calibrated_calorimeter_hit::collection_type& the_calorimeters =
              a_particle.get_associated_calorimeter_hits();
          // Only take care of the first associated calorimeter
          const snemo::datamodel::calibrated_calorimeter_hit& a_calo_hit =
              the_calorimeters.front().get();
          const double particle_time = a_calo_hit.get_time();
          const double particle_sigma_time = a_calo_hit.get_sigma_time();
          // Get block position and label
          geomtools::vector_3d a_block_position;
          std::string a_label;
          const geomtools::geom_id& a_gid = a_calo_hit.get_geom_id();
          if (get_calo_locator().is_calo_block_in_current_module(a_gid)) {
            get_calo_locator().get_block_position(a_gid, a_block_position);
            a_label = snemo::datamodel::particle_track::vertex_on_main_calorimeter_label();
          } else if (get_xcalo_locator().is_calo_block_in_current_module(a_gid)) {
            get_xcalo_locator().get_block_position(a_gid, a_block_position);
            a_label = snemo::datamodel::particle_track::vertex_on_x_calorimeter_label();
          } else if (get_gveto_locator().is_calo_block_in_current_module(a_gid)) {
            get_gveto_locator().get_block_position(a_gid, a_block_position);
            a_label = snemo::datamodel::particle_track::vertex_on_gamma_veto_label();
          } else {
            DT_THROW_IF(
                true, std::logic_error,
                "Current geom id '" << a_gid << "' does not match any scintillator block !");
          }

          const double track_length = (a_block_position - a_spot.get_position()).mag();
          const double dt_exp = std::abs(particle_time - gamma_time);
          const double dt_th = track_length / CLHEP::c_light;
          const double dt_int = dt_exp - dt_th;
          const double sigma = std::pow(particle_sigma_time, 2) + std::pow(gamma_sigma_time, 2);
          const double chi2_int = std::pow(dt_int, 2) / sigma;
          const double int_prob = gsl_cdf_chisq_Q(chi2_int, 1) * 100. * CLHEP::perCent;
          const double int_prob_limit = _add_gamma_from_annihilation_minimal_probability_;
          if (int_prob > int_prob_limit) {
            DT_LOG_DEBUG(
                get_logging_priority(),
                "Found gamma from annihilation (with internal probability = " << int_prob << ")");
            // Do not add the calorimeter hit from the charged particle to
            // the list of calorimeter hits of the gamma particle
            // snemo::datamodel::calibrated_calorimeter_hit::collection_type & hits =
            // a_gamma.grab_associated_calorimeter_hits(); hits.insert(hits.begin(),
            // the_calorimeters.front());
            snemo::datamodel::particle_track::handle_spot hBSv(new geomtools::blur_spot);
            a_gamma.grab_vertices().insert(a_gamma.grab_vertices().begin(), hBSv);
            geomtools::blur_spot& spot_v = hBSv.grab();
            spot_v.set_hit_id(a_calo_hit.get_hit_id());
            spot_v.set_geom_id(a_calo_hit.get_geom_id());
            spot_v.grab_auxiliaries().store(snemo::datamodel::particle_track::vertex_type_key(),
                                            a_label);
            spot_v.set_blur_dimension(geomtools::blur_spot::dimension_three);
            spot_v.set_position(a_block_position);
            a_gamma.grab_auxiliaries().update_flag("__gamma_from_annihilation");
            break;
          }
        }
      }
    }

    // Only take care of one particle (to be improved later)
    break;
  }

  return 0;
}

// static
void base_gamma_builder::ocd_support(datatools::object_configuration_description& ocd_,
                                     const std::string& prefix_) {
  datatools::logger::declare_ocd_logging_configuration(ocd_, "fatal", prefix_ + "BGB.");

  {
    // Description of the 'locator_plugin_name' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.locator_plugin_name")
        .set_terse_description("The name of the geometry locator plugin to be used")
        .set_from("snemo::processing::base_gamma_builder")
        .set_traits(datatools::TYPE_STRING)
        .set_long_description("Empty value means automatic search   \n")
        .set_mandatory(false)
        .add_example(
            "Set a specific value::                                   \n"
            "                                                         \n"
            "  BGB.locator_plugin_name : string = \"locators_driver\" \n"
            "                                                         \n");
  }

  {
    // Description of the 'BGB.add_foil_vertex_extrapolation' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.add_foil_vertex_extrapolation")
        .set_terse_description("Allow the extrapolation of gamma tracks to the source foil")
        .set_long_description(
            "Given the presence of charged particles with vertex on the source foil          \n"
            "a Time-Of-Flight computation is done to check if gamma are internals            \n"
            "with the charged particle. Then, the charged particle vertex on the source foil \n"
            "is also associated to the gamma.                                                \n")
        .set_from("snemo::processing::base_gamma_builder")
        .set_traits(datatools::TYPE_BOOLEAN)
        .set_default_value_boolean(true)
        .set_mandatory(false)
        .add_example(
            "Set the default value::                              \n"
            "                                                     \n"
            "  BGB.add_foil_vertex_extrapolation : boolean = true \n"
            "                                                     \n");
  }

  {
    // Description of the 'BGB.add_foil_vertex.minimal_probability' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.add_foil_vertex_extrapolation.minimal_probability")
        .set_terse_description("Set the minimal internal TOF probability")
        .set_from("snemo::processing::base_gamma_builder")
        .set_triggered_by_flag("BGB.add_foil_vertex_extrapolation")
        .set_traits(datatools::TYPE_REAL)
        .set_default_value_real(1 * CLHEP::perCent, "%")
        .set_mandatory(false)
        .add_example(
            "Set the default value::                                           \n"
            "                                                                  \n"
            "  BGB.add_foil_vertex_extrapolation.minimal_probability : real as fraction = 1% \n"
            "                                                                  \n");
  }

  {
    // Description of the 'BGB.add_gamma_from_annihilation' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.add_gamma_from_annihilation")
        .set_terse_description("Allow the search and the tagging of gamma from e+/e- annihilation")
        .set_long_description(
            "Given the presence of charged particles, a Time-Of-Flight computation       \n"
            "is done to check if gamma comes from an e+/e- annihilation.                 \n"
            "The TOF calculation is performed on calorimeter hit of the charged particle \n"
            "and the first calorimeter hits associated to gamma.                         \n")
        .set_from("snemo::processing::base_gamma_builder")
        .set_traits(datatools::TYPE_BOOLEAN)
        .set_default_value_boolean(true)
        .set_mandatory(false)
        .add_example(
            "Set the default value::                            \n"
            "                                                   \n"
            "  BGB.add_gamma_from_annihilation : boolean = true \n"
            "                                                   \n");
  }

  {
    // Description of the 'BGB.add_gamma_from_annihilation.minimal_probability' configuration
    // property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.add_gamma_from_annihilation.minimal_probability")
        .set_terse_description("Set the minimal TOF probability for gamma from e+/e- annihilation")
        .set_from("snemo::processing::base_gamma_builder")
        .set_triggered_by_flag("BGB.add_gamma_from_annihilation")
        .set_traits(datatools::TYPE_REAL)
        .set_default_value_real(1 * CLHEP::perCent, "%")
        .set_mandatory(false)
        .add_example(
            "Set the default value::                                                       \n"
            "                                                                              \n"
            "  BGB.add_gamma_from_annihilation.minimal_probability : real as fraction = 1% \n"
            "                                                                              \n");
  }

  {
    // Description of the 'BGB.select_calorimeter_hits' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.select_calorimeter_hits")
        .set_terse_description("Allow the pre-selection of calorimeter hits")
        .set_long_description(
            "The calorimeter hits can be pre-selected before sending them  \n"
            "to the gamma tracking/clustering algorithm based on a list of \n"
            "tags given by the used                                        \n")
        .set_from("snemo::processing::base_gamma_builder")
        .set_traits(datatools::TYPE_BOOLEAN)
        .set_default_value_boolean(false)
        .set_mandatory(false)
        .add_example(
            "Activate the selection::                       \n"
            "                                               \n"
            "  BGB.select_calorimeter_hits : boolean = true \n"
            "                                               \n");
  }

  {
    // Description of the 'BGB.select_calorimeter_hits.tags' configuration property :
    datatools::configuration_property_description& cpd = ocd_.add_property_info();
    cpd.set_name_pattern("BGB.select_calorimeter_hits.tags")
        .set_terse_description("Set the list of tags for calorimeter selection")
        .set_triggered_by_flag("BGB.select_calorimeter_hits")
        .set_from("snemo::processing::base_gamma_builder")
        .set_traits(datatools::TYPE_STRING, datatools::configuration_property_description::ARRAY)
        .set_mandatory(false)
        .add_example(
            "Set values::                                                    \n"
            "                                                                \n"
            "  BGB.select_calorimeter_hits.tags : string[1] = \"__isolated\" \n"
            "                                                                \n");
  }

  return;
}

}  // end of namespace processing

}  // end of namespace snemo
