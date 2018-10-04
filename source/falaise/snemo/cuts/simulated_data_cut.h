/** \file falaise/snemo/cuts/simulated_data_cut.h
 * Author(s)     : Francois Mauger <mauger@lpccaen.in2p3.fr>
 * Creation date : 2011-09-18
 * Last modified : 2015-06-20
 *
 * Copyright (C) 2011-2015 Francois Mauger <mauger@lpccaen.in2p3.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 *
 *
 * Description:
 *
 *   Simulated data cut.
 *
 * History:
 *
 */

#ifndef FALAISE_SNEMO_CUT_SIMULATED_DATA_CUT_H
#define FALAISE_SNEMO_CUT_SIMULATED_DATA_CUT_H 1

// Standard library:
#include <string>

// Third party:
// - Boost:
#include <boost/cstdint.hpp>
// - Bayeux/datatools:
#include <datatools/bit_mask.h>
// - Bayeux/mctools:
#include <bayeux/mctools/simulated_data.h>


// - Bayeux/cuts:
#include <cuts/i_cut.h>
// - Bayeux/geomtools:
#include <bayeux/geomtools/clhep.h>

namespace datatools {
class service_manager;
class properties;
}  // namespace datatools

namespace snemo {

namespace cut {

class simulated_data_cut : public cuts::i_cut {
 public:
  /// \brief The cut mode
  enum mode_type {
    MODE_UNDEFINED = 0,
    MODE_FLAG = datatools::bit_mask::bit00,
    MODE_HAS_HIT_CATEGORY = datatools::bit_mask::bit01,  // simulated_data::has_step_hits
    MODE_RANGE_HIT_CATEGORY =
        datatools::bit_mask::bit02,                     // simulated_data::get_number_of_step_hits
    MODE_HAS_HIT_PROPERTY = datatools::bit_mask::bit03,  //
    MODE_GAMMA_HAS_GONE_STRAIGHT = datatools::bit_mask::bit04,
    MODE_HAS_HIT_ON_WALL = datatools::bit_mask::bit05
  };

  /// Set the SD bank key
  void set_SD_label(const std::string& SD_label_);

  /// Return the SD bank key
  const std::string& get_SD_label() const;

  /// Return the cut mode
  uint32_t get_mode() const;

  /// Check mode MODE_FLAG:
  bool is_mode_flag() const;

  /// Check mode MODE_HAS_HIT_CATEGORY:
  bool is_mode_has_hit_category() const;

  /// Check mode MODE_RANGE_HIT_CATEGORY:
  bool is_mode_range_hit_category() const;

  /// Check mode MODE_HAS_HIT_PROPERTY:
  bool is_mode_has_hit_property() const;

  /// Check mode _MODE_GAMMA_HAS_GONE_STRAIGHT:
  bool is_mode_gamma_has_gone_straight() const;

  /// Check mode _MODE_HAS_HIT_ON_WALL:
  bool is_mode_has_hit_on_wall() const;
  
  /// Set the name of cut mode MODE_FLAG
  void set_flag_name(const std::string& flag_name_);

  /// Return the name of cut mode MODE_FLAG
  const std::string& get_flag_name() const;

  /// Constructor
  simulated_data_cut(datatools::logger::priority logging_priority_ = datatools::logger::PRIO_FATAL);

  /// Destructor
  virtual ~simulated_data_cut();

  /// Initilization
  virtual void initialize(const datatools::properties& configuration_,
                          datatools::service_manager& service_manager_,
                          cuts::cut_handle_dict_type& cut_dict_);

  //  Check if gamma_has_gone_straight
  bool is_inside_volume(const geomtools::vector_3d& Position);
  bool is_colinear(mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_min, mctools::simulated_data::hit_handle_collection_type::const_iterator it_t_max);
  bool is_colinear(const mctools::simulated_data::hit_handle_collection_type& step_handle_track, const mctools::simulated_data::hit_handle_collection_type& step_handle_calo, int track_id);
  
  bool has_hit_on_target_wall(const mctools::simulated_data::hit_handle_collection_type& step_handle, double wall);
  bool is_on_wall(mctools::simulated_data::hit_handle_collection_type::const_iterator it_step, double wall);

  /// Reset
  virtual void reset();

 protected:
  /// Default values
  void _set_defaults();

  /// Selection
  virtual int _accept();

 private:
  std::string _SD_label_;  //!< Name of the "Simulated data" bank
  uint32_t _mode_;         //!< Mode of the cut

  std::string _flag_name_;  //!< Name of the boolean property in the simulated data

  std::string _hit_category_;    //!< Name of the hit category to be checked
  int _hit_category_range_min_;  //!< Minimal number of hits in a category
  int _hit_category_range_max_;  //!< Maximal number of hits in a category
  double _gamma_epsilonY_;
  double _gamma_epsilonZ_;
  int _first_;
  int _last_;
  double _wall_;
  std::string _hit_property_logic_;  //!< Logic operation between property selection
  typedef std::map<std::string, std::vector<std::string> > property_values_dict_type;
  property_values_dict_type
      _hit_property_values_;  //!< Values of the 'step_hit' property to look for

  // Macro to automate the registration of the cut :
  CUT_REGISTRATION_INTERFACE(simulated_data_cut)
};

}  // end of namespace cut

}  // end of namespace snemo

// OCD support::
#include <datatools/ocd_macros.h>

// @arg snemo::cut::simulated_data_cut the name the registered class in the OCD system
DOCD_CLASS_DECLARATION(snemo::cut::simulated_data_cut)

#endif  // FALAISE_SNEMO_CUT_SIMULATED_DATA_CUT_H

/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** End: --
*/
