// -*- mode: c++ ; -*-
/// \file falaise/snemo/datamodels/mock_raw_tracker_hit.h
/* Author(s) :    Francois Mauger <mauger@lpccaen.in2p3.fr>
 * Creation date: 2010-03-15
 * Last modified: 2014-01-30
 *
 *
 * License:
 *
 * Description:
 *
 *   Raw tracker hit
 *
 * History:
 *  - was the 'snemo::core::model::raw_tracker_hit' class in legacy Falaise
 */

#ifndef FALAISE_SNEMO_DATAMODEL_MOCK_RAW_TRACKER_HIT_H
#define FALAISE_SNEMO_DATAMODEL_MOCK_RAW_TRACKER_HIT_H 1

// Third party:
// - Bayeux/geomtools:
#include <geomtools/base_hit.h>

namespace snemo {

namespace datamodel {

/// \brief A mock class to represent SuperNEMO raw tracker hit
class mock_raw_tracker_hit : public geomtools::base_hit {
 public:
  static const double INVALID_VALUE;
  static const std::string NOISY_FLAG;

  enum store_mask_type { STORE_REF_TIME = 0x8, STORE_TIMES = 0x10 };

  bool is_ref_time_missing() const;

  bool is_drift_time_missing() const;

  bool is_top_time_missing() const;

  bool is_bottom_time_missing() const;

  double get_ref_time() const;

  void set_ref_time(double);

  double get_sigma_ref_time() const;

  void set_sigma_ref_time(double);

  double get_drift_time() const;

  void set_drift_time(double);

  double get_sigma_drift_time() const;

  void set_sigma_drift_time(double);

  double get_top_time() const;

  void set_top_time(double);

  double get_sigma_top_time() const;

  void set_sigma_top_time(double);

  double get_bottom_time() const;

  void set_bottom_time(double);

  double get_sigma_bottom_time() const;

  void set_sigma_bottom_time(double);

  void invalidate_ref_time();

  void invalidate_times();

  /// Default constructor
  mock_raw_tracker_hit();

  /// Destructor
  virtual ~mock_raw_tracker_hit();

  /// Check if the hit is valid
  bool is_valid() const;

  /// Invalidate the hit
  void invalidate();

  /// Clear all attributes
  virtual void clear();

  /// Smart print
  virtual void tree_dump(std::ostream& out_ = std::clog, const std::string& title_ = "",
                         const std::string& indent_ = "", bool inherit_ = false) const;

  void dump() const;

 private:
  double _ref_time_;  //!< Relative time of the hit in relation to the absolute timestamp of the
                      //!< current event
  double _sigma_ref_time_;  //!< Relative time error

  // Geiger regime anode/cathodes drift times:
  double _drift_time_;         //!< Anode drift time
  double _sigma_drift_time_;   //!< Anode drift time error
  double _top_time_;           //!< Top cathode signal drift time
  double _sigma_top_time_;     //!< Top cathode signal drift time error
  double _bottom_time_;        //!< Bottom cathode signal drift time
  double _sigma_bottom_time_;  //!< Bottom cathode signal drift time error

  // DATATOOLS_SERIALIZATION_DECLARATION();
};

}  // end of namespace datamodel

}  // end of namespace snemo

//#include <boost/serialization/export.hpp>
// BOOST_CLASS_EXPORT_KEY2(snemo::datamodel::mock_raw_tracker_hit,
// "snemo::datamodel::mock_raw_tracker_hit")

#endif  // FALAISE_SNEMO_DATAMODEL_MOCK_RAW_TRACKER_HIT_H
