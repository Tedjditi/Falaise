// -*- mode: c++ ; -*-
/// \file falaise/snemo/datamodels/tracker_clustering_solution.ipp

#ifndef FALAISE_SNEMO_DATAMODELS_TRACKER_CLUSTERING_SOLUTION_IPP
#define FALAISE_SNEMO_DATAMODELS_TRACKER_CLUSTERING_SOLUTION_IPP 1

// Ourselves:
#include <falaise/snemo/datamodels/tracker_clustering_solution.h>

// Third party:
// - Boost
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
// - Bayeux/datatools:
#include <datatools/i_serializable.ipp>
#include <datatools/properties.ipp>

// This project:
#include <falaise/snemo/datamodels/calibrated_tracker_hit.h>
#include <falaise/snemo/datamodels/tracker_cluster.h>

namespace snemo {

namespace datamodel {

/// Serialization
template <class Archive>
void tracker_clustering_solution::serialize(Archive& ar_, const unsigned int /* version_ */) {
  ar_& DATATOOLS_SERIALIZATION_I_SERIALIZABLE_BASE_OBJECT_NVP;
  ar_& boost::serialization::make_nvp("solution_id", _solution_id_);
  ar_& boost::serialization::make_nvp("clusters", _clusters_);
  // if (version_ > 0) {
  //   ar_ & boost::serialization::make_nvp("delayed_clusters", _delayed__clusters_);
  // }
  ar_& boost::serialization::make_nvp("unclustered_hits", _unclustered_hits_);
  ar_& boost::serialization::make_nvp("auxiliaries", _auxiliaries_);
  return;
}

}  // end of namespace datamodel

}  // end of namespace snemo

#endif  // FALAISE_SNEMO_DATAMODELS_TRACKER_CLUSTERING_SOLUTION_IPP
