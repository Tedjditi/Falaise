#ifndef GT_TOF_COMPUTING_H
#define GT_TOF_COMPUTING_H 1

#include "event.h"

namespace gt {

  struct Position {
    double x;
    double y;
    double z;
  };
  
  class tof_computing {

  public:

    static double Get_Distance(Position R1,Position R2);
    static double get_theta(Position Phit2_,Position Phit1_);
    static double get_theta(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);
    static double get_sigma_tot(double energy, double theta = 45);
    static double get_sigma_tot(const event::calorimeter_hit& hit1_);
    static double get_sigma_tot_of_2nd_calo(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);

    static double get_mean_interaction_shift(double energy, double theta = 45);
    
    static double get_mean_interaction_shift(const event::calorimeter_hit& hit1_);
    static double get_mean_interaction_shift_of_2nd_calo(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);

    static double beta(double energy_, double mass_);

    static double get_theoritical_time(double energy_, double mass_, double track_length_);

    static double get_track_length(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);

    static double get_delta_time(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);

    static double get_chi2(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);

    static double get_internal_probability(double chi2_, size_t ndf_ = 1);

    static double get_sigma_length(double x, double p0=0.0301547, double p1=0.105565, double p2 = 0.0105275, double p3=0.066, double intersection=0.931309);
    
    static double get_sigma_length(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_);
  };

}  // namespace gt
#endif  // GT_TOF_COMPUTING_H
