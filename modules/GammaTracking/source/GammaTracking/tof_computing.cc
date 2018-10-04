// Ourselves:
#include "tof_computing.h"

// Standard library
#include <cmath>

// - Third party:
#include <fstream>
#include <iostream>

// - GSL:
#include <gsl/gsl_cdf.h>
// - Bayeux/datatools:
#include <datatools/clhep_units.h>


namespace gt {


  double tof_computing::Get_Distance(Position R1, Position R2) {

    return sqrt( (R2.x-R1.x)*(R2.x-R1.x) + (R2.y - R1.y)*(R2.y - R1.y) + (R2.z - R1.z)*(R2.z - R1.z) );

  }
  
  double tof_computing::get_theta(Position Phit2_,Position Phit1_){

    double hyp = tof_computing::Get_Distance(Phit2_,Phit1_);

    Phit1_.x = Phit2_.x;
    double cad = tof_computing::Get_Distance(Phit2_,Phit1_);

    return (180/3.14159265)*acos(cad/hyp);
      
  }
  
  double tof_computing::get_theta(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_){

    Position Phit2_;
    Phit2_.x = hit2_.position.getX();
    Phit2_.y = hit2_.position.getY();
    Phit2_.z = hit2_.position.getZ();
    
    Position Phit1_;
    Phit1_.x = hit1_.position.getX();
    Phit1_.y = hit1_.position.getY();
    Phit1_.z = hit1_.position.getZ();
    
    return tof_computing::get_theta(Phit2_,Phit1_);

  }
  
  double tof_computing::get_sigma_tot(double energy, double theta){
    
    bool check_theta = false;
    bool check_energy = false;
    double _theta_, _energy_ , _sigmatot_, _mean_t_, stock = 0;
    
    std::ifstream Sfile("/data/soft/Falaise.git/modules/GammaTracking/source/GammaTracking/SigmaToT_table");
    while(!Sfile.eof()){
      Sfile >> _theta_;
      Sfile >> _energy_;
      Sfile >> _sigmatot_;
      Sfile >> _mean_t_;

      for (int it_f=0; it_f < 23; it_f++){
	
	double theta_max = 90 - (180/3.14159265)*(atan( (it_f*259.0 ) /901));
	double theta_min = 90 - (180/3.14159265)*(atan( ((it_f+1)*259.0 )/901));

	if ( ( theta <= theta_max && _theta_<=theta_max ) && (theta > theta_min && _theta_ > theta_min) ){

	  check_theta = true;
	  break;
	}
      }

      if (!check_theta) continue;
    
      int S_int_E = 200; //keV
      
      for (int it_E=0; it_E < 25; it_E++){

      	double E_min = it_E*(S_int_E*1.0/1000);
	double E_max = (it_E+1)*(S_int_E*1.0/1000);

	if ( (energy < E_max && _energy_ < E_max) && (energy >= E_min && _energy_ >= E_min)){

	  check_energy = true;
	  break;
	}
      }
      
      if (check_theta && check_energy) {
	stock = _sigmatot_;
	Sfile.close();
	return stock;
      }
    }
    
    std::cout << "Error 404: SigmaTOT not FOUND" << std::endl;
    std::cout << "Theta Recherche = " << theta << " Energy Recherche = " <<energy<<std::endl<<std::endl;
    
    Sfile.close();
    
    return stock;
  }

  double tof_computing::get_mean_interaction_shift(double energy, double theta){
    
    bool check_theta = false;
    bool check_energy = false;
    double _theta_, _energy_ , _sigmatot_, _mean_t_, stock = 0;
    
    std::ifstream Sfile("/data/soft/Falaise.git/modules/GammaTracking/source/GammaTracking/SigmaToT_table");
    while(!Sfile.eof()){
      Sfile >> _theta_;
      Sfile >> _energy_;
      Sfile >> _sigmatot_;
      Sfile >> _mean_t_;

      for (int it_f=0; it_f < 23; it_f++){
	
	double theta_max = 90 - (180/3.14159265)*(atan( (it_f*259.0 ) /901));
	double theta_min = 90 - (180/3.14159265)*(atan( ((it_f+1)*259.0 )/901));

	if ( ( theta <= theta_max && _theta_<=theta_max ) && (theta > theta_min && _theta_ > theta_min) ){

	  check_theta = true;
	  break;
	}
      }

      if (!check_theta) continue;
    
      int S_int_E = 200; //keV
      
      for (int it_E=0; it_E < 25; it_E++){

      	double E_min = it_E*(S_int_E*1.0/1000);
	double E_max = (it_E+1)*(S_int_E*1.0/1000);

	if ( (energy < E_max && _energy_ < E_max) && (energy >= E_min && _energy_ >= E_min)){

	  check_energy = true;
	  break;
	}
      }
      
      if (check_theta && check_energy) {
	stock = _mean_t_;
	Sfile.close();
	return stock;
      }
    }
    
    std::cout << "Error 404: mean_t not FOUND" << std::endl;
    std::cout << "Theta Recherche = " << theta << " Energy Recherche = " <<energy<<std::endl<<std::endl;
    
    Sfile.close();
    
    return stock;
  }

  double tof_computing::get_mean_interaction_shift_of_2nd_calo(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_){

    double theta_ = abs(tof_computing::get_theta(hit1_,hit2_));
    if (theta_ == 0) return 0;
    return tof_computing::get_sigma_tot(hit2_.energy,theta_);

  }

  double tof_computing::get_mean_interaction_shift(const event::calorimeter_hit& hit1_){
    return tof_computing::get_sigma_tot(hit1_.energy);
  }

  double tof_computing::get_sigma_tot_of_2nd_calo(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_){

    double theta_ = abs(tof_computing::get_theta(hit1_,hit2_));
    if (theta_ == 0) return 1000;
    return tof_computing::get_sigma_tot(hit2_.energy,theta_);

  }

  double tof_computing::get_sigma_tot(const event::calorimeter_hit& hit1_){
    return tof_computing::get_sigma_tot(hit1_.energy);
  }
    
  double tof_computing::beta(double energy_, double mass_) {
    return std::sqrt(energy_ * (energy_ + 2. * mass_)) / (energy_ + mass_);
  }

  double tof_computing::get_theoritical_time(double energy_, double mass_, double track_length_) {
    return track_length_ / (tof_computing::beta(energy_, mass_) * CLHEP::c_light);
  }

  double tof_computing::get_track_length(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_) {
    return (hit1_.position - hit2_.position).mag();
  }

  double tof_computing::get_delta_time(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_) {

    const double track_length = tof_computing::get_track_length(hit1_, hit2_);
    const double t1 = hit1_.time;// + tof_computing::get_mean_interaction_shift(hit1_);
    const double t2 = hit2_.time;// + tof_computing::get_mean_interaction_shift_of_2nd_calo(hit1_,hit2_);
    const double th_time = tof_computing::get_theoritical_time(hit1_.energy, 0, track_length);
    // the mass being 0, it just gets t=l/c (since beta=1)
    return (t2 - t1) - th_time;
  }

  double tof_computing::get_sigma_length(double x, double p0, double p1, double p2, double p3, double intersection){

    if (x >= intersection)
      return (-1*p0*sqrt(x)/(x*x) + p1);

    if (x < intersection)
      return (p2*x+p3);

    if (x<0)
      std::cout<<"X Negatif"<<std::endl;

    return 0;

  }

  double tof_computing::get_sigma_length(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_){

    double track_length = tof_computing::get_track_length(hit1_, hit2_);
    double sigma_l = tof_computing::get_sigma_length(track_length);
    return sigma_l;
  }
  
  double tof_computing::get_chi2(const event::calorimeter_hit& hit1_, const event::calorimeter_hit& hit2_) {

   

    double sigma_exp = 0;
    
    const double sigma_length = tof_computing::get_sigma_length(hit1_, hit2_); 

    //const double sigma_exp = pow(hit1_.sigma_time, 2) + pow(hit2_.sigma_time, 2) + pow(sigma_length, 2); 

    
    const double sigma_tot_1 = tof_computing::get_sigma_tot(hit1_);
    const double sigma_tot_2 = tof_computing::get_sigma_tot_of_2nd_calo(hit1_,hit2_);

    if (sigma_tot_2 == 1000) {
      sigma_exp = pow(hit2_.sigma_time, 2) + pow(sigma_tot_1, 2) + pow(sigma_length, 2); 
    }
    else {
      sigma_exp = pow(sigma_tot_1, 2) + pow(sigma_tot_2, 2);
    }
    
    return pow(tof_computing::get_delta_time(hit1_, hit2_), 2) / sigma_exp;
  }


  double tof_computing::get_internal_probability(double chi2_, size_t ndf_) {
    return gsl_cdf_chisq_Q(chi2_, ndf_);
  }

}  // namespace gt
