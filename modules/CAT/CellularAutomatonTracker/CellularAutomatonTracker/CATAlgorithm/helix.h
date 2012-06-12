/* -*- mode: c++ -*- */
#ifndef __CATAlgorithm__IHELIX
#define __CATAlgorithm__IHELIX
#include <iostream>
#include <cmath>
#include <mybhep/error.h>
#include <mybhep/utilities.h>
#include <mybhep/point.h>
#include <CATAlgorithm/experimental_point.h>
#include <CATAlgorithm/experimental_vector.h>
#include <CATAlgorithm/circle.h>


namespace CAT {
  namespace topology{

    using namespace std;
    using namespace mybhep;

    class helix : public tracking_object {

      // a helix is identified by origin,
      // radius
      // and pitch

    private:
      string appname_;

      // experimental point
      experimental_point center_; 

      // radius
      experimental_double radius_;

      // pitch
      experimental_double pitch_;


      // points in the helix are given by:
      // x(phi) = center_.x() + radius*cos(phi)
      // y(phi) = center_.y() + pitch*phi
      // z(phi) = center_.z() + radius*sin(phi)


    public:   

      //!Default constructor 
      helix(prlevel level=mybhep::NORMAL, double nsigma=10.)
      {
        appname_= "helix: ";
        center_ = experimental_point();
        radius_ = experimental_double(small_neg, small_neg);
        pitch_ = experimental_double(small_neg, small_neg);
        set_print_level(level);
        set_nsigma(nsigma);
      }

      //!Default destructor
      virtual ~helix(){};

      //! constructor
      helix(experimental_point center, experimental_double radius, experimental_double  pitch, prlevel level=mybhep::NORMAL, double nsigma=10.){
        set_print_level(level);
        set_nsigma(nsigma);
        appname_= "helix: ";
        center_ = center;
        radius_ = radius;
        pitch_ = pitch;
      }

      //! constructor from a circle and a pitch
      helix(circle c, experimental_double  pitch, prlevel level=mybhep::NORMAL, double nsigma=10.){
        set_print_level(level);
        set_nsigma(nsigma);
        appname_= "helix: ";
        center_ = c.center();
        radius_ = c.radius();
        pitch_ = pitch;
      }

      /*** dump ***/
      virtual void dump (ostream & a_out         = clog,
                         const string & a_title  = "",
                         const string & a_indent = "",
                         bool a_inherit          = false)const{
        {
          string indent;
          if (! a_indent.empty ()) indent = a_indent;
          if (! a_title.empty ())
            {
              a_out << indent << a_title << endl;
            }

          a_out << indent << appname_ << " -------------- " << endl;
          a_out << indent << " center " << endl;
          this->center().dump(a_out, "", indent + "    ");
          a_out << indent << " radius: "; radius().dump(); a_out << " " << endl;
          a_out << indent << " pitch: "; pitch().dump(); a_out << " " << endl;
          a_out << indent << " one round: " << endl;
          for(size_t i=0; i<100; i++){
            experimental_double theta(i*3.1417/100., 0.);
            a_out << indent << " .. theta " << theta.value()*180./M_PI << " x " << position(theta).x().value() << " , z " << position(theta).z().value() << endl;
          }
          a_out << indent << " -------------- " << endl;

          return;
        }
      }



      //! set 
      void set(experimental_point center, experimental_double radius, experimental_double pitch)
      {
        center_ = center;
        radius_ = radius;
        pitch_ = pitch;
      }


      //! set center
      void set_center(experimental_point center)
      {
        center_ = center;
      }

      //! set radius
      void set_radius(experimental_double radius)
      {
        radius_ = radius;
      }

      //! set pitch
      void set_pitch(experimental_double pitch)
      {
        pitch_ = pitch;
      }


      //! get center
      const experimental_point& center()const
      {
        return center_;
      }      

      //! get radius
      const experimental_double& radius()const
      {
        return radius_;
      }      

      //! get pitch
      const experimental_double& pitch()const
      {
        return pitch_;
      }      

      //! get curvature
      experimental_double curvature()const
      {
        return 1./radius_;
      }      

      //! get circle
      circle get_circle()const
      {
        return circle(center(), radius(), print_level(), nsigma());
      }      

      // get the phi of a point
      experimental_double phi_of_point(experimental_point ep)const{

        experimental_double phi = experimental_vector(center_, ep).phi();

        if( phi.value() < 0. )
          phi.set_value(phi.value() + 2.*M_PI);

        return phi;

      }

      // get delta-phi of two points
      double delta_phi(experimental_point epa, experimental_point epb)const{

        double phi1 = experimental_vector(center_, epa).phi().value();
        double phi2 = experimental_vector(center_, epb).phi().value();

        fix_angles(&phi1, &phi2);

        return fabs(phi1 - phi2);

      }

      // get the position at parameter phi
      experimental_point position(experimental_double phi)const{
        experimental_double deltax = experimental_cos(phi)*radius();
        experimental_double deltay = phi*pitch();
        experimental_double deltaz = experimental_sin(phi)*radius();

        return (experimental_vector(center()) + experimental_vector(deltax, deltay, deltaz)).point_from_vector();
      }

      // get the position at the theta of point p
      experimental_point position(experimental_point ep)const{
        return position(phi_of_point(ep));
      }

      // get the chi2 with point p
      double chi2(experimental_point ep)const{
        experimental_point predicted = position(phi_of_point(ep));
        experimental_vector residual(ep , predicted);
        double res2 = square(residual.x().value()/residual.x().error()) +
          square(residual.y().value()/residual.y().error()) +
          square(residual.z().value()/residual.z().error());

        if( print_level() >= mybhep::VVERBOSE ){
          clog << " input point: ( "; ep.x().dump(); clog << " , "; ep.y().dump(); clog << " , "; ep.z().dump(); clog << " ) helix: ("; predicted.x().dump(); clog << " , "; predicted.y().dump(); clog << " , ";predicted.z().dump(); clog << " ) local chi2: " << res2  << " " << endl; 
        }

        return res2;
      }

      // compare with set of points and get total chi2
      double chi2(std::vector<experimental_point> ps)const{

        double _chi2 = 0.;
        for(std::vector<experimental_point>::const_iterator ip = ps.begin(); ip != ps.end(); ++ip){

          _chi2 += chi2(*ip);
        }

        return _chi2;

      }

      // compare with set of points and get set of chi2
      std::vector<double> chi2s(std::vector<experimental_point> ps)const{

        std::vector<double> _chi2s;
        for(std::vector<experimental_point>::const_iterator ip = ps.begin(); ip != ps.end(); ++ip){

          _chi2s.push_back(chi2(*ip));
        }

        return _chi2s;

      }


      bool intersect_plane(plane pl, experimental_point * ep, experimental_double _phi)const{
      
        if( pl.view() == "x" || pl.view() == "z" ){
          bool result = get_circle().intersect_plane(pl, ep, _phi);
          ep->set_y(position(*ep).y());
        
          return result;
        
        }
      
        if( pitch().value() == 0.)
          return false;
      
        *ep = position((pl.center().y() - center().y())/pitch());
        ep->set_y(pl.center().y());
        return true;
      
      
      }

      bool intersect_circle(circle c, experimental_point * ep, experimental_double _phi)const{
      
	bool result = get_circle().intersect_circle(c, ep, _phi);
	ep->set_y(position(*ep).y());
        
	return result;
        
      }



    };

    // average
    inline helix average (const std::vector<helix> vs)
    {
      helix mean;
    
      std::vector<experimental_double> radii;
      std::vector<experimental_double> pitches;
      std::vector<experimental_point> centers;
      for(std::vector<helix>::const_iterator iv=vs.begin(); iv!=vs.end(); ++iv){
        radii.push_back(iv->radius());
        pitches.push_back(iv->pitch());
        centers.push_back(iv->center());
      }

      return helix(average(centers), average(radii), average(pitches));
    }


    // get helix through three points
    inline helix three_points_helix(experimental_point epa, experimental_point epb, experimental_point epc){
    
      ////////////////////////////////////////////////////////////////////////
      //                                                                    //
      //  see http://local.wasp.uwa.edu.au/~pbourke/geometry/circlefrom3/   //
      //                                                                    //
      ////////////////////////////////////////////////////////////////////////
    
      experimental_double ma = experimental_vector(epa, epb).tan_phi();
      experimental_double mb = experimental_vector(epb, epc).tan_phi();
    
      experimental_double Xc = (ma*mb*(epa.z() - epc.z()) + mb*(epa.x() + epb.x()) - ma*(epb.x() + epc.x()))/((mb - ma)*2);
      experimental_double Zc;
      if( ma.value() != 0. )
        Zc = (epa.z() + epb.z())/2. - (Xc - (epa.x() + epb.x())/2.)/ma;
      else
        Zc = (epb.z() + epc.z())/2. - (Xc - (epb.x() + epc.x())/2.)/mb;
    
    
      experimental_double _radius = experimental_sqrt(experimental_square(Xc - epa.x()) + experimental_square(Zc - epa.z()));
    
      if( isnan(_radius.value()) || isinf(_radius.value()) )
        _radius.set_value(small_neg);
    
      experimental_double dist = epc.distance(epa);
    
      experimental_double deviation;
      if( dist.value() < 2.*_radius.value() )
        deviation = experimental_asin(dist/(_radius*2.))*2.;
      else
        deviation = experimental_asin(experimental_double(1.,0.))*2.;
    
      if( experimental_vector(epa, epb).phi().value() > experimental_vector(epb, epc).phi().value() )
        deviation.set_value(-deviation.value());

      experimental_double _pitch = (epc.y() - epa.y())/deviation;
    
      experimental_point _center(Xc, experimental_double(0.,0.), Zc);
      experimental_double Yc = epb.y() - _pitch*experimental_vector(_center, epb).phi();
      _center.set_y(Yc);

      helix h(_center,_radius,_pitch);

      return h;
    
    }


  }

}

#endif
