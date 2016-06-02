//////////////////////////////////////////////////////////////
// Name:      MCTools.h
// Date:      25 March 2016
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef MCTOOLS_H
#define MCTOOLS_H

// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larsim/Simulation/sim.h"
#include "larsim/Simulation/SimChannel.h"
#include "SimulationBase/MCParticle.h"

// ROOT includes
#include "TVector3.h"

// C++ includes
#include <utility>
#include <vector>

namespace EMShowerContainment {

  class MCTools {

   public:

    // constructor
    MCTools(fhicl::ParameterSet const& pset);

    // destructor
    ~MCTools();

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset);

    // returns dimensions of TPC
    double const detectorLength_();
    double const detectorWidth_();
    double const detectorHeight_();
    double const detectorDiagonal_();

    // returns true if point (x, y, z) is inside TPC
    bool isInsideTPC_(double const& x,
                      double const& y,
                      double const& z);

    // set fiducial cuboid volume
    void setFiducialCuboid_(double const (&fiducialCuboidBoundaryLower)[3],
                            double const (&fiducialCuboidBoundaryUpper)[3]);

    // returns true if point (x, y, z) is inside fiducial cuboid volume
    bool isInsideFiducialCuboid_(double const& x,
                                 double const& y,
                                 double const& z);

    // returns distance between two 3D points (x0, y0, z0) and (x1, y1, z1)
    double distance_(double const& x0,
                     double const& y0,
                     double const& z0,
                     double const& x1,
                     double const& y1,
                     double const& z1);

    // returns 3D spatial distance between two four-vectors
    double distance3D_(double const (& aXYZT)[4],
                       double const (& bXYZT)[4]);

    // returns 4D distance between two four-vectors
    double distance4D_(double const (& aXYZT)[4],
                       double const (& bXYZT)[4]);

    // returns 3D spatial angle between two four-vectors
    double angle3D_(double const (& aXYZT)[4],
                    double const (& bXYZT)[4]);

    // get four-vectors of start/end positions and momenta from particle
    void getFourVectors_(simb::MCParticle const& particle,
                         double               (& startXYZT)[4],
                         double               (& endXYZT)[4],
                         double               (& startPE)[4],
                         double               (& endPE)[4]);

    // get length of a particle's trajectory
    double getTrajectoryLength_(simb::MCParticle const& particle);

    // get daughter particles of showers
    void getShowerDaughters_(int                                      const& trackID,
                             std::map< int, const simb::MCParticle* > const& particleMap,
                             std::map< int, const simb::MCParticle* >      & showerParticleMap);

    // get distance from \vec{a} to TPC edge in the direction of
    // (\vec{a} - \vec{b})
    double distanceToTPCEdge_(double const (& aXYZT)[4],
                              double const (& bXYZT)[4]);

    // get distance from \vec{a} to TPC boundary in the direction
    // of (\vec{a} - \vec{b})
    void distanceToTPCBoundary_(double const (& aXYZT)[4],
                                double const (& bXYZT)[4],
                                double        & distanceX,
                                double        & distanceY,
                                double        & distanceZ,
                                double        & distance);

    // get longitudinal shower profile
    void longitudinalShowerProfile_(std::map< int, const simb::MCParticle* >           const& showerParticleMap,
                                    art::ValidHandle< std::vector< sim::SimChannel > > const& simChannelHandle,
                                    double                                             const& dz,
                                    //double                                             const& maxDepth,
                                    TVector3                                           const& startPoint,
                                    TVector3                                           const& longitudinalDirection,
                                    std::vector< double >                                   & longitudinalProfile);

    // get radial shower profile
    void radialShowerProfile_(std::map< int, const simb::MCParticle* >           const& showerParticleMap,
                              art::ValidHandle< std::vector< sim::SimChannel > > const& simChannelHandle,
                              double                                             const& dr,
                              //double                                             const& maxRadius,
                              TVector3                                           const& startPoint,
                              std::vector< double >                                   & radialProfile);

    // get shower variables
    void getShowerVariables_(int                                                const& trackID,
                             std::map< int, const simb::MCParticle* >           const& particleMap,
                             art::ValidHandle< std::vector< sim::SimChannel > > const& simChannelHandle,
                             std::map< int, const simb::MCParticle* >                & showerParticleMap,
                             double                                                 (& startXYZT)[4],
                             double                                                 (& endXYZT)[4],
                             double                                                 (& startPE)[4],
                             double                                                 (& endPE)[4],
                             double                                                  & photonAngleZ,
                             bool                                                    & photonConvert,
                             double                                                  & photonConversionLength,
                             double                                                  & photonEnergy,
                             double                                                  & photonEnergyDeposited,
                             double                                                  & photonTPCContainment,
                             double                                                  & photonTPCDistanceX,
                             double                                                  & photonTPCDistanceY,
                             double                                                  & photonTPCDistanceZ,
                             double                                                  & photonTPCDistance,
                             double                                                  & photonSphereEnergy,
                             double                                                  & photonSphereContainment,
                             double                                                  & photonSphereRadius);

    // this method is used for testing porpoises
    void hello_world();

   private:

    // this method is used for testing dolphins
    void hello_kitty();

    // pointer to geometry provider
    geo::GeometryCore const* fGeometry;

    // dimensions of TPC
    double fDetectorLength;
    double fDetectorWidth;
    double fDetectorHeight;
    double fDetectorDiagonal;

    // TPC boundaries
    double fTPCBoundaryLower[3];
    double fTPCBoundaryUpper[3];

    // fiducial volume boundaries
    double fFiducialCuboidBoundaryLower[3];
    double fFiducialCuboidBoundaryUpper[3];
    //double fFiducialCylinderCenter;
    //double fFiducialCylinderRadius;
    //double fFiducialCylinderHeight;

    // maximum distance to boundary 1D
    double fMaxDistanceToBoundary1D;

    bool fDebug;  // print messages out for debugging

  };

} // namespace EMShowerContainment

#endif
