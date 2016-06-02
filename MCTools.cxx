//////////////////////////////////////////////////////////////
// Name:      MCTools.cxx
// Date:      25 March 2016
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

// Class include
#include "uboone/EMShowerContainment/MCTools.h"

// Framework includes
#include "art/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "SimulationBase/MCTruth.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TVector3.h"

// C++ includes
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>

namespace EMShowerContainment {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  MCTools::MCTools(fhicl::ParameterSet const& pset)
  {
    // get a pointer to the geometry service provider
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

    // get dimensions of TPC
    fDetectorLength =      fGeometry->DetLength();
    fDetectorWidth  = 2. * fGeometry->DetHalfWidth();
    fDetectorHeight = 2. * fGeometry->DetHalfHeight();
    fDetectorDiagonal = std::sqrt(fDetectorLength*fDetectorLength +
                                  fDetectorWidth*fDetectorWidth +
                                  fDetectorHeight*fDetectorHeight);

    // TPC boundaries
    fTPCBoundaryLower[0] = 0;                  // x
    fTPCBoundaryLower[1] = -fDetectorHeight/2; // y
    fTPCBoundaryLower[2] = 0;                  // z
    fTPCBoundaryUpper[0] = fDetectorWidth;     // x
    fTPCBoundaryUpper[1] = fDetectorHeight/2;  // y
    fTPCBoundaryUpper[2] = fDetectorLength;    // z

    if (fDebug)
    {
      mf::LogVerbatim("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "  DetectorLength:      " << fDetectorLength      << std::endl
        << "  DetectorWidth:       " << fDetectorWidth       << std::endl
        << "  DetectorHeight:      " << fDetectorHeight      << std::endl
        << "  DetectorDiagonal:    " << fDetectorDiagonal    << std::endl
        << "  TPCBoundaryLower[0]: " << fTPCBoundaryLower[0] << std::endl
        << "  TPCBoundaryLower[1]: " << fTPCBoundaryLower[1] << std::endl
        << "  TPCBoundaryLower[2]: " << fTPCBoundaryLower[2] << std::endl
        << "  TPCBoundaryUpper[0]: " << fTPCBoundaryUpper[0] << std::endl
        << "  TPCBoundaryUpper[1]: " << fTPCBoundaryUpper[1] << std::endl
        << "  TPCBoundaryUpper[2]: " << fTPCBoundaryUpper[2] << std::endl
        << "//////////////////////////////////////////////////////////////\n";
    }

    // read in parameters from .fcl files
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  MCTools::~MCTools() {}

  //-----------------------------------------------------------------------
  void MCTools::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMaxDistanceToBoundary1D = pset.get< double >("MaxDistanceToBoundary1D", 1.0 /* cm */);
    fDebug = pset.get< bool >("Debug", true);

    return;
  }

  //-----------------------------------------------------------------------
  void MCTools::hello_world()
  {
    mf::LogVerbatim("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "Hello, World!\n"
        << "//////////////////////////////////////////////////////////////\n";

    this->hello_kitty();

    return;
  }

  //-----------------------------------------------------------------------
  void MCTools::hello_kitty()
  {
    return;
  }


  //-----------------------------------------------------------------------
  double const MCTools::detectorLength_()
  {
    return fDetectorLength;
  }


  //-----------------------------------------------------------------------
  double const MCTools::detectorWidth_()
  {
    return fDetectorWidth;
  }


  //-----------------------------------------------------------------------
  double const MCTools::detectorHeight_()
  {
    return fDetectorHeight;
  }


  //-----------------------------------------------------------------------
  double const MCTools::detectorDiagonal_()
  {
    return fDetectorDiagonal;
  }


  //-----------------------------------------------------------------------
  bool MCTools::isInsideTPC_(double const& x,
                             double const& y,
                             double const& z)
  {

    if (x >=  0.                  &&
        x <=  fDetectorWidth      &&
        y >= -fDetectorHeight/2.  &&
        y <=  fDetectorHeight/2.  &&
        z >=  0.                  &&
        z <=  fDetectorLength)
    {
      return true;
    }
    else
    {
      return false;
    }

    return false;

  } // MCTools::isInsideTPC_()


  //-----------------------------------------------------------------------
  void MCTools::setFiducialCuboid_(double const (&fiducialCuboidBoundaryLower)[3],
                                   double const (&fiducialCuboidBoundaryUpper)[3])
  {
    for (size_t i = 0; i < 3; ++i)
    {
      fFiducialCuboidBoundaryLower[i] = fiducialCuboidBoundaryLower[i];
      fFiducialCuboidBoundaryUpper[i] = fiducialCuboidBoundaryUpper[i];
    }
  } // MCTools::setFiducialCuboid_()


  //-----------------------------------------------------------------------
  bool MCTools::isInsideFiducialCuboid_(double const& x,
                                        double const& y,
                                        double const& z)
  {

    if (x >= fFiducialCuboidBoundaryLower[0] &&
        x <= fFiducialCuboidBoundaryUpper[0] &&
        y >= fFiducialCuboidBoundaryLower[1] &&
        y <= fFiducialCuboidBoundaryUpper[1] &&
        z >= fFiducialCuboidBoundaryLower[2] &&
        z <= fFiducialCuboidBoundaryUpper[2])
    {
      return true;
    }
    else
    {
      return false;
    }

    return false;

  } // MCTools::isInsideFiducialCuboid_()


  //-----------------------------------------------------------------------
  double MCTools::distance_(double const& x0,
                            double const& y0,
                            double const& z0,
                            double const& x1,
                            double const& y1,
                            double const& z1)
  {

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0));

  } // MCTools::distance_()


  //-----------------------------------------------------------------------
  double MCTools::distance3D_(double const (& aXYZT)[4],
                              double const (& bXYZT)[4])
  {

    double const& x0 = aXYZT[0];
    double const& y0 = aXYZT[1];
    double const& z0 = aXYZT[2];

    double const& x1 = bXYZT[0];
    double const& y1 = bXYZT[1];
    double const& z1 = bXYZT[2];

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0));

  } // MCTools::distance3D_()


  //-----------------------------------------------------------------------
  double MCTools::distance4D_(double const (& aXYZT)[4],
                              double const (& bXYZT)[4])
  {

    double const& x0 = aXYZT[0];
    double const& y0 = aXYZT[1];
    double const& z0 = aXYZT[2];
    double const& t0 = aXYZT[3];

    double const& x1 = bXYZT[0];
    double const& y1 = bXYZT[1];
    double const& z1 = bXYZT[2];
    double const& t1 = bXYZT[3];

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0) +
                     (t1 - t0)*(t1 - t0));

  } // MCTools::distance4D_()


  //-----------------------------------------------------------------------
  double MCTools::angle3D_(double const (& aXYZT)[4],
                           double const (& bXYZT)[4])
  {

    double const& x0 = aXYZT[0];
    double const& y0 = aXYZT[1];
    double const& z0 = aXYZT[2];

    double const& x1 = bXYZT[0];
    double const& y1 = bXYZT[1];
    double const& z1 = bXYZT[2];

    double const aLength = std::sqrt(x0*x0 + y0*y0 + z0*z0);
    double const bLength = std::sqrt(x1*x1 + y1*y1 + z1*z1);
    double const dotProduct = x0*x1 + y0*y1 + z0*z1;

    double angle = std::acos(dotProduct / (aLength * bLength));

    if (fDebug)
    {
      mf::LogVerbatim("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "  aXYZT: (" << aXYZT[0] << ", " << aXYZT[1] << ", " << aXYZT[2] << ")\n"
        << "  bXYZT: (" << bXYZT[0] << ", " << bXYZT[1] << ", " << bXYZT[2] << ")\n"
        << "  aLength: " << aLength << "\n"
        << "  bLength: " << bLength << "\n"
        << "  dotProduct: " << dotProduct << "\n"
        << "  angle: " << angle << "\n"
        << "//////////////////////////////////////////////////////////////\n";
    }

    return angle;

  } // MCTools::angle3D_()


  //-----------------------------------------------------------------------
  void MCTools::getFourVectors_(simb::MCParticle const& particle,
                                double               (& startXYZT)[4],
                                double               (& endXYZT)[4],
                                double               (& startPE)[4],
                                double               (& endPE)[4])
  {

    // get the number of trajectory points of particle
    size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

    // get the index of the last trajectory point
    int last = numberTrajectoryPoints - 1;

    // get start/end positions and momenta of particle
    const TLorentzVector & positionStart = particle.Position(0);
    const TLorentzVector & positionEnd   = particle.Position(last);
    const TLorentzVector & momentumStart = particle.Momentum(0);
    const TLorentzVector & momentumEnd   = particle.Momentum(last);

    // fill start/end positions and momenta arrays
    positionStart.GetXYZT(startXYZT);
    positionEnd.GetXYZT(endXYZT);
    momentumStart.GetXYZT(startPE);
    momentumEnd.GetXYZT(endPE);

    return;

  } // MCTools::getFourVectors_()


  //-----------------------------------------------------------------------
  double MCTools::getTrajectoryLength_(simb::MCParticle const& particle)
  {

    // initialize trajectory length
    double length = 0;

    // get the number of trajectory points of particle
    size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

    // get the index of the last trajectory point
    size_t last = numberTrajectoryPoints - 1;

    for (size_t i = 0; i < last; ++i)
    {

      // initialize position arrays
      double aXYZT[4];
      double bXYZT[4];

      // get this trajectory point and the next one
      const TLorentzVector & positionA = particle.Position(i);
      const TLorentzVector & positionB = particle.Position(i+1);

      // fill position arrays
      positionA.GetXYZT(aXYZT);
      positionB.GetXYZT(bXYZT);

      // add segment length to total length
      length += this->distance3D_(aXYZT, bXYZT);

    } // loop over trajectory points

    return length;

  } // MCTools::getTrajectoryLength_()


  //-----------------------------------------------------------------------
  void MCTools::getShowerDaughters_(int                                      const& trackID,
                                    std::map< int, const simb::MCParticle* > const& particleMap,
                                    std::map< int, const simb::MCParticle* >      & showerParticleMap)
  {

    // create iterator for particle map
    std::map< int, const simb::MCParticle* >::const_iterator particleIterator
      = particleMap.find(trackID);

    // return if the end of the particle map is reached
    if (particleIterator == particleMap.end()) return;

    // pointer to MCParticle object
    const simb::MCParticle* particle = (*particleIterator).second;

    // get number of daughters of particle
    int numberDaughters = particle->NumberDaughters();

    // add pointer of MCParticle object to the shower particle map
    showerParticleMap[trackID] = particle;

    // recur
    for (int i = 0; i < numberDaughters; ++i)
    {
      int daughterTrackID = particle->Daughter(i);
      this->getShowerDaughters_(daughterTrackID, particleMap, showerParticleMap);
    }

  } // MCTools::getShowerDaughters_()


  //-----------------------------------------------------------------------
  double MCTools::distanceToTPCEdge_(double const (& aXYZT)[4],
                                     double const (& bXYZT)[4])
  {
    // initialize directional vector array, \vec{c} = \vec{a} - \vec{b}
    double cXYZT[4];

    // fill \vec{c}
    for (size_t i = 0; i < 4; ++i)
    {
      cXYZT[i] = aXYZT[i] - bXYZT[i];
    }

    // get spatial magnitude of \vec{c}
    double cMagnitude = this->distance3D_(aXYZT, bXYZT);

    // initialize unit vector \vec{v} in the direction of \vec{c}
    double vXYZ[3];

    // fill \vec{v}
    for (size_t i = 0; i < 3; ++i)
    {
      vXYZ[i] = cXYZT[i] / cMagnitude;
    }

    // parameter for parametric line
    std::vector< double > parameters;

    // get parameters
    for (size_t i = 0; i < 3; ++i) {

      // initialize parameters
      double parameterA;
      double parameterB;

      // attempt to compute parameter with fTPCBoundaryLower
      // assign negative value if computation fails
      try
      {
        parameterA = (fTPCBoundaryLower[i] - bXYZT[i]) / vXYZ[i];
      }
      catch (cet::exception &e)
      {
        parameterA = -1;
      }

      // attempt to compute parameter with fTPCBoundaryUpper
      // assign negative value if computation fails
      try
      {
        parameterB = (fTPCBoundaryUpper[i] - bXYZT[i]) / vXYZ[i];
      }
      catch (cet::exception &e)
      {
        parameterB = -1;
      }

      // add parameter to parameter vector if it is positive
      // to ensure that we are only looking at points in the
      // direction of the vector (starting from \vec{b})
      if (parameterA > 0)
        parameters.push_back(parameterA);
      if (parameterB > 0)
        parameters.push_back(parameterB);
    }

    // initialize four-vector array for the TPC exit point
    double exitXYZT[4] = { 0, 0, 0, 0 };

    // flag
    bool exitPointFound = false;

    // loop over parameters
    for (auto const& parameter : parameters)
    {
      double const x = bXYZT[0] + parameter * vXYZ[0];
      double const y = bXYZT[1] + parameter * vXYZ[1];
      double const z = bXYZT[2] + parameter * vXYZ[2];

      // if point is inside TPC
      if (this->isInsideTPC_(x, y, z))
      {
        exitXYZT[0] = x;
        exitXYZT[1] = y;
        exitXYZT[2] = z;
        exitPointFound = true;
        break;
      } // if inside TPC

    } // loop over parameters

    // return negative distance if the exit point has not
    // been found
    if (!exitPointFound)
      return -1;

    // get distance from \vec{a} to TPC edge
    double distance = this->distance3D_(exitXYZT, aXYZT);

    if (fDebug)
    {
      mf::LogVerbatim("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "  aXYZT: (" << aXYZT[0] << ", " << aXYZT[1] << ", " << aXYZT[2] << ", " << aXYZT[3] << ")\n"
        << "  bXYZT: (" << bXYZT[0] << ", " << bXYZT[1] << ", " << bXYZT[2] << ", " << bXYZT[3] << ")\n"
        << "  cXYZT: (" << cXYZT[0] << ", " << cXYZT[1] << ", " << cXYZT[2] << ", " << cXYZT[3] << ")\n"
        << "  cMagnitude: " << cMagnitude << "\n"
        << "  vXYZT: (" << vXYZ[0]  << ", " << vXYZ[1]  << ", " << vXYZ[2]  << ")\n"
        << "  exitXYZT: (" << exitXYZT[0] << ", " << exitXYZT[1] << ", " << exitXYZT[2] << ")\n"
        << "//////////////////////////////////////////////////////////////\n";
    }

    return distance;

  } // MCTools::distanceToTPCEdge_()


  //-----------------------------------------------------------------------
  void MCTools::distanceToTPCBoundary_(double const (& aXYZT)[4],
                                       double const (& bXYZT)[4],
                                       double        & distanceX,
                                       double        & distanceY,
                                       double        & distanceZ,
                                       double        & distance)
  {
    // initialize directional vector array, \vec{c} = \vec{a} - \vec{b}
    double cXYZT[4];

    // fill \vec{c}
    for (size_t i = 0; i < 4; ++i)
    {
      cXYZT[i] = aXYZT[i] - bXYZT[i];
    }

    // get spatial magnitude of \vec{c}
    double cMagnitude = this->distance3D_(aXYZT, bXYZT);

    // initialize unit vector \vec{v} in the direction of \vec{c}
    double vXYZ[3];

    // fill \vec{v}
    for (size_t i = 0; i < 3; ++i)
    {
      vXYZ[i] = cXYZT[i] / cMagnitude;
    }

    // parameter for parametric line
    std::vector< double > parameters;

    // get parameters
    for (size_t i = 0; i < 3; ++i) {

      // initialize parameters
      double parameterA;
      double parameterB;

      // attempt to compute parameter with fTPCBoundaryLower
      // assign negative value if computation fails
      try
      {
        parameterA = (fTPCBoundaryLower[i] - bXYZT[i]) / vXYZ[i];
      }
      catch (cet::exception &e)
      {
        parameterA = -1;
      }

      // attempt to compute parameter with fTPCBoundaryUpper
      // assign negative value if computation fails
      try
      {
        parameterB = (fTPCBoundaryUpper[i] - bXYZT[i]) / vXYZ[i];
      }
      catch (cet::exception &e)
      {
        parameterB = -1;
      }

      // add parameter to parameter vector if it is positive
      // to ensure that we are only looking at points in the
      // direction of the vector (starting from \vec{b})
      if (parameterA > 0)
        parameters.push_back(parameterA);
      if (parameterB > 0)
        parameters.push_back(parameterB);
    }

    // initialize four-vector array for the TPC exit point
    double exitXYZT[4] = { 0, 0, 0, 0 };

    // flag
    bool exitPointFound = false;

    // loop over parameters
    for (auto const& parameter : parameters)
    {
      double const x = bXYZT[0] + parameter * vXYZ[0];
      double const y = bXYZT[1] + parameter * vXYZ[1];
      double const z = bXYZT[2] + parameter * vXYZ[2];

      // if point is inside TPC
      if (this->isInsideTPC_(x, y, z))
      {
        exitXYZT[0] = x;
        exitXYZT[1] = y;
        exitXYZT[2] = z;
        exitPointFound = true;
        break;
      } // if inside TPC

    } // loop over parameters

    // return negative distance if the exit point has not
    // been found
    if (!exitPointFound)
    {
      distanceX = -1;
      distanceY = -1;
      distanceZ = -1;
      distance = -1;
      return;
    }

    // get 1D distance to TPC boundary
    distanceX = std::abs(exitXYZT[0] - aXYZT[0]);
    distanceY = std::abs(exitXYZT[1] - aXYZT[1]);
    distanceZ = std::abs(exitXYZT[2] - aXYZT[2]);

    // get distance from \vec{a} to TPC boundary
    distance = this->distance3D_(exitXYZT, aXYZT);

    if (fDebug)
    {
      mf::LogVerbatim("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "  aXYZT: (" << aXYZT[0] << ", " << aXYZT[1] << ", " << aXYZT[2] << ", " << aXYZT[3] << ")\n"
        << "  bXYZT: (" << bXYZT[0] << ", " << bXYZT[1] << ", " << bXYZT[2] << ", " << bXYZT[3] << ")\n"
        << "  cXYZT: (" << cXYZT[0] << ", " << cXYZT[1] << ", " << cXYZT[2] << ", " << cXYZT[3] << ")\n"
        << "  cMagnitude: " << cMagnitude << "\n"
        << "  vXYZT: (" << vXYZ[0]  << ", " << vXYZ[1]  << ", " << vXYZ[2]  << ")\n"
        << "  exitXYZT: (" << exitXYZT[0] << ", " << exitXYZT[1] << ", " << exitXYZT[2] << ")\n"
        << "  distanceX: " << distanceX << "\n"
        << "  distanceY: " << distanceY << "\n"
        << "  distanceZ: " << distanceZ << "\n"
        << "  distance: "  << distance  << "\n"
        << "//////////////////////////////////////////////////////////////\n";
    }

    return;

  } // MCTools::distanceToTPCBoundary_()


  //-----------------------------------------------------------------------
  void MCTools::longitudinalShowerProfile_(std::map< int, const simb::MCParticle* >           const& showerParticleMap,
                                           art::ValidHandle< std::vector< sim::SimChannel > > const& simChannelHandle,
                                           double                                             const& dz,
                                           TVector3                                           const& startPoint,
                                           TVector3                                           const& longitudinalDirection,
                                           std::vector< double >                                   & longitudinalProfile)
  {
    longitudinalProfile.clear();

    // rotatation matrix for aligning with longitudinal direction
    TRotation m;
    m.RotateZ(longitudinalDirection.Phi());
    m.RotateX(longitudinalDirection.Theta());

    // rotate to align with longitudinal direction
    TVector3 startPoint_            = m * startPoint;
    TVector3 longitudinalDirection_ = m * longitudinalDirection;

    for (auto const& channel : (*simChannelHandle)) // sim::SimChannel
    {
      // get the numeric ID associated with this channel
      auto channelNumber = channel.Channel(); // raw::ChannelID_t

      // only look at the energy on the collection plane
      if (fGeometry->SignalType(channelNumber) != geo::kCollection)
        continue;

      // each channel has a map inside it that connects
      // a time slice to energy deposits in the detector
      auto const& timeSlices = channel.TDCIDEMap(); // std::map< unsigned short, std::vector< sim::IDE > >

      // for every time slice in this channel:
      for (auto const& timeSlice : timeSlices)
      {
        // get energy deposits from time slice
        auto const& energyDeposits = timeSlice.second; // std::vector< sim::IDE >

        // loop through energy deposits
        for (auto const& energyDeposit : energyDeposits) // sim::IDE
        {
          if (showerParticleMap.find(energyDeposit.trackID) != showerParticleMap.end() &&
              energyDeposit.trackID != sim::NoParticleId)
          {

            // position of energy deposit
            TVector3 point(energyDeposit.x, energyDeposit.y, energyDeposit.z);

            // rotate to align with longitudinal direction
            TVector3 point_ = m * point;

            // get distance between plane and energy deposit
            // plane:
            //   - normal vector is parallel to the longitudinal direction
            //   - intersects with the start point
            double distance = point_.Z() - startPoint_.Z();

            if (distance < 0) continue;

            // get bin where the energy deposit goes into
            const unsigned int bin = (unsigned int) (distance / dz);

//std::cout << " distance: " << distance << std::endl;
//std::cout << " dz:       " << dz       << std::endl;
//std::cout << " bin:      " << bin      << std::endl;

            // is the profile array large enough?
            if (longitudinalProfile.size() < bin+1)
            {
              // increase size of profile array and pad it with zeros
              longitudinalProfile.resize(bin+1, 0.);
            }

            // add energy deposit to that bin
            longitudinalProfile[bin] += energyDeposit.energy;

          } // if energy deposited by particle in particle map
        } // for each energy deposit
      } // for each time slice
    } // for each SimChannel
  } // MCTools::longitudinalShowerProfile_()


  //-----------------------------------------------------------------------
  void MCTools::radialShowerProfile_(std::map< int, const simb::MCParticle* >           const& showerParticleMap,
                                     art::ValidHandle< std::vector< sim::SimChannel > > const& simChannelHandle,
                                     double                                             const& dr,
                                     TVector3                                           const& startPoint,
                                     std::vector< double >                                   & radialProfile)
  {
    radialProfile.clear();

    for (auto const& channel : (*simChannelHandle)) // sim::SimChannel
    {
      // get the numeric ID associated with this channel
      auto channelNumber = channel.Channel(); // raw::ChannelID_t

      // only look at the energy on the collection plane
      if (fGeometry->SignalType(channelNumber) != geo::kCollection)
        continue;

      // each channel has a map inside it that connects
      // a time slice to energy deposits in the detector
      auto const& timeSlices = channel.TDCIDEMap(); // std::map< unsigned short, std::vector< sim::IDE > >

      // for every time slice in this channel:
      for (auto const& timeSlice : timeSlices)
      {
        // get energy deposits from time slice
        auto const& energyDeposits = timeSlice.second; // std::vector< sim::IDE >

        // loop through energy deposits
        for (auto const& energyDeposit : energyDeposits) // sim::IDE
        {
          if (showerParticleMap.find(energyDeposit.trackID) != showerParticleMap.end() &&
              energyDeposit.trackID != sim::NoParticleId)
          {

            // get distance between start point of shower and energy deposit
            double distance = this->distance_(startPoint.X(),
                                              startPoint.Y(),
                                              startPoint.z(),
                                              energyDeposit.x,
                                              energyDeposit.y,
                                              energyDeposit.z);

            // get bin where the energy deposit goes into
            const unsigned int bin = (unsigned int) (distance / dr);

            // is the profile array large enough?
            if (radialProfile.size() < bin+1)
            {
              // increase size of profile array and pad it with zeros
              radialProfile.resize(bin+1, 0.);
            }

            // add energy deposit to that bin
            radialProfile[bin] += energyDeposit.energy;

          } // if energy deposited by particle in particle map
        } // for each energy deposit
      } // for each time slice
    } // for each SimChannel

  } // MCTools::radialShowerProfile_()


  //-----------------------------------------------------------------------
  void MCTools::getShowerVariables_(int                                                const& trackID,
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
                                    double                                                  & photonSphereRadius)
  {
    // reset variables
    photonConvert = false;
    photonConversionLength  = fDetectorDiagonal;
    photonAngleZ            = 0;
    photonEnergy            = 0;
    photonEnergyDeposited   = 0;
    photonTPCContainment    = 0;
    photonTPCDistanceX      = fDetectorDiagonal;
    photonTPCDistanceY      = fDetectorDiagonal;
    photonTPCDistanceZ      = fDetectorDiagonal;
    photonTPCDistance       = fDetectorDiagonal;
    photonSphereEnergy      = 0;
    photonSphereContainment = 0;
    photonSphereRadius      = fDetectorDiagonal;

    // get particle
    const simb::MCParticle* particle = particleMap.at(trackID);

    // get start/end 4-positions and 4-momenta
    this->getFourVectors_(*particle,
                          startXYZT,
                          endXYZT,
                          startPE,
                          endPE);

    // get angle between photon and z-axis
    double zUnitVector[4] = { 0, 0, 1, 0 };
    photonAngleZ = this->angle3D_(startPE, zUnitVector);

    // get photon energy
    photonEnergy = startPE[3] * 1000.0; // MeV

    // check to see if photon converts inside TPC
    photonConvert = this->isInsideTPC_(endXYZT[0],
                                       endXYZT[1],
                                       endXYZT[2]);

    if (photonConvert)
    {
      // get photon conversion length
      photonConversionLength = this->distance3D_(startXYZT, endXYZT);

      // get distance from shower start point to TPC edge along
      // the forward direction of the shower
      //photonTPCDistance = this->distanceToTPCEdge_(endXYZT, startXYZT);
      this->distanceToTPCBoundary_(endXYZT,
                                   startXYZT,
                                   photonTPCDistanceX,
                                   photonTPCDistanceY,
                                   photonTPCDistanceZ,
                                   photonTPCDistance);
    }

    // fill shower particle map
    this->getShowerDaughters_(trackID, particleMap, showerParticleMap);

    // warn if there are less than 2 particles in the shower particle map
    if (showerParticleMap.size() < 2) 
    {
      mf::LogWarning("MCTools")
        << "//////////////////////////////////////////////////////////////\n"
        << "Size of showerParticleMap is " << showerParticleMap.size() << ".\n"
        << "KeepEMShowerDaughters may have been set to false!\n"
        << "//////////////////////////////////////////////////////////////";
    }

    for (auto const& channel : (*simChannelHandle)) // sim::SimChannel
    {
      // get the numeric ID associated with this channel
      auto channelNumber = channel.Channel(); // raw::ChannelID_t

      // only look at the energy on the collection plane
      if (fGeometry->SignalType(channelNumber) != geo::kCollection)
        continue;

      // each channel has a map inside it that connects
      // a time slice to energy deposits in the detector
      auto const& timeSlices = channel.TDCIDEMap(); // std::map< unsigned short, std::vector< sim::IDE > >

      // for every time slice in this channel:
      for (auto const& timeSlice : timeSlices)
      {
        // get energy deposits from time slice
        auto const& energyDeposits = timeSlice.second; // std::vector< sim::IDE >

        // loop through energy deposits
        for (auto const& energyDeposit : energyDeposits) // sim::IDE
        {
          if (showerParticleMap.find(energyDeposit.trackID) != showerParticleMap.end() &&
              energyDeposit.trackID != sim::NoParticleId)
          {
            // add energy to energy deposited by photon-initiated shower
            photonEnergyDeposited += energyDeposit.energy;

            if (photonConvert)
            {
              // get minimum distance to TPC edge
              // if distance is close enough, compute the distance to shower start
              // compare with current minimum value

              std::vector<double> distanceToBoundary1D;

              distanceToBoundary1D.push_back(std::abs(energyDeposit.x - fTPCBoundaryLower[0]));
              distanceToBoundary1D.push_back(std::abs(energyDeposit.x - fTPCBoundaryUpper[0]));
              distanceToBoundary1D.push_back(std::abs(energyDeposit.y - fTPCBoundaryLower[1]));
              distanceToBoundary1D.push_back(std::abs(energyDeposit.y - fTPCBoundaryUpper[1]));
              distanceToBoundary1D.push_back(std::abs(energyDeposit.z - fTPCBoundaryLower[2]));
              distanceToBoundary1D.push_back(std::abs(energyDeposit.z - fTPCBoundaryUpper[2]));

              double minDistanceToBoundary1D = *std::min_element(distanceToBoundary1D.begin(),
                                                                 distanceToBoundary1D.end());

              if (minDistanceToBoundary1D < fMaxDistanceToBoundary1D)
              {
                // array for energy deposit position
                double energyDepositXYZT[4] = { energyDeposit.x,
                                                energyDeposit.y,
                                                energyDeposit.z,
                                                0 };

                // get distance from energy deposit to shower start point
                double distanceToShowerStart = this->distance3D_(energyDepositXYZT,
                                                                 endXYZT);

                // get the minimum radius along TPC boundary
                if (distanceToShowerStart < photonSphereRadius)
                  photonSphereRadius = distanceToShowerStart;
              } // if minimum distance to boundary is within range
            } // if photon converts inside TPC

          } // if energy deposited by particle in particle map
        } // for each energy deposit
      } // for each time slice
    } // for each SimChannel

    if (photonConvert)
    {
      for (auto const& channel : (*simChannelHandle)) // sim::SimChannel
      {
        // get the numeric ID associated with this channel
        auto channelNumber = channel.Channel(); // raw::ChannelID_t

        // only look at the energy on the collection plane
        if (fGeometry->SignalType(channelNumber) != geo::kCollection)
          continue;

        // each channel has a map inside it that connects
        // a time slice to energy deposits in the detector
        auto const& timeSlices = channel.TDCIDEMap(); // std::map< unsigned short, std::vector< sim::IDE > >

        // for every time slice in this channel:
        for (auto const& timeSlice : timeSlices)
        {
          // get energy deposits from time slice
          auto const& energyDeposits = timeSlice.second; // std::vector< sim::IDE >

          // loop through energy deposits
          for (auto const& energyDeposit : energyDeposits) // sim::IDE
          {
            if (showerParticleMap.find(energyDeposit.trackID) != showerParticleMap.end() &&
                energyDeposit.trackID != sim::NoParticleId)
            {

              // array for energy deposit position
              double energyDepositXYZT[4] = { energyDeposit.x,
                                              energyDeposit.y,
                                              energyDeposit.z,
                                              0 };

              // get distance from energy deposit to shower start point
              double distanceToShowerStart = this->distance3D_(energyDepositXYZT,
                                                               endXYZT);

              // add energy deposit to sphere energy if it is within
              // radius if sphere
              if (distanceToShowerStart < photonSphereRadius)
                photonSphereEnergy += energyDeposit.energy;

            } // if energy deposited by particle in particle map
          } // for each energy deposit
        } // for each time slice
      } // for each SimChannel
    } // if photon converts inside TPC

    photonTPCContainment = photonEnergyDeposited / photonEnergy;
    photonSphereContainment = photonSphereEnergy / photonEnergy;

  } // MCTools::getShowerVariables_()


} // namespace EMShowerContainment
