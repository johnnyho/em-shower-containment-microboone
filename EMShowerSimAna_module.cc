//////////////////////////////////////////////////////////////
// Name:      EMShowerSimAna
// Date:      27 April 2016
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

// EMShowerSimAna_module.cc

#ifndef EMShowerSimAna_Module
#define EMShowerSimAna_Module

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"

// Framework includes
#include "art/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "lardata/MCBase/MCTrack.h"
#include "lardata/MCBase/MCShower.h"
#include "lardata/MCBase/MCStep.h"
#include "larsim/Simulation/sim.h"
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

// uboonecode includes
#include "uboone/EMShowerContainment/MCTools.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ includes
#include <algorithm>
#include <cmath>
#include <map>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>


namespace EMShowerSimAna {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class EMShowerSimAna : public art::EDAnalyzer
  {
   public:

    // standard constructor and destructor
    explicit EMShowerSimAna(fhicl::ParameterSet const& parameterSet);
    virtual ~EMShowerSimAna();

    // this method is called once, at the start of the job
    virtual void beginJob() override;

    // this method is called once, at the start of each run
    virtual void beginRun(const art::Run& run) override;

    // this method reads in any parameters from the .fcl files
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;

    // the analysis routine, called once per event
    virtual void analyze (const art::Event& event) override;

   private:

    // parameters read from the .fcl file
    std::string fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    std::string fMCTrackProducerLabel;    ///< The name of the producer that created MCTracks
    std::string fMCShowerProducerLabel;   ///< The name of the producer that created MCShowers

    // MCTools
    EMShowerContainment::MCTools fMCTools;

    // pointers to n-tuple
    TTree * fSimulationTree;

    // variables that will go into the n-tuple
    int fEvent;     ///< number of the event being processed
    int fRun;       ///< number of the run being processed
    int fSubRun;    ///< number of the sub-run being processed

    // PDG code and track ID for primary particle
    int fPrimaryPDGCode;
    int fPrimaryTrackID;

    // Arrays for 4-vectors: (x, y, z, t) and (Px, Py, Pz, E).
    double fStartXYZT[4]; ///< (x, y, z, t) of the true start of the particle
    double fEndXYZT[4];   ///< (x, y, z, t) of the true end of the particle
    double fStartPE[4];   ///< (Px, Py, Pz, E) at the true start of the particle
    double fEndPE[4];     ///< (Px, Py, Pz, E) at the true end of the particle

    // shower variables
    bool   fPhotonConvert;           // flag, true of photon converted inside TPC, false otherwise
    double fPhotonConversionLength;  // photon conversion length
    double fPhotonAngleZ;            // initial angle of photon w.r.t. z-axis
    double fPhotonEnergy;            // initial energy of photon
    double fPhotonEnergyDeposited;   // total energy deposited by photon-initiated shower
    double fEnergyDeposited;         // total energy deposited on TPC wires
    double fPhotonTPCContainment;    // fraction of energy contained within the TPC
    double fPhotonTPCDistance;       // distance from start point of shower to TPC boundary along the direction of the shower
    double fPhotonTPCDistanceX;      // x-component of PhotonTPCDistance
    double fPhotonTPCDistanceY;      // y-component of PhotonTPCDistance
    double fPhotonTPCDistanceZ;      // z-component of PhotonTPCDistance
    double fPhotonSphereEnergy;      // energy deposited by photon-initiated shower within a sphere of radius PhotonSphereRadius centered at start point of shower
    double fPhotonSphereContainment; // fraction energy contained by photon-initiated shower within sphere
    double fPhotonSphereRadius;      // radius of sphere centered at start point of shower

    // dE/dx of EM shower
    double fdEdx;

    // bin sizes of longitudinal and radial shower profiles
    double fLongitudinalBinSize;
    double fRadialBinSize;

    // number of longitudinal and radial shower profile bins
    int fLongitudinalNumberBins;
    int fRadialNumberBins;

    // longitudinal and radial shower profiles
    std::vector< double > fLongitudinalShowerProfile;
    std::vector< double > fRadialShowerProfile;

    // pointer to geometry provider
    geo::GeometryCore const* fGeometry;

  }; // class EMShowerSimAna


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  EMShowerSimAna::EMShowerSimAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
    , fMCTools(parameterSet.get<fhicl::ParameterSet>("MCTools"))
  {
    // get a pointer to the geometry service provider
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

    // read in parameters from the .fcl file
    reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------
  // destructor
  EMShowerSimAna::~EMShowerSimAna() {}


  //-----------------------------------------------------------------------
  void EMShowerSimAna::beginJob()
  {
    // access art's TFileService
    art::ServiceHandle<art::TFileService> tfs;

    fSimulationTree = tfs->make<TTree>("SimulationTree", "SimulationTree");

    fSimulationTree->Branch("Event",     &fEvent,    "Event/I");
    fSimulationTree->Branch("SubRun",    &fSubRun,   "SubRun/I");
    fSimulationTree->Branch("Run",       &fRun,      "Run/I");

    fSimulationTree->Branch("StartXYZT", fStartXYZT, "StartXYZT[4]/D");
    fSimulationTree->Branch("EndXYZT",   fEndXYZT,   "EndXYZT[4]/D");
    fSimulationTree->Branch("StartPE",   fStartPE,   "StartPE[4]/D");
    fSimulationTree->Branch("EndPE",     fEndPE,     "EndPE[4]/D");

    fSimulationTree->Branch("PrimaryPDGCode", &fPrimaryPDGCode, "PrimaryPDGCode/I");
    fSimulationTree->Branch("PrimaryTrackID", &fPrimaryTrackID, "PrimaryTrackID/I");

    fSimulationTree->Branch("PhotonAngleZ",            &fPhotonAngleZ,            "PhotonAngleZ/D");
    fSimulationTree->Branch("PhotonConvert",           &fPhotonConvert,           "PhotonConvert/O");
    fSimulationTree->Branch("PhotonConversionLength",  &fPhotonConversionLength,  "PhotonConversionLength/D");
    fSimulationTree->Branch("PhotonEnergy",            &fPhotonEnergy,            "PhotonEnergy/D");
    fSimulationTree->Branch("PhotonEnergyDeposited",   &fPhotonEnergyDeposited,   "PhotonEnergyDeposited/D");
    fSimulationTree->Branch("PhotonTPCContainment",    &fPhotonTPCContainment,    "PhotonTPCContainment/D");
    fSimulationTree->Branch("PhotonTPCDistanceX",      &fPhotonTPCDistanceX,      "PhotonTPCDistanceX/D");
    fSimulationTree->Branch("PhotonTPCDistanceY",      &fPhotonTPCDistanceY,      "PhotonTPCDistanceY/D");
    fSimulationTree->Branch("PhotonTPCDistanceZ",      &fPhotonTPCDistanceZ,      "PhotonTPCDistanceZ/D");
    fSimulationTree->Branch("PhotonTPCDistance",       &fPhotonTPCDistance,       "PhotonTPCDistance/D");
    fSimulationTree->Branch("PhotonSphereEnergy",      &fPhotonSphereEnergy,      "PhotonSphereEnergy/D");
    fSimulationTree->Branch("PhotonSphereContainment", &fPhotonSphereContainment, "PhotonSphereContainment/D");
    fSimulationTree->Branch("PhotonSphereRadius",      &fPhotonSphereRadius,      "PhotonSphereRadius/D");

    fSimulationTree->Branch("dEdx", &fdEdx, "dEdx/D");

    fSimulationTree->Branch("LongitudinalBinSize",    &fLongitudinalBinSize,    "LongitudinalBinSize/D");
    fSimulationTree->Branch("RadialBinSize",          &fRadialBinSize,          "RadialBinSize/D");
    fSimulationTree->Branch("LongitudinalNumberBins", &fLongitudinalNumberBins, "LongitudinalNumberBins/I");
    fSimulationTree->Branch("RadialNumberBins",       &fRadialNumberBins,       "RadialNumberBins/I");

    fSimulationTree->Branch("LongitudinalShowerProfile", &fLongitudinalShowerProfile);
    fSimulationTree->Branch("RadialShowerProfile",       &fRadialShowerProfile);
  }


  //-----------------------------------------------------------------------
  void EMShowerSimAna::beginRun(const art::Run& /*run*/)
  {}


  //-----------------------------------------------------------------------
  void EMShowerSimAna::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // read in parameters from the .fcl file
    fSimulationProducerLabel = parameterSet.get< std::string >("SimulationLabel",     "largeant");
    fMCTrackProducerLabel    = parameterSet.get< std::string >("MCTrackLabel",        "mcreco");
    fMCShowerProducerLabel   = parameterSet.get< std::string >("MCShowerLabel",       "mcreco");
    fLongitudinalBinSize     = parameterSet.get< double      >("LongitudinalBinSize", 1.0);
    fRadialBinSize           = parameterSet.get< double      >("RadialBinSize",       1.0);
  }


  //-----------------------------------------------------------------------
  void EMShowerSimAna::analyze(const art::Event& event) 
  {
    // get event, run, and subrun numbers
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    // get all the simulated particles for the event
    art::ValidHandle< std::vector< simb::MCParticle > > particleHandle = event.getValidHandle< std::vector< simb::MCParticle > >(fSimulationProducerLabel);

    // get all the simulated channels for the event
    art::ValidHandle< std::vector< sim::SimChannel > > simChannelHandle = event.getValidHandle< std::vector< sim::SimChannel > >(fSimulationProducerLabel);

    // read in the MCTracks
    art::ValidHandle< std::vector< sim::MCTrack > > mcTrackHandle = event.getValidHandle< std::vector< sim::MCTrack > >(fMCTrackProducerLabel);

    // read in the MCShowers
    art::ValidHandle< std::vector< sim::MCShower > > mcShowerHandle = event.getValidHandle< std::vector< sim::MCShower > >(fMCShowerProducerLabel);

    // map of pointers to MCParticle objects
    std::map< int, const simb::MCParticle* > particleMap;

    // map of pointers to shower daughters particles
    std::map< int, const simb::MCParticle* > showerParticleMap;

    // map of pointers to MCTrack objects
    std::map< int, const sim::MCTrack* > mcTrackMap;

    // map of pointers to MCShower objects
    std::map< int, const sim::MCShower* > mcShowerMap;

    // flags for checking if the primary particle is a photon that converts in the TPC
    bool primaryPhoton              = false;
    bool primaryPhotonConvertsInTPC = false;

    // reset variables for new event
    fPhotonEnergy          = -1;
    fEnergyDeposited       = -1;
    fPhotonTPCDistance = -1;

    fPhotonConvert           = false;
    fPhotonConversionLength  = fMCTools.detectorDiagonal_();
    fPhotonAngleZ            = 0;
    fPhotonEnergy            = 0;
    fPhotonEnergyDeposited   = 0;
    fPhotonTPCContainment    = 0;
    fPhotonTPCDistanceX      = fMCTools.detectorDiagonal_();
    fPhotonTPCDistanceY      = fMCTools.detectorDiagonal_();
    fPhotonTPCDistanceZ      = fMCTools.detectorDiagonal_();
    fPhotonTPCDistance       = fMCTools.detectorDiagonal_();
    fPhotonSphereEnergy      = 0;
    fPhotonSphereContainment = 0;
    fPhotonSphereRadius      = fMCTools.detectorDiagonal_();

    fdEdx = 0;

    fLongitudinalNumberBins = 0;
    fRadialNumberBins       = 0;

    fLongitudinalShowerProfile.clear();
    fRadialShowerProfile.clear();

    // loop through MCTrack objects and fill MCTrack map
    for (auto const& mcTrack : *mcTrackHandle) {
      if (!mcTrack.size()) continue;
      int trackID = (int) mcTrack.TrackID();
      mcTrackMap[trackID] = &mcTrack; 
    }

    // loop through MCShower objects and fill MCShower map
    for (auto const& mcShower : *mcShowerHandle) {
      int trackID = (int) mcShower.TrackID();
      mcShowerMap[trackID] = &mcShower; 
    }

    // loop through MCParticle objects
    for (auto const& particle : (*particleHandle)) // simb::MCParticle
    {
      // get track ID of particle
      int trackID = particle.TrackId();

      // get PDG code of particle
      int pdgCode = particle.PdgCode();

      // add the address of the MCParticle to the map, with the track ID as the key
      particleMap[trackID] = &particle; 

      // if particle is primary and photon
      if (particle.Process() == "primary" && std::abs(pdgCode) == 22)
      {
        // change flag
        primaryPhoton = true;

        // get primary track ID
        fPrimaryTrackID = trackID;

        // get primary PDG code
        fPrimaryPDGCode = pdgCode;

        // get primary start and end 4-positions and 4-momenta
        fMCTools.getFourVectors_(particle, fStartXYZT, fEndXYZT, fStartPE, fEndPE);

        // check to see if primary photon converts inside TPC
        primaryPhotonConvertsInTPC = fMCTools.isInsideTPC_(fEndXYZT[0], fEndXYZT[1], fEndXYZT[2]);

        // get 3D distance between the start and end positions of the primary
        const double trajectoryLength = fMCTools.getTrajectoryLength_(particle);

        LOG_DEBUG("EMShowerSimAna") << "Trajectory length: " << trajectoryLength << " cm";

        fPhotonEnergy = fStartPE[3] * 1000.0; // MeV

        // return if primary photon does not convert inside TPC
        if (!primaryPhotonConvertsInTPC) return;

      } // if primary and photon
    } // loop over all particles in the event

    // get shower daughters of primary photon and fill shower particle map
    fMCTools.getShowerDaughters_(fPrimaryTrackID, particleMap, showerParticleMap);

    // get various variables for shower
    fMCTools.getShowerVariables_(fPrimaryTrackID,
                                 particleMap,
                                 simChannelHandle,
                                 showerParticleMap,
                                 fStartXYZT,
                                 fEndXYZT,
                                 fStartPE,
                                 fEndPE,
                                 fPhotonAngleZ,
                                 fPhotonConvert,
                                 fPhotonConversionLength,
                                 fPhotonEnergy,
                                 fPhotonEnergyDeposited,
                                 fPhotonTPCContainment,
                                 fPhotonTPCDistanceX,
                                 fPhotonTPCDistanceY,
                                 fPhotonTPCDistanceZ,
                                 fPhotonTPCDistance,
                                 fPhotonSphereEnergy,
                                 fPhotonSphereContainment,
                                 fPhotonSphereRadius);

    // set start point and direction of shower TVector3 objects
    TVector3 startPoint(fEndXYZT[0], fEndXYZT[1], fEndXYZT[2]);
    TVector3 direction(fEndXYZT[0] - fStartXYZT[0], fEndXYZT[1] - fStartXYZT[1], fEndXYZT[2] - fStartXYZT[2]);

    // get radial shower profile
    fMCTools.radialShowerProfile_(showerParticleMap, simChannelHandle, fRadialBinSize, startPoint, fRadialShowerProfile);
    fRadialNumberBins = fRadialShowerProfile.size();

    // get longitudinal shower profile
    fMCTools.longitudinalShowerProfile_(showerParticleMap, simChannelHandle, fLongitudinalBinSize, startPoint, direction, fLongitudinalShowerProfile);
    fLongitudinalNumberBins = fLongitudinalShowerProfile.size();

    // get MCShower associated with primary photon
    if (mcShowerMap.count(fPrimaryTrackID))
    {
      const sim::MCShower * mcShower = mcShowerMap[fPrimaryTrackID];

      // get dE/dx of EM shower
      fdEdx = mcShower->dEdx();
    }

    // look at all the energy deposited on the TPC wires
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

          if (particleMap.find(energyDeposit.trackID) != particleMap.end() &&
              energyDeposit.trackID != sim::NoParticleId)
          {
            // add to total energy deposited on TPC wire
            fEnergyDeposited += energyDeposit.energy;
          } // if energy deposited by particle in particle map

        } // for each energy deposit
      } // for each time slice
    } // for each SimChannel

    fSimulationTree->Fill();

  } // EMShowerSimAna::analyze()


  DEFINE_ART_MODULE(EMShowerSimAna)

} // namespace EMShowerSimAna

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // EMShowerSimAna_Module
