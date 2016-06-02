//////////////////////////////////////////////////////////////
// Name:      NeutralPionSimAna
// Date:      19 October 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

// NeutralPionSimAna_module.cc

#ifndef NeutralPionSimAna_Module
#define NeutralPionSimAna_Module

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

// C++ Includes
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>


namespace NeutralPionSimAna {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class NeutralPionSimAna : public art::EDAnalyzer
  {
   public:

    // standard constructor and destructor
    explicit NeutralPionSimAna(fhicl::ParameterSet const& parameterSet);
    virtual ~NeutralPionSimAna();

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

    // total energy deposited on TPC wires
    double fEnergyDeposited;

    // Arrays for 4-vectors: (x, y, z, t) and (Px, Py, Pz, E).
    double fStartXYZT[4]; ///< (x, y, z, t) of the true start of the particle
    double fEndXYZT[4];   ///< (x, y, z, t) of the true end of the particle
    double fStartPE[4];   ///< (Px, Py, Pz, E) at the true start of the particle
    double fEndPE[4];     ///< (Px, Py, Pz, E) at the true end of the particle

    // counters!
    int fNumberNeutralPionsInEvent; // counter for neutral pions in event
    int fNumberPhotonConversions;   // counter for photons from neutral pion decay that convert inside TPC

    double fPhotonAStartXYZT[4];
    double fPhotonAEndXYZT[4];
    double fPhotonAStartPE[4];
    double fPhotonAEndPE[4];

    double fPhotonBStartXYZT[4];
    double fPhotonBEndXYZT[4];
    double fPhotonBStartPE[4];
    double fPhotonBEndPE[4];

    double fPhotonAAngleZ;
    bool   fPhotonAConvert;
    double fPhotonAConversionLength;
    double fPhotonAEnergy;
    double fPhotonAEnergyDeposited;
    double fPhotonATPCContainment;
    double fPhotonATPCDistanceX;
    double fPhotonATPCDistanceY;
    double fPhotonATPCDistanceZ;
    double fPhotonATPCDistance;
    double fPhotonASphereEnergy;
    double fPhotonASphereContainment;
    double fPhotonASphereRadius;

    double fPhotonBAngleZ;
    bool   fPhotonBConvert;
    double fPhotonBConversionLength;
    double fPhotonBEnergy;
    double fPhotonBEnergyDeposited;
    double fPhotonBTPCContainment;
    double fPhotonBTPCDistanceX;
    double fPhotonBTPCDistanceY;
    double fPhotonBTPCDistanceZ;
    double fPhotonBTPCDistance;
    double fPhotonBSphereEnergy;
    double fPhotonBSphereContainment;
    double fPhotonBSphereRadius;

    double fPhotonAngle; // opening angle between the 2 photons from the neutral pion decay

    std::vector< double > fPhotonEnergy;
    std::vector< double > fNeutralPionEnergy;
    std::vector< double > fNeutralPionMomentum;

    // pointer to geometry provider
    geo::GeometryCore const* fGeometry;

  }; // class NeutralPionSimAna


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  NeutralPionSimAna::NeutralPionSimAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
    , fMCTools(parameterSet.get<fhicl::ParameterSet>("MCTools"))
  {
    // get a pointer to the geometry service provider
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

    // read in parameters from the .fcl file
    reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------
  // Destructor
  NeutralPionSimAna::~NeutralPionSimAna() {}


  //-----------------------------------------------------------------------
  void NeutralPionSimAna::beginJob()
  {
    // access art's TFileService
    art::ServiceHandle<art::TFileService> tfs;

    fSimulationTree = tfs->make<TTree>("SimulationTree", "SimulationTree");

    fSimulationTree->Branch("Event",       &fEvent,          "Event/I");
    fSimulationTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
    fSimulationTree->Branch("Run",         &fRun,            "Run/I");

    fSimulationTree->Branch("PhotonEnergy",        &fPhotonEnergy);
    fSimulationTree->Branch("NeutralPionEnergy",   &fNeutralPionEnergy);
    fSimulationTree->Branch("NeutralPionMomentum", &fNeutralPionMomentum);

    fSimulationTree->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fSimulationTree->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fSimulationTree->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fSimulationTree->Branch("EndPE",       fEndPE,           "EndPE[4]/D");

    fSimulationTree->Branch("PhotonAStartXYZT", fPhotonAStartXYZT, "PhotonAStartXYZT[4]/D");
    fSimulationTree->Branch("PhotonAEndXYZT",   fPhotonAEndXYZT,   "PhotonAEndXYZT[4]/D");
    fSimulationTree->Branch("PhotonAStartPE",   fPhotonAStartPE,   "PhotonAStartPE[4]/D");
    fSimulationTree->Branch("PhotonAEndPE",     fPhotonAEndPE,     "PhotonAEndPE[4]/D");

    fSimulationTree->Branch("PhotonBStartXYZT", fPhotonBStartXYZT, "PhotonBStartXYZT[4]/D");
    fSimulationTree->Branch("PhotonBEndXYZT",   fPhotonBEndXYZT,   "PhotonBEndXYZT[4]/D");
    fSimulationTree->Branch("PhotonBStartPE",   fPhotonBStartPE,   "PhotonBStartPE[4]/D");
    fSimulationTree->Branch("PhotonBEndPE",     fPhotonBEndPE,     "PhotonBEndPE[4]/D");

    fSimulationTree->Branch("PhotonAngle",              &fPhotonAngle,              "PhotonAngle/D");

    fSimulationTree->Branch("PhotonAAngleZ",            &fPhotonAAngleZ,            "PhotonAAngleZ/D");
    fSimulationTree->Branch("PhotonAConvert",           &fPhotonAConvert,           "PhotonAConvert/O");
    fSimulationTree->Branch("PhotonAConversionLength",  &fPhotonAConversionLength,  "PhotonAConversionLength/D");
    fSimulationTree->Branch("PhotonAEnergy",            &fPhotonAEnergy,            "PhotonAEnergy/D");
    fSimulationTree->Branch("PhotonAEnergyDeposited",   &fPhotonAEnergyDeposited,   "PhotonAEnergyDeposited/D");
    fSimulationTree->Branch("PhotonATPCContainment",    &fPhotonATPCContainment,    "PhotonATPCContainment/D");
    fSimulationTree->Branch("PhotonATPCDistanceX",      &fPhotonATPCDistanceX,      "PhotonATPCDistanceX/D");
    fSimulationTree->Branch("PhotonATPCDistanceY",      &fPhotonATPCDistanceY,      "PhotonATPCDistanceY/D");
    fSimulationTree->Branch("PhotonATPCDistanceZ",      &fPhotonATPCDistanceZ,      "PhotonATPCDistanceZ/D");
    fSimulationTree->Branch("PhotonATPCDistance",       &fPhotonATPCDistance,       "PhotonATPCDistance/D");
    fSimulationTree->Branch("PhotonASphereEnergy",      &fPhotonASphereEnergy,      "PhotonASphereEnergy/D");
    fSimulationTree->Branch("PhotonASphereContainment", &fPhotonASphereContainment, "PhotonASphereContainment/D");
    fSimulationTree->Branch("PhotonASphereRadius",      &fPhotonASphereRadius,      "PhotonASphereRadius/D");

    fSimulationTree->Branch("PhotonBAngleZ",            &fPhotonBAngleZ,            "PhotonBAngleZ/D");
    fSimulationTree->Branch("PhotonBConvert",           &fPhotonBConvert,           "PhotonBConvert/O");
    fSimulationTree->Branch("PhotonBConversionLength",  &fPhotonBConversionLength,  "PhotonBConversionLength/D");
    fSimulationTree->Branch("PhotonBEnergy",            &fPhotonBEnergy,            "PhotonBEnergy/D");
    fSimulationTree->Branch("PhotonBEnergyDeposited",   &fPhotonBEnergyDeposited,   "PhotonBEnergyDeposited/D");
    fSimulationTree->Branch("PhotonBTPCContainment",    &fPhotonBTPCContainment,    "PhotonBTPCContainment/D");
    fSimulationTree->Branch("PhotonBTPCDistanceX",      &fPhotonBTPCDistanceX,      "PhotonBTPCDistanceX/D");
    fSimulationTree->Branch("PhotonBTPCDistanceY",      &fPhotonBTPCDistanceY,      "PhotonBTPCDistanceY/D");
    fSimulationTree->Branch("PhotonBTPCDistanceZ",      &fPhotonBTPCDistanceZ,      "PhotonBTPCDistanceZ/D");
    fSimulationTree->Branch("PhotonBTPCDistance",       &fPhotonBTPCDistance,       "PhotonBTPCDistance/D");
    fSimulationTree->Branch("PhotonBSphereEnergy",      &fPhotonBSphereEnergy,      "PhotonBSphereEnergy/D");
    fSimulationTree->Branch("PhotonBSphereContainment", &fPhotonBSphereContainment, "PhotonBSphereContainment/D");
    fSimulationTree->Branch("PhotonBSphereRadius",      &fPhotonBSphereRadius,      "PhotonBSphereRadius/D");
  }


  //-----------------------------------------------------------------------
  void NeutralPionSimAna::beginRun(const art::Run& /*run*/)
  {}


  //-----------------------------------------------------------------------
  void NeutralPionSimAna::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // read in parameters from the .fcl file
    fSimulationProducerLabel = parameterSet.get< std::string >("SimulationLabel");
  }


  //-----------------------------------------------------------------------
  void NeutralPionSimAna::analyze(const art::Event& event) 
  {
    // get event, run, and subrun numbers
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    // get all the simulated particles for the event
    art::ValidHandle< std::vector< simb::MCParticle > > particleHandle
      = event.getValidHandle< std::vector< simb::MCParticle > >
      (fSimulationProducerLabel);

    // get all the simulated channels for the event
    art::ValidHandle< std::vector< sim::SimChannel > > simChannelHandle
      = event.getValidHandle< std::vector< sim::SimChannel > >
      (fSimulationProducerLabel);

    // map of pointers to MCParticle objects
    std::map< int, const simb::MCParticle* > particleMap;

    // flags for checking if the primary particle is a neutral
    // pion that stops in the TPC
    bool primaryNeutralPion = false;
    bool primaryStopsInTPC  = false;

    // reset variables for new event
    fEnergyDeposited            = 0;
    fNumberNeutralPionsInEvent  = 0;
    fNumberPhotonConversions    = 0;

    fPhotonAConvert = false;
    fPhotonAConversionLength  = fMCTools.detectorDiagonal_();
    fPhotonAAngleZ            = 0;
    fPhotonAEnergy            = 0;
    fPhotonAEnergyDeposited   = 0;
    fPhotonATPCContainment    = 0;
    fPhotonATPCDistanceX      = fMCTools.detectorDiagonal_();
    fPhotonATPCDistanceY      = fMCTools.detectorDiagonal_();
    fPhotonATPCDistanceZ      = fMCTools.detectorDiagonal_();
    fPhotonATPCDistance       = fMCTools.detectorDiagonal_();
    fPhotonASphereEnergy      = 0;
    fPhotonASphereContainment = 0;
    fPhotonASphereRadius      = fMCTools.detectorDiagonal_();

    fPhotonBConvert = false;
    fPhotonBConversionLength  = fMCTools.detectorDiagonal_();
    fPhotonBAngleZ            = 0;
    fPhotonBEnergy            = 0;
    fPhotonBEnergyDeposited   = 0;
    fPhotonBTPCContainment    = 0;
    fPhotonATPCDistanceX      = fMCTools.detectorDiagonal_();
    fPhotonATPCDistanceY      = fMCTools.detectorDiagonal_();
    fPhotonATPCDistanceZ      = fMCTools.detectorDiagonal_();
    fPhotonBTPCDistance       = fMCTools.detectorDiagonal_();
    fPhotonBSphereEnergy      = 0;
    fPhotonBSphereContainment = 0;
    fPhotonBSphereRadius      = fMCTools.detectorDiagonal_();

    fPhotonEnergy.clear();
    fNeutralPionEnergy.clear();
    fNeutralPionMomentum.clear();

    // declare vector of neutral pion track IDs
    std::vector<int> neutralPionTrackIDs;

    // declare vector of photon track IDs
    //std::vector<int> photonTrackIDs;
    std::vector< std::pair< double, int > > photonEnergyTrackID;

    // declare map of shower particle map where the key is the
    // track ID of the photon
    std::map< int, std::map< int, const simb::MCParticle* > > showerParticleMaps;

    // loop through MCParticle objects
    for (auto const& particle : (*particleHandle)) // simb::MCParticle
    {
      // get track ID of particle
      int trackID = particle.TrackId();

      // get PDG code of particle
      int pdgCode = particle.PdgCode();

      // add the address of the MCParticle to the map, with 
      // the track ID as the key
      particleMap[trackID] = &particle; 

      // if particle is primary and neutral pion
      if (particle.Process() == "primary" && std::abs(pdgCode) == 111)
      {
        // change flag
        primaryNeutralPion = true;

        // get primary track ID
        fPrimaryTrackID = trackID;

        // get primary PDG code
        fPrimaryPDGCode = pdgCode;

        // get primary start and end 4-positions and 4-momenta
        fMCTools.getFourVectors_(particle,
                                 fStartXYZT,
                                 fEndXYZT,
                                 fStartPE,
                                 fEndPE);

        // check to see if primary neutral pion stops inside TPC
        primaryStopsInTPC = fMCTools.isInsideTPC_(fEndXYZT[0],
                                                  fEndXYZT[1],
                                                  fEndXYZT[2]);

        // get 3D distance between the start and end
        // positions of the primary
        const double trajectoryLength = fMCTools.getTrajectoryLength_(particle);
        LOG_DEBUG("NeutralPionSimAna")
          << "Trajectory length: " << trajectoryLength << " cm";

        // return if primary neutral pion does not stop inside TPC
        if (!primaryStopsInTPC)
          return;

        fNeutralPionEnergy.push_back(fStartPE[3] * 1000);
        fNeutralPionMomentum.push_back(std::sqrt(fStartPE[0]*fStartPE[0] + fStartPE[1]*fStartPE[1] + fStartPE[2]*fStartPE[2]) * 1000.0);

        neutralPionTrackIDs.push_back(trackID);

      } // if primary and neutral pion

      // if neutral pion, increment neutral pion counter and keep track of track ID
      if (particle.PdgCode() == 111)
        ++fNumberNeutralPionsInEvent;

    } // loop over all particles in the event

    // return if primary is not a neutral pion
    if (!primaryNeutralPion)
      return;

    // loop over neutral pions in the event
    for (auto const& trackID : neutralPionTrackIDs)
    {
      // get particle from map of particles using track ID as the key
      auto const& particle = particleMap[trackID];

      // get number of daughter particles
      int numberDaughters = particle->NumberDaughters();

      // for each daughter particle
      for (int i = 0; i < numberDaughters; ++i)
      {
        // get daughter track ID
        int daughterTrackID = particle->Daughter(i);

        // get daughter particle from map of particles using daughter
        // track ID as the key
        auto const& daughterParticle = particleMap[daughterTrackID];

        // if daughter particle is a photon from the decay of
        // the neutral pion
        if (daughterParticle->PdgCode() == 22 &&
            daughterParticle->Process() == "Decay")
        {
          // declare start/end positions and momenta arrays
          double startXYZT[4];
          double endXYZT[4];
          double startPE[4];
          double endPE[4];

          // get start and end 4-positions and 4-momenta
          fMCTools.getFourVectors_(*daughterParticle,
                                   startXYZT,
                                   endXYZT,
                                   startPE,
                                   endPE);

          fPhotonEnergy.push_back(startPE[3] * 1000.0);

          // keep track of photons that decay from the
          // neutral pion
          photonEnergyTrackID.push_back(std::pair<double, int>(startPE[3] * 1000.0, // MeV
                                                               daughterTrackID));

          // check to see if photon converts inside TPC
          bool photonConvertsInTPC = fMCTools.isInsideTPC_(endXYZT[0],
                                                           endXYZT[1],
                                                           endXYZT[2]);

          // if photon converts inside TPC
          if (photonConvertsInTPC)
          {
            // increment counter for photons from neutral pion
            // decay that convert inside TPC
            ++fNumberPhotonConversions;
          }

        } // if particle is a photon from neutral pion decay
      } // loop over daughter particles of neutral pion
    } // loop over neutral pions in the event

    // return if the number of photons is not 2
    if (photonEnergyTrackID.size() != 2)
      return;

    // sort photons by energy
    std::sort(photonEnergyTrackID.begin(), photonEnergyTrackID.end());
    std::reverse(photonEnergyTrackID.begin(), photonEnergyTrackID.end());

    // set photon A and B track IDs
    int photonATrackID = photonEnergyTrackID.at(0).second;
    int photonBTrackID = photonEnergyTrackID.at(1).second;

    // get shower variables for photon A
    fMCTools.getShowerVariables_(photonATrackID,
                                 particleMap,
                                 simChannelHandle,
                                 showerParticleMaps[photonATrackID],
                                 fPhotonAStartXYZT,
                                 fPhotonAEndXYZT,
                                 fPhotonAStartPE,
                                 fPhotonAEndPE,
                                 fPhotonAAngleZ,
                                 fPhotonAConvert,
                                 fPhotonAConversionLength,
                                 fPhotonAEnergy,
                                 fPhotonAEnergyDeposited,
                                 fPhotonATPCContainment,
                                 fPhotonATPCDistanceX,
                                 fPhotonATPCDistanceY,
                                 fPhotonATPCDistanceZ,
                                 fPhotonATPCDistance,
                                 fPhotonASphereEnergy,
                                 fPhotonASphereContainment,
                                 fPhotonASphereRadius);

    // get shower variables for photon B
    fMCTools.getShowerVariables_(photonBTrackID,
                                 particleMap,
                                 simChannelHandle,
                                 showerParticleMaps[photonBTrackID],
                                 fPhotonBStartXYZT,
                                 fPhotonBEndXYZT,
                                 fPhotonBStartPE,
                                 fPhotonBEndPE,
                                 fPhotonBAngleZ,
                                 fPhotonBConvert,
                                 fPhotonBConversionLength,
                                 fPhotonBEnergy,
                                 fPhotonBEnergyDeposited,
                                 fPhotonBTPCContainment,
                                 fPhotonBTPCDistanceX,
                                 fPhotonBTPCDistanceY,
                                 fPhotonBTPCDistanceZ,
                                 fPhotonBTPCDistance,
                                 fPhotonBSphereEnergy,
                                 fPhotonBSphereContainment,
                                 fPhotonBSphereRadius);

    // get opening angle between photon A and photon B
    fPhotonAngle = fMCTools.angle3D_(fPhotonAStartPE, fPhotonBStartPE);

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

  } // NeutralPionSimAna::analyze()


  DEFINE_ART_MODULE(NeutralPionSimAna)

} // namespace NeutralPionSimAna


#endif // NeutralPionSimAna_Module
