#include <string>
#include <numeric>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// Adding header files for muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// Adding header files for generating Ntuples
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

// eff vs PU test libraries
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimTracker/VertexAssociation/interface/calculateVertexSharedTracks.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

// for test
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

class MtdMuonIsoValidation : public DQMEDAnalyzer {
public:
  explicit MtdMuonIsoValidation(const edm::ParameterSet&);
  ~MtdMuonIsoValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMinEta_;
  const float trackMaxEta_;
  const double rel_iso_cut_;

  const bool track_match_PV_;
  const bool dt_sig_track_;
  const bool optionalPlots_;

  const float min_dR_cut;
  const float max_dR_cut;
  const float min_pt_cut_EB;
  const float min_pt_cut_EE;
  const float max_dz_cut_EB;
  const float max_dz_cut_EE;
  const float max_dz_vtx_cut;
  const float max_dxy_vtx_cut;
///  const float min_strip_cut;
  const float min_track_mtd_mva_cut;
  const std::vector<double> max_dt_vtx_cut{0.30, 0.27, 0.24, 0.21, 0.18, 0.15, 0.12};
  const std::vector<double> max_dt_track_cut{0.30, 0.27, 0.24, 0.21, 0.18, 0.15, 0.12};
  const std::vector<double> max_dt_significance_cut{4.0, 3.0, 2.0};
  const std::vector<double> pT_bins_dt_distrb{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const std::vector<double> eta_bins_dt_distrib{0.0, 0.5, 1.0, 1.5, 2.0, 2.4, 2.7, 3};
  static constexpr double avg_sim_sigTrk_t_err = 0.03668;
  static constexpr double avg_sim_PUtrack_t_err = 0.03461;

  // test
  static constexpr double c_ = 2.99792458e1;  // c in cm/ns
  static constexpr double maxTry_ = 10.;
  static constexpr double zWosMatchMax_ = 1.;
  static constexpr double maxRank_ = 8.;
  const reco::RecoToSimCollection* r2s__;
  const reco::SimToRecoCollection* s2r_;
  edm::EDGetTokenT<reco::BeamSpot> RecBeamSpotToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> Rec4DVerToken_;
  // test end

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<reco::MuonCollection> MuonToken_;
  edm::EDGetTokenT<reco::TrackCollection> GlobalMuonTrk_; // Adding token for global muon track collection
							  // It is unnecessary. I checked GlobalMuonTrk is the same as muon.combinedMuon().
							  // Now tracker muon is used
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;

  // TFile Services
  edm::Service<TFileService> fs_;
  TTree* tree_;
  int run_;
  int event_;
  std::vector<float> track_type_v1_, track_type_v2_;
  std::vector<float> muon_pt_, track_pt_;
  std::vector<float> muon_time_, vtx_time_, track_time_;
  std::vector<float> dtsig_muon_track_, dtsig_vtx_track_;
  std::vector<float> muon_PVweight_, track_PVweight_;
  std::vector<float> muon_time_err_, vtx_time_err_, track_time_err_;
  std::vector<bool>  muon_prompt_, muon_isBarrel_;
  std::vector<int>   track_bx_, track_evtId_;

  int vtx_index_;
  int recovtx_sim_, simvtx_reco_, simvtx_bx_, simvtx_evtId_;
  std::vector<bool> selectedVtxMatching_, selectedLV_, match_vtx_reco2sim_, match_vtx_sim2reco_;
  float simvtx_pt_, simvtx_ptsq_, recovtx_pt_, recovtx_ptsq_;
  int simvtx_nGenTrk_, simvtx_num_matched_reco_tracks_, recovtx_nRecoTrk_, recovtx_num_matched_sim_tracks_;


  // test
  typedef math::XYZTLorentzVector LorentzVector;
  static constexpr unsigned int NOT_MATCHED = 66666;
  // auxiliary class holding simulated vertices
  struct simPrimaryVertex {
    simPrimaryVertex(double x1, double y1, double z1, double t1)
        : x(x1),
          y(y1),
          z(z1),
          t(t1),
          pt(0),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          nGenTrk(0),
          num_matched_reco_tracks(0),
          average_match_quality(0.0) {
      ptot.setPx(0);
      ptot.setPy(0);
      ptot.setPz(0);
      ptot.setE(0);
      p4 = LorentzVector(0, 0, 0, 0);
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r, t;
    HepMC::FourVector ptot;
    LorentzVector p4;
    double pt;
    double ptsq;
    double closest_vertex_distance_z;
    int nGenTrk;
    int num_matched_reco_tracks;
    float average_match_quality;
    EncodedEventId eventId;
    TrackingVertexRef sim_vertex;
    int OriginalIndex = -1;

    unsigned int nwosmatch = 0;                    // number of recvertices dominated by this simevt (by wos)
    unsigned int nwntmatch = 0;                    // number of recvertices dominated by this simevt  (by nt)
    std::vector<unsigned int> wos_dominated_recv;  // list of dominated recv (by wos, size==nwosmatch)

    std::map<unsigned int, double> wnt;  // weighted number of tracks in recvtx (by index)
    std::map<unsigned int, double> wos;  // sum of wos in recvtx (by index) // oops -> this was int before 04-22
    double sumwos = 0;                   // sum of wos in any recvtx
    double sumwnt = 0;                   // sum of weighted tracks
    unsigned int rec = NOT_MATCHED;      // best match (NO_MATCH if not matched)
    unsigned int matchQuality = 0;       // quality flag
    void addTrack(unsigned int irecv, double twos, double twt) {
      sumwnt += twt;
      if (wnt.find(irecv) == wnt.end()) {
        wnt[irecv] = twt;
      } else {
        wnt[irecv] += twt;
      }

      sumwos += twos;
      if (wos.find(irecv) == wos.end()) {
        wos[irecv] = twos;
      } else {
        wos[irecv] += twos;
      }
    }
  };
  // auxiliary class holding reconstructed vertices
  struct recoPrimaryVertex {
    recoPrimaryVertex(double x1, double y1, double z1)
        : x(x1),
          y(y1),
          z(z1),
          pt(0),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          nRecoTrk(0),
          num_matched_sim_tracks(0),
          ndof(0.),
          recVtx(nullptr) {
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r;
    double pt;
    double ptsq;
    double closest_vertex_distance_z;
    int nRecoTrk;
    int num_matched_sim_tracks;
    double ndof;
    const reco::Vertex* recVtx;
    reco::VertexBaseRef recVtxRef;
    int OriginalIndex = -1;

    std::map<unsigned int, double> wos;  // simevent -> wos
    std::map<unsigned int, double> wnt;  // simevent -> weighted number of truth matched tracks
    unsigned int wosmatch;               // index of the simevent providing the largest contribution to wos
    unsigned int wntmatch;               // index of the simevent providing the highest number of tracks
    double sumwos = 0;                   // total sum of wos of all truth matched tracks
    double sumwnt = 0;                   // total weighted number of truth matchted tracks
    double maxwos = 0;                   // largest wos sum from one sim event (wosmatch)
    double maxwnt = 0;                   // largest weighted number of tracks from one sim event (ntmatch)
    int maxwosnt = 0;                    // number of tracks from the simevt with highest wos
    unsigned int sim = NOT_MATCHED;      // best match (NO_MATCH if not matched)
    unsigned int matchQuality = 0;       // quality flag

    bool is_real() { return (matchQuality > 0) && (matchQuality < 99); }

    bool is_fake() { return (matchQuality <= 0) || (matchQuality >= 99); }

    bool is_signal() { return (sim == 0); }

    int split_from() {
      if (is_real())
        return -1;
      if ((maxwos > 0) && (maxwos > 0.3 * sumwos))
        return wosmatch;
      return -1;
    }
    bool other_fake() { return (is_fake() && (split_from() < 0)); }

    void addTrack(unsigned int iev, double twos, double wt) {
      sumwnt += wt;
      if (wnt.find(iev) == wnt.end()) {
        wnt[iev] = wt;
      } else {
        wnt[iev] += wt;
      }

      sumwos += twos;
      if (wos.find(iev) == wos.end()) {
        wos[iev] = twos;
      } else {
        wos[iev] += twos;
      }
    }
  };

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;

  void matchReco2Sim(std::vector<recoPrimaryVertex>&,
                     std::vector<simPrimaryVertex>&,
                     const edm::ValueMap<float>&,
                     const edm::ValueMap<float>&,
                     const edm::Handle<reco::BeamSpot>&);
  void getWosWnt(const reco::Vertex&,
                 const reco::TrackBaseRef&,
                 const edm::ValueMap<float>&,
                 const edm::Handle<reco::BeamSpot>&,
                 double&,
                 double&);
  std::pair<const edm::Ref<std::vector<TrackingParticle>>*, int> getMatchedTP(const reco::TrackBaseRef&,
                                                                              const TrackingVertexRef&);
  std::vector<MtdMuonIsoValidation::simPrimaryVertex> getSimPVs(const edm::Handle<TrackingVertexCollection>&);
  std::vector<MtdMuonIsoValidation::recoPrimaryVertex> getRecoPVs(const edm::Handle<edm::View<reco::Vertex>>&);

  void printSimVtxRecoVtxInfo(const struct MtdMuonIsoValidation::simPrimaryVertex&,
                              const struct MtdMuonIsoValidation::recoPrimaryVertex&);

  // test end


  // histograms for testing  FIXME: Need to be cleaned
  
  // test for GEN CASE
  std::vector<MonitorElement*> Muon_pT_gen_MTD_EB_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_pT_gen_MTD_EE_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_pT_gen_MTD_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_pT_gen_MTD_EE_list_Significance_Bkg;

  MonitorElement* meMuon_pt_gen_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuon_pt_gen_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuon_pt_gen_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuon_pt_gen_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuon_pt_gen_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuon_pt_gen_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuon_pt_gen_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_gen_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_gen_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_gen_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuon_pt_gen_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuon_pt_gen_MTD_4sigma_Bkg_EE_;

  std::vector<MonitorElement*> ch_iso_gen_EB_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_gen_EE_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_gen_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> ch_iso_gen_EE_list_Significance_Bkg;

  MonitorElement* meMuonISO_chIso_MTD_gen_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_gen_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_gen_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_gen_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_gen_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_gen_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_gen_4sigma_Bkg_EE_;
  
  std::vector<MonitorElement*> rel_ch_iso_gen_EB_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_gen_EE_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_gen_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_gen_EE_list_Significance_Bkg;
  
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EE_;
  
  // Several tests
  
  // test to checking the type of tracks
  MonitorElement* meMuonISO_trk_type_Sig_EB_;
  MonitorElement* meMuonISO_trk_type_Sig_EE_;
  MonitorElement* meMuonISO_trk_type_Bkg_EB_;
  MonitorElement* meMuonISO_trk_type_Bkg_EE_;
  MonitorElement* meMuonISO_trk_type_v2_Sig_EB_;
  MonitorElement* meMuonISO_trk_type_v2_Sig_EE_;
  MonitorElement* meMuonISO_trk_type_v2_Bkg_EB_;
  MonitorElement* meMuonISO_trk_type_v2_Bkg_EE_;

  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EE_;

  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EE_;

  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EE_;

  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EE_;
  
  //////////
  
  MonitorElement* meMuonISO_trk_genMatched_Sig_;
  MonitorElement* meMuonISO_trk_genMatched_Sig_EB_;
  MonitorElement* meMuonISO_trk_genMatched_Sig_EE_;
  MonitorElement* meMuonISO_trk_genMatched_Bkg_;
  MonitorElement* meMuonISO_trk_genMatched_Bkg_EB_;
  MonitorElement* meMuonISO_trk_genMatched_Bkg_EE_;

  MonitorElement* meMuonISO_pTdiff_reco_tracker_muon_Sig_;
  MonitorElement* meMuonISO_pTdiff_reco_tracker_muon_Bkg_;
  
  MonitorElement* meMuonISO_pT_muon_sim_;
  MonitorElement* meMuonISO_pT_trk_sim_;
  
  MonitorElement* meMuonISO_time_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_time_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_time_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_time_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_tErr_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_tErr_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_tErr_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_tErr_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_time_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_time_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_time_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_time_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_tErr_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_tErr_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_tErr_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_tErr_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_time_muon_sim_Sig_EB_;
  MonitorElement* meMuonISO_time_muon_sim_Sig_EE_;
  MonitorElement* meMuonISO_time_muon_sim_Bkg_EB_;
  MonitorElement* meMuonISO_time_muon_sim_Bkg_EE_;
  MonitorElement* meMuonISO_time_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_time_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_time_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_time_trk_sim_Bkg_EE_;

  MonitorElement* meMuonISO_has_time_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_time_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_time_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_time_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_time_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_time_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_time_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_time_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_tErr_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_tErr_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_tErr_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_tErr_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_tErr_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_tErr_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_tErr_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_tErr_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_mva_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_mva_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_mva_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_mva_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_mva_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_has_mva_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_has_mva_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_has_mva_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_has_time_muon_sim_Sig_EB_;
  MonitorElement* meMuonISO_has_time_muon_sim_Sig_EE_;
  MonitorElement* meMuonISO_has_time_muon_sim_Bkg_EB_;
  MonitorElement* meMuonISO_has_time_muon_sim_Bkg_EE_;
  MonitorElement* meMuonISO_has_time_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_has_time_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_has_time_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_has_time_trk_sim_Bkg_EE_;

  MonitorElement* meMuonISO_dt_muon_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_dt_muon_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_dt_muon_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_dt_muon_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_dt_muon_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_dt_muon_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_dt_muon_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_dt_muon_trk_sim_Bkg_EE_;
  MonitorElement* meMuonISO_dt_PV_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_dt_PV_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_dt_PV_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_dt_PV_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_dt_PV_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_dt_PV_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_dt_PV_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_dt_PV_trk_sim_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EE_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EB_;
  MonitorElement* meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EE_;

  MonitorElement* meMuonISO_mva_muon_reco_Sig_;
  MonitorElement* meMuonISO_mva_muon_reco_Sig_EB_;
  MonitorElement* meMuonISO_mva_muon_reco_Sig_EE_;
  MonitorElement* meMuonISO_mva_muon_reco_Bkg_;
  MonitorElement* meMuonISO_mva_muon_reco_Bkg_EB_;
  MonitorElement* meMuonISO_mva_muon_reco_Bkg_EE_;
  MonitorElement* meMuonISO_mva_trk_reco_Sig_;
  MonitorElement* meMuonISO_mva_trk_reco_Sig_EB_;
  MonitorElement* meMuonISO_mva_trk_reco_Sig_EE_;
  MonitorElement* meMuonISO_mva_trk_reco_Bkg_;
  MonitorElement* meMuonISO_mva_trk_reco_Bkg_EB_;
  MonitorElement* meMuonISO_mva_trk_reco_Bkg_EE_;

  MonitorElement* meMuonISO_dz_muon_;
  MonitorElement* meMuonISO_dz_muon_EB_;
  MonitorElement* meMuonISO_dz_muon_EE_;
  MonitorElement* meMuonISO_dxy_muon_;
  MonitorElement* meMuonISO_dxy_muon_EB_;
  MonitorElement* meMuonISO_dxy_muon_EE_;

  MonitorElement* meMuonISO_Nmuons_Sig_;
  MonitorElement* meMuonISO_Nmuons_Bkg_;
  MonitorElement* meMuonISO_Nmuons_Sig_EB_;
  MonitorElement* meMuonISO_Nmuons_Sig_EE_;
  MonitorElement* meMuonISO_Nmuons_Bkg_EB_;
  MonitorElement* meMuonISO_Nmuons_Bkg_EE_;

  // Signal histograms

  MonitorElement* meMuon_no_dt_check_;
  MonitorElement* meTrk_genMatch_check_;

  MonitorElement* meMuon_avg_error_SigTrk_check_;
  MonitorElement* meMuon_avg_error_PUTrk_check_;
  MonitorElement* meMuon_avg_error_vtx_check_;

  // Adding histograms for barrel muons
  MonitorElement* meMuonISO_Ntracks_Sig_EB_;
  MonitorElement* meMuonISO_chIso_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_Sig_EB_;
  MonitorElement* meMuonISO_Ntracks_MTD_1_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_1_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_1_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_2_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_2_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_3_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_3_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_4_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_4_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_5_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_5_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_5_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_6_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_6_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_6_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_7_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_7_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_7_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_1_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_1_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_1_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_5_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_5_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_5_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_6_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_6_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_6_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_7_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_7_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_7_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_gen_Sig_EB_;
  MonitorElement* meMuonISO_chIso_gen_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_gen_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2sigma_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3sigma_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4sigma_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4sigma_Sig_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EB_;

  MonitorElement* meMuon_pt_tot_Sig_EB_;
  MonitorElement* meMuon_pt_sim_tot_Sig_EB_;  // for GEN case is the same 
  MonitorElement* meMuon_eta_tot_Sig_EB_;
  MonitorElement* meMuon_phi_tot_Sig_EB_;
  MonitorElement* meMuon_pt_MTD_1_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_1_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_1_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_1_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_2_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_2_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_2_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_2_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_3_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_3_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_3_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_3_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_4_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_4_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_4_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_4_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_5_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_5_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_5_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_5_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_6_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_6_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_6_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_6_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_7_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_7_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_7_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_7_Sig_EB_;

  MonitorElement* meMuon_pt_noMTD_Sig_EB_;
  MonitorElement* meMuon_eta_noMTD_Sig_EB_;
  MonitorElement* meMuon_phi_noMTD_Sig_EB_;

  MonitorElement* meMuon_pt_gen_Sig_EB_;
  MonitorElement* meMuon_eta_gen_Sig_EB_;
  MonitorElement* meMuon_phi_gen_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_2sigma_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_2sigma_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_3sigma_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_3sigma_Sig_EB_;

  MonitorElement* meMuon_pt_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuon_eta_MTD_4sigma_Sig_EB_;
  MonitorElement* meMuon_phi_MTD_4sigma_Sig_EB_;

  // Adding histograms for endcap muons
  MonitorElement* meMuonISO_Ntracks_Sig_EE_;
  MonitorElement* meMuonISO_chIso_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_Sig_EE_;
  MonitorElement* meMuonISO_Ntracks_MTD_1_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_1_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_1_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_2_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_2_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_3_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_3_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_4_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_4_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_5_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_5_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_5_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_6_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_6_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_6_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_7_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_7_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_7_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_1_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_1_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_1_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_5_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_5_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_5_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_6_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_6_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_6_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_7_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_7_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_7_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_gen_Sig_EE_;
  MonitorElement* meMuonISO_chIso_gen_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_gen_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2sigma_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3sigma_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4sigma_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4sigma_Sig_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EE_;

  MonitorElement* meMuon_pt_tot_Sig_EE_;
  MonitorElement* meMuon_pt_sim_tot_Sig_EE_;
  MonitorElement* meMuon_eta_tot_Sig_EE_;
  MonitorElement* meMuon_phi_tot_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_1_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_1_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_1_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_1_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_2_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_2_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_2_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_2_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_3_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_3_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_3_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_3_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_4_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_4_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_4_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_4_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_5_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_5_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_5_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_5_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_6_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_6_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_6_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_6_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_7_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_7_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_7_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_7_Sig_EE_;

  MonitorElement* meMuon_pt_noMTD_Sig_EE_;
  MonitorElement* meMuon_eta_noMTD_Sig_EE_;
  MonitorElement* meMuon_phi_noMTD_Sig_EE_;

  MonitorElement* meMuon_pt_gen_Sig_EE_;
  MonitorElement* meMuon_eta_gen_Sig_EE_;
  MonitorElement* meMuon_phi_gen_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_2sigma_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_2sigma_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_3sigma_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_3sigma_Sig_EE_;

  MonitorElement* meMuon_pt_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuon_eta_MTD_4sigma_Sig_EE_;
  MonitorElement* meMuon_phi_MTD_4sigma_Sig_EE_;

  // Signal histograms end

  // Background histograms
  // Adding histograms for barrel muons
  MonitorElement* meMuonISO_Ntracks_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_Bkg_EB_;
  MonitorElement* meMuonISO_Ntracks_MTD_1_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_1_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_1_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_2_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_2_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_3_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_3_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_4_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_4_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_5_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_5_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_5_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_6_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_6_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_6_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_7_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_7_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_7_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_1_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_1_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_1_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_5_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_5_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_5_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_6_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_6_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_6_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_7_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_7_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_7_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_gen_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_gen_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_gen_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2sigma_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3sigma_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4sigma_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EB_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4sigma_Bkg_EB_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EB_;

  MonitorElement* meMuon_pt_tot_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_tot_Bkg_EB_;
  MonitorElement* meMuon_eta_tot_Bkg_EB_;
  MonitorElement* meMuon_phi_tot_Bkg_EB_;
  MonitorElement* meMuon_pt_MTD_1_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_1_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_1_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_1_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_2_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_2_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_2_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_2_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_3_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_3_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_3_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_3_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_4_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_4_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_4_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_4_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_5_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_5_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_5_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_5_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_6_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_6_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_6_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_6_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_7_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_7_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_7_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_7_Bkg_EB_;

  MonitorElement* meMuon_pt_noMTD_Bkg_EB_;
  MonitorElement* meMuon_eta_noMTD_Bkg_EB_;
  MonitorElement* meMuon_phi_noMTD_Bkg_EB_;

  MonitorElement* meMuon_pt_gen_Bkg_EB_;
  MonitorElement* meMuon_eta_gen_Bkg_EB_;
  MonitorElement* meMuon_phi_gen_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_2sigma_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_2sigma_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_3sigma_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_3sigma_Bkg_EB_;

  MonitorElement* meMuon_pt_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuon_eta_MTD_4sigma_Bkg_EB_;
  MonitorElement* meMuon_phi_MTD_4sigma_Bkg_EB_;

  // Adding histograms for endcap muons
  MonitorElement* meMuonISO_Ntracks_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_Bkg_EE_;
  MonitorElement* meMuonISO_Ntracks_MTD_1_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_1_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_1_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_2_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_2_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_3_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_3_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_4_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_4_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_5_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_5_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_5_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_6_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_6_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_6_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_7_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_7_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_7_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_1_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_1_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_1_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_5_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_5_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_5_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_6_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_6_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_6_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_7_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_7_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_7_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_gen_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_gen_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_gen_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_2sigma_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_3sigma_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_4sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_4sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_4sigma_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_2sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_3sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EE_;

  MonitorElement* meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EE_;
  MonitorElement* meMuonISO_chIso_MTD_sim_4sigma_Bkg_EE_;
  MonitorElement* meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EE_;

  MonitorElement* meMuon_pt_tot_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_tot_Bkg_EE_;
  MonitorElement* meMuon_eta_tot_Bkg_EE_;
  MonitorElement* meMuon_phi_tot_Bkg_EE_;
  MonitorElement* meMuon_pt_MTD_1_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_1_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_1_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_1_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_2_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_2_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_2_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_2_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_3_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_3_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_3_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_3_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_4_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_4_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_4_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_4_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_5_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_5_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_5_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_5_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_6_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_6_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_6_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_6_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_7_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_7_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_7_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_7_Bkg_EE_;

  MonitorElement* meMuon_pt_noMTD_Bkg_EE_;
  MonitorElement* meMuon_eta_noMTD_Bkg_EE_;
  MonitorElement* meMuon_phi_noMTD_Bkg_EE_;

  MonitorElement* meMuon_pt_gen_Bkg_EE_;
  MonitorElement* meMuon_eta_gen_Bkg_EE_;
  MonitorElement* meMuon_phi_gen_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_2sigma_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_2sigma_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_3sigma_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_3sigma_Bkg_EE_;

  MonitorElement* meMuon_pt_MTD_4sigma_Bkg_EE_;
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Bkg_EE_;
  MonitorElement* meMuon_eta_MTD_4sigma_Bkg_EE_;
  MonitorElement* meMuon_phi_MTD_4sigma_Bkg_EE_;
  // Background histograms end

  // promt part for histogram vectors
  std::vector<MonitorElement*> Ntracks_EB_list_Sig;
  std::vector<MonitorElement*> ch_iso_EB_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Sig;

  std::vector<MonitorElement*> Ntracks_EE_list_Sig;
  std::vector<MonitorElement*> ch_iso_EE_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Sig;

  std::vector<MonitorElement*> Ntracks_sim_EB_list_Sig;
  std::vector<MonitorElement*> ch_iso_sim_EB_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_sim_EB_list_Sig;

  std::vector<MonitorElement*> Ntracks_sim_EE_list_Sig;
  std::vector<MonitorElement*> ch_iso_sim_EE_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_sim_EE_list_Sig;

  std::vector<MonitorElement*> Muon_pT_MTD_EB_list_Sig;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EB_list_Sig;
  std::vector<MonitorElement*> Muon_eta_MTD_EB_list_Sig;
  std::vector<MonitorElement*> Muon_phi_MTD_EB_list_Sig;

  std::vector<MonitorElement*> Muon_pT_MTD_EE_list_Sig;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EE_list_Sig;
  std::vector<MonitorElement*> Muon_eta_MTD_EE_list_Sig;
  std::vector<MonitorElement*> Muon_phi_MTD_EE_list_Sig;

  std::vector<MonitorElement*> Ntracks_EB_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_EB_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Significance_Sig;

  std::vector<MonitorElement*> Ntracks_EE_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_EE_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Significance_Sig;

  std::vector<MonitorElement*> Ntracks_sim_EB_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_sim_EB_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_sim_EB_list_Significance_Sig;

  std::vector<MonitorElement*> Ntracks_sim_EE_list_Significance_Sig;
  std::vector<MonitorElement*> ch_iso_sim_EE_list_Significance_Sig;
  std::vector<MonitorElement*> rel_ch_iso_sim_EE_list_Significance_Sig;

  std::vector<MonitorElement*> Muon_pT_MTD_EB_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EB_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_eta_MTD_EB_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_phi_MTD_EB_list_Significance_Sig;

  std::vector<MonitorElement*> Muon_pT_MTD_EE_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EE_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_eta_MTD_EE_list_Significance_Sig;
  std::vector<MonitorElement*> Muon_phi_MTD_EE_list_Significance_Sig;

  // Non-promt part for histogram vectors
  std::vector<MonitorElement*> Ntracks_EB_list_Bkg;
  std::vector<MonitorElement*> ch_iso_EB_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Bkg;

  std::vector<MonitorElement*> Ntracks_EE_list_Bkg;
  std::vector<MonitorElement*> ch_iso_EE_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Bkg;

  std::vector<MonitorElement*> Ntracks_sim_EB_list_Bkg;
  std::vector<MonitorElement*> ch_iso_sim_EB_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_sim_EB_list_Bkg;

  std::vector<MonitorElement*> Ntracks_sim_EE_list_Bkg;
  std::vector<MonitorElement*> ch_iso_sim_EE_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_sim_EE_list_Bkg;

  std::vector<MonitorElement*> Muon_pT_MTD_EB_list_Bkg;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EB_list_Bkg;
  std::vector<MonitorElement*> Muon_eta_MTD_EB_list_Bkg;
  std::vector<MonitorElement*> Muon_phi_MTD_EB_list_Bkg;

  std::vector<MonitorElement*> Muon_pT_MTD_EE_list_Bkg;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EE_list_Bkg;
  std::vector<MonitorElement*> Muon_eta_MTD_EE_list_Bkg;
  std::vector<MonitorElement*> Muon_phi_MTD_EE_list_Bkg;

  std::vector<MonitorElement*> Ntracks_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> ch_iso_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Significance_Bkg;

  std::vector<MonitorElement*> Ntracks_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> ch_iso_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Significance_Bkg;

  std::vector<MonitorElement*> Ntracks_sim_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> ch_iso_sim_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_sim_EB_list_Significance_Bkg;

  std::vector<MonitorElement*> Ntracks_sim_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> ch_iso_sim_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_sim_EE_list_Significance_Bkg;

  std::vector<MonitorElement*> Muon_pT_MTD_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_eta_MTD_EB_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_phi_MTD_EB_list_Significance_Bkg;

  std::vector<MonitorElement*> Muon_pT_MTD_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_pT_sim_MTD_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_eta_MTD_EE_list_Significance_Bkg;
  std::vector<MonitorElement*> Muon_phi_MTD_EE_list_Significance_Bkg;

  // dt distribution part for histogram vectors
  std::vector<MonitorElement*> general_pT_list;
  std::vector<MonitorElement*> general_eta_list;

  std::vector<MonitorElement*> general_pT_Signif_list;
  std::vector<MonitorElement*> general_eta_Signif_list;

///
  // histograms for dt distributions in pT/eta bins

  MonitorElement* meMuon_dt_general_pT_1;
  MonitorElement* meMuon_dt_general_pT_2;
  MonitorElement* meMuon_dt_general_pT_3;
  MonitorElement* meMuon_dt_general_pT_4;
  MonitorElement* meMuon_dt_general_pT_5;
  MonitorElement* meMuon_dt_general_pT_6;
  MonitorElement* meMuon_dt_general_pT_7;
  MonitorElement* meMuon_dt_general_pT_8;
  MonitorElement* meMuon_dt_general_pT_9;

  MonitorElement* meMuon_dtSignif_general_pT_1;
  MonitorElement* meMuon_dtSignif_general_pT_2;
  MonitorElement* meMuon_dtSignif_general_pT_3;
  MonitorElement* meMuon_dtSignif_general_pT_4;
  MonitorElement* meMuon_dtSignif_general_pT_5;
  MonitorElement* meMuon_dtSignif_general_pT_6;
  MonitorElement* meMuon_dtSignif_general_pT_7;
  MonitorElement* meMuon_dtSignif_general_pT_8;
  MonitorElement* meMuon_dtSignif_general_pT_9;

  MonitorElement* meMuon_dt_general_eta_1;
  MonitorElement* meMuon_dt_general_eta_2;
  MonitorElement* meMuon_dt_general_eta_3;
  MonitorElement* meMuon_dt_general_eta_4;
  MonitorElement* meMuon_dt_general_eta_5;
  MonitorElement* meMuon_dt_general_eta_6;
  MonitorElement* meMuon_dt_general_eta_7;

  MonitorElement* meMuon_dtSignif_general_eta_1;
  MonitorElement* meMuon_dtSignif_general_eta_2;
  MonitorElement* meMuon_dtSignif_general_eta_3;
  MonitorElement* meMuon_dtSignif_general_eta_4;
  MonitorElement* meMuon_dtSignif_general_eta_5;
  MonitorElement* meMuon_dtSignif_general_eta_6;
  MonitorElement* meMuon_dtSignif_general_eta_7;
};

// ------------ constructor and destructor --------------
MtdMuonIsoValidation::MtdMuonIsoValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMinEta_(iConfig.getParameter<double>("trackMinimumEta")),
      trackMaxEta_(iConfig.getParameter<double>("trackMaximumEta")),
      rel_iso_cut_(iConfig.getParameter<double>("rel_iso_cut")),
      track_match_PV_(iConfig.getParameter<bool>("optionTrackMatchToPV")),
      dt_sig_track_(iConfig.getParameter<bool>("option_dtToTrack")),
      optionalPlots_(iConfig.getParameter<bool>("option_plots")),
      min_dR_cut(iConfig.getParameter<double>("min_dR_cut")),
      max_dR_cut(iConfig.getParameter<double>("max_dR_cut")),
      min_pt_cut_EB(iConfig.getParameter<double>("min_pt_cut_EB")),
      min_pt_cut_EE(iConfig.getParameter<double>("min_pt_cut_EE")),
      max_dz_cut_EB(iConfig.getParameter<double>("max_dz_cut_EB")),
      max_dz_cut_EE(iConfig.getParameter<double>("max_dz_cut_EE")),
      max_dz_vtx_cut(iConfig.getParameter<double>("max_dz_vtx_cut")),
      max_dxy_vtx_cut(iConfig.getParameter<double>("max_dxy_vtx_cut")),
//      min_strip_cut(iConfig.getParameter<double>("min_strip_cut")),
      min_track_mtd_mva_cut(iConfig.getParameter<double>("min_track_mtd_mva_cut")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecVertexToken_ =
      consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTag_vtx"));  // Vtx 4D collection

  MuonToken_ = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("inputMuon"));
  GlobalMuonTrk_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("globalMuonTrk"));  // It is unnecessary

  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));

  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));

  // test
  Rec4DVerToken_ = consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTag_vtx"));
  RecBeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBS"));
  sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  trackingVertexCollectionToken_ = consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  // test end

  // for ntuple
  tree_ = fs_->make<TTree>("muonIso", "muonIso");
  tree_->Branch("run_",			&run_);
  tree_->Branch("event_",		&event_);
  tree_->Branch("track_type_v1_",	&track_type_v1_);
  tree_->Branch("track_type_v2_",	&track_type_v2_);
  tree_->Branch("muon_pt_",		&muon_pt_);
  tree_->Branch("track_pt_",		&track_pt_);
  tree_->Branch("muon_time_",		&muon_time_);
  tree_->Branch("vtx_time_",		&vtx_time_);
  tree_->Branch("track_time_",		&track_time_);
  tree_->Branch("dtsig_muon_track_",	&dtsig_muon_track_);
  tree_->Branch("dtsig_vtx_track_",	&dtsig_vtx_track_);
  tree_->Branch("muon_PVweight_",	&muon_PVweight_);
  tree_->Branch("track_PVweight_",	&track_PVweight_);
  tree_->Branch("muon_time_err_",	&muon_time_err_);
  tree_->Branch("vtx_time_err_",	&vtx_time_err_);
  tree_->Branch("track_time_err_",	&track_time_err_);
  tree_->Branch("muon_prompt_",		&muon_prompt_);
  tree_->Branch("muon_isBarrel_",	&muon_isBarrel_);
  tree_->Branch("track_bx_",		&track_bx_);
  tree_->Branch("track_evtId_",		&track_evtId_);

  tree_->Branch("vtx_index_",		&vtx_index_);
  tree_->Branch("recovtx_sim_",		&recovtx_sim_);
  tree_->Branch("simvtx_reco_",		&simvtx_reco_);
  tree_->Branch("simvtx_bx_",		&simvtx_bx_);
  tree_->Branch("simvtx_evtId_",	&simvtx_evtId_);
  tree_->Branch("selectedVtxMatching_",	&selectedVtxMatching_);
  tree_->Branch("selectedLV_",		&selectedLV_);
  tree_->Branch("match_vtx_reco2sim_",	&match_vtx_reco2sim_);
  tree_->Branch("match_vtx_sim2reco_",	&match_vtx_sim2reco_);
  tree_->Branch("simvtx_pt_",		&simvtx_pt_);
  tree_->Branch("simvtx_ptsq_",		&simvtx_ptsq_);
  tree_->Branch("simvtx_nGenTrk_",	&simvtx_nGenTrk_);
  tree_->Branch("simvtx_num_matched_reco_tracks_", &simvtx_num_matched_reco_tracks_);
  tree_->Branch("recovtx_pt_",		&recovtx_pt_);
  tree_->Branch("recovtx_ptsq_",	&recovtx_ptsq_);
  tree_->Branch("recovtx_nRecoTrk_",	&recovtx_nRecoTrk_);
  tree_->Branch("recovtx_num_matched_sim_tracks_", &recovtx_num_matched_sim_tracks_);
}

MtdMuonIsoValidation::~MtdMuonIsoValidation() {}

// test
void MtdMuonIsoValidation::matchReco2Sim(std::vector<recoPrimaryVertex>& recopv,
                                              std::vector<simPrimaryVertex>& simpv,
                                              const edm::ValueMap<float>& sigmat0,
                                              const edm::ValueMap<float>& MVA,
                                              const edm::Handle<reco::BeamSpot>& BS) {
  for (auto vv : simpv) {
    vv.wnt.clear();
    vv.wos.clear();
  }
  for (auto rv : recopv) {
    rv.wnt.clear();
    rv.wos.clear();
  }

  for (unsigned int iv = 0; iv < recopv.size(); iv++) {
    const reco::Vertex* vertex = recopv.at(iv).recVtx;

    for (unsigned int iev = 0; iev < simpv.size(); iev++) {
      double wnt = 0;
      double wos = 0;
      double evwnt = 0;
      double evwos = 0;
      double evnt = 0;

      for (auto iTrack = vertex->tracks_begin(); iTrack != vertex->tracks_end(); ++iTrack) {
        if (vertex->trackWeight(*iTrack) < 0.5)
          continue;
        if (MVA[(*iTrack)] < 0.1)
          continue;

        auto tp_info = getMatchedTP(*iTrack, simpv.at(iev).sim_vertex).first;
        int matchCategory = getMatchedTP(*iTrack, simpv.at(iev).sim_vertex).second;
        // matched TP equal to any TP of a given sim vertex
        if (tp_info != nullptr && matchCategory == 0) {
          getWosWnt(*vertex, *iTrack, sigmat0, BS, wos, wnt);
          simpv.at(iev).addTrack(iv, wos, wnt);
          recopv.at(iv).addTrack(iev, wos, wnt);
          evwos += wos;
          evwnt += wnt;
          evnt++;
        }
      }  // RecoTracks loop

      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > recopv.at(iv).maxwos) && (evnt > 1)) {
        recopv.at(iv).wosmatch = iev;
        recopv.at(iv).maxwos = evwos;
        recopv.at(iv).maxwosnt = evnt;

        simpv.at(iev).wos_dominated_recv.push_back(iv);
        simpv.at(iev).nwosmatch++;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > recopv.at(iv).maxwnt)) {
        recopv.at(iv).wntmatch = iev;
        recopv.at(iv).maxwnt = evwnt;
      }
    }  // TrackingVertex loop

  }  // RecoPrimaryVertex

  // after filling infos, goes for the sim-reco match
  for (auto& vrec : recopv) {
    vrec.sim = NOT_MATCHED;
    vrec.matchQuality = 0;
  }
  unsigned int iev = 0;
  for (auto& vv : simpv) {
    if (0) {
      edm::LogPrint("MtdMuonIsoValidation") << "iev: " << iev;
      edm::LogPrint("MtdMuonIsoValidation") << "wos_dominated_recv.size: " << vv.wos_dominated_recv.size();
    }
    for (unsigned int i = 0; i < vv.wos_dominated_recv.size(); i++) {
      auto recov = vv.wos_dominated_recv.at(i);
      if (0) {
        edm::LogPrint("MtdMuonIsoValidation")
            << "index of reco vertex: " << recov << " that has a wos: " << vv.wos.at(recov) << " at position " << i;
      }
    }
    vv.rec = NOT_MATCHED;
    vv.matchQuality = 0;
    iev++;
  }
  // this tries a one-to-one match, taking simPV with highest wos if there are > 1 simPV candidates
  for (unsigned int rank = 1; rank < maxRank_; rank++) {
    for (unsigned int iev = 0; iev < simpv.size(); iev++) {  //loop on SimPV
      if (simpv.at(iev).rec != NOT_MATCHED)
        continue;
      if (simpv.at(iev).nwosmatch == 0)
        continue;
      if (simpv.at(iev).nwosmatch > rank)
        continue;
      unsigned int iv = NOT_MATCHED;
      for (unsigned int k = 0; k < simpv.at(iev).wos_dominated_recv.size(); k++) {
        unsigned int rec = simpv.at(iev).wos_dominated_recv.at(k);
        auto vrec = recopv.at(rec);
        if (vrec.sim != NOT_MATCHED)
          continue;  // already matched
        if (std::abs(simpv.at(iev).z - vrec.z) > zWosMatchMax_)
          continue;  // insanely far away
        if ((iv == NOT_MATCHED) || simpv.at(iev).wos.at(rec) > simpv.at(iev).wos.at(iv)) {
          iv = rec;
        }
      }
      if (iv !=
          NOT_MATCHED) {  // if the rec vertex has already been associated is possible that iv remains NOT_MATCHED at this point
        recopv.at(iv).sim = iev;
        simpv.at(iev).rec = iv;
        recopv.at(iv).matchQuality = rank;
        simpv.at(iev).matchQuality = rank;
      }
    }
  }
  //give vertices a chance that have a lot of overlap, but are still recognizably
  //caused by a specific simvertex (without being classified as dominating)
  //like a small peak sitting on the flank of a larger nearby peak
  unsigned int ntry = 0;
  while (ntry++ < maxTry_) {
    unsigned nmatch = 0;
    for (unsigned int iev = 0; iev < simpv.size(); iev++) {
      if ((simpv.at(iev).rec != NOT_MATCHED) || (simpv.at(iev).wos.empty()))
        continue;
      // find a rec vertex for the NOT_MATCHED sim vertex
      unsigned int rec = NOT_MATCHED;
      for (auto rv : simpv.at(iev).wos) {
        if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wos.at(rec))) {
          rec = rv.first;
        }
      }

      if (rec == NOT_MATCHED) {  // try with wnt match
        for (auto rv : simpv.at(iev).wnt) {
          if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wnt.at(rec))) {
            rec = rv.first;
          }
        }
      }
      if (rec == NOT_MATCHED)
        continue;
      if (recopv.at(rec).sim != NOT_MATCHED)
        continue;  // already gone

      // check if the recvertex can be  matched
      unsigned int rec2sim = NOT_MATCHED;
      for (auto sv : recopv.at(rec).wos) {
        if (simpv.at(sv.first).rec != NOT_MATCHED)
          continue;  // already used
        if ((rec2sim == NOT_MATCHED) || (sv.second > recopv.at(rec).wos.at(rec2sim))) {
          rec2sim = sv.first;
        }
      }
      if (iev == rec2sim) {
        // do the match and assign lowest quality (i.e. max rank)
        recopv.at(rec).sim = iev;
        recopv.at(rec).matchQuality = maxRank_;
        simpv.at(iev).rec = rec;
        simpv.at(iev).matchQuality = maxRank_;
        nmatch++;
      }
    }  // sim loop
    if (nmatch == 0) {
      break;
    }
  }  // ntry
}

void MtdMuonIsoValidation::getWosWnt(const reco::Vertex& recoVtx,
                                          const reco::TrackBaseRef& recoTrk,
                                          const edm::ValueMap<float>& sigmat0,
                                          const edm::Handle<reco::BeamSpot>& BS,
                                          double& wos,
                                          double& wnt) {
  double dz2_beam = pow((*BS).BeamWidthX() * cos(recoTrk->phi()) / tan(recoTrk->theta()), 2) +
                    pow((*BS).BeamWidthY() * sin(recoTrk->phi()) / tan(recoTrk->theta()), 2);
  double dz2 =
      pow(recoTrk->dzError(), 2) + dz2_beam + pow(0.0020, 2);  // added 20 um, some tracks have crazy small resolutions
  wos = recoVtx.trackWeight(recoTrk) / dz2;
  wnt = recoVtx.trackWeight(recoTrk) * std::min(recoTrk->pt(), 1.0);

  if (sigmat0[recoTrk] > 0) {
    double sigmaZ = (*BS).sigmaZ();
    double sigmaT = sigmaZ / c_;  // c in cm/ns
    wos = wos / erf(sigmat0[recoTrk] / sigmaT);
  }
}

std::pair<const edm::Ref<std::vector<TrackingParticle>>*, int> MtdMuonIsoValidation::getMatchedTP(
    const reco::TrackBaseRef& recoTrack, const TrackingVertexRef& vsim) {
//std::cout << "DEBUG 1" << std::endl;
  auto found = r2s__->find(recoTrack);

//std::cout << "DEBUG 2" << std::endl;
  // reco track not matched to any TP (fake tracks)
  if (found == r2s__->end())
    return std::make_pair(nullptr, -1);
//std::cout << "DEBUG 3" << std::endl;

  // matched TP equal to any TP of a given sim vertex
  for (const auto& tp : found->val) {
//std::cout << "DEBUG 4" << std::endl;
    if (std::find_if(vsim->daughterTracks_begin(), vsim->daughterTracks_end(), [&](const TrackingParticleRef& vtp) {
          return tp.first == vtp;
        }) != vsim->daughterTracks_end())
      return std::make_pair(&tp.first, 0);
    // matched TP not associated to any daughter track of a given sim vertex but having the same eventID (track from secondary vtx)
    else if (tp.first->eventId().bunchCrossing() == vsim->eventId().bunchCrossing() &&
             tp.first->eventId().event() == vsim->eventId().event()) {
      return std::make_pair(&tp.first, 1);
    }
    // matched TP not associated to any sim vertex of a given simulated event (PU track)
    else {
      return std::make_pair(&tp.first, 2);
    }
  }
//std::cout << "DEBUG 5" << std::endl;

  // reco track not matched to any TP from vertex
  return std::make_pair(nullptr, -1);
}
/* Extract information form TrackingParticles/TrackingVertex and fill
 * the helper class simPrimaryVertex with proper generation-level
 * information */
std::vector<MtdMuonIsoValidation::simPrimaryVertex> MtdMuonIsoValidation::getSimPVs(
    const edm::Handle<TrackingVertexCollection>& tVC) {
  std::vector<MtdMuonIsoValidation::simPrimaryVertex> simpv;
  int current_event = -1;
  int s = -1;
  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    // We keep only the first vertex from all the events at BX=0.
    if (v->eventId().bunchCrossing() != 0)
      continue;
    if (v->eventId().event() != current_event) {
      current_event = v->eventId().event();
    } else {
      continue;
    }
    s++;
    if (std::abs(v->position().z()) > 1000)
      continue;  // skip junk vertices

    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    simPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z(), v->position().t());
    sv.eventId = v->eventId();
    sv.sim_vertex = TrackingVertexRef(tVC, std::distance(tVC->begin(), v));
    sv.OriginalIndex = s;

    for (TrackingParticleRefVector::iterator iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end();
         ++iTrack) {
      assert((**iTrack).eventId().bunchCrossing() == 0);
    }
    simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) && (std::abs(sv.x - v0->x) < 1e-5) && (std::abs(sv.y - v0->y) < 1e-5) &&
          (std::abs(sv.z - v0->z) < 1e-5)) {
        vp = &(*v0);
        break;
      }
    }
    if (!vp) {
      // this is a new vertex, add it to the list of sim-vertices
      simpv.push_back(sv);
      vp = &simpv.back();
    }
// Loop over daughter track(s) as Tracking Particles
    for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
      auto momentum = (*(*iTP)).momentum();
      const reco::Track* matched_best_reco_track = nullptr;
      double match_quality = -1;
      if ((**iTP).charge() == 0)
        continue;
      if (s2r_->find(*iTP) != s2r_->end()) {
        matched_best_reco_track = (*s2r_)[*iTP][0].first.get();
        match_quality = (*s2r_)[*iTP][0].second;
      }

      vp->ptot.setPx(vp->ptot.x() + momentum.x());
      vp->ptot.setPy(vp->ptot.y() + momentum.y());
      vp->ptot.setPz(vp->ptot.z() + momentum.z());
      vp->ptot.setE(vp->ptot.e() + (**iTP).energy());
      vp->pt += (**iTP).pt();
      vp->ptsq += ((**iTP).pt() * (**iTP).pt());
      vp->nGenTrk++;

      if (matched_best_reco_track) {
        vp->num_matched_reco_tracks++;
        vp->average_match_quality += match_quality;
      }
    }  // End of for loop on daughters sim-particles
    if (vp->num_matched_reco_tracks)
      vp->average_match_quality /= static_cast<float>(vp->num_matched_reco_tracks);
    if (false) {
      edm::LogPrint("MtdMuonIsoValidation")
          << "average number of associated tracks: " << vp->num_matched_reco_tracks / static_cast<float>(vp->nGenTrk)
          << " with average quality: " << vp->average_match_quality;
    }
  }  // End of for loop on tracking vertices

  // In case of no simulated vertices, break here
  if (simpv.empty())
    return simpv;

  // Now compute the closest distance in z between all simulated vertex
  // first initialize
  auto prev_z = simpv.back().z;
  for (simPrimaryVertex& vsim : simpv) {
    vsim.closest_vertex_distance_z = std::abs(vsim.z - prev_z);
    prev_z = vsim.z;
  }
  // then calculate
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
    std::vector<simPrimaryVertex>::iterator vsim2 = vsim;
    vsim2++;
    for (; vsim2 != simpv.end(); vsim2++) {
      double distance = std::abs(vsim->z - vsim2->z);
      // need both to be complete
      vsim->closest_vertex_distance_z = std::min(vsim->closest_vertex_distance_z, distance);
      vsim2->closest_vertex_distance_z = std::min(vsim2->closest_vertex_distance_z, distance);
    }
  }
  return simpv;
}
/* Extract information form recoVertex and fill the helper class
 * recoPrimaryVertex with proper reco-level information */
std::vector<MtdMuonIsoValidation::recoPrimaryVertex> MtdMuonIsoValidation::getRecoPVs(
    const edm::Handle<edm::View<reco::Vertex>>& tVC) {
  std::vector<MtdMuonIsoValidation::recoPrimaryVertex> recopv;
  int r = -1;
  for (auto v = tVC->begin(); v != tVC->end(); ++v) {
    r++;
    // Skip junk vertices
    if (std::abs(v->z()) > 1000)
      continue;
    if (v->isFake() || !v->isValid())
      continue;

    recoPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z());
    sv.recVtx = &(*v);
    sv.recVtxRef = reco::VertexBaseRef(tVC, std::distance(tVC->begin(), v));

    sv.OriginalIndex = r;
    sv.ndof = v->ndof();
    // this is a new vertex, add it to the list of reco-vertices
    recopv.push_back(sv);
    MtdMuonIsoValidation::recoPrimaryVertex* vp = &recopv.back();

    // Loop over daughter track(s)
    for (auto iTrack = v->tracks_begin(); iTrack != v->tracks_end(); ++iTrack) {
      auto momentum = (*(*iTrack)).innerMomentum();
      if (momentum.mag2() == 0)
        momentum = (*(*iTrack)).momentum();
      vp->pt += std::sqrt(momentum.perp2());
      vp->ptsq += (momentum.perp2());
      vp->nRecoTrk++;

      auto matched = r2s__->find(*iTrack);
      if (matched != r2s__->end()) {
        vp->num_matched_sim_tracks++;
      }

    }  // End of for loop on daughters reconstructed tracks
  }    // End of for loop on tracking vertices

  // In case of no reco vertices, break here
  if (recopv.empty())
    return recopv;

  // Now compute the closest distance in z between all reconstructed vertex
  // first initialize
  auto prev_z = recopv.back().z;
  for (recoPrimaryVertex& vreco : recopv) {
    vreco.closest_vertex_distance_z = std::abs(vreco.z - prev_z);
    prev_z = vreco.z;
  }
  for (std::vector<recoPrimaryVertex>::iterator vreco = recopv.begin(); vreco != recopv.end(); vreco++) {
    std::vector<recoPrimaryVertex>::iterator vreco2 = vreco;
    vreco2++;
    for (; vreco2 != recopv.end(); vreco2++) {
      double distance = std::abs(vreco->z - vreco2->z);
      // need both to be complete
      vreco->closest_vertex_distance_z = std::min(vreco->closest_vertex_distance_z, distance);
      vreco2->closest_vertex_distance_z = std::min(vreco2->closest_vertex_distance_z, distance);
    }
  }
  return recopv;
}

// test end


// ------------ method called for each event  ------------
void MtdMuonIsoValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  // for ntuple
  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
//  track_type_v1_.clear(), track_type_v2_.clear();
//  muon_pt_.clear(), track_pt_.clear();
//  muon_time_.clear(), vtx_time_.clear(), track_time_.clear();
//  dtsig_muon_track_.clear(), dtsig_vtx_track_.clear();
//  muon_PVweight_.clear(), track_PVweight_.clear();
//  muon_time_err_.clear(), vtx_time_err_.clear(), track_time_err_.clear();
//  muon_prompt_.clear(), muon_isBarrel_.clear();
//  track_bx_.clear(), track_evtId_.clear();
//  selectedVtxMatching_.clear(), selectedLV_.clear(), match_vtx_reco2sim_.clear(), match_vtx_sim2reco_.clear();

  vtx_index_=-999, recovtx_sim_=-999, simvtx_reco_=-999, simvtx_bx_=-999, simvtx_evtId_=-99;
  simvtx_pt_=-999., simvtx_ptsq_=-999., recovtx_pt_=-999., recovtx_ptsq_=-999.;
  simvtx_nGenTrk_=-999, simvtx_num_matched_reco_tracks_=-999, recovtx_nRecoTrk_=-999, recovtx_num_matched_sim_tracks_=-999;


  auto GenRecTrackHandle = iEvent.getHandle(GenRecTrackToken_);

  auto VertexHandle = iEvent.getHandle(RecVertexToken_);
  std::vector<reco::Vertex> vertices = *VertexHandle;

  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);

  auto muHandle = makeValid(iEvent.getHandle(MuonToken_));
  reco::MuonCollection muColl = *(muHandle.product());
  auto GloMuonTrackHandle = iEvent.getHandle(GlobalMuonTrk_);

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  const reco::RecoToSimCollection* r2s_ = recoToSimH.product();

  // Creating muon collection
  std::vector<reco::Muon> localMuonCollection;
  for (const auto& mu_ : muColl) {
    if(mu_.passed(reco::Muon::CutBasedIdLoose)) {     // loose ID
      localMuonCollection.emplace_back(mu_);
    }
  }
  localMuonCollection.shrink_to_fit();

  // test
  // Check the fraction of loose_cut
  std::vector<reco::Muon> localMuonCollection_loose;
  for (const auto& mu_ : muColl) {
    if(mu_.passed(reco::Muon::CutBasedIdLoose)) localMuonCollection_loose.emplace_back(mu_);
  }
  localMuonCollection_loose.shrink_to_fit();
  // test end

  reco::Vertex Vtx_chosen;
//  reco::Vertex* aa;
  unsigned int vtx_index=0;
  // This part has to be included, because in ~1% of the events, the "good" vertex is the 1st one not the 0th one in the collection
  for (int iVtx = 0; iVtx < (int)vertices.size(); iVtx++) {
    const reco::Vertex& vertex = vertices.at(iVtx);
    if (!vertex.isFake() && vertex.ndof() >= 4) {
      Vtx_chosen = vertex;
      vtx_index = iVtx;
//      cout << endl;
//      cout << "iVtx: " << iVtx << endl;
      break;
    }
  }
  vtx_index_ = vtx_index;
//  aa = &Vtx_chosen;
//  cout << "vtx:   " << vertices.size() << endl;

  // test
  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  edm::Handle<reco::SimToRecoCollection> simToRecoH;
  iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
  if (simToRecoH.isValid())
    s2r_ = simToRecoH.product();
  else
    edm::LogWarning("MtdMuonIsoValidation") << "simToRecoH is not valid";
  r2s__ = recoToSimH.product();
  edm::Handle<reco::BeamSpot> BeamSpotH;
  iEvent.getByToken(RecBeamSpotToken_, BeamSpotH);
  if (!BeamSpotH.isValid())
    edm::LogWarning("MtdMuonIsoValidation") << "BeamSpotH is not valid";
//  if (recoToSimH.isValid())
//    r2s__ = recoToSimH.product();
//  else
//    edm::LogWarning("MtdMuonIsoValidation") << "recoToSimH is not valid";

  std::vector<simPrimaryVertex> simpv; // a list of simulated MC PVs
  simpv = getSimPVs(TVCollectionH);

  std::vector<recoPrimaryVertex> recopv;  // a list of reconstructed primary MC vertices
  edm::Handle<edm::View<reco::Vertex>> recVtxs;
  //iEvent.getByToken(RecVertexToken_, recVtxs);
  iEvent.getByToken(Rec4DVerToken_, recVtxs);
  recopv = getRecoPVs(recVtxs);

//  cout << "simpv: " << simpv.size() << endl;
//  cout << "recpv: " << vertices.size() << endl;

  const auto& sigmat0Safe = iEvent.get(sigmat0SafePidToken_);
  matchReco2Sim(recopv, simpv, sigmat0Safe, mtdQualMVA, BeamSpotH);

  // need to be check whether vtx has the highest pT
//  cout << "simpv.size(): " << simpv.size() << endl;
  
//  cout << "recopv.at(" << vtx_index << ").sim: " << recopv.at(vtx_index).sim << endl;
  recovtx_sim_ = recopv.at(vtx_index).sim;
  recovtx_pt_ = recopv.at(vtx_index).pt;
  recovtx_ptsq_ = recopv.at(vtx_index).ptsq;
  recovtx_nRecoTrk_ = recopv.at(vtx_index).nRecoTrk;
  recovtx_num_matched_sim_tracks_ = recopv.at(vtx_index).num_matched_sim_tracks;

//  simpv_reco_ = simpv.at(recopv.at(vtx_index).sim).rec;
  for (unsigned int iev=0; iev<simpv.size(); iev++) {
    if (recopv.at(vtx_index).sim == iev) {
//      cout << "simpv.at(" << iev << ").rec: " << simpv.at(iev).rec << endl;
//      cout << "recopv.at(" << simpv.at(iev).rec << ").sim: " << recopv.at(simpv.at(iev).rec).sim << endl;
//      cout << "[simpv] BX: " << simpv.at(iev).eventId.bunchCrossing() << ", evtId: " << simpv.at(iev).eventId.event() << endl;
//      printSimVtxRecoVtxInfo(simpv.at(iev), recopv.at(vtx_index));
      simvtx_reco_ = simpv.at(iev).rec;
      simvtx_bx_ = simpv.at(iev).eventId.bunchCrossing();
      simvtx_evtId_ = simpv.at(iev).eventId.event();
      simvtx_pt_ = simpv.at(iev).pt;
      simvtx_ptsq_ = simpv.at(iev).ptsq;
      simvtx_nGenTrk_ = simpv.at(iev).nGenTrk;
      simvtx_num_matched_reco_tracks_ = simpv.at(iev).num_matched_reco_tracks;
    }
//    if (recopv.at(vtx_index).sim == iev && simpv.at(iev).rec == vtx_index) {
//      cout << "Vtx Matching" << endl;
//      cout << "recopv.at(" << simpv.at(iev).rec << ").sim: " << recopv.at(simpv.at(iev).rec).sim << endl;
//      printSimVtxRecoVtxInfo(simpv.at(iev), recopv.at(vtx_index));
//    }

//    if (simpv.at(iev).eventId.bunchCrossing()==0 && simpv.at(iev).eventId.event()==0) {
//      cout << "BX/evtId matching for simpv" << endl;
//      cout << "iev: " << iev << endl;
//      cout << "simpv.at(" << iev << ").rec: " << simpv.at(iev).rec << endl;
//      if(simpv.at(iev).rec!=66666) cout << "recopv.at(" << simpv.at(iev).rec << ").sim: " << recopv.at(simpv.at(iev).rec).sim << endl;
//      printSimVtxRecoVtxInfo(simpv.at(iev), recopv.at(vtx_index));
//    }
  }
/*
  for (unsigned int iv=0; iv<recopv.size(); iv++) {
    if (recopv.at(iv).recVtx->t() == Vtx_chosen.t()) cout << "iv" << iv << " vtx is found." << endl;
    if (recopv.at(iv).recVtx == aa) cout << "iv: " << iv << " vtx is found2." << endl;
    //if (recopv.at(iv) == aa) cout << "[" << iv << ", " << iev << "] vtx is found3." << endl;
    for (unsigned int iev=0; iev<simpv.size(); iev++) {
      bool selectedVtxMatching = recopv.at(iv).sim == iev && simpv.at(iev).rec == iv;
      bool selectedLV = simpv.at(iev).eventId.bunchCrossing() == 0 && simpv.at(iev).eventId.event() == 0 &&
                          recopv.at(iv).OriginalIndex == 0;
      bool selectedLVMatching = selectedVtxMatching && selectedLV;
      if (selectedVtxMatching == true) cout << "v1" << endl;
      if (selectedLV == true) cout << "v2" << endl;
      if (selectedLVMatching == true) cout << "v3" << endl;
//      cout << "[" << iv << ", " << iev << "] recopv.at(iv).sim: " << recopv.at(iv).sim << endl;
//      cout << "[" << iv << ", " << iev << "] simpv.at(iev).rec: " << simpv.at(iev).rec << endl;
    }
  }
*/
 
  // test end

  auto pdgCheck = [](int pdg) {
    pdg = std::abs(pdg);
    return (pdg == 23 or pdg == 24 or pdg == 15 or pdg == 13);
  };

  int nmuons_Sig=0, nmuons_Sig_EB=0, nmuons_Sig_EE=0;
  int nmuons_Bkg=0, nmuons_Bkg_EB=0, nmuons_Bkg_EE=0;

  for (const auto& muon : localMuonCollection) {
  // for ntuple
  track_type_v1_.clear(), track_type_v2_.clear();
  muon_pt_.clear(), track_pt_.clear();
  muon_time_.clear(), vtx_time_.clear(), track_time_.clear();
  dtsig_muon_track_.clear(), dtsig_vtx_track_.clear();
  muon_PVweight_.clear(), track_PVweight_.clear();
  muon_time_err_.clear(), vtx_time_err_.clear(), track_time_err_.clear();
  muon_prompt_.clear(), muon_isBarrel_.clear();
  track_bx_.clear(), track_evtId_.clear();
  selectedVtxMatching_.clear(), selectedLV_.clear(), match_vtx_reco2sim_.clear(), match_vtx_sim2reco_.clear();

    bool muon_Prompt = false;
    float muon_track_source_dz = std::abs(muon.track()->dz(Vtx_chosen.position()));
    float muon_track_source_dxy = std::abs(muon.track()->dxy(Vtx_chosen.position()));

    // selecting "good" RECO muons
    // PARAM

    // test
    if (muon.track()->pt() >= 10 || std::abs(muon.track()->eta()) <= 2.4) {
      meMuonISO_dz_muon_->Fill(muon_track_source_dz);
      meMuonISO_dxy_muon_->Fill(muon_track_source_dxy);
      if (std::abs(muon.track()->eta()) < 1.5) {
        meMuonISO_dz_muon_EB_->Fill(muon_track_source_dz);
        meMuonISO_dxy_muon_EB_->Fill(muon_track_source_dxy);
      }
      else {
        meMuonISO_dz_muon_EE_->Fill(muon_track_source_dz);
        meMuonISO_dxy_muon_EE_->Fill(muon_track_source_dxy);
      }
    }
    // test end
    
    if (muon.track()->pt() < 10 || std::abs(muon.track()->eta()) > 2.4 || muon_track_source_dz > max_dz_vtx_cut || muon_track_source_dxy > max_dxy_vtx_cut)
      continue;

    bool Barrel_muon = 0;
    if(std::abs(muon.track()->eta()) < 1.5) Barrel_muon = 1;
    const reco::TrackRef muon_SigTrkRef = muon.track();

    double tsim_muon = -1.;
    double muon_sim_pt = -1.;
    double muon_sim_phi = -1.;
    double muon_sim_eta = -1.;

    // association with tracking particle to have sim info
    const reco::TrackBaseRef trkrefb(muon_SigTrkRef);
    auto found = r2s_->find(trkrefb);
    if (found != r2s_->end()) {
      const auto& tp = (found->val)[0];
      tsim_muon = (tp.first)->parentVertex()->position().t() * 1e9;
      muon_sim_pt = (tp.first)->pt();
      muon_sim_phi = (tp.first)->phi();
      muon_sim_eta = (tp.first)->eta();
      // check that the genParticle vector is not empty
      if (tp.first->status() != -99) {
        const auto genParticle = *(tp.first->genParticles()[0]);
        // check if prompt (not from hadron, muon, or tau decay) and final state
        // or if is a direct decay product of a prompt tau and is final state
        if ((genParticle.isPromptFinalState() or genParticle.isDirectPromptTauDecayProductFinalState()) and
            pdgCheck(genParticle.mother()->pdgId())) {
          muon_Prompt = true;
          // TODO get simtrackster from mtd, simtrack to tp and check that a recocluster was there
        }
      }
    }
    // test
      // Prompt
    if(muon_Prompt) {
      meMuonISO_mva_muon_reco_Sig_->Fill(mtdQualMVA[muon_SigTrkRef]);
      if (Barrel_muon) {
        if(t0Pid[muon_SigTrkRef]!=0) meMuonISO_time_muon_reco_Sig_EB_->Fill(t0Pid[muon_SigTrkRef]);
        meMuonISO_tErr_muon_reco_Sig_EB_->Fill(Sigmat0Pid[muon_SigTrkRef]);
        meMuonISO_mva_muon_reco_Sig_EB_->Fill(mtdQualMVA[muon_SigTrkRef]);
        // Check whether muon has observables below or not
          // time_reco
        if(t0Pid[muon_SigTrkRef]==0) meMuonISO_has_time_muon_reco_Sig_EB_->Fill(0);
        else meMuonISO_has_time_muon_reco_Sig_EB_->Fill(1);
          // tErr_reco
        if(Sigmat0Pid[muon_SigTrkRef]==-1) meMuonISO_has_tErr_muon_reco_Sig_EB_->Fill(0);
        else meMuonISO_has_tErr_muon_reco_Sig_EB_->Fill(1);
          // mva_reco
        if(mtdQualMVA[muon_SigTrkRef]==-1) meMuonISO_has_mva_muon_reco_Sig_EB_->Fill(0);
        else meMuonISO_has_mva_muon_reco_Sig_EB_->Fill(1);
          // time_sim
        if(tsim_muon!=-1) meMuonISO_time_muon_sim_Sig_EB_->Fill(tsim_muon);
        if(tsim_muon==-1) meMuonISO_has_time_muon_sim_Sig_EB_->Fill(0);
        else meMuonISO_has_time_muon_sim_Sig_EB_->Fill(1);
      }
      else {
        if(t0Pid[muon_SigTrkRef]!=0) meMuonISO_time_muon_reco_Sig_EE_->Fill(t0Pid[muon_SigTrkRef]);
        meMuonISO_tErr_muon_reco_Sig_EE_->Fill(Sigmat0Pid[muon_SigTrkRef]);
        meMuonISO_mva_muon_reco_Sig_EE_->Fill(mtdQualMVA[muon_SigTrkRef]);
        // Check whether muon has observables below or not
          // time_reco
        if(t0Pid[muon_SigTrkRef]==0) meMuonISO_has_time_muon_reco_Sig_EE_->Fill(0);
        else meMuonISO_has_time_muon_reco_Sig_EE_->Fill(1);
          // tErr_reco
        if(Sigmat0Pid[muon_SigTrkRef]==-1) meMuonISO_has_tErr_muon_reco_Sig_EE_->Fill(0);
        else meMuonISO_has_tErr_muon_reco_Sig_EE_->Fill(1);
          // mva_reco
        if(mtdQualMVA[muon_SigTrkRef]==-1) meMuonISO_has_mva_muon_reco_Sig_EE_->Fill(0);
        else meMuonISO_has_mva_muon_reco_Sig_EE_->Fill(1);
          // time_sim
        if(tsim_muon!=-1) meMuonISO_time_muon_sim_Sig_EE_->Fill(tsim_muon);
        if(tsim_muon==-1) meMuonISO_has_time_muon_sim_Sig_EE_->Fill(0);
        else meMuonISO_has_time_muon_sim_Sig_EE_->Fill(1);
      }
      // pTdiff between reco muon and tracker muon
      meMuonISO_pTdiff_reco_tracker_muon_Sig_->Fill(std::abs(muon.pt() - muon.track()->pt()));
    }
    else {
      meMuonISO_mva_muon_reco_Bkg_->Fill(mtdQualMVA[muon_SigTrkRef]);
      if (Barrel_muon) {
        if(t0Pid[muon_SigTrkRef]!=0) meMuonISO_time_muon_reco_Bkg_EB_->Fill(t0Pid[muon_SigTrkRef]);
        meMuonISO_tErr_muon_reco_Bkg_EB_->Fill(Sigmat0Pid[muon_SigTrkRef]);
        meMuonISO_mva_muon_reco_Bkg_EB_->Fill(mtdQualMVA[muon_SigTrkRef]);
          // time_reco
        if(t0Pid[muon_SigTrkRef]==0) meMuonISO_has_time_muon_reco_Bkg_EB_->Fill(0);
        else meMuonISO_has_time_muon_reco_Bkg_EB_->Fill(1);
          // tErr_reco
        if(Sigmat0Pid[muon_SigTrkRef]==-1) meMuonISO_has_tErr_muon_reco_Bkg_EB_->Fill(0);
        else meMuonISO_has_tErr_muon_reco_Bkg_EB_->Fill(1);
          // mva_reco
        if(mtdQualMVA[muon_SigTrkRef]==-1) meMuonISO_has_mva_muon_reco_Bkg_EB_->Fill(0);
        else meMuonISO_has_mva_muon_reco_Bkg_EB_->Fill(1);
          // time_sim
        if(tsim_muon!=-1) meMuonISO_time_muon_sim_Bkg_EB_->Fill(tsim_muon);
        if(tsim_muon==-1) meMuonISO_has_time_muon_sim_Bkg_EB_->Fill(0);
        else meMuonISO_has_time_muon_sim_Bkg_EB_->Fill(1);
      }
      else {
        if(t0Pid[muon_SigTrkRef]!=0) meMuonISO_time_muon_reco_Bkg_EE_->Fill(t0Pid[muon_SigTrkRef]);
        meMuonISO_tErr_muon_reco_Bkg_EE_->Fill(Sigmat0Pid[muon_SigTrkRef]);
        meMuonISO_mva_muon_reco_Bkg_EE_->Fill(mtdQualMVA[muon_SigTrkRef]);
          // time_reco
        if(t0Pid[muon_SigTrkRef]==0) meMuonISO_has_time_muon_reco_Bkg_EE_->Fill(0);
        else meMuonISO_has_time_muon_reco_Bkg_EE_->Fill(1);
          // tErr_reco
        if(Sigmat0Pid[muon_SigTrkRef]==-1) meMuonISO_has_tErr_muon_reco_Bkg_EE_->Fill(0);
        else meMuonISO_has_tErr_muon_reco_Bkg_EE_->Fill(1);
          // mva_reco
        if(mtdQualMVA[muon_SigTrkRef]==-1) meMuonISO_has_mva_muon_reco_Bkg_EE_->Fill(0);
        else meMuonISO_has_mva_muon_reco_Bkg_EE_->Fill(1);
          // time_sim
        if(tsim_muon!=-1) meMuonISO_time_muon_sim_Bkg_EE_->Fill(tsim_muon);
        if(tsim_muon==-1) meMuonISO_has_time_muon_sim_Bkg_EE_->Fill(0);
        else meMuonISO_has_time_muon_sim_Bkg_EE_->Fill(1);
      }
      // pTdiff between reco muon and tracker muon
      meMuonISO_pTdiff_reco_tracker_muon_Bkg_->Fill(std::abs(muon.pt() - muon.track()->pt()));
    }
    // pT of sim muon
    meMuonISO_pT_muon_sim_->Fill(muon_sim_pt);


    math::XYZVector MuonSigTrackMomentumAtVtx = muon.track()->momentum();

    double muon_sigTrkTime = -1;
    double muon_sigTrkTimeErr = -1;
    double muon_sigTrkMtdMva = -1;

    // if we found a track-matching, we add MTD timing information for it
    if (muon_SigTrkRef.isNonnull()) {
      // track pT/dz cuts
      float min_pt_cut = Barrel_muon ? min_pt_cut_EB : min_pt_cut_EE;
      float max_dz_cut = Barrel_muon ? max_dz_cut_EB : max_dz_cut_EE;

      muon_sigTrkTime = t0Pid[muon_SigTrkRef];
      muon_sigTrkMtdMva = mtdQualMVA[muon_SigTrkRef];
      muon_sigTrkTimeErr = (muon_sigTrkMtdMva > min_track_mtd_mva_cut) ? Sigmat0Pid[muon_SigTrkRef] : -1;

      meMuon_avg_error_SigTrk_check_->Fill(muon_sigTrkTimeErr);

      if (muon_Prompt) {
	nmuons_Sig++;
        // For signal (prompt)
        if (Barrel_muon) {
          // All selected muon information for efficiency plots later
          meMuon_pt_tot_Sig_EB_->Fill(muon.track()->pt());
          if(muon_sim_pt!=-1) meMuon_pt_sim_tot_Sig_EB_->Fill(muon_sim_pt);
          meMuon_eta_tot_Sig_EB_->Fill(std::abs(muon.track()->eta()));
          meMuon_phi_tot_Sig_EB_->Fill(muon.track()->phi());
          nmuons_Sig_EB++;
        } else {
          // All selected muon information for efficiency plots later
          meMuon_pt_tot_Sig_EE_->Fill(muon.track()->pt());
          if(muon_sim_pt!=-1) meMuon_pt_sim_tot_Sig_EE_->Fill(muon_sim_pt);
          meMuon_eta_tot_Sig_EE_->Fill(std::abs(muon.track()->eta()));
          meMuon_phi_tot_Sig_EE_->Fill(muon.track()->phi());
          nmuons_Sig_EE++;
        }
      } else {
        // For background (non-promt)
	nmuons_Bkg++;
        if (Barrel_muon) {
          meMuon_pt_tot_Bkg_EB_->Fill(muon.track()->pt());
          if(muon_sim_pt!=-1) meMuon_pt_sim_tot_Bkg_EB_->Fill(muon_sim_pt);
          meMuon_eta_tot_Bkg_EB_->Fill(std::abs(muon.track()->eta()));
          meMuon_phi_tot_Bkg_EB_->Fill(muon.track()->phi());
          nmuons_Bkg_EB++;
        } else {
          meMuon_pt_tot_Bkg_EE_->Fill(muon.track()->pt());
          if(muon_sim_pt!=-1) meMuon_pt_sim_tot_Bkg_EE_->Fill(muon_sim_pt);
          meMuon_eta_tot_Bkg_EE_->Fill(std::abs(muon.track()->eta()));
          meMuon_phi_tot_Bkg_EE_->Fill(muon.track()->phi());
          nmuons_Bkg_EE++;
        }
      }

      int N_tracks_noMTD = 0;
      double pT_sum_noMTD = 0;
      double rel_pT_sum_noMTD = 0;
      std::vector<int> N_tracks_MTD{0, 0, 0, 0, 0, 0, 0};
      std::vector<double> pT_sum_MTD{0, 0, 0, 0, 0, 0, 0};
      std::vector<double> rel_pT_sum_MTD{0, 0, 0, 0, 0, 0, 0};

      std::vector<int> N_tracks_sim_MTD{0, 0, 0, 0, 0, 0, 0};
      std::vector<double> pT_sum_sim_MTD{0, 0, 0, 0, 0, 0, 0};
      std::vector<double> rel_pT_sum_sim_MTD{0, 0, 0, 0, 0, 0, 0};
      int N_tracks_gen = 0;
      double pT_sum_gen = 0;
      double rel_pT_sum_gen = 0;

      std::vector<int> N_tracks_MTD_significance{0, 0, 0};
      std::vector<double> pT_sum_MTD_significance{0, 0, 0};
      std::vector<double> rel_pT_sum_MTD_significance{0, 0, 0};

      std::vector<int> N_tracks_sim_MTD_significance{0, 0, 0};
      std::vector<double> pT_sum_sim_MTD_significance{0, 0, 0};
      std::vector<double> rel_pT_sum_sim_MTD_significance{0, 0, 0};
      // test
      std::vector<double> pT_sum_gen_MTD_significance{0, 0, 0};
      std::vector<double> rel_pT_sum_gen_MTD_significance{0, 0, 0};

      // for ntuple
//      muon_isBarrel_.emplace_back(Barrel_muon);
//      muon_prompt_.emplace_back(muon_Prompt);
//      muon_pt_.emplace_back(muon_SigTrkRef->pt());
//      muon_time_.emplace_back(muon_sigTrkTime);
//      muon_time_err_.emplace_back(muon_sigTrkTimeErr);
//      muon_PVweight_.emplace_back(Vtx_chosen.trackWeight(muon_SigTrkRef));
//      vtx_time_.emplace_back(Vtx_chosen.t());
//      vtx_time_err_.emplace_back(Vtx_chosen.tError());


      int general_index = 0;

      for (const auto& trackGen : *GenRecTrackHandle) {
        const reco::TrackRef trackref_general(GenRecTrackHandle, general_index);
        general_index++;

        // Skip muon track
        if (trackref_general == muon_SigTrkRef) {
          continue;
        }

        if (trackGen.pt() < min_pt_cut) {
          continue;
        }

        if (std::abs(trackGen.vz() - muon.track()->vz()) > max_dz_cut) {
          continue;
        }

        // cut for general track matching to PV
        if (track_match_PV_) {
          if (Vtx_chosen.trackWeight(trackref_general) < 0.5) {
            continue;
          }
        }

        double dR = reco::deltaR(trackGen.momentum(), MuonSigTrackMomentumAtVtx);

        // restrict to tracks in the isolation cone
        if (dR < min_dR_cut || dR > max_dR_cut) {
          continue;
        }

         // for ntuple
        muon_isBarrel_.emplace_back(Barrel_muon);
        muon_prompt_.emplace_back(muon_Prompt);
        muon_pt_.emplace_back(muon_SigTrkRef->pt());
        muon_time_.emplace_back(muon_sigTrkTime);
        muon_time_err_.emplace_back(muon_sigTrkTimeErr);
        muon_PVweight_.emplace_back(Vtx_chosen.trackWeight(muon_SigTrkRef));
        vtx_time_.emplace_back(Vtx_chosen.t());
        vtx_time_err_.emplace_back(Vtx_chosen.tError());

        // no MTD case
        ++N_tracks_noMTD;
        pT_sum_noMTD += trackGen.pt();

        // MTD case
        const reco::TrackBaseRef trkrefBase(trackref_general);
        auto TPmatched = r2s_->find(trkrefBase);
        double tsim_trk = -1.;
        double trk_ptSim = -1.;
        bool genMatched = false;
	// test to checking the type of tracks
	// 0: tracks from PV, 1: PU tracks, 2: fake tracks
	if (TPmatched == r2s_->end()) {
	  track_type_v1_.emplace_back(2);
	  if (muon_Prompt) {
	    if (Barrel_muon) meMuonISO_trk_type_Sig_EB_->Fill(2);
	    else meMuonISO_trk_type_Sig_EE_->Fill(2);
	  }
	  else {
	    if (Barrel_muon) meMuonISO_trk_type_Bkg_EB_->Fill(2);
	    else meMuonISO_trk_type_Bkg_EE_->Fill(2);
	  }
	}
	else if (TPmatched != r2s_->end()) {
          // reco track matched to a TP
          const auto& tp = (TPmatched->val)[0];
          tsim_trk = (tp.first)->parentVertex()->position().t() * 1e9;
          trk_ptSim = (tp.first)->pt();
          // check that the genParticle vector is not empty
          if (tp.first->status() != -99) {
            genMatched = true;
            meTrk_genMatch_check_->Fill(1);
	    // test to checking the type of tracks
	    if (muon_Prompt) {
	      if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
	        track_type_v1_.emplace_back(0);
	        if (Barrel_muon) meMuonISO_trk_type_Sig_EB_->Fill(0);
		else meMuonISO_trk_type_Sig_EE_->Fill(0);
	      }
	      // there are no else cases
	      else {
	        track_type_v1_.emplace_back(1);
		if (Barrel_muon) meMuonISO_trk_type_Sig_EB_->Fill(1);
		else meMuonISO_trk_type_Sig_EE_->Fill(1);
	      }
	    }
	    else {
	      if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
	        track_type_v1_.emplace_back(0);
                if (Barrel_muon) meMuonISO_trk_type_Bkg_EB_->Fill(0);
		else meMuonISO_trk_type_Bkg_EE_->Fill(0);
              }
	      // there are no else cases
              else {
	        track_type_v1_.emplace_back(1);
                if (Barrel_muon) meMuonISO_trk_type_Bkg_EB_->Fill(1);
                else meMuonISO_trk_type_Bkg_EE_->Fill(1);
              }
	    }
          } else {
	    track_type_v1_.emplace_back(1);
            meTrk_genMatch_check_->Fill(0);
	    if (muon_Prompt) {
	      if (Barrel_muon) meMuonISO_trk_type_Sig_EB_->Fill(1);
	      else meMuonISO_trk_type_Sig_EE_->Fill(1);
	    }
	    else {
	      if (Barrel_muon) meMuonISO_trk_type_Bkg_EB_->Fill(1);
	      else meMuonISO_trk_type_Bkg_EE_->Fill(1);
	    }
          }
        }

        // test
        if(muon_Prompt) {
    	  meMuonISO_mva_trk_reco_Sig_->Fill(mtdQualMVA[trackref_general]);
	  if (Barrel_muon) {
    	    if(t0Pid[trackref_general]!=0) meMuonISO_time_trk_reco_Sig_EB_->Fill(t0Pid[trackref_general]);
    	    meMuonISO_tErr_trk_reco_Sig_EB_->Fill(Sigmat0Pid[trackref_general]);
    	    meMuonISO_mva_trk_reco_Sig_EB_->Fill(mtdQualMVA[trackref_general]);
            // Check whether track has observables below or not
  	      // time_reco
  	    if(t0Pid[trackref_general]==0) meMuonISO_has_time_trk_reco_Sig_EB_->Fill(0);
  	    else meMuonISO_has_time_trk_reco_Sig_EB_->Fill(1);
              // tErr_reco
  	    if(Sigmat0Pid[trackref_general]==-1) meMuonISO_has_tErr_trk_reco_Sig_EB_->Fill(0);
  	    else meMuonISO_has_tErr_trk_reco_Sig_EB_->Fill(1);
              // mva_reco
            if(mtdQualMVA[trackref_general]==-1) meMuonISO_has_mva_trk_reco_Sig_EB_->Fill(0);
            else meMuonISO_has_mva_trk_reco_Sig_EB_->Fill(1);
              // time_sim
            if(tsim_trk!=-1) meMuonISO_time_trk_sim_Sig_EB_->Fill(tsim_trk);
            if(tsim_trk==-1) meMuonISO_has_time_trk_sim_Sig_EB_->Fill(0);
            else meMuonISO_has_time_trk_sim_Sig_EB_->Fill(1);

            // Check on tracks matched with GenParticles
            meMuonISO_trk_genMatched_Sig_EB_->Fill(genMatched);
	  }
	  else {
    	    if(t0Pid[trackref_general]!=0) meMuonISO_time_trk_reco_Sig_EE_->Fill(t0Pid[trackref_general]);
    	    meMuonISO_tErr_trk_reco_Sig_EE_->Fill(Sigmat0Pid[trackref_general]);
    	    meMuonISO_mva_trk_reco_Sig_EE_->Fill(mtdQualMVA[trackref_general]);
            // Check whether track has observables below or not
  	      // time_reco
  	    if(t0Pid[trackref_general]==0) meMuonISO_has_time_trk_reco_Sig_EE_->Fill(0);
  	    else meMuonISO_has_time_trk_reco_Sig_EE_->Fill(1);
              // tErr_reco
  	    if(Sigmat0Pid[trackref_general]==-1) meMuonISO_has_tErr_trk_reco_Sig_EE_->Fill(0);
  	    else meMuonISO_has_tErr_trk_reco_Sig_EE_->Fill(1);
              // mva_reco
            if(mtdQualMVA[trackref_general]==-1) meMuonISO_has_mva_trk_reco_Sig_EE_->Fill(0);
            else meMuonISO_has_mva_trk_reco_Sig_EE_->Fill(1);
              // time_sim
            if(tsim_trk!=-1) meMuonISO_time_trk_sim_Sig_EE_->Fill(tsim_trk);
            if(tsim_trk==-1) meMuonISO_has_time_trk_sim_Sig_EE_->Fill(0);
            else meMuonISO_has_time_trk_sim_Sig_EE_->Fill(1);

            // Check on tracks matched with GenParticles
	    meMuonISO_trk_genMatched_Sig_EE_->Fill(genMatched);
	  }
	  meMuonISO_trk_genMatched_Sig_->Fill(genMatched);
	}
	else {
	  meMuonISO_mva_trk_reco_Bkg_->Fill(mtdQualMVA[trackref_general]);
	  if (Barrel_muon) {
	    if(t0Pid[trackref_general]!=0) meMuonISO_time_trk_reco_Bkg_EB_->Fill(t0Pid[trackref_general]);
	    meMuonISO_tErr_trk_reco_Bkg_EB_->Fill(Sigmat0Pid[trackref_general]);
	    meMuonISO_mva_trk_reco_Bkg_EB_->Fill(mtdQualMVA[trackref_general]);
	      // time_reco
	    if(t0Pid[trackref_general]==0) meMuonISO_has_time_trk_reco_Bkg_EB_->Fill(0);
	    else meMuonISO_has_time_trk_reco_Bkg_EB_->Fill(1);
              // tErr_reco
	    if(Sigmat0Pid[trackref_general]==-1) meMuonISO_has_tErr_trk_reco_Bkg_EB_->Fill(0);
	    else meMuonISO_has_tErr_trk_reco_Bkg_EB_->Fill(1);
              // mva_reco
            if(mtdQualMVA[trackref_general]==-1) meMuonISO_has_mva_trk_reco_Bkg_EB_->Fill(0);
            else meMuonISO_has_mva_trk_reco_Bkg_EB_->Fill(1);
              // time_sim
            if(tsim_trk!=-1) meMuonISO_time_trk_sim_Bkg_EB_->Fill(tsim_trk);
            if(tsim_trk==-1) meMuonISO_has_time_trk_sim_Bkg_EB_->Fill(0);
            else meMuonISO_has_time_trk_sim_Bkg_EB_->Fill(1);

            // Check on tracks matched with GenParticles
            meMuonISO_trk_genMatched_Bkg_EB_->Fill(genMatched);
	  }
	  else {
	    if(t0Pid[trackref_general]!=0) meMuonISO_time_trk_reco_Bkg_EE_->Fill(t0Pid[trackref_general]);
	    meMuonISO_tErr_trk_reco_Bkg_EE_->Fill(Sigmat0Pid[trackref_general]);
	    meMuonISO_mva_trk_reco_Bkg_EE_->Fill(mtdQualMVA[trackref_general]);
	      // time_reco
	    if(t0Pid[trackref_general]==0) meMuonISO_has_time_trk_reco_Bkg_EE_->Fill(0);
	    else meMuonISO_has_time_trk_reco_Bkg_EE_->Fill(1);
              // tErr_reco
	    if(Sigmat0Pid[trackref_general]==-1) meMuonISO_has_tErr_trk_reco_Bkg_EE_->Fill(0);
	    else meMuonISO_has_tErr_trk_reco_Bkg_EE_->Fill(1);
              // mva_reco
            if(mtdQualMVA[trackref_general]==-1) meMuonISO_has_mva_trk_reco_Bkg_EE_->Fill(0);
            else meMuonISO_has_mva_trk_reco_Bkg_EE_->Fill(1);
              // time_sim
            if(tsim_trk!=-1) meMuonISO_time_trk_sim_Bkg_EE_->Fill(tsim_trk);
            if(tsim_trk==-1) meMuonISO_has_time_trk_sim_Bkg_EE_->Fill(0);
            else meMuonISO_has_time_trk_sim_Bkg_EE_->Fill(1);

            // Check on tracks matched with GenParticles
	    meMuonISO_trk_genMatched_Bkg_EE_->Fill(genMatched);
	  }
	  meMuonISO_trk_genMatched_Bkg_->Fill(genMatched);
	}
        // pT of sim track
        meMuonISO_pT_trk_sim_->Fill(trk_ptSim);

        double TrkMTDTime = t0Pid[trackref_general];
        double TrkMTDMva = mtdQualMVA[trackref_general];
        double TrkMTDTimeErr = (TrkMTDMva > min_track_mtd_mva_cut) ? Sigmat0Pid[trackref_general] : -1;

        meMuon_avg_error_PUTrk_check_->Fill(TrkMTDTimeErr);

	// test
	bool match_vtx_reco2sim = false;
	bool match_vtx_sim2reco = false;
	bool selectedVtxMatching = false;
	bool selectedLV = false;

	bool check_match_vtx_reco2sim = false;
	bool check_match_vtx_sim2reco = false;
	bool check_selectedVtxMatching = false;
	bool check_selectedLV = false;
	const reco::TrackBaseRef trkrefBase2(trackref_general);
	cout << "AAA" << endl;
  	for (unsigned int iev=0; iev < simpv.size(); iev++) {
    	  auto vsim = simpv.at(iev).sim_vertex;
	  match_vtx_reco2sim = (recopv.at(vtx_index).sim == iev);
	  match_vtx_sim2reco = (simpv.at(iev).rec == vtx_index);
	  selectedVtxMatching = match_vtx_reco2sim && match_vtx_sim2reco;
	  //selectedVtxMatching = (recopv.at(vtx_index).sim == iev) && (simpv.at(iev).rec == vtx_index);
	  selectedLV = (simpv.at(iev).eventId.bunchCrossing() == 0) && (simpv.at(iev).eventId.event() == 0) && (recopv.at(vtx_index).OriginalIndex == 0);

	  if (match_vtx_reco2sim) check_match_vtx_reco2sim = true;
	  if (match_vtx_sim2reco) check_match_vtx_sim2reco = true;
	  if (selectedVtxMatching) check_selectedVtxMatching = true;
	  if (selectedLV) check_selectedLV = true;

	  //if (selectedVtxMatching == true) 
	  if (match_vtx_reco2sim) {
	    cout << "BBB" << endl;
//	    match_vtx_sim2reco_.emplace_back(match_vtx_sim2reco);
//	    selectedVtxMatching_.emplace_back(selectedVtxMatching);
//	    selectedLV_.emplace_back(selectedLV);
//            cout << "reco-sim vtx matching" << endl;
      	    auto tp_info = getMatchedTP(trkrefBase2, vsim).first;
  	    if (tp_info != nullptr) {
//  	      cout << "(*tp_info)->eventId().bunchCrossing(): "<< (*tp_info)->eventId().bunchCrossing() << endl;
//  	      cout << "(*tp_info)->eventId().event(): "<< (*tp_info)->eventId().event() << endl;
  	      track_bx_.emplace_back((*tp_info)->eventId().bunchCrossing());
  	      track_evtId_.emplace_back((*tp_info)->eventId().event());
  	    }
  	    else if (tp_info == nullptr) {
  	      track_bx_.emplace_back(-999);
  	      track_evtId_.emplace_back(-999);
  	    }
  	    // test to checking the type of tracks
  	    // 0: tracks from PV, 1: PU tracks, 2: fake tracks, 3: tracks from secondary vertices
            int matchCategory = getMatchedTP(trkrefBase2, vsim).second;
  	    if (matchCategory == -1) {
  	      track_type_v2_.emplace_back(2);
  	      if (muon_Prompt) {
  	        if (Barrel_muon) meMuonISO_trk_type_v2_Sig_EB_->Fill(2);
  	        else meMuonISO_trk_type_v2_Sig_EE_->Fill(2);
  	      }
  	      else {
  	        if (Barrel_muon) meMuonISO_trk_type_v2_Bkg_EB_->Fill(2);
                else meMuonISO_trk_type_v2_Bkg_EE_->Fill(2);
  	      }
  	    }
  	    else if(matchCategory == 0) {
  	      track_type_v2_.emplace_back(0);
  	      if (muon_Prompt) {
                if (Barrel_muon) meMuonISO_trk_type_v2_Sig_EB_->Fill(0);
                else meMuonISO_trk_type_v2_Sig_EE_->Fill(0);
              }
              else {
                if (Barrel_muon) meMuonISO_trk_type_v2_Bkg_EB_->Fill(0);
                else meMuonISO_trk_type_v2_Bkg_EE_->Fill(0);
              }
  	    }
  	    else if(matchCategory == 1) {
  	      track_type_v2_.emplace_back(3);
  	      if (muon_Prompt) {
                if (Barrel_muon) meMuonISO_trk_type_v2_Sig_EB_->Fill(3);
                else meMuonISO_trk_type_v2_Sig_EE_->Fill(3);
              }
              else {
                if (Barrel_muon) meMuonISO_trk_type_v2_Bkg_EB_->Fill(3);
                else meMuonISO_trk_type_v2_Bkg_EE_->Fill(3);
              }
  	    }
  	    else if(matchCategory == 2) {
  	      track_type_v2_.emplace_back(1);
  	      if (muon_Prompt) {
                if (Barrel_muon) meMuonISO_trk_type_v2_Sig_EB_->Fill(1);
                else meMuonISO_trk_type_v2_Sig_EE_->Fill(1);
              }
              else {
                if (Barrel_muon) meMuonISO_trk_type_v2_Bkg_EB_->Fill(1);
                else meMuonISO_trk_type_v2_Bkg_EE_->Fill(1);
              }
  	    }
//	  cout << "matchCategory: "<< matchCategory << endl;
	    else {
  	      track_bx_.emplace_back(-999);
  	      track_evtId_.emplace_back(-999);
	      track_type_v2_.emplace_back(-999);
	    }
	  }
//	  else { // not reco-sim vtx matching
//	    cout << "not reco-sim vtx matching" << endl;
//	  }
  	}
	match_vtx_reco2sim_.emplace_back(check_match_vtx_reco2sim);
	match_vtx_sim2reco_.emplace_back(check_match_vtx_sim2reco);
	selectedVtxMatching_.emplace_back(check_selectedVtxMatching);
	selectedLV_.emplace_back(check_selectedLV);
        // test end
	
	// for ntuple
	track_PVweight_.emplace_back(Vtx_chosen.trackWeight(trackref_general));
	track_pt_.emplace_back(trackGen.pt());
	track_time_.emplace_back(TrkMTDTime);
	track_time_err_.emplace_back(TrkMTDTimeErr);
	if (TrkMTDTimeErr > 0 && muon_sigTrkTimeErr > 0) {
	  dtsig_muon_track_.emplace_back(std::abs(TrkMTDTime - muon_sigTrkTime) / std::sqrt(TrkMTDTimeErr * TrkMTDTimeErr + muon_sigTrkTimeErr * muon_sigTrkTimeErr));
	}
	else dtsig_muon_track_.emplace_back(-999);
	if (TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0) {
	  dtsig_vtx_track_.emplace_back(std::abs(TrkMTDTime - Vtx_chosen.t()) / std::sqrt(TrkMTDTimeErr * TrkMTDTimeErr + Vtx_chosen.tError() * Vtx_chosen.tError()));
	}
	else dtsig_vtx_track_.emplace_back(-999);


        // MTD GEN case
        if (genMatched) {
          N_tracks_gen++;
          pT_sum_gen += trk_ptSim;
        }

        //const reco::TrackBaseRef trkrefBase(trackref_general);
        //auto TPmatched = r2s_->find(trkrefBase);
	
        // dt with the track
        if (dt_sig_track_) {
          double dt_sigTrk = 0;
          double dt_sigTrk_signif = 0;
          double dt_sim_sigTrk = 0;
          double dt_sim_sigTrk_signif = 0;

          // MTD SIM CASE
          if (tsim_trk != -1 && tsim_muon != -1 && trk_ptSim > 0) {
            dt_sim_sigTrk = std::abs(tsim_trk - tsim_muon);
            dt_sim_sigTrk_signif = dt_sim_sigTrk / std::sqrt(avg_sim_PUtrack_t_err * avg_sim_PUtrack_t_err +
                                                             avg_sim_sigTrk_t_err * avg_sim_sigTrk_t_err);
	    if (muon_Prompt) {
	      if (Barrel_muon) { // Barrel region
	        meMuonISO_dt_muon_trk_sim_Sig_EB_->Fill(dt_sim_sigTrk);
	        meMuonISO_dtSig_muon_trk_sim_Sig_EB_->Fill(dt_sim_sigTrk_signif);
	        if(genMatched) meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EB_->Fill(dt_sim_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EB_->Fill(dt_sim_sigTrk_signif);
		  }
	          else meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EB_->Fill(dt_sim_sigTrk_signif);
	        }
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EB_->Fill(1);
		}
	      }
	      else { // Endcap region
	        meMuonISO_dt_muon_trk_sim_Sig_EE_->Fill(dt_sim_sigTrk);
	        meMuonISO_dtSig_muon_trk_sim_Sig_EE_->Fill(dt_sim_sigTrk_signif);
	        if(genMatched) meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EE_->Fill(dt_sim_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EE_->Fill(dt_sim_sigTrk_signif);
		  }
	          else meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EE_->Fill(dt_sim_sigTrk_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EE_->Fill(1);
		}
	      }
	    }
	    else {
	      if (Barrel_muon) { // Barrel region
	        meMuonISO_dt_muon_trk_sim_Bkg_EB_->Fill(dt_sim_sigTrk);
	        meMuonISO_dtSig_muon_trk_sim_Bkg_EB_->Fill(dt_sim_sigTrk_signif);
	        if(genMatched) meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EB_->Fill(dt_sim_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EB_->Fill(dt_sim_sigTrk_signif);
		  }
	          else meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EB_->Fill(dt_sim_sigTrk_signif);
	        }
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EB_->Fill(1);
		}
	      }
	      else { // Endcap region
	        meMuonISO_dt_muon_trk_sim_Bkg_EE_->Fill(dt_sim_sigTrk);
	        meMuonISO_dtSig_muon_trk_sim_Bkg_EE_->Fill(dt_sim_sigTrk_signif);
	        if(genMatched) meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EE_->Fill(dt_sim_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EE_->Fill(dt_sim_sigTrk_signif);
		  }
	          else meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EE_->Fill(dt_sim_sigTrk_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EE_->Fill(1);
		}
	      }
	    }

            if (optionalPlots_) {
              // absolute timing cuts
              for (long unsigned int i = 0; i < N_tracks_sim_MTD.size(); i++) {
                if (dt_sim_sigTrk < max_dt_track_cut[i]) {
                  N_tracks_sim_MTD[i] = N_tracks_sim_MTD[i] + 1;
                  pT_sum_sim_MTD[i] = pT_sum_sim_MTD[i] + trk_ptSim;
                }
              }
            }
            // significance cuts
            for (long unsigned int i = 0; i < N_tracks_sim_MTD_significance.size(); i++) {
              if (dt_sim_sigTrk_signif < max_dt_significance_cut[i]) {
                N_tracks_sim_MTD_significance[i]++;
                pT_sum_sim_MTD_significance[i] += trk_ptSim;
		// test
		if(genMatched) pT_sum_gen_MTD_significance[i] += trk_ptSim;
              }
            }

          } else if (trk_ptSim > 0) {
            // if there is no error for MTD information, we count the MTD isolation case same as noMTD
            if (optionalPlots_) {
              for (long unsigned int i = 0; i < N_tracks_sim_MTD.size(); i++) {
                N_tracks_sim_MTD[i] = N_tracks_sim_MTD[i] + 1;
                pT_sum_sim_MTD[i] = pT_sum_sim_MTD[i] + trk_ptSim;
		
		// test for checking the type of tracks
		if (muon_Prompt) {
		  if (Barrel_muon) meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EB_->Fill(1);
		  else meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EE_->Fill(1);
		}
		else {
		  if (Barrel_muon) meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EB_->Fill(1);
		  else meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EE_->Fill(1);
		}
              }
            }
            for (long unsigned int i = 0; i < N_tracks_sim_MTD_significance.size(); i++) {
              N_tracks_sim_MTD_significance[i]++;
              pT_sum_sim_MTD_significance[i] += trk_ptSim;
	      // test
	      if(genMatched) pT_sum_gen_MTD_significance[i] += trk_ptSim;
            }
          }

          // MTD reco case
          if (TrkMTDTimeErr > 0 && muon_sigTrkTimeErr > 0) {    //FIXME For tracks, there exist cases where track has time or error of time even without an MVA score. It seems MVA score is needed. // FIXME It is already considered by defining tErr of track above
            dt_sigTrk = std::abs(TrkMTDTime - muon_sigTrkTime);
            dt_sigTrk_signif =
                dt_sigTrk / std::sqrt(TrkMTDTimeErr * TrkMTDTimeErr + muon_sigTrkTimeErr * muon_sigTrkTimeErr);

	    if (muon_Prompt) {
	      if (Barrel_muon) {
	        meMuonISO_dt_muon_trk_reco_Sig_EB_->Fill(dt_sigTrk);
	        meMuonISO_dtSig_muon_trk_reco_Sig_EB_->Fill(dt_sigTrk_signif);
	        if (genMatched) meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EB_->Fill(dt_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EB_->Fill(dt_sigTrk_signif);
  		  }
  	          else meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EB_->Fill(dt_sigTrk_signif);
	        }
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EB_->Fill(1);
		}
	      }
	      else {
	        meMuonISO_dt_muon_trk_reco_Sig_EE_->Fill(dt_sigTrk);
	        meMuonISO_dtSig_muon_trk_reco_Sig_EE_->Fill(dt_sigTrk_signif);
	        if (genMatched) meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EE_->Fill(dt_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EE_->Fill(dt_sigTrk_signif);
  		  }
  	          else meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EE_->Fill(dt_sigTrk_signif);
		}
		else {
		  meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EE_->Fill(1);
		}
	      }
	    }
	    else { // Non-prompt
	      if (Barrel_muon) {
	        meMuonISO_dt_muon_trk_reco_Bkg_EB_->Fill(dt_sigTrk);
	        meMuonISO_dtSig_muon_trk_reco_Bkg_EB_->Fill(dt_sigTrk_signif);
	        if (genMatched) meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EB_->Fill(dt_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EB_->Fill(dt_sigTrk_signif);
  		  }
  	          else meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EB_->Fill(dt_sigTrk_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EB_->Fill(1);
		}
	      }
	      else {
	        meMuonISO_dt_muon_trk_reco_Bkg_EE_->Fill(dt_sigTrk);
	        meMuonISO_dtSig_muon_trk_reco_Bkg_EE_->Fill(dt_sigTrk_signif);
	        if (genMatched) meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EE_->Fill(dt_sigTrk_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EE_->Fill(dt_sigTrk_signif);
  		  }
  	          else meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EE_->Fill(dt_sigTrk_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EE_->Fill(1);
		}
	      }
	    }
            meMuon_no_dt_check_->Fill(1);

            if (optionalPlots_) {
              // absolute timing cuts
              for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
                if (dt_sigTrk < max_dt_track_cut[i]) {
                  N_tracks_MTD[i] = N_tracks_MTD[i] + 1;
                  pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();
                }
              }
            }
            // significance cuts
            for (long unsigned int i = 0; i < N_tracks_MTD_significance.size(); i++) {
              if (dt_sigTrk_signif < max_dt_significance_cut[i]) {
                N_tracks_MTD_significance[i]++;
                pT_sum_MTD_significance[i] += trackGen.pt();
              }
            }

          } else {
            // if there is no error for MTD information, we count the MTD isolation case same as noMTD
            if (optionalPlots_) {
              for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
                N_tracks_MTD[i] = N_tracks_MTD[i] + 1;
                pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();
              }
	      
	      // test for checking the type of tracks
	      if (muon_Prompt) {
		if (Barrel_muon) meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EB_->Fill(1);
		else meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EE_->Fill(1);
	      }
	      else {
		if (Barrel_muon) meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EB_->Fill(1);
		else meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EE_->Fill(1);
	      }
            }
            for (long unsigned int i = 0; i < N_tracks_MTD_significance.size(); i++) {
              N_tracks_MTD_significance[i]++;
              pT_sum_MTD_significance[i] += trackGen.pt();
            }
            meMuon_no_dt_check_->Fill(0);
          }

          if (optionalPlots_) {
            for (long unsigned int i = 0; i < (pT_bins_dt_distrb.size() - 1); i++) {
              //stuff general pT
              if (muon.track()->pt() > pT_bins_dt_distrb[i] && muon.track()->pt() < pT_bins_dt_distrb[i + 1]) {
                general_pT_list[i]->Fill(dt_sigTrk);
                general_pT_Signif_list[i]->Fill(dt_sigTrk_signif);
              }
            }

            for (long unsigned int i = 0; i < (eta_bins_dt_distrib.size() - 1); i++) {
              //stuff general eta
              if (std::abs(muon.track()->eta()) > eta_bins_dt_distrib[i] && std::abs(muon.track()->eta()) < eta_bins_dt_distrib[i + 1]) {
                general_eta_list[i]->Fill(dt_sigTrk);
                general_eta_Signif_list[i]->Fill(dt_sigTrk_signif);
              }
            }
          }  // End of optional dt distributions plots

          // dt with the vertex
        } else {
          double dt_vtx = 0;  // dt regular track vs vtx
          double dt_vtx_signif = 0;

          double dt_sim_vtx = 0;  // dt regular track vs vtx
          double dt_sim_vtx_signif = 0;

          // MTD SIM case
          if (tsim_trk != -1 && Vtx_chosen.tError() > 0 && trk_ptSim > 0) {
            dt_sim_vtx = std::abs(tsim_trk - Vtx_chosen.t());
            dt_sim_vtx_signif = dt_sim_vtx / std::sqrt(avg_sim_PUtrack_t_err * avg_sim_PUtrack_t_err +
                                                       Vtx_chosen.tError() * Vtx_chosen.tError());
	    if (muon_Prompt) {
	      if (Barrel_muon) { // Barrel region
	        meMuonISO_dt_PV_trk_sim_Sig_EB_->Fill(dt_sim_vtx);
	        meMuonISO_dtSig_PV_trk_sim_Sig_EB_->Fill(dt_sim_vtx_signif);
	        if(genMatched) meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EB_->Fill(dt_sim_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EB_->Fill(dt_sim_vtx_signif);
		  }
	          else meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EB_->Fill(dt_sim_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EB_->Fill(1);
		}
	      }
	      else { // Endcap region
	        meMuonISO_dt_PV_trk_sim_Sig_EE_->Fill(dt_sim_vtx);
	        meMuonISO_dtSig_PV_trk_sim_Sig_EE_->Fill(dt_sim_vtx_signif);
	        if(genMatched) meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EE_->Fill(dt_sim_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EE_->Fill(dt_sim_vtx_signif);
		  }
	          else meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EE_->Fill(dt_sim_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EE_->Fill(1);
		}
	      }
	    }
	    else {
	      if (Barrel_muon) { // Barrel region
	        meMuonISO_dt_PV_trk_sim_Bkg_EB_->Fill(dt_sim_vtx);
	        meMuonISO_dtSig_PV_trk_sim_Bkg_EB_->Fill(dt_sim_vtx_signif);
	        if(genMatched) meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EB_->Fill(dt_sim_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EB_->Fill(dt_sim_vtx_signif);
		  }
	          else meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EB_->Fill(dt_sim_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EB_->Fill(1);
		}
	      }
	      else { // Endcap region
	        meMuonISO_dt_PV_trk_sim_Bkg_EE_->Fill(dt_sim_vtx);
	        meMuonISO_dtSig_PV_trk_sim_Bkg_EE_->Fill(dt_sim_vtx_signif);
	        if(genMatched) meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EE_->Fill(dt_sim_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
	          if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
		    meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EE_->Fill(dt_sim_vtx_signif);
		  }
	          else meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EE_->Fill(dt_sim_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EE_->Fill(1);
		}
	      }
	    }

            if (optionalPlots_) {
              // absolute timing cuts
              for (long unsigned int i = 0; i < N_tracks_sim_MTD.size(); i++) {
                if (dt_sim_vtx < max_dt_vtx_cut[i]) {
                  N_tracks_sim_MTD[i] = N_tracks_sim_MTD[i] + 1;
                  pT_sum_sim_MTD[i] = pT_sum_sim_MTD[i] + trk_ptSim;
                }
              }
            }
            // significance timing cuts
            for (long unsigned int i = 0; i < N_tracks_sim_MTD_significance.size(); i++) {
              if (dt_sim_vtx_signif < max_dt_significance_cut[i]) {
                N_tracks_sim_MTD_significance[i]++;
                pT_sum_sim_MTD_significance[i] += trk_ptSim;
		// test
		if(genMatched) pT_sum_gen_MTD_significance[i] += trk_ptSim;
              }
            }
          } else if (trk_ptSim > 0) {
            if (optionalPlots_) {
              for (long unsigned int i = 0; i < N_tracks_sim_MTD.size(); i++) {
                N_tracks_sim_MTD[i] = N_tracks_sim_MTD[i] + 1;      // N_tracks_noMTD
                pT_sum_sim_MTD[i] = pT_sum_sim_MTD[i] + trk_ptSim;  // pT_sum_noMTD

		// test for checking the type of tracks
                if (muon_Prompt) {
		  if (Barrel_muon) meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EB_->Fill(1);
		  else meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EE_->Fill(1);
		}
                else {
		  if (Barrel_muon) meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EB_->Fill(1);
		  else meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EE_->Fill(1);
		}
              }
            }
            for (long unsigned int i = 0; i < N_tracks_sim_MTD_significance.size(); i++) {
              N_tracks_sim_MTD_significance[i]++;
              pT_sum_sim_MTD_significance[i] += trk_ptSim;
	      // test
	      if(genMatched) pT_sum_gen_MTD_significance[i] += trk_ptSim;
            }
          }

          // MTD RECO case
          if (TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0) {    //FIXME For tracks, there exist cases where track has time or error of time even without an MVA score. It seems MVA score is needed.
            dt_vtx = std::abs(TrkMTDTime - Vtx_chosen.t());
            dt_vtx_signif =
                dt_vtx / std::sqrt(TrkMTDTimeErr * TrkMTDTimeErr + Vtx_chosen.tError() * Vtx_chosen.tError());

	    if (muon_Prompt) {
	      if (Barrel_muon) {
	        meMuonISO_dt_PV_trk_reco_Sig_EB_->Fill(dt_vtx);
	        meMuonISO_dtSig_PV_trk_reco_Sig_EB_->Fill(dt_vtx_signif);
	        if (genMatched) meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EB_->Fill(dt_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EB_->Fill(dt_vtx_signif);
  		  }
  	          else meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EB_->Fill(dt_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EB_->Fill(1);
		}
	      }
	      else {
	        meMuonISO_dt_PV_trk_reco_Sig_EE_->Fill(dt_vtx);
	        meMuonISO_dtSig_PV_trk_reco_Sig_EE_->Fill(dt_vtx_signif);
	        if (genMatched) meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EE_->Fill(dt_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EE_->Fill(dt_vtx_signif);
  		  }
  	          else meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EE_->Fill(dt_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EE_->Fill(1);
		}
	      }
	    }
	    else { // Non-prompt
	      if (Barrel_muon) {
	        meMuonISO_dt_PV_trk_reco_Bkg_EB_->Fill(dt_vtx);
	        meMuonISO_dtSig_PV_trk_reco_Bkg_EB_->Fill(dt_vtx_signif);
	        if (genMatched) meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EB_->Fill(dt_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EB_->Fill(dt_vtx_signif);
  		  }
  	          else meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EB_->Fill(dt_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EB_->Fill(1);
		}
	      }
	      else {
	        meMuonISO_dt_PV_trk_reco_Bkg_EE_->Fill(dt_vtx);
	        meMuonISO_dtSig_PV_trk_reco_Bkg_EE_->Fill(dt_vtx_signif);
	        if (genMatched) meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EE_->Fill(dt_vtx_signif);

	        // test for checking the type of tracks
		if (TPmatched != r2s_->end()) {
		  const auto& tp = (TPmatched->val)[0];
  		  if ((tp.first)->eventId().bunchCrossing()==0 && (tp.first)->eventId().event()==0) {
  		    meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EE_->Fill(dt_vtx_signif);
  		  }
  	          else meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EE_->Fill(dt_vtx_signif);
		}
		else { // tracks are not TPmatched but have t, tErr, MVA -> maybe fake track?
		  meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EE_->Fill(1);
		}
	      }
	    }
            meMuon_no_dt_check_->Fill(1);
            meMuon_avg_error_vtx_check_->Fill(Vtx_chosen.tError());

            if (optionalPlots_) {
              // absolute timing cuts
              for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
                if (dt_vtx < max_dt_vtx_cut[i]) {
                  N_tracks_MTD[i] = N_tracks_MTD[i] + 1;
                  pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();
                }
              }
            }
            // significance timing cuts
            for (long unsigned int i = 0; i < N_tracks_MTD_significance.size(); i++) {
              if (dt_vtx_signif < max_dt_significance_cut[i]) {
                N_tracks_MTD_significance[i]++;
                pT_sum_MTD_significance[i] += trackGen.pt();
              }
            }
          } else {
            if (optionalPlots_) {
              for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
                N_tracks_MTD[i] = N_tracks_MTD[i] + 1;          // N_tracks_noMTD
                pT_sum_MTD[i] = pT_sum_MTD[i] + trackGen.pt();  // pT_sum_noMTD
	      }

	      // test for checking the type of tracks
              if (muon_Prompt) {
		if (Barrel_muon) meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EB_->Fill(1);
		else meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EE_->Fill(1);
	      }
              else {
		if (Barrel_muon) meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EB_->Fill(1);
		else meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EE_->Fill(1);
              }
            }
            for (long unsigned int i = 0; i < N_tracks_MTD_significance.size(); i++) {
              N_tracks_MTD_significance[i]++;
              pT_sum_MTD_significance[i] += trackGen.pt();
            }
            meMuon_no_dt_check_->Fill(0);
          }

          // Optional dt distribution plots
          if (optionalPlots_) {
            for (long unsigned int i = 0; i < (pT_bins_dt_distrb.size() - 1); i++) {
              //stuff general pT
              if (muon.track()->pt() > pT_bins_dt_distrb[i] && muon.track()->pt() < pT_bins_dt_distrb[i + 1]) {
                general_pT_list[i]->Fill(dt_vtx);
                general_pT_Signif_list[i]->Fill(dt_vtx_signif);
              }
            }

            for (long unsigned int i = 0; i < (eta_bins_dt_distrib.size() - 1); i++) {
              //stuff general eta
              if (std::abs(muon.track()->eta()) > eta_bins_dt_distrib[i] && std::abs(muon.track()->eta()) < eta_bins_dt_distrib[i + 1]) {
                general_eta_list[i]->Fill(dt_vtx);
                general_eta_Signif_list[i]->Fill(dt_vtx_signif);
              }
            }
          }  // End of optional dt distributions plots
        }
      }

      rel_pT_sum_noMTD = pT_sum_noMTD / muon.track()->pt();  // rel_ch_iso calculation

      if (optionalPlots_) {
        for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
          rel_pT_sum_MTD[i] = pT_sum_MTD[i] / muon.track()->pt();
          if(muon_sim_pt!=-1) rel_pT_sum_sim_MTD[i] = pT_sum_sim_MTD[i] / muon_sim_pt;
        }
        // now compute the isolation
        rel_pT_sum_noMTD = pT_sum_noMTD / muon.track()->pt();

        if(muon_sim_pt!=-1) rel_pT_sum_gen = pT_sum_gen / muon_sim_pt;
      }

      for (long unsigned int i = 0; i < N_tracks_MTD_significance.size(); i++) {
        rel_pT_sum_MTD_significance[i] = pT_sum_MTD_significance[i] / muon.track()->pt();
        if(muon_sim_pt!=-1) rel_pT_sum_sim_MTD_significance[i] = pT_sum_sim_MTD_significance[i] / muon_sim_pt;
        if(muon_sim_pt!=-1) rel_pT_sum_gen_MTD_significance[i] = pT_sum_gen_MTD_significance[i] / muon_sim_pt;
      }

      if (muon_Prompt) {  // promt part
        if (Barrel_muon) {
          meMuonISO_Ntracks_Sig_EB_->Fill(N_tracks_noMTD);
          meMuonISO_chIso_Sig_EB_->Fill(pT_sum_noMTD);
          meMuonISO_rel_chIso_Sig_EB_->Fill(rel_pT_sum_noMTD);
          if (optionalPlots_) {
            for (long unsigned int j = 0; j < Ntracks_EB_list_Sig.size(); j++) {
              Ntracks_EB_list_Sig[j]->Fill(N_tracks_MTD[j]);
              ch_iso_EB_list_Sig[j]->Fill(pT_sum_MTD[j]);
              rel_ch_iso_EB_list_Sig[j]->Fill(rel_pT_sum_MTD[j]);

              Ntracks_sim_EB_list_Sig[j]->Fill(N_tracks_sim_MTD[j]);
              ch_iso_sim_EB_list_Sig[j]->Fill(pT_sum_sim_MTD[j]);
              rel_ch_iso_sim_EB_list_Sig[j]->Fill(rel_pT_sum_sim_MTD[j]);
            }
            meMuonISO_Ntracks_gen_Sig_EB_->Fill(N_tracks_gen);
            meMuonISO_chIso_gen_Sig_EB_->Fill(pT_sum_gen);
            meMuonISO_rel_chIso_gen_Sig_EB_->Fill(rel_pT_sum_gen);
          }

          for (long unsigned int j = 0; j < Ntracks_EB_list_Significance_Sig.size(); j++) {
            Ntracks_EB_list_Significance_Sig[j]->Fill(N_tracks_MTD_significance[j]);
            ch_iso_EB_list_Significance_Sig[j]->Fill(pT_sum_MTD_significance[j]);
            rel_ch_iso_EB_list_Significance_Sig[j]->Fill(rel_pT_sum_MTD_significance[j]);

            if (optionalPlots_) {
              Ntracks_sim_EB_list_Significance_Sig[j]->Fill(N_tracks_sim_MTD_significance[j]);
              ch_iso_sim_EB_list_Significance_Sig[j]->Fill(pT_sum_sim_MTD_significance[j]);
              rel_ch_iso_sim_EB_list_Significance_Sig[j]->Fill(rel_pT_sum_sim_MTD_significance[j]);
              ch_iso_gen_EB_list_Significance_Sig[j]->Fill(pT_sum_gen_MTD_significance[j]);
	      rel_ch_iso_gen_EB_list_Significance_Sig[j]->Fill(rel_pT_sum_gen_MTD_significance[j]);
            }
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meMuon_pt_noMTD_Sig_EB_->Fill(muon.track()->pt());
            meMuon_eta_noMTD_Sig_EB_->Fill(std::abs(muon.track()->eta()));
            meMuon_phi_noMTD_Sig_EB_->Fill(muon.track()->phi());
          }
          if (optionalPlots_) {
            for (long unsigned int k = 0; k < Ntracks_EB_list_Sig.size(); k++) {
              if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
                Muon_pT_MTD_EB_list_Sig[k]->Fill(muon.track()->pt());
                Muon_eta_MTD_EB_list_Sig[k]->Fill(std::abs(muon.track()->eta()));
                Muon_phi_MTD_EB_list_Sig[k]->Fill(muon.track()->phi());

                if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EB_list_Sig[k]->Fill(muon_sim_pt);
              }
            }
            if (rel_pT_sum_gen < rel_iso_cut_) {
              if(muon_sim_pt!=-1) meMuon_pt_gen_Sig_EB_->Fill(muon_sim_pt);
              meMuon_eta_gen_Sig_EB_->Fill(std::abs(muon_sim_eta));
              meMuon_phi_gen_Sig_EB_->Fill(muon_sim_phi);
            }
          }

          for (long unsigned int k = 0; k < Ntracks_EB_list_Significance_Sig.size(); k++) {
            if (rel_pT_sum_MTD_significance[k] < rel_iso_cut_) {
              Muon_pT_MTD_EB_list_Significance_Sig[k]->Fill(muon.track()->pt());
              Muon_eta_MTD_EB_list_Significance_Sig[k]->Fill(std::abs(muon.track()->eta()));
              Muon_phi_MTD_EB_list_Significance_Sig[k]->Fill(muon.track()->phi());
            }
            if (optionalPlots_ and rel_pT_sum_sim_MTD_significance[k] < rel_iso_cut_) {
              if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EB_list_Significance_Sig[k]->Fill(muon_sim_pt);
	    }
	    if (optionalPlots_ and rel_pT_sum_gen_MTD_significance[k] < rel_iso_cut_) {
	      if(muon_sim_pt!=-1) Muon_pT_gen_MTD_EB_list_Significance_Sig[k]->Fill(muon_sim_pt);
	    }
          }

        } else {  // for endcap
          meMuonISO_Ntracks_Sig_EE_->Fill(N_tracks_noMTD);
          meMuonISO_chIso_Sig_EE_->Fill(pT_sum_noMTD);
          meMuonISO_rel_chIso_Sig_EE_->Fill(rel_pT_sum_noMTD);
          if (optionalPlots_) {
            for (long unsigned int j = 0; j < Ntracks_EE_list_Sig.size(); j++) {
              Ntracks_EE_list_Sig[j]->Fill(N_tracks_MTD[j]);
              ch_iso_EE_list_Sig[j]->Fill(pT_sum_MTD[j]);
              rel_ch_iso_EE_list_Sig[j]->Fill(rel_pT_sum_MTD[j]);

              Ntracks_sim_EE_list_Sig[j]->Fill(N_tracks_sim_MTD[j]);
              ch_iso_sim_EE_list_Sig[j]->Fill(pT_sum_sim_MTD[j]);
              rel_ch_iso_sim_EE_list_Sig[j]->Fill(rel_pT_sum_sim_MTD[j]);
            }
            meMuonISO_Ntracks_gen_Sig_EE_->Fill(N_tracks_gen);
            meMuonISO_chIso_gen_Sig_EE_->Fill(pT_sum_gen);
            meMuonISO_rel_chIso_gen_Sig_EE_->Fill(rel_pT_sum_gen);
          }

          for (long unsigned int j = 0; j < Ntracks_EE_list_Significance_Sig.size(); j++) {
            Ntracks_EE_list_Significance_Sig[j]->Fill(N_tracks_MTD_significance[j]);
            ch_iso_EE_list_Significance_Sig[j]->Fill(pT_sum_MTD_significance[j]);
            rel_ch_iso_EE_list_Significance_Sig[j]->Fill(rel_pT_sum_MTD_significance[j]);

            if (optionalPlots_) {
              Ntracks_sim_EE_list_Significance_Sig[j]->Fill(N_tracks_sim_MTD_significance[j]);
              ch_iso_sim_EE_list_Significance_Sig[j]->Fill(pT_sum_sim_MTD_significance[j]);
              rel_ch_iso_sim_EE_list_Significance_Sig[j]->Fill(rel_pT_sum_sim_MTD_significance[j]);
              ch_iso_gen_EE_list_Significance_Sig[j]->Fill(pT_sum_gen_MTD_significance[j]);
              rel_ch_iso_gen_EE_list_Significance_Sig[j]->Fill(rel_pT_sum_gen_MTD_significance[j]);
            }
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meMuon_pt_noMTD_Sig_EE_->Fill(muon.track()->pt());
            meMuon_eta_noMTD_Sig_EE_->Fill(std::abs(muon.track()->eta()));
            meMuon_phi_noMTD_Sig_EE_->Fill(muon.track()->phi());
          }
          if (optionalPlots_) {
            for (long unsigned int k = 0; k < Ntracks_EE_list_Sig.size(); k++) {
              if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
                Muon_pT_MTD_EE_list_Sig[k]->Fill(muon.track()->pt());
                Muon_eta_MTD_EE_list_Sig[k]->Fill(std::abs(muon.track()->eta()));
                Muon_phi_MTD_EE_list_Sig[k]->Fill(muon.track()->phi());

                if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EE_list_Sig[k]->Fill(muon_sim_pt);
              }
            }
            if (rel_pT_sum_gen < rel_iso_cut_) {
              if(muon_sim_pt!=-1) meMuon_pt_gen_Sig_EE_->Fill(muon_sim_pt);
              meMuon_eta_gen_Sig_EE_->Fill(std::abs(muon_sim_eta));
              meMuon_phi_gen_Sig_EE_->Fill(muon_sim_phi);
            }
          }
          for (long unsigned int k = 0; k < Ntracks_EE_list_Significance_Sig.size(); k++) {
            if (rel_pT_sum_MTD_significance[k] < rel_iso_cut_) {
              Muon_pT_MTD_EE_list_Significance_Sig[k]->Fill(muon.track()->pt());
              Muon_eta_MTD_EE_list_Significance_Sig[k]->Fill(std::abs(muon.track()->eta()));
              Muon_phi_MTD_EE_list_Significance_Sig[k]->Fill(muon.track()->phi());
	    }
            if (optionalPlots_ and rel_pT_sum_sim_MTD_significance[k] < rel_iso_cut_) {
              if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EE_list_Significance_Sig[k]->Fill(muon_sim_pt);
	    }
	    if (optionalPlots_ and rel_pT_sum_gen_MTD_significance[k] < rel_iso_cut_) {
	      if(muon_sim_pt!=-1) Muon_pT_gen_MTD_EE_list_Significance_Sig[k]->Fill(muon_sim_pt);
	    }
          }
        }
      } else {  // non-promt part
        if (Barrel_muon) {
          meMuonISO_Ntracks_Bkg_EB_->Fill(N_tracks_noMTD);
          meMuonISO_chIso_Bkg_EB_->Fill(pT_sum_noMTD);
          meMuonISO_rel_chIso_Bkg_EB_->Fill(rel_pT_sum_noMTD);
          if (optionalPlots_) {
            for (long unsigned int j = 0; j < Ntracks_EB_list_Bkg.size(); j++) {
              Ntracks_EB_list_Bkg[j]->Fill(N_tracks_MTD[j]);
              ch_iso_EB_list_Bkg[j]->Fill(pT_sum_MTD[j]);
              rel_ch_iso_EB_list_Bkg[j]->Fill(rel_pT_sum_MTD[j]);

              Ntracks_sim_EB_list_Bkg[j]->Fill(N_tracks_sim_MTD[j]);
              ch_iso_sim_EB_list_Bkg[j]->Fill(pT_sum_sim_MTD[j]);
              rel_ch_iso_sim_EB_list_Bkg[j]->Fill(rel_pT_sum_sim_MTD[j]);
            }
            meMuonISO_Ntracks_gen_Bkg_EB_->Fill(N_tracks_gen);
            meMuonISO_chIso_gen_Bkg_EB_->Fill(pT_sum_gen);
            meMuonISO_rel_chIso_gen_Bkg_EB_->Fill(rel_pT_sum_gen);
          }

          for (long unsigned int j = 0; j < Ntracks_EB_list_Significance_Bkg.size(); j++) {
            Ntracks_EB_list_Significance_Bkg[j]->Fill(N_tracks_MTD_significance[j]);
            ch_iso_EB_list_Significance_Bkg[j]->Fill(pT_sum_MTD_significance[j]);
            rel_ch_iso_EB_list_Significance_Bkg[j]->Fill(rel_pT_sum_MTD_significance[j]);

            if (optionalPlots_) {
              Ntracks_sim_EB_list_Significance_Bkg[j]->Fill(N_tracks_sim_MTD_significance[j]);
              ch_iso_sim_EB_list_Significance_Bkg[j]->Fill(pT_sum_sim_MTD_significance[j]);
              rel_ch_iso_sim_EB_list_Significance_Bkg[j]->Fill(rel_pT_sum_sim_MTD_significance[j]);
              ch_iso_gen_EB_list_Significance_Bkg[j]->Fill(pT_sum_gen_MTD_significance[j]);
              rel_ch_iso_gen_EB_list_Significance_Bkg[j]->Fill(rel_pT_sum_gen_MTD_significance[j]);
            }
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meMuon_pt_noMTD_Bkg_EB_->Fill(muon.track()->pt());
            meMuon_eta_noMTD_Bkg_EB_->Fill(std::abs(muon.track()->eta()));
            meMuon_phi_noMTD_Bkg_EB_->Fill(muon.track()->phi());
          }
          if (optionalPlots_) {
            for (long unsigned int k = 0; k < Ntracks_EB_list_Bkg.size(); k++) {
              if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
                Muon_pT_MTD_EB_list_Bkg[k]->Fill(muon.track()->pt());
                Muon_eta_MTD_EB_list_Bkg[k]->Fill(std::abs(muon.track()->eta()));
                Muon_phi_MTD_EB_list_Bkg[k]->Fill(muon.track()->phi());

                if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EB_list_Bkg[k]->Fill(muon_sim_pt);
              }
            }
            if (rel_pT_sum_gen < rel_iso_cut_) {
              if(muon_sim_pt!=-1) meMuon_pt_gen_Bkg_EB_->Fill(muon_sim_pt);
              meMuon_eta_gen_Bkg_EB_->Fill(std::abs(muon_sim_eta));
              meMuon_phi_gen_Bkg_EB_->Fill(muon_sim_phi);
            }
          }
          for (long unsigned int k = 0; k < Ntracks_EB_list_Significance_Bkg.size(); k++) {
            if (rel_pT_sum_MTD_significance[k] < rel_iso_cut_) {
              Muon_pT_MTD_EB_list_Significance_Bkg[k]->Fill(muon.track()->pt());
              Muon_eta_MTD_EB_list_Significance_Bkg[k]->Fill(std::abs(muon.track()->eta()));
              Muon_phi_MTD_EB_list_Significance_Bkg[k]->Fill(muon.track()->phi());
	    }

            if (optionalPlots_ and rel_pT_sum_sim_MTD_significance[k] < rel_iso_cut_) {
              if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EB_list_Significance_Bkg[k]->Fill(muon_sim_pt);
	    }
	    if (optionalPlots_ and rel_pT_sum_gen_MTD_significance[k] < rel_iso_cut_) {
	      if(muon_sim_pt!=-1) Muon_pT_gen_MTD_EB_list_Significance_Bkg[k]->Fill(muon_sim_pt);
	    }
          }

        } else {  // for endcap
          meMuonISO_Ntracks_Bkg_EE_->Fill(N_tracks_noMTD);
          meMuonISO_chIso_Bkg_EE_->Fill(pT_sum_noMTD);
          meMuonISO_rel_chIso_Bkg_EE_->Fill(rel_pT_sum_noMTD);
          if (optionalPlots_) {
            for (long unsigned int j = 0; j < Ntracks_EE_list_Bkg.size(); j++) {
              Ntracks_EE_list_Bkg[j]->Fill(N_tracks_MTD[j]);
              ch_iso_EE_list_Bkg[j]->Fill(pT_sum_MTD[j]);
              rel_ch_iso_EE_list_Bkg[j]->Fill(rel_pT_sum_MTD[j]);

              Ntracks_sim_EE_list_Bkg[j]->Fill(N_tracks_sim_MTD[j]);
              ch_iso_sim_EE_list_Bkg[j]->Fill(pT_sum_sim_MTD[j]);
              rel_ch_iso_sim_EE_list_Bkg[j]->Fill(rel_pT_sum_sim_MTD[j]);
            }
            meMuonISO_Ntracks_gen_Bkg_EE_->Fill(N_tracks_gen);
            meMuonISO_chIso_gen_Bkg_EE_->Fill(pT_sum_gen);
            meMuonISO_rel_chIso_gen_Bkg_EE_->Fill(rel_pT_sum_gen);
          }

          for (long unsigned int j = 0; j < Ntracks_EE_list_Significance_Bkg.size(); j++) {
            Ntracks_EE_list_Significance_Bkg[j]->Fill(N_tracks_MTD_significance[j]);
            ch_iso_EE_list_Significance_Bkg[j]->Fill(pT_sum_MTD_significance[j]);
            rel_ch_iso_EE_list_Significance_Bkg[j]->Fill(rel_pT_sum_MTD_significance[j]);

            if (optionalPlots_) {
              Ntracks_sim_EE_list_Significance_Bkg[j]->Fill(N_tracks_sim_MTD_significance[j]);
              ch_iso_sim_EE_list_Significance_Bkg[j]->Fill(pT_sum_sim_MTD_significance[j]);
              rel_ch_iso_sim_EE_list_Significance_Bkg[j]->Fill(rel_pT_sum_sim_MTD_significance[j]);
              ch_iso_gen_EE_list_Significance_Bkg[j]->Fill(pT_sum_gen_MTD_significance[j]);
              rel_ch_iso_gen_EE_list_Significance_Bkg[j]->Fill(rel_pT_sum_gen_MTD_significance[j]);
            }
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meMuon_pt_noMTD_Bkg_EE_->Fill(muon.track()->pt());
            meMuon_eta_noMTD_Bkg_EE_->Fill(std::abs(muon.track()->eta()));
            meMuon_phi_noMTD_Bkg_EE_->Fill(muon.track()->phi());
          }
          if (optionalPlots_) {
            for (long unsigned int k = 0; k < Ntracks_EE_list_Bkg.size(); k++) {
              if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
                Muon_pT_MTD_EE_list_Bkg[k]->Fill(muon.track()->pt());
                Muon_eta_MTD_EE_list_Bkg[k]->Fill(std::abs(muon.track()->eta()));
                Muon_phi_MTD_EE_list_Bkg[k]->Fill(muon.track()->phi());

                if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EE_list_Bkg[k]->Fill(muon_sim_pt);
              }
            }
            if (rel_pT_sum_gen < rel_iso_cut_) {
              if(muon_sim_pt!=-1) meMuon_pt_gen_Bkg_EE_->Fill(muon_sim_pt);
              meMuon_eta_gen_Bkg_EE_->Fill(std::abs(muon_sim_eta));
              meMuon_phi_gen_Bkg_EE_->Fill(muon_sim_phi);
            }
          }

          for (long unsigned int k = 0; k < Ntracks_EE_list_Significance_Bkg.size(); k++) {
            if (rel_pT_sum_MTD_significance[k] < rel_iso_cut_) {
              Muon_pT_MTD_EE_list_Significance_Bkg[k]->Fill(muon.track()->pt());
              Muon_eta_MTD_EE_list_Significance_Bkg[k]->Fill(std::abs(muon.track()->eta()));
              Muon_phi_MTD_EE_list_Significance_Bkg[k]->Fill(muon.track()->phi());
	    }
            if (optionalPlots_ and rel_pT_sum_sim_MTD_significance[k] < rel_iso_cut_) {
              if(muon_sim_pt!=-1) Muon_pT_sim_MTD_EE_list_Significance_Bkg[k]->Fill(muon_sim_pt);
	    }
	    if (optionalPlots_ and rel_pT_sum_gen_MTD_significance[k] < rel_iso_cut_) {
	      if(muon_sim_pt!=-1) Muon_pT_gen_MTD_EE_list_Significance_Bkg[k]->Fill(muon_sim_pt);
	    }
          }
        }
      }
    }  // muon matched to a track
  tree_->Fill();
  }    // muon collection inside single event
//  tree_->Fill();
  meMuonISO_Nmuons_Sig_->Fill(nmuons_Sig);
  meMuonISO_Nmuons_Bkg_->Fill(nmuons_Bkg);
  if(nmuons_Sig_EB!=0) meMuonISO_Nmuons_Sig_EB_->Fill(nmuons_Sig_EB);
  if(nmuons_Sig_EE!=0) meMuonISO_Nmuons_Sig_EE_->Fill(nmuons_Sig_EE);
  if(nmuons_Bkg_EB!=0) meMuonISO_Nmuons_Bkg_EB_->Fill(nmuons_Bkg_EB);
  if(nmuons_Bkg_EE!=0) meMuonISO_Nmuons_Bkg_EE_->Fill(nmuons_Bkg_EE);
}

// ------------ method for histogram booking ------------
void MtdMuonIsoValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // for regular Validation use a reduced binning, for detailed analysis and ROC curves use the larger one
  //int nbin_1 = 40;
  //int nbin_2 = 40;
  int nbin_1 = 1000;
  int nbin_2 = 2000;
  if (optionalPlots_) {
    nbin_1 = 1000;
    nbin_2 = 2000;
  }

  // histogram booking

  meMuon_avg_error_SigTrk_check_ =
      ibook.book1D("SigTrk_avg_timing_err",
                   "Average signal muon track MTD timing uncertainty;Time Error (ns);Counts",
                   200,
                   0,
                   0.1);
  meMuon_avg_error_PUTrk_check_ = ibook.book1D(
      "PUTrk_avg_timing_err", "Average PU track MTD timing uncertainty;Time Error (ns);Counts", 200, 0, 0.1);
  meMuon_avg_error_vtx_check_ =
      ibook.book1D("Vtx_avg_timing_err", "Average vertex timing uncertainty;Time Error (ns);Counts", 200, 0, 0.1);

  meMuon_no_dt_check_ =
      ibook.book1D("Track_dt_info_check",
                   "Tracks dt check - ratio between tracks with (value 1) and without (value 0) timing info",
                   2,
                   0,
                   2);

  meTrk_genMatch_check_ = ibook.book1D(
      "Track_genMatch_info_check", "Check on tracks matched with a GenParticle (matched=1, non matched=0)", 2, 0, 2);

  // signal
  meMuonISO_Ntracks_Sig_EB_ = ibook.book1D("Muon_Iso_Ntracks_Sig_EB",
                                          "Number of tracks in isolation cone around muon track after basic cuts - "
                                          "Signal Barrel;Number of tracks;Counts",
                                          20,
                                          0,
                                          20);

  meMuonISO_chIso_Sig_EB_ = ibook.book1D(
      "Muon_chIso_sum_Sig_EB",
      "Track pT sum in isolation cone around muon track after basic cuts - Signal Barrel;p_{T} (GeV);Counts",
      nbin_2,
      0,
      20);

  meMuonISO_rel_chIso_Sig_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_Sig_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_1_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_1_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);

    meMuonISO_chIso_MTD_1_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_1_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_1_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_1_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
    // gen
    meMuonISO_Ntracks_gen_Sig_EB_ = ibook.book1D("Muon_Iso_Ntracks_gen_Sig_EB",
                                                "Number of tracks in isolation cone around muon track after basic "
                                                "cuts using genInfo - Signal Barrel;Number of tracks;Counts",
                                                20,
                                                0,
                                                20);

    meMuonISO_chIso_gen_Sig_EB_ = ibook.book1D("Muon_chIso_sum_gen_Sig_EB",
                                              "Track pT sum in isolation cone around muon track after basic cuts "
                                              "using genInfo - Signal Barrel;p_{T} (GeV);Counts",
                                              nbin_2,
                                              0,
                                              20);

    meMuonISO_rel_chIso_gen_Sig_EB_ = ibook.book1D("Muon_rel_chIso_sum_gen_Sig_EB",
                                                  "Track relative pT sum in isolation cone around muon track after "
                                                  "basic cuts using genInfo - Signal Barrel;Isolation;Counts",
                                                  nbin_1,
                                                  0,
                                                  4);

    meMuonISO_Ntracks_MTD_2_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_2_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);

    meMuonISO_chIso_MTD_2_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_2_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_2_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_2_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_3_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_3_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_3_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_3_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_3_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_3_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_4_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_4_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_4_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_4_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_4_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_4_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_5_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_5_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_5_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_5_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_5_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_5_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_6_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_6_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_6_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_6_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_6_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_6_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_7_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_7_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_7_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_7_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_7_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_7_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_1_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_1_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);

    meMuonISO_chIso_MTD_sim_1_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_1_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_1_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_1_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_2_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);

    meMuonISO_chIso_MTD_sim_2_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_2_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_2_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_2_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_3_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_3_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_3_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_3_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_4_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_4_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_4_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_4_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_5_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_5_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_5_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_5_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_5_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_5_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_6_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_6_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_6_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_6_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_6_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_6_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_7_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_7_Sig_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_7_Sig_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_7_Sig_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_7_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_7_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
  }
  meMuonISO_Ntracks_MTD_4sigma_Sig_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_4sigma_Sig_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 4 sigma compatibiliy - "
                   "Signal Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_4sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_4sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 4 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_4sigma_Sig_EB_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_4sigma_Sig_EB",
                   "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
                   "compatibiliy - Signal Barrel;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_3sigma_Sig_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_3sigma_Sig_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 3 sigma compatibiliy - "
                   "Signal Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_3sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_3sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 3 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_3sigma_Sig_EB_ = ibook.book1D("Muon_rel_chIso_sum_MTD_3sigma_Sig_EB",
                                                       "Track relative pT sum in isolation cone around muon track "
                                                       "after basic cuts with MTD - 3 sigma;Isolation;Counts"
                                                       "compatibiliy - Signal Barrel",
                                                       nbin_1,
                                                       0,
                                                       4);

  meMuonISO_Ntracks_MTD_2sigma_Sig_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_2sigma_Sig_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 2 sigma compatibiliy - "
                   "Signal Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_2sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_2sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 2 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_2sigma_Sig_EB_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_2sigma_Sig_EB",
                   "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
                   "compatibiliy - Signal Barrel;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuon_pt_tot_Sig_EB_ =
      ibook.book1D("Muon_pT_tot_Sig_EB", "Muon pT tot - Signal Barrel;p_{T} (GeV);Counts", 9, 10, 100);
      //ibook.book1D("Muon_pT_tot_Sig_EB", "Muon pT tot - Signal Barrel;p_{T} (GeV);Counts", 500, 0, 5000); //FIXME
  meMuon_pt_noMTD_Sig_EB_ =
      ibook.book1D("Muon_pT_noMTD_Sig_EB", "Muon pT noMTD - Signal Barrel;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_pt_sim_tot_Sig_EB_ =
      ibook.book1D("Muon_pT_sim_tot_Sig_EB", "Muon SIM pT tot - Signal Barrel;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_eta_tot_Sig_EB_ =
      ibook.book1D("Muon_eta_tot_Sig_EB", "Muon eta tot - Signal Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_eta_noMTD_Sig_EB_ =
      ibook.book1D("Muon_eta_noMTD_Sig_EB", "Muon eta noMTD - Signal Barrel;#eta;Counts", 64, 0., 3.2);

  meMuon_phi_tot_Sig_EB_ =
      ibook.book1D("Muon_phi_tot_Sig_EB", "Muon phi tot - Signal Barrel;#phi;Counts", 64, -3.2, 3.2);
  meMuon_phi_noMTD_Sig_EB_ =
      ibook.book1D("Muon_phi_noMTD_Sig_EB", "Muon phi noMTD - Signal Barrel;#phi;Counts", 64, -3.2, 3.2);

  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4sigma_Sig_EB",
                     "Number of tracks in isolation cone around muon track after basic cuts with MTD - 4 sigma "
                     "compatibiliy - Signal Barrel;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4sigma_Sig_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_4sigma_Sig_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD - 4 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_4sigma_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
        "compatibiliy - Signal Barrel;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3sigma_Sig_EB",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD  - 3 sigma compatibiliy - Signal Barrel;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3sigma_Sig_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_3sigma_Sig_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD - 3 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_3sigma_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
        "compatibiliy - Signal Barrel;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EB_ = ibook.book1D(
        "Muon_Iso_Ntracks_MTD_sim_2sigma_Sig_EB",
        "Tracks in isolation cone around muon track after basic cuts with MTD - 2 sigma compatibiliy - "
        "Signal Barrel;Number of tracks;Counts",
        20,
        0,
        20);
    meMuonISO_chIso_MTD_sim_2sigma_Sig_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_2sigma_Sig_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD - 2 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_2sigma_Sig_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
        "compatibiliy - Signal Barrel;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuon_pt_gen_Sig_EB_ =
        ibook.book1D("Muon_pT_gen_Sig_EB", "Muon pT genInfo - Signal Barrel;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_gen_Sig_EB_ =
        ibook.book1D("Muon_eta_gen_Sig_EB", "Muon eta genInfo - Signal Barrel;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_gen_Sig_EB_ =
        ibook.book1D("Muon_phi_gen_Sig_EB", "Muon phi genInfo - Signal Barrel;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_1_Sig_EB_ = ibook.book1D("Muon_pT_MTD_1_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_1_Sig_EB_ = ibook.book1D("Muon_eta_MTD_1_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_1_Sig_EB_ = ibook.book1D("Muon_phi_MTD_1_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_2_Sig_EB_ = ibook.book1D("Muon_pT_MTD_2_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_2_Sig_EB_ = ibook.book1D("Muon_eta_MTD_2_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_2_Sig_EB_ = ibook.book1D("Muon_phi_MTD_2_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_3_Sig_EB_ = ibook.book1D("Muon_pT_MTD_3_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_3_Sig_EB_ = ibook.book1D("Muon_eta_MTD_3_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_3_Sig_EB_ = ibook.book1D("Muon_phi_MTD_3_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_4_Sig_EB_ = ibook.book1D("Muon_pT_MTD_4_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_4_Sig_EB_ = ibook.book1D("Muon_eta_MTD_4_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_4_Sig_EB_ = ibook.book1D("Muon_phi_MTD_4_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_5_Sig_EB_ = ibook.book1D("Muon_pT_MTD_5_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_5_Sig_EB_ = ibook.book1D("Muon_eta_MTD_5_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_5_Sig_EB_ = ibook.book1D("Muon_phi_MTD_5_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_6_Sig_EB_ = ibook.book1D("Muon_pT_MTD_6_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_6_Sig_EB_ = ibook.book1D("Muon_eta_MTD_6_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_6_Sig_EB_ = ibook.book1D("Muon_phi_MTD_6_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_7_Sig_EB_ = ibook.book1D("Muon_pT_MTD_7_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_7_Sig_EB_ = ibook.book1D("Muon_eta_MTD_7_Sig_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_7_Sig_EB_ = ibook.book1D("Muon_phi_MTD_7_Sig_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_sim_MTD_1_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_1_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_2_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_2_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_3_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_3_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_4_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_4_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_5_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_5_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_6_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_6_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_7_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_7_Sig_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
  }

  meMuon_pt_MTD_4sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_MTD_4sigma_Sig_EB",
                   "Muon pT MTD - 4 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_4sigma_Sig_EB_ = ibook.book1D(
      "Muon_eta_MTD_4sigma_Sig_EB", "Muon eta MTD - 4 sigma compatibility - Signal Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_4sigma_Sig_EB_ = ibook.book1D("Muon_phi_MTD_4sigma_Sig_EB",
                                              "Muon phi MTD - 4 sigma compatibility - Signal Barrel;#phi;Counts",
                                              64,
                                              -3.2,
                                              3.2);

  meMuon_pt_MTD_3sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_MTD_3sigma_Sig_EB",
                   "Muon pT MTD - 3 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_3sigma_Sig_EB_ = ibook.book1D(
      "Muon_eta_MTD_3sigma_Sig_EB", "Muon eta MTD - 3 sigma compatibility - Signal Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_3sigma_Sig_EB_ = ibook.book1D("Muon_phi_MTD_3sigma_Sig_EB",
                                              "Muon phi MTD - 3 sigma compatibility - Signal Barrel;#phi;Counts",
                                              64,
                                              -3.2,
                                              3.2);

  meMuon_pt_MTD_2sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_MTD_2sigma_Sig_EB",
                   "Muon pT MTD - 2 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_2sigma_Sig_EB_ = ibook.book1D(
      "Muon_eta_MTD_2sigma_Sig_EB", "Muon eta MTD - 2 sigma compatibility - Signal Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_2sigma_Sig_EB_ = ibook.book1D("Muon_phi_MTD_2sigma_Sig_EB",
                                              "Muon phi MTD - 2 sigma compatibility - Signal Barrel;#phi;Counts",
                                              64,
                                              -3.2,
                                              3.2);

  meMuonISO_Ntracks_Sig_EE_ = ibook.book1D("Muon_Iso_Ntracks_Sig_EE",
                                          "Number of tracks in isolation cone around muon track after basic cuts - "
                                          "Signal Endcap;Number of tracks;Counts",
                                          20,
                                          0,
                                          20);
  meMuonISO_chIso_Sig_EE_ = ibook.book1D(
      "Muon_chIso_sum_Sig_EE",
      "Track pT sum in isolation cone around muon track after basic cuts - Signal Endcap;p_{T} (GeV);Counts",
      nbin_2,
      0,
      20);
  meMuonISO_rel_chIso_Sig_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_Sig_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts - Signal Endcap;Isolation;Counts",
      nbin_1,
      0,
      4);

  if (optionalPlots_) {
    meMuon_pt_sim_MTD_4sigma_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_4sigma_Sig_EB",
                     "Muon pT MTD SIM - 4 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_3sigma_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_3sigma_Sig_EB",
                     "Muon pT MTD SIM - 3 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_2sigma_Sig_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_2sigma_Sig_EB",
                     "Muon pT MTD SIM - 2 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);

    meMuonISO_Ntracks_MTD_1_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_1_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_1_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_1_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_1_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_1_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_2_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_2_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_2_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_2_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_2_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_2_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_gen_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_gen_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts using genInfo - Signal Endcap",
                     20,
                     0,
                     20);
    meMuonISO_chIso_gen_Sig_EE_ =
        ibook.book1D("Muon_chIso_sum_gen_Sig_EE",
                     "Track pT sum in isolation cone around muon track after basic cuts - Signal Endcap",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_gen_Sig_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_gen_Sig_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts - Signal Endcap",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_3_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_3_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_3_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_3_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_3_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_3_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_4_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_4_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_4_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_4_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_4_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_4_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_5_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_5_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_5_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_5_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_5_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_5_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_6_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_6_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_6_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_6_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_6_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_6_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_7_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_7_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_7_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_7_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_7_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_7_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_1_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_1_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_1_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_1_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_1_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_1_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_2_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_2_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_2_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_2_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_3_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_3_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_3_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_3_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_4_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_4_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_4_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_4_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_5_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_5_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_5_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_5_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_5_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_5_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_6_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_6_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_6_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_6_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_6_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_6_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_7_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_7_Sig_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_7_Sig_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_7_Sig_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_7_Sig_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_7_Sig_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
  }
  meMuonISO_Ntracks_MTD_4sigma_Sig_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_4sigma_Sig_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 4 sigma significance - "
                   "Signal Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_4sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_4sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 4 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_4sigma_Sig_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_4sigma_Sig_EE",
                   "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
                   "significance - Signal Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_3sigma_Sig_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_3sigma_Sig_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 3 sigma significance - "
                   "Signal Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_3sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_3sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 3 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_3sigma_Sig_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_3sigma_Sig_EE",
                   "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
                   "significance - Signal Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_2sigma_Sig_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_2sigma_Sig_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 2 sigma significance - "
                   "Signal Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_2sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_2sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 2 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_2sigma_Sig_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_2sigma_Sig_EE",
                   "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
                   "significance - Signal Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuon_pt_tot_Sig_EE_ =
      ibook.book1D("Muon_pT_tot_Sig_EE", "Muon pT tot - Signal Endcap;p_{T} (GeV);Counts", 9, 10, 100);
  meMuon_pt_noMTD_Sig_EE_ =
      ibook.book1D("Muon_pT_noMTD_Sig_EE", "Muon pT noMTD - Signal Endcap;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_pt_sim_tot_Sig_EE_ =
      ibook.book1D("Muon_pT_sim_tot_Sig_EE", "Muon pT tot - Signal Endcap;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_eta_tot_Sig_EE_ =
      ibook.book1D("Muon_eta_tot_Sig_EE", "Muon eta tot - Signal Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_eta_noMTD_Sig_EE_ =
      ibook.book1D("Muon_eta_noMTD_Sig_EE", "Muon eta noMTD - Signal Endcap;#eta;Counts", 64, 0., 3.2);

  meMuon_phi_tot_Sig_EE_ =
      ibook.book1D("Muon_phi_tot_Sig_EE", "Muon phi tot - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);
  meMuon_phi_noMTD_Sig_EE_ =
      ibook.book1D("Muon_phi_noMTD_Sig_EE", "Muon phi noMTD - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);

  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4sigma_Sig_EE",
                     "Number of tracks in isolation cone around muon track after basic cuts with MTD SIM - 4 sigma "
                     "significance - Signal Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4sigma_Sig_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_4sigma_Sig_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 4 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_4sigma_Sig_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 4 "
                     "sigma significance - Signal Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3sigma_Sig_EE",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 3 sigma significance - Signal Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3sigma_Sig_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_3sigma_Sig_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 3 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_3sigma_Sig_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 3 "
                     "sigma significance - Signal Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2sigma_Sig_EE",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 2 sigma significance - Signal Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2sigma_Sig_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_2sigma_Sig_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 2 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_2sigma_Sig_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 2 "
                     "sigma significance - Signal Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuon_pt_MTD_1_Sig_EE_ = ibook.book1D("Muon_pT_MTD_1_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_1_Sig_EE_ = ibook.book1D("Muon_eta_MTD_1_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_1_Sig_EE_ = ibook.book1D("Muon_phi_MTD_1_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);
    meMuon_pt_gen_Sig_EE_ =
        ibook.book1D("Muon_pT_gen_Sig_EE", "Muon pT genInfo - Signal Endcap;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_gen_Sig_EE_ =
        ibook.book1D("Muon_eta_gen_Sig_EE", "Muon eta genInfo - Signal Endcap;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_gen_Sig_EE_ =
        ibook.book1D("Muon_phi_gen_Sig_EE", "Muon phi genInfo - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_2_Sig_EE_ = ibook.book1D("Muon_pT_MTD_2_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_2_Sig_EE_ = ibook.book1D("Muon_eta_MTD_2_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_2_Sig_EE_ = ibook.book1D("Muon_phi_MTD_2_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_3_Sig_EE_ = ibook.book1D("Muon_pT_MTD_3_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_3_Sig_EE_ = ibook.book1D("Muon_eta_MTD_3_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_3_Sig_EE_ = ibook.book1D("Muon_phi_MTD_3_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_4_Sig_EE_ = ibook.book1D("Muon_pT_MTD_4_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_4_Sig_EE_ = ibook.book1D("Muon_eta_MTD_4_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_4_Sig_EE_ = ibook.book1D("Muon_phi_MTD_4_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_5_Sig_EE_ = ibook.book1D("Muon_pT_MTD_5_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_5_Sig_EE_ = ibook.book1D("Muon_eta_MTD_5_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_5_Sig_EE_ = ibook.book1D("Muon_phi_MTD_5_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_6_Sig_EE_ = ibook.book1D("Muon_pT_MTD_6_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_6_Sig_EE_ = ibook.book1D("Muon_eta_MTD_6_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_6_Sig_EE_ = ibook.book1D("Muon_phi_MTD_6_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_7_Sig_EE_ = ibook.book1D("Muon_pT_MTD_7_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_7_Sig_EE_ = ibook.book1D("Muon_eta_MTD_7_Sig_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_7_Sig_EE_ = ibook.book1D("Muon_phi_MTD_7_Sig_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_sim_MTD_1_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_1_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_2_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_2_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_3_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_3_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_4_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_4_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_5_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_5_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_6_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_6_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_7_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_7_Sig_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);

    meMuon_pt_sim_MTD_4sigma_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_4sigma_Sig_EE",
                     "Muon pT MTD SIM - 4 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_3sigma_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_3sigma_Sig_EE",
                     "Muon pT MTD SIM - 3 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_2sigma_Sig_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_2sigma_Sig_EE",
                     "Muon pT MTD SIM - 2 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
  }

  meMuon_pt_MTD_4sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_MTD_4sigma_Sig_EE",
                   "Muon pT MTD - 4 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_4sigma_Sig_EE_ = ibook.book1D(
      "Muon_eta_MTD_4sigma_Sig_EE", "Muon eta MTD - 4 sigma significance - Signal Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_4sigma_Sig_EE_ = ibook.book1D(
      "Muon_phi_MTD_4sigma_Sig_EE", "Muon phi MTD - 4 sigma significance - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_3sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_MTD_3sigma_Sig_EE",
                   "Muon pT MTD - 3 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_3sigma_Sig_EE_ = ibook.book1D(
      "Muon_eta_MTD_3sigma_Sig_EE", "Muon eta MTD - 3 sigma significance - Signal Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_3sigma_Sig_EE_ = ibook.book1D(
      "Muon_phi_MTD_3sigma_Sig_EE", "Muon phi MTD - 3 sigma significance - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_2sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_MTD_2sigma_Sig_EE",
                   "Muon pT MTD - 2 sigma significance - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_eta_MTD_2sigma_Sig_EE_ = ibook.book1D(
      "Muon_eta_MTD_2sigma_Sig_EE", "Muon eta MTD - 2 sigma significance - Signal Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_2sigma_Sig_EE_ = ibook.book1D(
      "Muon_phi_MTD_2sigma_Sig_EE", "Muon phi MTD - 2 sigma significance - Signal Endcap;#phi;Counts", 64, -3.2, 3.2);

  // background
  meMuonISO_Ntracks_Bkg_EB_ = ibook.book1D(
      "Muon_Iso_Ntracks_Bkg_EB",
      "Number of tracks in isolation cone around muon track after basic cuts - Bkg Barrel;Number of tracks;Counts",
      20,
      0,
      20);
  meMuonISO_chIso_Bkg_EB_ = ibook.book1D(
      "Muon_chIso_sum_Bkg_EB",
      "Track pT sum in isolation cone around muon track after basic cuts - Bkg Barrel;p_{T} (GeV);Counts",
      nbin_2,
      0,
      20);
  meMuonISO_rel_chIso_Bkg_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_Bkg_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts - Bkg Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_1_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_1_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_1_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_1_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_1_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_1_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_2_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_2_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_2_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_2_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_2_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_2_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
    meMuonISO_Ntracks_gen_Bkg_EB_ = ibook.book1D("Muon_Iso_Ntracks_gen_Bkg_EB",
                                                "Tracks in isolation cone around muon track after basic cuts using "
                                                "genInfo - Bkg Barrel;Number of tracks;Counts",
                                                20,
                                                0,
                                                20);
    meMuonISO_chIso_gen_Bkg_EB_ = ibook.book1D("Muon_chIso_sum_gen_Bkg_EB",
                                              "Track pT sum in isolation cone around muon track after basic cuts "
                                              "using genInfo - Bkg Barrel;p_{T} (GeV);Counts",
                                              nbin_2,
                                              0,
                                              20);
    meMuonISO_rel_chIso_gen_Bkg_EB_ = ibook.book1D("Muon_rel_chIso_sum_gen_Bkg_EB",
                                                  "Track relative pT sum in isolation cone around muon track after "
                                                  "basic cuts using genInfo - Bkg Barrel;Isolation;Counts",
                                                  nbin_1,
                                                  0,
                                                  4);
    meMuonISO_Ntracks_MTD_3_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_3_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_3_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_3_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_3_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_3_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_4_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_4_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_4_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_4_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_4_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_4_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_5_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_5_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_5_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_5_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_5_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_5_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_6_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_6_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_6_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_6_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_6_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_6_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_7_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_7_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_7_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_7_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_7_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_7_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_1_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_1_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_1_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_1_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_1_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_1_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_2_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_2_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_2_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_2_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_3_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_3_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_3_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_3_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_4_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_4_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_4_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_4_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_5_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_5_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_5_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_5_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_5_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_5_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_6_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_6_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_6_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_6_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_6_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_6_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_7_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_7_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_7_Bkg_EB_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_7_Bkg_EB",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_7_Bkg_EB_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_7_Bkg_EB",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
  }
  meMuonISO_Ntracks_MTD_4sigma_Bkg_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_4sigma_Bkg_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 4 sigma significance - "
                   "Bkg Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_4sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_4sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 4 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_4sigma_Bkg_EB_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 4 sigma significance - Bkg Barrel;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_3sigma_Bkg_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_3sigma_Bkg_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 3 sigma significance - "
                   "Bkg Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_3sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_3sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 3 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_3sigma_Bkg_EB_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 3 sigma significance - Bkg Barrel;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_2sigma_Bkg_EB_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_2sigma_Bkg_EB",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 2 sigma significance - "
                   "Bkg Barrel;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_2sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_2sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 2 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_2sigma_Bkg_EB_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 2 sigma significance - Bkg Barrel;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuon_pt_tot_Bkg_EB_ =
      ibook.book1D("Muon_pT_tot_Bkg_EB", "Muon pT tot - Bkg Barrel;p_{T} (GeV);Counts", 9, 10, 100);
  meMuon_pt_noMTD_Bkg_EB_ =
      ibook.book1D("Muon_pT_noMTD_Bkg_EB", "Muon pT noMTD - Bkg Barrel;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_pt_sim_tot_Bkg_EB_ =
      ibook.book1D("Muon_pT_sim_tot_Bkg_EB", "Muon pT tot - Bkg Barrel;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_eta_tot_Bkg_EB_ =
      ibook.book1D("Muon_eta_tot_Bkg_EB", "Muon eta tot - Bkg Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_eta_noMTD_Bkg_EB_ =
      ibook.book1D("Muon_eta_noMTD_Bkg_EB", "Muon eta noMTD - Bkg Barrel;#eta;Counts", 64, 0., 3.2);

  meMuon_phi_tot_Bkg_EB_ =
      ibook.book1D("Muon_phi_tot_Bkg_EB", "Muon phi tot - Bkg Barrel;#phi;#Counts", 64, -3.2, 3.2);
  meMuon_phi_noMTD_Bkg_EB_ =
      ibook.book1D("Muon_phi_noMTD_Bkg_EB", "Muon phi noMTD - Bkg Barrel;#phi;#Counts", 64, -3.2, 3.2);

  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4sigma_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 4 sigma significance - Bkg Barrel;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4sigma_Bkg_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_4sigma_Bkg_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 4 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EB_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_4sigma_Bkg_EB",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 4 "
                     "sigma significance - Bkg Barrel;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3sigma_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 3 sigma significance - Bkg Barrel;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3sigma_Bkg_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_3sigma_Bkg_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 3 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EB_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_3sigma_Bkg_EB",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 3 "
                     "sigma significance - Bkg Barrel;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EB_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2sigma_Bkg_EB",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 3 sigma significance - Bkg Barrel;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2sigma_Bkg_EB_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_2sigma_Bkg_EB",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 2 sigma significance - Bkg Barrel;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EB_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_2sigma_Bkg_EB",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 2 "
                     "sigma significance - Bkg Barrel;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuon_pt_MTD_1_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_1_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_1_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_1_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_1_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_1_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);
    meMuon_pt_gen_Bkg_EB_ =
        ibook.book1D("Muon_pT_gen_Bkg_EB", "Muon pT genInfo - Bkg Barrel;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_gen_Bkg_EB_ =
        ibook.book1D("Muon_eta_gen_Bkg_EB", "Muon eta genInfo - Bkg Barrel;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_gen_Bkg_EB_ =
        ibook.book1D("Muon_phi_gen_Bkg_EB", "Muon phi genInfo - Bkg Barrel;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_2_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_2_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_2_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_2_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_2_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_2_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_3_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_3_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_3_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_3_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_3_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_3_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_4_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_4_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_4_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_4_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_4_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_4_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_5_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_5_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_5_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_5_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_5_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_5_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_6_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_6_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_6_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_6_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_6_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_6_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_7_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_7_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_7_Bkg_EB_ = ibook.book1D("Muon_eta_MTD_7_Bkg_EB", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_7_Bkg_EB_ = ibook.book1D("Muon_phi_MTD_7_Bkg_EB", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_sim_MTD_1_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_1_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_2_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_2_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_3_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_3_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_4_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_4_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_5_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_5_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_6_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_6_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_7_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_7_Bkg_EB", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
  }
  meMuon_pt_MTD_4sigma_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_4sigma_Bkg_EB",
                                             "Muon pT MTD - 4 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_4sigma_Bkg_EB_ = ibook.book1D(
      "Muon_eta_MTD_4sigma_Bkg_EB", "Muon eta MTD - 4 sigma compatibility - Bkg Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_4sigma_Bkg_EB_ = ibook.book1D(
      "Muon_phi_MTD_4sigma_Bkg_EB", "Muon phi MTD - 4 sigma compatibility - Bkg Barrel;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_3sigma_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_3sigma_Bkg_EB",
                                             "Muon pT MTD - 3 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_3sigma_Bkg_EB_ = ibook.book1D(
      "Muon_eta_MTD_3sigma_Bkg_EB", "Muon eta MTD - 3 sigma compatibility - Bkg Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_3sigma_Bkg_EB_ = ibook.book1D(
      "Muon_phi_MTD_3sigma_Bkg_EB", "Muon phi MTD - 3 sigma compatibility - Bkg Barrel;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_2sigma_Bkg_EB_ = ibook.book1D("Muon_pT_MTD_2sigma_Bkg_EB",
                                             "Muon pT MTD - 2 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_2sigma_Bkg_EB_ = ibook.book1D(
      "Muon_eta_MTD_2sigma_Bkg_EB", "Muon eta MTD - 2 sigma compatibility - Bkg Barrel;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_2sigma_Bkg_EB_ = ibook.book1D(
      "Muon_phi_MTD_2sigma_Bkg_EB", "Muon phi MTD - 2 sigma compatibility - Bkg Barrel;#phi;Counts", 64, -3.2, 3.2);

  meMuonISO_Ntracks_Bkg_EE_ = ibook.book1D(
      "Muon_Iso_Ntracks_Bkg_EE",
      "Number of tracks in isolation cone around muon track after basic cuts - Bkg Endcap;Number of tracks;Counts",
      20,
      0,
      20);
  meMuonISO_chIso_Bkg_EE_ = ibook.book1D(
      "Muon_chIso_sum_Bkg_EE",
      "Track pT sum in isolation cone around muon track after basic cuts - Bkg Endcap;p_{T} (GeV);Counts",
      nbin_2,
      0,
      20);
  meMuonISO_rel_chIso_Bkg_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_Bkg_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts - Bkg Endcap;Isolation;Counts",
      nbin_1,
      0,
      4);
  if (optionalPlots_) {
    meMuon_pt_sim_MTD_4sigma_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_4sigma_Bkg_EB",
                     "Muon pT MTD SIM - 4 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_3sigma_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_3sigma_Bkg_EB",
                     "Muon pT MTD SIM - 3 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_2sigma_Bkg_EB_ =
        ibook.book1D("Muon_pT_sim_MTD_2sigma_Bkg_EB",
                     "Muon pT MTD SIM - 2 sigma compatibility - Bkg Barrel;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);

    meMuonISO_Ntracks_MTD_1_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_1_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_1_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_1_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_1_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_1_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_2_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_2_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_2_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_2_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_2_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_2_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
    meMuonISO_Ntracks_gen_Bkg_EE_ = ibook.book1D("Muon_Iso_Ntracks_gen_Bkg_EE",
                                                "Tracks in isolation cone around muon track after basic cuts using "
                                                "genInfo - Bkg Endcap;Number of tracks;Counts",
                                                20,
                                                0,
                                                20);
    meMuonISO_chIso_gen_Bkg_EE_ = ibook.book1D("Muon_chIso_sum_gen_Bkg_EE",
                                              "Track pT sum in isolation cone around muon track after basic cuts "
                                              "using genInfo - Bkg Endcap;p_{T} (GeV);Counts",
                                              nbin_2,
                                              0,
                                              20);
    meMuonISO_rel_chIso_gen_Bkg_EE_ = ibook.book1D("Muon_rel_chIso_sum_gen_Bkg_EE",
                                                  "Track relative pT sum in isolation cone around muon track after "
                                                  "basic cuts using genInfo - Bkg Endcap;Isolation;Counts",
                                                  nbin_1,
                                                  0,
                                                  4);

    meMuonISO_Ntracks_MTD_3_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_3_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_3_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_3_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_3_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_3_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_4_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_4_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_4_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_4_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_4_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_4_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_5_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_5_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_5_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_5_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_5_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_5_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_6_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_6_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_6_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_6_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_6_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_6_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_7_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_7_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_7_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_7_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_7_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_7_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_1_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_1_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_1_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_1_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_1_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_1_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_2_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_2_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_2_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_2_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_3_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_3_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_3_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_3_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_4_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_4_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_4_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_4_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_5_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_5_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_5_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_5_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_5_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_5_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_6_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_6_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_6_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_6_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_6_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_6_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);

    meMuonISO_Ntracks_MTD_sim_7_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_7_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic cuts with MTD;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_7_Bkg_EE_ = ibook.book1D(
        "Muon_chIso_sum_MTD_sim_7_Bkg_EE",
        "Track pT sum in isolation cone around muon track after basic cuts with MTD;p_{T} (GeV);Counts",
        nbin_2,
        0,
        20);
    meMuonISO_rel_chIso_MTD_sim_7_Bkg_EE_ = ibook.book1D(
        "Muon_rel_chIso_sum_MTD_sim_7_Bkg_EE",
        "Track relative pT sum in isolation cone around muon track after basic cuts with MTD;Isolation;Counts",
        nbin_1,
        0,
        4);
  }
  meMuonISO_Ntracks_MTD_4sigma_Bkg_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_4sigma_Bkg_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 4 sigma compatibility - "
                   "Bkg Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_4sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_4sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 4 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_4sigma_Bkg_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 4 sigma compatibility - Bkg Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_3sigma_Bkg_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_3sigma_Bkg_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 3 sigma compatibility - "
                   "Bkg Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_3sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_3sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 3 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_3sigma_Bkg_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 3 sigma compatibility - Bkg Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuonISO_Ntracks_MTD_2sigma_Bkg_EE_ =
      ibook.book1D("Muon_Iso_Ntracks_MTD_2sigma_Bkg_EE",
                   "Tracks in isolation cone around muon track after basic cuts with MTD - 2 sigma compatibility - "
                   "Bkg Endcap;Number of tracks;Counts",
                   20,
                   0,
                   20);
  meMuonISO_chIso_MTD_2sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_2sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after basic "
                   "cuts with MTD - 2 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_2sigma_Bkg_EE_ =
      ibook.book1D("Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE",
                   "Track relative pT sum in isolation cone around muon track "
                   "after basic cuts with MTD - 2 sigma compatibility - Bkg Endcap;Isolation;Counts",
                   nbin_1,
                   0,
                   4);

  meMuon_pt_tot_Bkg_EE_ =
      ibook.book1D("Muon_pT_tot_Bkg_EE", "Muon pT tot - Bkg Endcap;p_{T} (GeV);Counts", 9, 10, 100);
  meMuon_pt_noMTD_Bkg_EE_ =
      ibook.book1D("Muon_pT_noMTD_Bkg_EE", "Muon pT noMTD - Bkg Endcap;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_pt_sim_tot_Bkg_EE_ =
      ibook.book1D("Muon_pT_sim_tot_Bkg_EE", "Muon pT tot - Bkg Endcap;p_{T} (GeV);Counts", 9, 10, 100);

  meMuon_eta_tot_Bkg_EE_ =
      ibook.book1D("Muon_eta_tot_Bkg_EE", "Muon eta tot - Bkg Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_eta_noMTD_Bkg_EE_ =
      ibook.book1D("Muon_eta_noMTD_Bkg_EE", "Muon eta noMTD - Bkg Endcap;#eta;Counts", 64, 0., 3.2);

  meMuon_phi_tot_Bkg_EE_ =
      ibook.book1D("Muon_phi_tot_Bkg_EE", "Muon phi tot - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);
  meMuon_phi_noMTD_Bkg_EE_ =
      ibook.book1D("Muon_phi_noMTD_Bkg_EE", "Muon phi noMTD - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);
  if (optionalPlots_) {
    meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_4sigma_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 4 sigma compatibility - Bkg Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_4sigma_Bkg_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_4sigma_Bkg_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 4 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_4sigma_Bkg_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 4 "
                     "sigma compatibility - Bkg Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_3sigma_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 3 sigma compatibility - Bkg Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_3sigma_Bkg_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_3sigma_Bkg_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 3 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_3sigma_Bkg_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 3 "
                     "sigma compatibility - Bkg Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EE_ =
        ibook.book1D("Muon_Iso_Ntracks_MTD_sim_2sigma_Bkg_EE",
                     "Tracks in isolation cone around muon track after basic "
                     "cuts with MTD SIM - 2 sigma compatibility - Bkg Endcap;Number of tracks;Counts",
                     20,
                     0,
                     20);
    meMuonISO_chIso_MTD_sim_2sigma_Bkg_EE_ =
        ibook.book1D("Muon_chIso_sum_MTD_sim_2sigma_Bkg_EE",
                     "Track pT sum in isolation cone around muon track after "
                     "basic cuts with MTD SIM - 2 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                     nbin_2,
                     0,
                     20);
    meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EE_ =
        ibook.book1D("Muon_rel_chIso_sum_MTD_sim_2sigma_Bkg_EE",
                     "Track relative pT sum in isolation cone around muon track after basic cuts with MTD SIM - 2 "
                     "sigma compatibility - Bkg Endcap;Isolation;Counts",
                     nbin_1,
                     0,
                     4);

    meMuon_pt_MTD_1_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_1_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_1_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_1_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_1_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_1_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);
    meMuon_pt_gen_Bkg_EE_ =
        ibook.book1D("Muon_pT_gen_Bkg_EE", "Muon pT genInfo - Bkg Endcap;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_gen_Bkg_EE_ =
        ibook.book1D("Muon_eta_gen_Bkg_EE", "Muon eta genInfo - Bkg Endcap;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_gen_Bkg_EE_ =
        ibook.book1D("Muon_phi_gen_Bkg_EE", "Muon phi genInfo - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_2_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_2_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_2_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_2_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_2_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_2_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_3_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_3_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_3_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_3_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_3_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_3_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_4_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_4_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_4_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_4_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_4_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_4_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_5_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_5_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_5_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_5_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_5_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_5_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_6_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_6_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_6_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_6_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_6_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_6_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_MTD_7_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_7_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_eta_MTD_7_Bkg_EE_ = ibook.book1D("Muon_eta_MTD_7_Bkg_EE", "Muon eta MTD;#eta;Counts", 64, 0., 3.2);
    meMuon_phi_MTD_7_Bkg_EE_ = ibook.book1D("Muon_phi_MTD_7_Bkg_EE", "Muon phi MTD;#phi;Counts", 64, -3.2, 3.2);

    meMuon_pt_sim_MTD_1_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_1_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_2_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_2_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_3_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_3_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_4_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_4_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_5_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_5_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_6_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_6_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
    meMuon_pt_sim_MTD_7_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_7_Bkg_EE", "Muon pT MTD;p_{T} (GeV);Counts", 9, 10, 100);
  }

  meMuon_pt_MTD_4sigma_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_4sigma_Bkg_EE",
                                             "Muon pT MTD - 4 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_4sigma_Bkg_EE_ = ibook.book1D(
      "Muon_eta_MTD_4sigma_Bkg_EE", "Muon eta MTD - 4 sigma compatibility - Bkg Endcapi;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_4sigma_Bkg_EE_ = ibook.book1D(
      "Muon_phi_MTD_4sigma_Bkg_EE", "Muon phi MTD - 4 sigma compatibility - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_3sigma_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_3sigma_Bkg_EE",
                                             "Muon pT MTD - 3 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_3sigma_Bkg_EE_ = ibook.book1D(
      "Muon_eta_MTD_3sigma_Bkg_EE", "Muon eta MTD - 3 sigma compatibility - Bkg Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_3sigma_Bkg_EE_ = ibook.book1D(
      "Muon_phi_MTD_3sigma_Bkg_EE", "Muon phi MTD - 3 sigma compatibility - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);

  meMuon_pt_MTD_2sigma_Bkg_EE_ = ibook.book1D("Muon_pT_MTD_2sigma_Bkg_EE",
                                             "Muon pT MTD - 2 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                                             9,
                                             10,
                                             100);
  meMuon_eta_MTD_2sigma_Bkg_EE_ = ibook.book1D(
      "Muon_eta_MTD_2sigma_Bkg_EE", "Muon eta MTD - 2 sigma compatibility - Bkg Endcap;#eta;Counts", 64, 0., 3.2);
  meMuon_phi_MTD_2sigma_Bkg_EE_ = ibook.book1D(
      "Muon_phi_MTD_2sigma_Bkg_EE", "Muon phi MTD - 2 sigma compatibility - Bkg Endcap;#phi;Counts", 64, -3.2, 3.2);

  if (optionalPlots_) {
    meMuon_pt_sim_MTD_4sigma_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_4sigma_Bkg_EE",
                     "Muon pT MTD SIM - 4 sigma compatibility - Bkg Endcap;p_{T} (GeV);Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_3sigma_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_3sigma_Bkg_EE",
                     "Muon pT MTD SIM - 3 sigma compatibility - Bkg Endcap;#eta;Counts",
                     9,
                     10,
                     100);
    meMuon_pt_sim_MTD_2sigma_Bkg_EE_ =
        ibook.book1D("Muon_pT_sim_MTD_2sigma_Bkg_EE",
                     "Muon pT MTD SIM - 2 sigma compatibility - Bkg Endcap;#phi;Counts",
                     9,
                     10,
                     100);
  }
  meMuon_dt_general_pT_1 = ibook.book1D("Iso_track_dt_general_pT_10_20",
                                       "Iso cone track dt distribution in pT bin 10-20 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_2 = ibook.book1D("Iso_track_dt_general_pT_20_30",
                                       "Iso cone track dt distribution in pT bin 20-30 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_3 = ibook.book1D("Iso_track_dt_general_pT_30_40",
                                       "Iso cone track dt distribution in pT bin 30-40 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_4 = ibook.book1D("Iso_track_dt_general_pT_40_50",
                                       "Iso cone track dt distribution in pT bin 40-50 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_5 = ibook.book1D("Iso_track_dt_general_pT_50_60",
                                       "Iso cone track dt distribution in pT bin 50-60 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_6 = ibook.book1D("Iso_track_dt_general_pT_60_70",
                                       "Iso cone track dt distribution in pT bin 60-70 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_7 = ibook.book1D("Iso_track_dt_general_pT_70_80",
                                       "Iso cone track dt distribution in pT bin 70-80 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_8 = ibook.book1D("Iso_track_dt_general_pT_80_90",
                                       "Iso cone track dt distribution in pT bin 80-90 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);
  meMuon_dt_general_pT_9 = ibook.book1D("Iso_track_dt_general_pT_90_100",
                                       "Iso cone track dt distribution in pT bin 90-100 GeV ;Time (ns);Counts",
                                       100,
                                       0,
                                       1);

  meMuon_dtSignif_general_pT_1 = ibook.book1D("Iso_track_dtSignif_general_pT_10_20",
                                             "Iso cone track dt distribution in pT bin 10-20 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_2 = ibook.book1D("Iso_track_dtSignif_general_pT_20_30",
                                             "Iso cone track dt distribution in pT bin 20-30 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_3 = ibook.book1D("Iso_track_dtSignif_general_pT_30_40",
                                             "Iso cone track dt distribution in pT bin 30-40 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_4 = ibook.book1D("Iso_track_dtSignif_general_pT_40_50",
                                             "Iso cone track dt distribution in pT bin 40-50 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_5 = ibook.book1D("Iso_track_dtSignif_general_pT_50_60",
                                             "Iso cone track dt distribution in pT bin 50-60 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_6 = ibook.book1D("Iso_track_dtSignif_general_pT_60_70",
                                             "Iso cone track dt distribution in pT bin 60-70 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_7 = ibook.book1D("Iso_track_dtSignif_general_pT_70_80",
                                             "Iso cone track dt distribution in pT bin 70-80 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_8 = ibook.book1D("Iso_track_dtSignif_general_pT_80_90",
                                             "Iso cone track dt distribution in pT bin 80-90 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);
  meMuon_dtSignif_general_pT_9 = ibook.book1D("Iso_track_dtSignif_general_pT_90_100",
                                             "Iso cone track dt distribution in pT bin 90-100 GeV ;Time (ns);Counts",
                                             nbin_1,
                                             0,
                                             10);

  meMuon_dt_general_eta_1 = ibook.book1D(
      "Iso_track_dt_general_eta_0_05", "Iso cone track dt distribution in eta bin 0.0-0.5 ;Time (ns);Counts", 100, 0, 1);
  meMuon_dt_general_eta_2 = ibook.book1D("Iso_track_dt_general_eta_05_10",
                                        "Iso cone track dt distribution in eta bin 0.5-1.0 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);
  meMuon_dt_general_eta_3 = ibook.book1D("Iso_track_dt_general_eta_10_15",
                                        "Iso cone track dt distribution in eta bin 1.0-1.5 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);
  meMuon_dt_general_eta_4 = ibook.book1D("Iso_track_dt_general_eta_15_20",
                                        "Iso cone track dt distribution in eta bin 1.5-2.0 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);
  meMuon_dt_general_eta_5 = ibook.book1D("Iso_track_dt_general_eta_20_24",
                                        "Iso cone track dt distribution in eta bin 2.0-2.4 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);
  meMuon_dt_general_eta_6 = ibook.book1D("Iso_track_dt_general_eta_24_27",
                                        "Iso cone track dt distribution in eta bin 2.4-2.7 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);
  meMuon_dt_general_eta_7 = ibook.book1D("Iso_track_dt_general_eta_27_30",
                                        "Iso cone track dt distribution in eta bin 2.7-3.0 ;Time (ns);Counts",
                                        100,
                                        0,
                                        1);

  meMuon_dtSignif_general_eta_1 = ibook.book1D("Iso_track_dtSignif_general_eta_0_05",
                                              "Iso cone track dt distribution in eta bin 0.0-0.5 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_2 = ibook.book1D("Iso_track_dtSignif_general_eta_05_10",
                                              "Iso cone track dt distribution in eta bin 0.5-1.0 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_3 = ibook.book1D("Iso_track_dtSignif_general_eta_10_15",
                                              "Iso cone track dt distribution in eta bin 1.0-1.5 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_4 = ibook.book1D("Iso_track_dtSignif_general_eta_15_20",
                                              "Iso cone track dt distribution in eta bin 1.5-2.0 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_5 = ibook.book1D("Iso_track_dtSignif_general_eta_20_24",
                                              "Iso cone track dt distribution in eta bin 2.0-2.4 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_6 = ibook.book1D("Iso_track_dtSignif_general_eta_24_27",
                                              "Iso cone track dt distribution in eta bin 2.4-2.7 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);
  meMuon_dtSignif_general_eta_7 = ibook.book1D("Iso_track_dtSignif_general_eta_27_30",
                                              "Iso cone track dt distribution in eta bin 2.7-3.0 ;Time (ns);Counts",
                                              nbin_1,
                                              0,
                                              10);

  // test
    // GEN CASE
  meMuon_pt_gen_MTD_2sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_2sigma_Sig_EB",
                   "Muon pT MTD GEN - 2 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_3sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_3sigma_Sig_EB",
                   "Muon pT MTD GEN - 3 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_4sigma_Sig_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_4sigma_Sig_EB",
                   "Muon pT MTD GEN - 4 sigma compatibility - Signal Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_2sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_2sigma_Sig_EE",
                   "Muon pT MTD GEN - 2 sigma compatibility - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_3sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_3sigma_Sig_EE",
                   "Muon pT MTD GEN - 3 sigma compatibility - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_4sigma_Sig_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_4sigma_Sig_EE",
                   "Muon pT MTD GEN - 4 sigma compatibility - Signal Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_2sigma_Bkg_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_2sigma_Bkg_EB",
                   "Muon pT MTD GEN - 2 sigma compatibility - Background Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_3sigma_Bkg_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_3sigma_Bkg_EB",
                   "Muon pT MTD GEN - 3 sigma compatibility - Background Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_4sigma_Bkg_EB_ =
      ibook.book1D("Muon_pT_gen_MTD_4sigma_Bkg_EB",
                   "Muon pT MTD GEN - 4 sigma compatibility - Background Barrel;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_2sigma_Bkg_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_2sigma_Bkg_EE",
                   "Muon pT MTD GEN - 2 sigma compatibility - Background Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_3sigma_Bkg_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_3sigma_Bkg_EE",
                   "Muon pT MTD GEN - 3 sigma compatibility - Background Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);
  meMuon_pt_gen_MTD_4sigma_Bkg_EE_ =
      ibook.book1D("Muon_pT_gen_MTD_4sigma_Bkg_EE",
                   "Muon pT MTD GEN - 4 sigma compatibility - Background Endcap;p_{T} (GeV);Counts",
                   9,
                   10,
                   100);

  meMuonISO_chIso_MTD_gen_2sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_2sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 2 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_3sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_3sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 3 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_4sigma_Sig_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_4sigma_Sig_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 4 sigma compatibiliy - Signal Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_2sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_2sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 2 sigma compatibiliy - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_3sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_3sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 3 sigma compatibiliy - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_4sigma_Sig_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_4sigma_Sig_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 4 sigma compatibiliy - Signal Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_2sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_2sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 2 sigma compatibiliy - Background Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_3sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_3sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 3 sigma compatibiliy - Background Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_4sigma_Bkg_EB_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_4sigma_Bkg_EB",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 4 sigma compatibiliy - Background Barrel;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_2sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_2sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 2 sigma compatibiliy - Background Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_3sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_3sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 3 sigma compatibiliy - Background Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_chIso_MTD_gen_4sigma_Bkg_EE_ =
      ibook.book1D("Muon_chIso_sum_MTD_gen_4sigma_Bkg_EE",
                   "Track pT sum in isolation cone around muon track after "
                   "basic cuts with MTD - 4 sigma compatibiliy - Background Endcap;p_{T} (GeV);Counts",
                   nbin_2,
                   0,
                   20);
  meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EB_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EB",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
      "compatibiliy - Signal Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 2 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 3 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);
  meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EE_ = ibook.book1D(
      "Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EE",
      "Track relative pT sum in isolation cone around muon track after basic cuts with MTD - 4 sigma "
      "compatibiliy - Background Barrel;Isolation;Counts",
      nbin_1,
      0,
      4);

  // Several tests
  meMuonISO_trk_genMatched_Sig_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Sig", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Signal", 2, 0, 2);
  meMuonISO_trk_genMatched_Sig_EB_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Sig_EB", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Signal Barrel", 2, 0, 2);
  meMuonISO_trk_genMatched_Sig_EE_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Sig_EE", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Signal Endcap", 2, 0, 2);
  meMuonISO_trk_genMatched_Bkg_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Bkg", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Background", 2, 0, 2);
  meMuonISO_trk_genMatched_Bkg_EB_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Bkg_EB", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Background Barrel", 2, 0, 2);
  meMuonISO_trk_genMatched_Bkg_EE_ = ibook.book1D(
      "Muon_Iso_track_genMatched_Bkg_EE", "Check on tracks matched with a GenParticle (matched=1, non matched=0) - Background Endcap", 2, 0, 2);

  meMuonISO_pTdiff_reco_tracker_muon_Sig_ = ibook.book1D("Muon_Iso_pTdiff_reco_tracker_muon_Sig",
                                                "pT difference between reco and tracker muons - Signal;pT (GeV);Counts",
                                                200,
                                                0,
                                                1000);
  meMuonISO_pTdiff_reco_tracker_muon_Bkg_ = ibook.book1D("Muon_Iso_pTdiff_reco_tracker_muon_Bkg",
                                                "pT difference between reco and tracker muons - Background;pT (GeV);Counts",
                                                200,
                                                0,
                                                1000);

  meMuonISO_pT_muon_sim_ = ibook.book1D("Muon_Iso_pT_muon_sim",
                                                "pT of sim muon ;pT (GeV);Counts",
                                                101,
                                                -1,
                                                100);
  meMuonISO_pT_trk_sim_ = ibook.book1D("Muon_Iso_pT_track_sim",
                                                "pT of sim track ;pT (GeV);Counts",
                                                101,
                                                -1,
                                                100);
  meMuonISO_dz_muon_ = ibook.book1D("Muon_Iso_dz_muon_reco",
                                                "dz (muon, PV) distribution for reco muons ;dz (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_dz_muon_EB_ = ibook.book1D("Muon_Iso_dz_muon_reco_EB",
                                                "dz (muon, PV) distribution for reco muons - Barrel;dz (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_dz_muon_EE_ = ibook.book1D("Muon_Iso_dz_muon_reco_EE",
                                                "dz (muon, PV) distribution for reco muons - Endcap;dz (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_dxy_muon_ = ibook.book1D("Muon_Iso_dxy_muon_reco",
                                                "dxy (muon, PV) distribution for reco muons ;dxy (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_dxy_muon_EB_ = ibook.book1D("Muon_Iso_dxy_muon_reco_EB",
                                                "dxy (muon, PV) distribution for reco muons - Barrel;dxy (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_dxy_muon_EE_ = ibook.book1D("Muon_Iso_dxy_muon_reco_EE",
                                                "dxy (muon, PV) distribution for reco muons - Endcap;dxy (cm);Counts",
                                                40,
                                                0,
                                                2);
  meMuonISO_time_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_time_muon_reco_Sig_EB",
                                                "Time distribution for reco muons - Signal Barrel;t^{reco}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_time_muon_reco_Sig_EE",
                                                "Time distribution for reco muons - Signal Endcap;t^{reco}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_time_muon_reco_Bkg_EB",
                                                "Time distribution for reco muons - Background Barrel;t^{reco}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_time_muon_reco_Bkg_EE",
                                                "Time distribution for reco muons - Background Endcap;t^{reco}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_time_track_reco_Sig_EB",
                                                "Time distribution for reco tracks - Signal Barrel;t^{reco}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_time_track_reco_Sig_EE",
                                                "Time distribution for reco tracks - Signal Endcap;t^{reco}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_time_track_reco_Bkg_EB",
                                                "Time distribution for reco tracks - Background Barrel;t^{reco}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_time_track_reco_Bkg_EE",
                                                "Time distribution for reco tracks - Background Endcap;t^{reco}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_tErr_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_tErr_muon_reco_Sig_EB",
                                                "Error of time distribution for reco muons - Signal Barrel;t^{reco}_{muon} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_tErr_muon_reco_Sig_EE",
                                                "Error of time distribution for reco muons - Signal Endcap;t^{reco}_{muon} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_tErr_muon_reco_Bkg_EB",
                                                "Error of time distribution for reco muons - Background Barrel;t^{reco}_{muon} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_tErr_muon_reco_Bkg_EE",
                                                "Error of time distribution for reco muons - Background Endcap;t^{reco}_{muon} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_tErr_track_reco_Sig_EB",
                                                "Error of time distribution for reco tracks - Signal Barrel;t^{reco}_{track} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_tErr_track_reco_Sig_EE",
                                                "Error of time distribution for reco tracks - Signal Endcap;t^{reco}_{track} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_tErr_track_reco_Bkg_EB",
                                                "Error of time distribution for reco tracks - Background Barrel;t^{reco}_{track} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_tErr_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_tErr_track_reco_Bkg_EE",
                                                "Error of time distribution for reco tracks - Background Endcap;t^{reco}_{track} (ns);Counts",
                                                60,
                                                -1,
                                                2);
  meMuonISO_time_muon_sim_Sig_EB_ = ibook.book1D("Muon_Iso_time_muon_sim_Sig_EB",
                                                "Time distribution for sim muons - Signal Barrel;t^{sim}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_sim_Sig_EE_ = ibook.book1D("Muon_Iso_time_muon_sim_Sig_EE",
                                                "Time distribution for sim muons - Signal Endcap;t^{sim}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_time_muon_sim_Bkg_EB",
                                                "Time distribution for sim muons - Background Barrel;t^{sim}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_muon_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_time_muon_sim_Bkg_EE",
                                                "Time distribution for sim muons - Background Endcap;t^{sim}_{muon} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_time_track_sim_Sig_EB",
                                                "Time distribution for sim tracks - Signal Barrel;t^{sim}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_time_track_sim_Sig_EE",
                                                "Time distribution for sim tracks - Signal Endcap;t^{sim}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_time_track_sim_Bkg_EB",
                                                "Time distribution for sim tracks - Background Barrel;t^{sim}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_time_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_time_track_sim_Bkg_EE",
                                                "Time distribution for sim tracks - Background Endcap;t^{sim}_{track} (ns);Counts",
                                                80,
                                                -2,
                                                2);
  meMuonISO_has_time_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_time_reco_Sig_EB",
                                                "Check whether reco muon has timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_time_reco_Sig_EE",
                                                "Check whether reco muon has timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_time_reco_Bkg_EB",
                                                "Check whether reco muon has timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_time_reco_Bkg_EE",
                                                "Check whether reco muon has timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_track_has_time_reco_Sig_EB",
                                                "Check whether reco track has timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_track_has_time_reco_Sig_EE",
                                                "Check whether reco track has timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_track_has_time_reco_Bkg_EB",
                                                "Check whether reco track has timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_track_has_time_reco_Bkg_EE",
                                                "Check whether reco track has timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_tErr_reco_Sig_EB",
                                                "Check whether reco muon has error of timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_tErr_reco_Sig_EE",
                                                "Check whether reco muon has error of timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_tErr_reco_Bkg_EB",
                                                "Check whether reco muon has error of timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_tErr_reco_Bkg_EE",
                                                "Check whether reco muon has error of timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_track_has_tErr_reco_Sig_EB",
                                                "Check whether reco track has error of timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_track_has_tErr_reco_Sig_EE",
                                                "Check whether reco track has error of timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_track_has_tErr_reco_Bkg_EB",
                                                "Check whether reco track has error of timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_tErr_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_track_has_tErr_reco_Bkg_EE",
                                                "Check whether reco track has error of timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_mva_reco_Sig_EB",
                                                "Check whether reco muon track has mva score - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_mva_reco_Sig_EE",
                                                "Check whether reco muon track has mva score - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_mva_reco_Bkg_EB",
                                                "Check whether reco muon track has mva score - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_mva_reco_Bkg_EE",
                                                "Check whether reco muon track has mva score - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_whether_track_has_mva_reco_Sig_EB",
                                                "Check whether reco track has mva score - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_whether_track_has_mva_reco_Sig_EE",
                                                "Check whether reco track has mva score - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_track_has_mva_reco_Bkg_EB",
                                                "Check whether reco track has mva score - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_mva_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_track_has_mva_reco_Bkg_EE",
                                                "Check whether reco track has mva score - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_sim_Sig_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_time_sim_Sig_EB",
                                                "Check whether sim muon has timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_sim_Sig_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_time_sim_Sig_EE",
                                                "Check whether sim muon has timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_muon_has_time_sim_Bkg_EB",
                                                "Check whether sim muon has timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_muon_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_muon_has_time_sim_Bkg_EE",
                                                "Check whether sim muon has timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_whether_track_has_time_sim_Sig_EB",
                                                "Check whether sim track has timing information - Signal Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_whether_track_has_time_sim_Sig_EE",
                                                "Check whether sim track has timing information - Signal Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_whether_track_has_time_sim_Bkg_EB",
                                                "Check whether sim track has timing information - Background Barrel",
                                                2,
                                                0,
                                                2);
  meMuonISO_has_time_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_whether_track_has_time_sim_Bkg_EE",
                                                "Check whether sim track has timing information - Background Endcap",
                                                2,
                                                0,
                                                2);
  meMuonISO_dt_muon_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_dt_muon_track_reco_Sig_EB",
                                                "dt distribution for reco track and reco muon - Signal Barrel;|t^{reco}_{track}-t^{reco}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_dt_muon_track_reco_Sig_EE",
                                                "dt distribution for reco track and reco muon - Signal Endcap;|t^{reco}_{track}-t^{reco}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_dt_muon_track_reco_Bkg_EB",
                                                "dt distribution for reco track and reco muon - Background Barrel;|t^{reco}_{track}-t^{reco}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_dt_muon_track_reco_Bkg_EE",
                                                "dt distribution for reco track and reco muon - Background Endcap;|t^{reco}_{track}-t^{reco}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_dt_muon_track_sim_Sig_EB",
                                                "dt distribution for sim track and sim muon - Signal Barrel;|t^{sim}_{track}-t^{sim}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_dt_muon_track_sim_Sig_EE",
                                                "dt distribution for sim track and sim muon - Signal Endcap;|t^{sim}_{track}-t^{sim}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_dt_muon_track_sim_Bkg_EB",
                                                "dt distribution for sim track and sim muon - Background Barrel;|t^{sim}_{track}-t^{sim}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_muon_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_dt_muon_track_sim_Bkg_EE",
                                                "dt distribution for sim track and sim muon - Background Endcap;|t^{sim}_{track}-t^{sim}_{muon}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_dt_PV_track_reco_Sig_EB",
                                                "dt distribution for reco track and PV - Signal Barrel;|t^{reco}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_dt_PV_track_reco_Sig_EE",
                                                "dt distribution for reco track and PV - Signal Endcap;|t^{reco}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_dt_PV_track_reco_Bkg_EB",
                                                "dt distribution for reco track and PV - Background Barrel;|t^{reco}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_dt_PV_track_reco_Bkg_EE",
                                                "dt distribution for reco track and PV - Background Endcap;|t^{reco}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_dt_PV_track_sim_Sig_EB",
                                                "dt distribution for sim track and PV - Signal Barrel;|t^{sim}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_dt_PV_track_sim_Sig_EE",
                                                "dt distribution for sim track and PV - Signal Endcap;|t^{sim}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_dt_PV_track_sim_Bkg_EB",
                                                "dt distribution for sim track and PV - Background Barrel;|t^{sim}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);
  meMuonISO_dt_PV_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_dt_PV_track_sim_Bkg_EE",
                                                "dt distribution for sim track and PV - Background Endcap;|t^{sim}_{track}-t_{PV}| (ns);Counts",
                                                30,
                                                0,
                                                3);

  meMuonISO_dtSig_muon_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_Sig_EB",
                                                "dtSig distribution for reco track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_Sig_EE",
                                                "dtSig distribution for reco track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_Bkg_EB",
                                                "dtSig distribution for reco track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_Bkg_EE",
                                                "dtSig distribution for reco track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_Sig_EB",
                                                "dtSig distribution for sim track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_Sig_EE",
                                                "dtSig distribution for sim track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_Bkg_EB",
                                                "dtSig distribution for sim track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_Bkg_EE",
                                                "dtSig distribution for sim track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_genMatched_Sig_EB",
                                                "dtSig distribution for genMatched reco track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_genMatched_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_genMatched_Sig_EE",
                                                "dtSig distribution for genMatched reco track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_genMatched_Bkg_EB",
                                                "dtSig distribution for genMatched reco track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_genMatched_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_genMatched_Bkg_EE",
                                                "dtSig distribution for genMatched reco track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_genMatched_Sig_EB",
                                                "dtSig distribution for genMatched sim track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_genMatched_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_genMatched_Sig_EE",
                                                "dtSig distribution for genMatched sim track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_genMatched_Bkg_EB",
                                                "dtSig distribution for genMatched sim track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_genMatched_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_genMatched_Bkg_EE",
                                                "dtSig distribution for genMatched sim track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_Sig_EB",
                                                "dtSig distribution for reco track and PV - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_Sig_EE",
                                                "dtSig distribution for reco track and PV - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_Bkg_EB",
                                                "dtSig distribution for reco track and PV - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_Bkg_EE",
                                                "dtSig distribution for reco track and PV - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_Sig_EB",
                                                "dtSig distribution for sim track and PV - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_Sig_EE",
                                                "dtSig distribution for sim track and PV - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_Bkg_EB",
                                                "dtSig distribution for sim track and PV - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_Bkg_EE",
                                                "dtSig distribution for sim track and PV - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_genMatched_Sig_EB",
                                                "dtSig distribution for genMatched reco track and PV - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_genMatched_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_genMatched_Sig_EE",
                                                "dtSig distribution for genMatched reco track and PV - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_genMatched_Bkg_EB",
                                                "dtSig distribution for genMatched reco track and PV - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_genMatched_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_genMatched_Bkg_EE",
                                                "dtSig distribution for genMatched reco track and PV - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_genMatched_Sig_EB",
                                                "dtSig distribution for genMatched sim track and PV - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_genMatched_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_genMatched_Sig_EE",
                                                "dtSig distribution for genMatched sim track and PV - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_genMatched_Bkg_EB",
                                                "dtSig distribution for genMatched sim track and PV - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_genMatched_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_genMatched_Bkg_EE",
                                                "dtSig distribution for genMatched sim track and PV - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);

  meMuonISO_mva_muon_reco_Sig_ = ibook.book1D("Muon_Iso_mva_muon_reco_Sig",
                                                "MVA score of reco muons ;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_muon_reco_Sig_EB_ = ibook.book1D("Muon_Iso_mva_muon_reco_Sig_EB",
                                                "MVA score of reco muons - Signal Barrel;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_muon_reco_Sig_EE_ = ibook.book1D("Muon_Iso_mva_muon_reco_Sig_EE",
                                                "MVA score of reco muons - Signal Endcal;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_muon_reco_Bkg_ = ibook.book1D("Muon_Iso_mva_muon_reco_Bkg",
                                                "MVA score of reco muons ;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_muon_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_mva_muon_reco_Bkg_EB",
                                                "MVA score of reco muons - Background Barrel;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_muon_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_mva_muon_reco_Bkg_EE",
                                                "MVA score of reco muons - Background Endcap;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Sig_ = ibook.book1D("Muon_Iso_mva_track_reco_Sig",
                                                "MVA score of reco tracks ;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Sig_EB_ = ibook.book1D("Muon_Iso_mva_track_reco_Sig_EB",
                                                "MVA score of reco tracks Signal - Barrel;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Sig_EE_ = ibook.book1D("Muon_Iso_mva_track_reco_Sig_EE",
                                                "MVA score of reco tracks Signal - Endcap;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Bkg_ = ibook.book1D("Muon_Iso_mva_track_reco_Bkg",
                                                "MVA score of reco tracks ;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Bkg_EB_ = ibook.book1D("Muon_Iso_mva_track_reco_Bkg_EB",
                                                "MVA score of reco tracks - Background Barrel;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
  meMuonISO_mva_trk_reco_Bkg_EE_ = ibook.book1D("Muon_Iso_mva_track_reco_Bkg_EE",
                                                "MVA score of reco tracks - Background Endcap;MVA Score;Counts",
                                                100,
                                                0,
                                                1);
    // test to checking the type of tracks
  meMuonISO_trk_type_Sig_EB_ = ibook.book1D(
      "Muon_Iso_track_type_Sig_EB", "Check the type of tracks (from PV=0, from PU=1, fake=2) - Signal Barrel", 4, 0, 4);
  meMuonISO_trk_type_Sig_EE_ = ibook.book1D(
      "Muon_Iso_track_type_Sig_EE", "Check the type of tracks (from PV=0, from PU=1, fake=2) - Signal Endcap", 4, 0, 4);
  meMuonISO_trk_type_Bkg_EB_ = ibook.book1D(
      "Muon_Iso_track_type_Bkg_EB", "Check the type of tracks (from PV=0, from PU=1, fake=2) - Background Barrel", 4, 0, 4);
  meMuonISO_trk_type_Bkg_EE_ = ibook.book1D(
      "Muon_Iso_track_type_Bkg_EE", "Check the type of tracks (from PV=0, from PU=1, fake=2) - Background Endcap", 4, 0, 4);
  meMuonISO_trk_type_v2_Sig_EB_ = ibook.book1D(
      "Muon_Iso_track_type_v2_Sig_EB", "Check the type of tracks (from PV=0, from PU=1, fake=2, from SV=3) - Signal Barrel", 4, 0, 4);
  meMuonISO_trk_type_v2_Sig_EE_ = ibook.book1D(
      "Muon_Iso_track_type_v2_Sig_EE", "Check the type of tracks (from PV=0, from PU=1, fake=2, from SV=3) - Signal Endcap", 4, 0, 4);
  meMuonISO_trk_type_v2_Bkg_EB_ = ibook.book1D(
      "Muon_Iso_track_type_v2_Bkg_EB", "Check the type of tracks (from PV=0, from PU=1, fake=2, from SV=3) - Background Barrel", 4, 0, 4);
  meMuonISO_trk_type_v2_Bkg_EE_ = ibook.book1D(
      "Muon_Iso_track_type_v2_Bkg_EE", "Check the type of tracks (from PV=0, from PU=1, fake=2, from SV=3) - Background Endcap", 4, 0, 4);

  meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EB",
                                                "dtSig distribution for reco PV track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PVtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EE",
                                                "dtSig distribution for reco PV track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EB",
                                                "dtSig distribution for reco PV track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PVtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EE",
                                                "dtSig distribution for reco PV track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EB",
                                                "dtSig distribution for reco PU track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PUtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EE",
                                                "dtSig distribution for reco PU track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EB",
                                                "dtSig distribution for reco PU track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_reco_PUtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EE",
                                                "dtSig distribution for reco PU track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  // FIXME the name of histograms below need to be changed (it is not about signif)
  meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EB",
                                                "check whether it is a fake track or not - 1: fake track - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_faketrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EE",
                                                "check whether it is a fake track or not - 1: fake track - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EB",
                                                "check whether it is a fake track or not - 1: fake track - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_faketrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EE",
                                                "check whether it is a fake track or not - 1: fake track - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_without_tErr_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_reco_without_tErr_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);

  meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PVtrk_Sig_EB",
                                                "dtSig distribution for sim PV track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PVtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PVtrk_Sig_EE",
                                                "dtSig distribution for sim PV track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PVtrk_Bkg_EB",
                                                "dtSig distribution for sim PV track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PVtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PVtrk_Bkg_EE",
                                                "dtSig distribution for sim PV track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PUtrk_Sig_EB",
                                                "dtSig distribution for sim PU track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PUtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PUtrk_Sig_EE",
                                                "dtSig distribution for sim PU track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PUtrk_Bkg_EB",
                                                "dtSig distribution for sim PU track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_muon_trk_sim_PUtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_PUtrk_Bkg_EE",
                                                "dtSig distribution for sim PU track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  // FIXME the name of histograms below need to be changed (it is not about signif)
  meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_faketrk_Sig_EB",
                                                "check whether it is a fake track or not - 1: fake track - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_faketrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_faketrk_Sig_EE",
                                                "check whether it is a fake track or not - 1: fake track - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_faketrk_Bkg_EB",
                                                "check whether it is a fake track or not - 1: fake track - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_faketrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_faketrk_Bkg_EE",
                                                "check whether it is a fake track or not - 1: fake track - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_without_tErr_Sig_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_without_tErr_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_without_tErr_Sig_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_without_tErr_Bkg_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_muon_trk_sim_without_tErr_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_muon_track_sim_without_tErr_Bkg_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);

  meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EB",
                                                "dtSig distribution for reco PV track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PVtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EE",
                                                "dtSig distribution for reco PV track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EB",
                                                "dtSig distribution for reco PV track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PVtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EE",
                                                "dtSig distribution for reco PV track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EB",
                                                "dtSig distribution for reco PU track and reco muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PUtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EE",
                                                "dtSig distribution for reco PU track and reco muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EB",
                                                "dtSig distribution for reco PU track and reco muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_reco_PUtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EE",
                                                "dtSig distribution for reco PU track and reco muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  // FIXME the name of histograms below need to be changed (it is not about signif)
  meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EB",
                                                "check whether it is a fake track or not - 1: fake track - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_faketrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EE",
                                                "check whether it is a fake track or not - 1: fake track - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EB",
                                                "check whether it is a fake track or not - 1: fake track - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_faketrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EE",
                                                "check whether it is a fake track or not - 1: fake track - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_without_tErr_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_reco_without_tErr_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);

  meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PVtrk_Sig_EB",
                                                "dtSig distribution for sim PV track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PVtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PVtrk_Sig_EE",
                                                "dtSig distribution for sim PV track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PVtrk_Bkg_EB",
                                                "dtSig distribution for sim PV track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PVtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PVtrk_Bkg_EE",
                                                "dtSig distribution for sim PV track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PUtrk_Sig_EB",
                                                "dtSig distribution for sim PU track and sim muon - Signal Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PUtrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PUtrk_Sig_EE",
                                                "dtSig distribution for sim PU track and sim muon - Signal Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PUtrk_Bkg_EB",
                                                "dtSig distribution for sim PU track and sim muon - Background Barrel;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  meMuonISO_dtSig_PV_trk_sim_PUtrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_PUtrk_Bkg_EE",
                                                "dtSig distribution for sim PU track and sim muon - Background Endcap;#sigma;Counts",
                                                10,
                                                0,
                                                10);
  // FIXME the name of histograms below need to be changed (it is not about signif)
  meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_faketrk_Sig_EB",
                                                "check whether it is a fake track or not - 1: fake track - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_faketrk_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_faketrk_Sig_EE",
                                                "check whether it is a fake track or not - 1: fake track - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_faketrk_Bkg_EB",
                                                "check whether it is a fake track or not - 1: fake track - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_faketrk_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_faketrk_Bkg_EE",
                                                "check whether it is a fake track or not - 1: fake track - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_without_tErr_Sig_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_without_tErr_Sig_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_without_tErr_Sig_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Signal Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EB_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_without_tErr_Bkg_EB",
                                                "check on tracks not having error of time - 1: not having error of time - Background Barrel;#sigma;Counts",
                                                2,
                                                0,
                                                2);
  meMuonISO_dtSig_PV_trk_sim_without_tErr_Bkg_EE_ = ibook.book1D("Muon_Iso_dtSig_PV_track_sim_without_tErr_Bkg_EE",
                                                "check on tracks not having error of time - 1: not having error of time - Background Endcap;#sigma;Counts",
                                                2,
                                                0,
                                                2);






  meMuonISO_Nmuons_Sig_ = ibook.book1D("Muon_Iso_Nmuons_Sig",
                                                "Number of muons after basic "
                                                "cuts - Signal;Number of muons;Counts",
                                                5,
                                                0,
                                                5);
  meMuonISO_Nmuons_Sig_EB_ = ibook.book1D("Muon_Iso_Nmuons_Sig_EB",
                                                "Number of muons after basic "
                                                "cuts - Signal Barrel;Number of muons;Counts",
                                                5,
                                                0,
                                                5);
  meMuonISO_Nmuons_Sig_EE_ = ibook.book1D("Muon_Iso_Nmuons_Sig_EE",
                                                "Number of muons after basic "
                                                "cuts - Signal Endcap;Number of muons;Counts",
                                                5,
                                                0,
                                                5);
  meMuonISO_Nmuons_Bkg_ = ibook.book1D("Muon_Iso_Nmuons_Bkg",
                                                "Number of muons after basic "
                                                "cuts - Background;Number of muons;Counts",
                                                5,
                                                0,
                                                5);
  meMuonISO_Nmuons_Bkg_EB_ = ibook.book1D("Muon_Iso_Nmuons_Bkg_EB",
                                                "Number of muons after basic "
                                                "cuts - Background Barrel;Number of muons;Counts",
                                                5,
                                                0,
                                                5);
  meMuonISO_Nmuons_Bkg_EE_ = ibook.book1D("Muon_Iso_Nmuons_Bkg_EE",
                                                "Number of muons after basic "
                                                "cuts - Background Endcap;Number of muons;Counts",
                                                5,
                                                0,
                                                5);

  // defining vectors for more efficient hist filling
  // Promt part
  if (optionalPlots_) {
    Ntracks_EB_list_Sig = {meMuonISO_Ntracks_MTD_1_Sig_EB_,
                           meMuonISO_Ntracks_MTD_2_Sig_EB_,
                           meMuonISO_Ntracks_MTD_3_Sig_EB_,
                           meMuonISO_Ntracks_MTD_4_Sig_EB_,
                           meMuonISO_Ntracks_MTD_5_Sig_EB_,
                           meMuonISO_Ntracks_MTD_6_Sig_EB_,
                           meMuonISO_Ntracks_MTD_7_Sig_EB_};
    ch_iso_EB_list_Sig = {meMuonISO_chIso_MTD_1_Sig_EB_,
                          meMuonISO_chIso_MTD_2_Sig_EB_,
                          meMuonISO_chIso_MTD_3_Sig_EB_,
                          meMuonISO_chIso_MTD_4_Sig_EB_,
                          meMuonISO_chIso_MTD_5_Sig_EB_,
                          meMuonISO_chIso_MTD_6_Sig_EB_,
                          meMuonISO_chIso_MTD_7_Sig_EB_};
    rel_ch_iso_EB_list_Sig = {meMuonISO_rel_chIso_MTD_1_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_2_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_3_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_4_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_5_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_6_Sig_EB_,
                              meMuonISO_rel_chIso_MTD_7_Sig_EB_};
  }
  Ntracks_EB_list_Significance_Sig = {
      meMuonISO_Ntracks_MTD_4sigma_Sig_EB_, meMuonISO_Ntracks_MTD_3sigma_Sig_EB_, meMuonISO_Ntracks_MTD_2sigma_Sig_EB_};
  ch_iso_EB_list_Significance_Sig = {
      meMuonISO_chIso_MTD_4sigma_Sig_EB_, meMuonISO_chIso_MTD_3sigma_Sig_EB_, meMuonISO_chIso_MTD_2sigma_Sig_EB_};
  rel_ch_iso_EB_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_4sigma_Sig_EB_,
                                         meMuonISO_rel_chIso_MTD_3sigma_Sig_EB_,
                                         meMuonISO_rel_chIso_MTD_2sigma_Sig_EB_};

  if (optionalPlots_) {
    Ntracks_EE_list_Sig = {meMuonISO_Ntracks_MTD_1_Sig_EE_,
                           meMuonISO_Ntracks_MTD_2_Sig_EE_,
                           meMuonISO_Ntracks_MTD_3_Sig_EE_,
                           meMuonISO_Ntracks_MTD_4_Sig_EE_,
                           meMuonISO_Ntracks_MTD_5_Sig_EE_,
                           meMuonISO_Ntracks_MTD_6_Sig_EE_,
                           meMuonISO_Ntracks_MTD_7_Sig_EE_};
    ch_iso_EE_list_Sig = {meMuonISO_chIso_MTD_1_Sig_EE_,
                          meMuonISO_chIso_MTD_2_Sig_EE_,
                          meMuonISO_chIso_MTD_3_Sig_EE_,
                          meMuonISO_chIso_MTD_4_Sig_EE_,
                          meMuonISO_chIso_MTD_5_Sig_EE_,
                          meMuonISO_chIso_MTD_6_Sig_EE_,
                          meMuonISO_chIso_MTD_7_Sig_EE_};
    rel_ch_iso_EE_list_Sig = {meMuonISO_rel_chIso_MTD_1_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_2_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_3_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_4_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_5_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_6_Sig_EE_,
                              meMuonISO_rel_chIso_MTD_7_Sig_EE_};
  }
  Ntracks_EE_list_Significance_Sig = {
      meMuonISO_Ntracks_MTD_4sigma_Sig_EE_, meMuonISO_Ntracks_MTD_3sigma_Sig_EE_, meMuonISO_Ntracks_MTD_2sigma_Sig_EE_};
  ch_iso_EE_list_Significance_Sig = {
      meMuonISO_chIso_MTD_4sigma_Sig_EE_, meMuonISO_chIso_MTD_3sigma_Sig_EE_, meMuonISO_chIso_MTD_2sigma_Sig_EE_};
  rel_ch_iso_EE_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_4sigma_Sig_EE_,
                                         meMuonISO_rel_chIso_MTD_3sigma_Sig_EE_,
                                         meMuonISO_rel_chIso_MTD_2sigma_Sig_EE_};

  if (optionalPlots_) {
    Muon_pT_MTD_EB_list_Sig = {meMuon_pt_MTD_1_Sig_EB_,
                              meMuon_pt_MTD_2_Sig_EB_,
                              meMuon_pt_MTD_3_Sig_EB_,
                              meMuon_pt_MTD_4_Sig_EB_,
                              meMuon_pt_MTD_5_Sig_EB_,
                              meMuon_pt_MTD_6_Sig_EB_,
                              meMuon_pt_MTD_7_Sig_EB_};
    Muon_eta_MTD_EB_list_Sig = {meMuon_eta_MTD_1_Sig_EB_,
                               meMuon_eta_MTD_2_Sig_EB_,
                               meMuon_eta_MTD_3_Sig_EB_,
                               meMuon_eta_MTD_4_Sig_EB_,
                               meMuon_eta_MTD_5_Sig_EB_,
                               meMuon_eta_MTD_6_Sig_EB_,
                               meMuon_eta_MTD_7_Sig_EB_};
    Muon_phi_MTD_EB_list_Sig = {meMuon_phi_MTD_1_Sig_EB_,
                               meMuon_phi_MTD_2_Sig_EB_,
                               meMuon_phi_MTD_3_Sig_EB_,
                               meMuon_phi_MTD_4_Sig_EB_,
                               meMuon_phi_MTD_5_Sig_EB_,
                               meMuon_phi_MTD_6_Sig_EB_,
                               meMuon_phi_MTD_7_Sig_EB_};
  }

  Muon_pT_MTD_EB_list_Significance_Sig = {
      meMuon_pt_MTD_4sigma_Sig_EB_, meMuon_pt_MTD_3sigma_Sig_EB_, meMuon_pt_MTD_2sigma_Sig_EB_};
  Muon_eta_MTD_EB_list_Significance_Sig = {
      meMuon_eta_MTD_4sigma_Sig_EB_, meMuon_eta_MTD_3sigma_Sig_EB_, meMuon_eta_MTD_2sigma_Sig_EB_};
  Muon_phi_MTD_EB_list_Significance_Sig = {
      meMuon_phi_MTD_4sigma_Sig_EB_, meMuon_phi_MTD_3sigma_Sig_EB_, meMuon_phi_MTD_2sigma_Sig_EB_};

  if (optionalPlots_) {
    Muon_pT_MTD_EE_list_Sig = {meMuon_pt_MTD_1_Sig_EE_,
                              meMuon_pt_MTD_2_Sig_EE_,
                              meMuon_pt_MTD_3_Sig_EE_,
                              meMuon_pt_MTD_4_Sig_EE_,
                              meMuon_pt_MTD_5_Sig_EE_,
                              meMuon_pt_MTD_6_Sig_EE_,
                              meMuon_pt_MTD_7_Sig_EE_};
    Muon_eta_MTD_EE_list_Sig = {meMuon_eta_MTD_1_Sig_EE_,
                               meMuon_eta_MTD_2_Sig_EE_,
                               meMuon_eta_MTD_3_Sig_EE_,
                               meMuon_eta_MTD_4_Sig_EE_,
                               meMuon_eta_MTD_5_Sig_EE_,
                               meMuon_eta_MTD_6_Sig_EE_,
                               meMuon_eta_MTD_7_Sig_EE_};
    Muon_phi_MTD_EE_list_Sig = {meMuon_phi_MTD_1_Sig_EE_,
                               meMuon_phi_MTD_2_Sig_EE_,
                               meMuon_phi_MTD_3_Sig_EE_,
                               meMuon_phi_MTD_4_Sig_EE_,
                               meMuon_phi_MTD_5_Sig_EE_,
                               meMuon_phi_MTD_6_Sig_EE_,
                               meMuon_phi_MTD_7_Sig_EE_};
  }
  Muon_pT_MTD_EE_list_Significance_Sig = {
      meMuon_pt_MTD_4sigma_Sig_EE_, meMuon_pt_MTD_3sigma_Sig_EE_, meMuon_pt_MTD_2sigma_Sig_EE_};
  Muon_eta_MTD_EE_list_Significance_Sig = {
      meMuon_eta_MTD_4sigma_Sig_EE_, meMuon_eta_MTD_3sigma_Sig_EE_, meMuon_eta_MTD_2sigma_Sig_EE_};
  Muon_phi_MTD_EE_list_Significance_Sig = {
      meMuon_phi_MTD_4sigma_Sig_EE_, meMuon_phi_MTD_3sigma_Sig_EE_, meMuon_phi_MTD_2sigma_Sig_EE_};

  // For SIM CASE
  if (optionalPlots_) {
    Ntracks_sim_EB_list_Sig = {meMuonISO_Ntracks_MTD_sim_1_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_2_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_3_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_4_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_5_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_6_Sig_EB_,
                               meMuonISO_Ntracks_MTD_sim_7_Sig_EB_};
    ch_iso_sim_EB_list_Sig = {meMuonISO_chIso_MTD_sim_1_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_2_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_3_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_4_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_5_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_6_Sig_EB_,
                              meMuonISO_chIso_MTD_sim_7_Sig_EB_};
    rel_ch_iso_sim_EB_list_Sig = {meMuonISO_rel_chIso_MTD_sim_1_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_2_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_3_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_4_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_5_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_6_Sig_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_7_Sig_EB_};

    Ntracks_sim_EB_list_Significance_Sig = {meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EB_,
                                            meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EB_,
                                            meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EB_};
    ch_iso_sim_EB_list_Significance_Sig = {meMuonISO_chIso_MTD_sim_4sigma_Sig_EB_,
                                           meMuonISO_chIso_MTD_sim_3sigma_Sig_EB_,
                                           meMuonISO_chIso_MTD_sim_2sigma_Sig_EB_};
    rel_ch_iso_sim_EB_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EB_,
                                               meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EB_,
                                               meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EB_};

    Ntracks_sim_EE_list_Sig = {meMuonISO_Ntracks_MTD_sim_1_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_2_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_3_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_4_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_5_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_6_Sig_EE_,
                               meMuonISO_Ntracks_MTD_sim_7_Sig_EE_};
    ch_iso_sim_EE_list_Sig = {meMuonISO_chIso_MTD_sim_1_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_2_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_3_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_4_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_5_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_6_Sig_EE_,
                              meMuonISO_chIso_MTD_sim_7_Sig_EE_};
    rel_ch_iso_sim_EE_list_Sig = {meMuonISO_rel_chIso_MTD_sim_1_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_2_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_3_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_4_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_5_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_6_Sig_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_7_Sig_EE_};

    Ntracks_sim_EE_list_Significance_Sig = {meMuonISO_Ntracks_MTD_sim_4sigma_Sig_EE_,
                                            meMuonISO_Ntracks_MTD_sim_3sigma_Sig_EE_,
                                            meMuonISO_Ntracks_MTD_sim_2sigma_Sig_EE_};
    ch_iso_sim_EE_list_Significance_Sig = {meMuonISO_chIso_MTD_sim_4sigma_Sig_EE_,
                                           meMuonISO_chIso_MTD_sim_3sigma_Sig_EE_,
                                           meMuonISO_chIso_MTD_sim_2sigma_Sig_EE_};
    rel_ch_iso_sim_EE_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_sim_4sigma_Sig_EE_,
                                               meMuonISO_rel_chIso_MTD_sim_3sigma_Sig_EE_,
                                               meMuonISO_rel_chIso_MTD_sim_2sigma_Sig_EE_};

    Muon_pT_sim_MTD_EB_list_Sig = {meMuon_pt_sim_MTD_1_Sig_EB_,
                                  meMuon_pt_sim_MTD_2_Sig_EB_,
                                  meMuon_pt_sim_MTD_3_Sig_EB_,
                                  meMuon_pt_sim_MTD_4_Sig_EB_,
                                  meMuon_pt_sim_MTD_5_Sig_EB_,
                                  meMuon_pt_sim_MTD_6_Sig_EB_,
                                  meMuon_pt_sim_MTD_7_Sig_EB_};

    Muon_pT_sim_MTD_EB_list_Significance_Sig = {
        meMuon_pt_sim_MTD_4sigma_Sig_EB_, meMuon_pt_sim_MTD_3sigma_Sig_EB_, meMuon_pt_sim_MTD_2sigma_Sig_EB_};

    Muon_pT_sim_MTD_EE_list_Sig = {meMuon_pt_sim_MTD_1_Sig_EE_,
                                  meMuon_pt_sim_MTD_2_Sig_EE_,
                                  meMuon_pt_sim_MTD_3_Sig_EE_,
                                  meMuon_pt_sim_MTD_4_Sig_EE_,
                                  meMuon_pt_sim_MTD_5_Sig_EE_,
                                  meMuon_pt_sim_MTD_6_Sig_EE_,
                                  meMuon_pt_sim_MTD_7_Sig_EE_};
    Muon_pT_sim_MTD_EE_list_Significance_Sig = {
        meMuon_pt_sim_MTD_4sigma_Sig_EE_, meMuon_pt_sim_MTD_3sigma_Sig_EE_, meMuon_pt_sim_MTD_2sigma_Sig_EE_};
  }
  // test
  // For GEN CASE
    Muon_pT_gen_MTD_EB_list_Significance_Sig = {
        meMuon_pt_gen_MTD_4sigma_Sig_EB_, meMuon_pt_gen_MTD_3sigma_Sig_EB_, meMuon_pt_gen_MTD_2sigma_Sig_EB_};
    Muon_pT_gen_MTD_EE_list_Significance_Sig = {
        meMuon_pt_gen_MTD_4sigma_Sig_EE_, meMuon_pt_gen_MTD_3sigma_Sig_EE_, meMuon_pt_gen_MTD_2sigma_Sig_EE_};
    Muon_pT_gen_MTD_EB_list_Significance_Bkg = {
        meMuon_pt_gen_MTD_4sigma_Bkg_EB_, meMuon_pt_gen_MTD_3sigma_Bkg_EB_, meMuon_pt_gen_MTD_2sigma_Bkg_EB_};
    Muon_pT_gen_MTD_EE_list_Significance_Bkg = {
        meMuon_pt_gen_MTD_4sigma_Bkg_EE_, meMuon_pt_gen_MTD_3sigma_Bkg_EE_, meMuon_pt_gen_MTD_2sigma_Bkg_EE_};

    ch_iso_gen_EB_list_Significance_Sig = {meMuonISO_chIso_MTD_gen_4sigma_Sig_EB_,
                                           meMuonISO_chIso_MTD_gen_3sigma_Sig_EB_,
                                           meMuonISO_chIso_MTD_gen_2sigma_Sig_EB_};
    rel_ch_iso_gen_EB_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EB_,
                                               meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EB_,
                                               meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EB_};
    ch_iso_gen_EE_list_Significance_Sig = {meMuonISO_chIso_MTD_gen_4sigma_Sig_EE_,
                                           meMuonISO_chIso_MTD_gen_3sigma_Sig_EE_,
                                           meMuonISO_chIso_MTD_gen_2sigma_Sig_EE_};
    rel_ch_iso_gen_EE_list_Significance_Sig = {meMuonISO_rel_chIso_MTD_gen_4sigma_Sig_EE_,
                                               meMuonISO_rel_chIso_MTD_gen_3sigma_Sig_EE_,
                                               meMuonISO_rel_chIso_MTD_gen_2sigma_Sig_EE_};
    ch_iso_gen_EB_list_Significance_Bkg = {meMuonISO_chIso_MTD_gen_4sigma_Bkg_EB_,
                                           meMuonISO_chIso_MTD_gen_3sigma_Bkg_EB_,
                                           meMuonISO_chIso_MTD_gen_2sigma_Bkg_EB_};
    rel_ch_iso_gen_EB_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EB_,
                                               meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EB_,
                                               meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EB_};
    ch_iso_gen_EE_list_Significance_Bkg = {meMuonISO_chIso_MTD_gen_4sigma_Bkg_EE_,
                                           meMuonISO_chIso_MTD_gen_3sigma_Bkg_EE_,
                                           meMuonISO_chIso_MTD_gen_2sigma_Bkg_EE_};
    rel_ch_iso_gen_EE_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_gen_4sigma_Bkg_EE_,
                                               meMuonISO_rel_chIso_MTD_gen_3sigma_Bkg_EE_,
                                               meMuonISO_rel_chIso_MTD_gen_2sigma_Bkg_EE_};
    // test end

  // Non-promt part
  if (optionalPlots_) {
    Ntracks_EB_list_Bkg = {meMuonISO_Ntracks_MTD_1_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_2_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_3_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_4_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_5_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_6_Bkg_EB_,
                           meMuonISO_Ntracks_MTD_7_Bkg_EB_};
    ch_iso_EB_list_Bkg = {meMuonISO_chIso_MTD_1_Bkg_EB_,
                          meMuonISO_chIso_MTD_2_Bkg_EB_,
                          meMuonISO_chIso_MTD_3_Bkg_EB_,
                          meMuonISO_chIso_MTD_4_Bkg_EB_,
                          meMuonISO_chIso_MTD_5_Bkg_EB_,
                          meMuonISO_chIso_MTD_6_Bkg_EB_,
                          meMuonISO_chIso_MTD_7_Bkg_EB_};
    rel_ch_iso_EB_list_Bkg = {meMuonISO_rel_chIso_MTD_1_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_2_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_3_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_4_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_5_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_6_Bkg_EB_,
                              meMuonISO_rel_chIso_MTD_7_Bkg_EB_};
  }
  Ntracks_EB_list_Significance_Bkg = {
      meMuonISO_Ntracks_MTD_4sigma_Bkg_EB_, meMuonISO_Ntracks_MTD_3sigma_Bkg_EB_, meMuonISO_Ntracks_MTD_2sigma_Bkg_EB_};
  ch_iso_EB_list_Significance_Bkg = {
      meMuonISO_chIso_MTD_4sigma_Bkg_EB_, meMuonISO_chIso_MTD_3sigma_Bkg_EB_, meMuonISO_chIso_MTD_2sigma_Bkg_EB_};
  rel_ch_iso_EB_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_4sigma_Bkg_EB_,
                                         meMuonISO_rel_chIso_MTD_3sigma_Bkg_EB_,
                                         meMuonISO_rel_chIso_MTD_2sigma_Bkg_EB_};

  if (optionalPlots_) {
    Ntracks_EE_list_Bkg = {meMuonISO_Ntracks_MTD_1_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_2_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_3_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_4_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_5_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_6_Bkg_EE_,
                           meMuonISO_Ntracks_MTD_7_Bkg_EE_};
    ch_iso_EE_list_Bkg = {meMuonISO_chIso_MTD_1_Bkg_EE_,
                          meMuonISO_chIso_MTD_2_Bkg_EE_,
                          meMuonISO_chIso_MTD_3_Bkg_EE_,
                          meMuonISO_chIso_MTD_4_Bkg_EE_,
                          meMuonISO_chIso_MTD_5_Bkg_EE_,
                          meMuonISO_chIso_MTD_6_Bkg_EE_,
                          meMuonISO_chIso_MTD_7_Bkg_EE_};
    rel_ch_iso_EE_list_Bkg = {meMuonISO_rel_chIso_MTD_1_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_2_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_3_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_4_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_5_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_6_Bkg_EE_,
                              meMuonISO_rel_chIso_MTD_7_Bkg_EE_};
  }
  Ntracks_EE_list_Significance_Bkg = {
      meMuonISO_Ntracks_MTD_4sigma_Bkg_EE_, meMuonISO_Ntracks_MTD_3sigma_Bkg_EE_, meMuonISO_Ntracks_MTD_2sigma_Bkg_EE_};
  ch_iso_EE_list_Significance_Bkg = {
      meMuonISO_chIso_MTD_4sigma_Bkg_EE_, meMuonISO_chIso_MTD_3sigma_Bkg_EE_, meMuonISO_chIso_MTD_2sigma_Bkg_EE_};
  rel_ch_iso_EE_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_4sigma_Bkg_EE_,
                                         meMuonISO_rel_chIso_MTD_3sigma_Bkg_EE_,
                                         meMuonISO_rel_chIso_MTD_2sigma_Bkg_EE_};
  if (optionalPlots_) {
    Muon_pT_MTD_EB_list_Bkg = {meMuon_pt_MTD_1_Bkg_EB_,
                              meMuon_pt_MTD_2_Bkg_EB_,
                              meMuon_pt_MTD_3_Bkg_EB_,
                              meMuon_pt_MTD_4_Bkg_EB_,
                              meMuon_pt_MTD_5_Bkg_EB_,
                              meMuon_pt_MTD_6_Bkg_EB_,
                              meMuon_pt_MTD_7_Bkg_EB_};
    Muon_eta_MTD_EB_list_Bkg = {meMuon_eta_MTD_1_Bkg_EB_,
                               meMuon_eta_MTD_2_Bkg_EB_,
                               meMuon_eta_MTD_3_Bkg_EB_,
                               meMuon_eta_MTD_4_Bkg_EB_,
                               meMuon_eta_MTD_5_Bkg_EB_,
                               meMuon_eta_MTD_6_Bkg_EB_,
                               meMuon_eta_MTD_7_Bkg_EB_};
    Muon_phi_MTD_EB_list_Bkg = {meMuon_phi_MTD_1_Bkg_EB_,
                               meMuon_phi_MTD_2_Bkg_EB_,
                               meMuon_phi_MTD_3_Bkg_EB_,
                               meMuon_phi_MTD_4_Bkg_EB_,
                               meMuon_phi_MTD_5_Bkg_EB_,
                               meMuon_phi_MTD_6_Bkg_EB_,
                               meMuon_phi_MTD_7_Bkg_EB_};
  }
  Muon_pT_MTD_EB_list_Significance_Bkg = {
      meMuon_pt_MTD_4sigma_Bkg_EB_, meMuon_pt_MTD_3sigma_Bkg_EB_, meMuon_pt_MTD_2sigma_Bkg_EB_};
  Muon_eta_MTD_EB_list_Significance_Bkg = {
      meMuon_eta_MTD_4sigma_Bkg_EB_, meMuon_eta_MTD_3sigma_Bkg_EB_, meMuon_eta_MTD_2sigma_Bkg_EB_};
  Muon_phi_MTD_EB_list_Significance_Bkg = {
      meMuon_phi_MTD_4sigma_Bkg_EB_, meMuon_phi_MTD_3sigma_Bkg_EB_, meMuon_phi_MTD_2sigma_Bkg_EB_};

  if (optionalPlots_) {
    Muon_pT_MTD_EE_list_Bkg = {meMuon_pt_MTD_1_Bkg_EE_,
                              meMuon_pt_MTD_2_Bkg_EE_,
                              meMuon_pt_MTD_3_Bkg_EE_,
                              meMuon_pt_MTD_4_Bkg_EE_,
                              meMuon_pt_MTD_5_Bkg_EE_,
                              meMuon_pt_MTD_6_Bkg_EE_,
                              meMuon_pt_MTD_7_Bkg_EE_};
    Muon_eta_MTD_EE_list_Bkg = {meMuon_eta_MTD_1_Bkg_EE_,
                               meMuon_eta_MTD_2_Bkg_EE_,
                               meMuon_eta_MTD_3_Bkg_EE_,
                               meMuon_eta_MTD_4_Bkg_EE_,
                               meMuon_eta_MTD_5_Bkg_EE_,
                               meMuon_eta_MTD_6_Bkg_EE_,
                               meMuon_eta_MTD_7_Bkg_EE_};
    Muon_phi_MTD_EE_list_Bkg = {meMuon_phi_MTD_1_Bkg_EE_,
                               meMuon_phi_MTD_2_Bkg_EE_,
                               meMuon_phi_MTD_3_Bkg_EE_,
                               meMuon_phi_MTD_4_Bkg_EE_,
                               meMuon_phi_MTD_5_Bkg_EE_,
                               meMuon_phi_MTD_6_Bkg_EE_,
                               meMuon_phi_MTD_7_Bkg_EE_};
  }
  Muon_pT_MTD_EE_list_Significance_Bkg = {
      meMuon_pt_MTD_4sigma_Bkg_EE_, meMuon_pt_MTD_3sigma_Bkg_EE_, meMuon_pt_MTD_2sigma_Bkg_EE_};
  Muon_eta_MTD_EE_list_Significance_Bkg = {
      meMuon_eta_MTD_4sigma_Bkg_EE_, meMuon_eta_MTD_3sigma_Bkg_EE_, meMuon_eta_MTD_2sigma_Bkg_EE_};
  Muon_phi_MTD_EE_list_Significance_Bkg = {
      meMuon_phi_MTD_4sigma_Bkg_EE_, meMuon_phi_MTD_3sigma_Bkg_EE_, meMuon_phi_MTD_2sigma_Bkg_EE_};

  // SIM CASE
  if (optionalPlots_) {
    Ntracks_sim_EB_list_Bkg = {meMuonISO_Ntracks_MTD_sim_1_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_2_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_3_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_4_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_5_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_6_Bkg_EB_,
                               meMuonISO_Ntracks_MTD_sim_7_Bkg_EB_};
    ch_iso_sim_EB_list_Bkg = {meMuonISO_chIso_MTD_sim_1_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_2_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_3_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_4_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_5_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_6_Bkg_EB_,
                              meMuonISO_chIso_MTD_sim_7_Bkg_EB_};
    rel_ch_iso_sim_EB_list_Bkg = {meMuonISO_rel_chIso_MTD_sim_1_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_2_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_3_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_4_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_5_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_6_Bkg_EB_,
                                  meMuonISO_rel_chIso_MTD_sim_7_Bkg_EB_};
    Ntracks_sim_EB_list_Significance_Bkg = {meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EB_,
                                            meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EB_,
                                            meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EB_};
    ch_iso_sim_EB_list_Significance_Bkg = {meMuonISO_chIso_MTD_sim_4sigma_Bkg_EB_,
                                           meMuonISO_chIso_MTD_sim_3sigma_Bkg_EB_,
                                           meMuonISO_chIso_MTD_sim_2sigma_Bkg_EB_};
    rel_ch_iso_sim_EB_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EB_,
                                               meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EB_,
                                               meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EB_};

    Ntracks_sim_EE_list_Bkg = {meMuonISO_Ntracks_MTD_sim_1_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_2_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_3_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_4_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_5_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_6_Bkg_EE_,
                               meMuonISO_Ntracks_MTD_sim_7_Bkg_EE_};
    ch_iso_sim_EE_list_Bkg = {meMuonISO_chIso_MTD_sim_1_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_2_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_3_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_4_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_5_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_6_Bkg_EE_,
                              meMuonISO_chIso_MTD_sim_7_Bkg_EE_};
    rel_ch_iso_sim_EE_list_Bkg = {meMuonISO_rel_chIso_MTD_sim_1_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_2_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_3_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_4_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_5_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_6_Bkg_EE_,
                                  meMuonISO_rel_chIso_MTD_sim_7_Bkg_EE_};

    Ntracks_sim_EE_list_Significance_Bkg = {meMuonISO_Ntracks_MTD_sim_4sigma_Bkg_EE_,
                                            meMuonISO_Ntracks_MTD_sim_3sigma_Bkg_EE_,
                                            meMuonISO_Ntracks_MTD_sim_2sigma_Bkg_EE_};
    ch_iso_sim_EE_list_Significance_Bkg = {meMuonISO_chIso_MTD_sim_4sigma_Bkg_EE_,
                                           meMuonISO_chIso_MTD_sim_3sigma_Bkg_EE_,
                                           meMuonISO_chIso_MTD_sim_2sigma_Bkg_EE_};
    rel_ch_iso_sim_EE_list_Significance_Bkg = {meMuonISO_rel_chIso_MTD_sim_4sigma_Bkg_EE_,
                                               meMuonISO_rel_chIso_MTD_sim_3sigma_Bkg_EE_,
                                               meMuonISO_rel_chIso_MTD_sim_2sigma_Bkg_EE_};

    Muon_pT_sim_MTD_EB_list_Bkg = {meMuon_pt_sim_MTD_1_Bkg_EB_,
                                  meMuon_pt_sim_MTD_2_Bkg_EB_,
                                  meMuon_pt_sim_MTD_3_Bkg_EB_,
                                  meMuon_pt_sim_MTD_4_Bkg_EB_,
                                  meMuon_pt_sim_MTD_5_Bkg_EB_,
                                  meMuon_pt_sim_MTD_6_Bkg_EB_,
                                  meMuon_pt_sim_MTD_7_Bkg_EB_};

    Muon_pT_sim_MTD_EB_list_Significance_Bkg = {
        meMuon_pt_sim_MTD_4sigma_Bkg_EB_, meMuon_pt_sim_MTD_3sigma_Bkg_EB_, meMuon_pt_sim_MTD_2sigma_Bkg_EB_};

    Muon_pT_sim_MTD_EE_list_Bkg = {meMuon_pt_sim_MTD_1_Bkg_EE_,
                                  meMuon_pt_sim_MTD_2_Bkg_EE_,
                                  meMuon_pt_sim_MTD_3_Bkg_EE_,
                                  meMuon_pt_sim_MTD_4_Bkg_EE_,
                                  meMuon_pt_sim_MTD_5_Bkg_EE_,
                                  meMuon_pt_sim_MTD_6_Bkg_EE_,
                                  meMuon_pt_sim_MTD_7_Bkg_EE_};

    Muon_pT_sim_MTD_EE_list_Significance_Bkg = {
        meMuon_pt_sim_MTD_4sigma_Bkg_EE_, meMuon_pt_sim_MTD_3sigma_Bkg_EE_, meMuon_pt_sim_MTD_2sigma_Bkg_EE_};
  }
  // dt distribution hist vecotrs

  general_pT_list = {meMuon_dt_general_pT_1,
                     meMuon_dt_general_pT_2,
                     meMuon_dt_general_pT_3,
                     meMuon_dt_general_pT_4,
                     meMuon_dt_general_pT_5,
                     meMuon_dt_general_pT_6,
                     meMuon_dt_general_pT_7,
                     meMuon_dt_general_pT_8,
                     meMuon_dt_general_pT_9};

  general_pT_Signif_list = {meMuon_dtSignif_general_pT_1,
                            meMuon_dtSignif_general_pT_2,
                            meMuon_dtSignif_general_pT_3,
                            meMuon_dtSignif_general_pT_4,
                            meMuon_dtSignif_general_pT_5,
                            meMuon_dtSignif_general_pT_6,
                            meMuon_dtSignif_general_pT_7,
                            meMuon_dtSignif_general_pT_8,
                            meMuon_dtSignif_general_pT_9};
  general_eta_list = {meMuon_dt_general_eta_1,
                      meMuon_dt_general_eta_2,
                      meMuon_dt_general_eta_3,
                      meMuon_dt_general_eta_4,
                      meMuon_dt_general_eta_5,
                      meMuon_dt_general_eta_6,
                      meMuon_dt_general_eta_7};

  general_eta_Signif_list = {meMuon_dtSignif_general_eta_1,
                             meMuon_dtSignif_general_eta_2,
                             meMuon_dtSignif_general_eta_3,
                             meMuon_dtSignif_general_eta_4,
                             meMuon_dtSignif_general_eta_5,
                             meMuon_dtSignif_general_eta_6,
                             meMuon_dtSignif_general_eta_7};
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdMuonIsoValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/MuonIso");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTag_vtx", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputMuon", edm::InputTag("muons"));
  desc.add<edm::InputTag>("globalMuonTrk", edm::InputTag("globalMuons"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("trackMinimumPt", 1.0);  // [GeV]
  desc.add<double>("trackMinimumEta", 1.5);
  desc.add<double>("trackMaximumEta", 3.2);
  desc.add<double>("rel_iso_cut", 0.08);
  //desc.add<bool>("optionTrackMatchToPV", true);
  desc.add<bool>("optionTrackMatchToPV", false);
  //desc.add<bool>("option_dtToTrack", false);  // default is dt with track, if false will do dt to vertex
  desc.add<bool>("option_dtToTrack", true);  // default is dt with track, if false will do dt to vertex
  desc.add<bool>("option_plots", true);
  desc.add<double>("min_dR_cut", 0.01);
  desc.add<double>("max_dR_cut", 0.3);
  desc.add<double>("min_pt_cut_EB", 0.7);
  desc.add<double>("min_pt_cut_EE", 0.4);
  desc.add<double>("max_dz_cut_EB", 0.1);  // PARAM
  desc.add<double>("max_dz_cut_EE", 0.2);  // PARAM
  desc.add<double>("max_dz_vtx_cut", 0.5);
  desc.add<double>("max_dxy_vtx_cut", 0.2);
//  desc.add<double>("min_strip_cut", 0.01);
  desc.add<double>("min_track_mtd_mva_cut", 0.5);

  // test
  desc.add<edm::InputTag>("offlineBS", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  // test end

  descriptions.add("mtdMuonIsoValid", desc);
}

// test
void MtdMuonIsoValidation::printSimVtxRecoVtxInfo(
    const struct MtdMuonIsoValidation::simPrimaryVertex& simpVtx,
    const struct MtdMuonIsoValidation::recoPrimaryVertex& recopVtx) {
  edm::LogPrint("MtdMuonIsoValidation")
      << "Sim vtx (x,y,z,t) = (" << std::setprecision(4) << simpVtx.x << "," << std::setprecision(4) << simpVtx.y << ","
      << std::setprecision(4) << simpVtx.z << "," << std::setprecision(4) << simpVtx.t * 1e9 << ")";
  edm::LogPrint("MtdMuonIsoValidation")
      << "Sim vtx: pt = " << std::setprecision(4) << simpVtx.pt << " ptsq = " << std::setprecision(6) << simpVtx.ptsq
      << " nGenTrk = " << simpVtx.nGenTrk << " nmatch recotrks = " << simpVtx.num_matched_reco_tracks;
  edm::LogPrint("MtdMuonIsoValidation")
      << "Reco vtx (x,y,z) = (" << std::setprecision(4) << recopVtx.x << "," << std::setprecision(4) << recopVtx.y
      << "," << std::setprecision(4) << recopVtx.z << ")";
  edm::LogPrint("MtdMuonIsoValidation")
      << "Reco vtx: pt = " << std::setprecision(4) << recopVtx.pt << " ptsq = " << std::setprecision(6) << recopVtx.ptsq
      << " nrecotrks = " << recopVtx.nRecoTrk << " nmatch simtrks = " << recopVtx.num_matched_sim_tracks;
  edm::LogPrint("MtdMuonIsoValidation") << "wnt " << recopVtx.sumwnt << " wos = " << recopVtx.sumwos;
  for (auto iTP = simpVtx.sim_vertex->daughterTracks_begin(); iTP != simpVtx.sim_vertex->daughterTracks_end(); ++iTP) {
    if ((**iTP).charge() == 0) {
      continue;
    }
    edm::LogPrint("MtdMuonIsoValidation")
        << "Daughter track of sim vertex: pt =" << std::setw(6) << std::setprecision(2) << (*iTP)->pt()
        << "  eta =" << std::setw(6) << std::setprecision(2) << (*iTP)->eta();
  }
}
// test end

DEFINE_FWK_MODULE(MtdMuonIsoValidation);

//*/
