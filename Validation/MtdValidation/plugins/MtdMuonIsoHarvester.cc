#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

class MtdMuonIsoHarvester : public DQMEDHarvester {
public:
  explicit MtdMuonIsoHarvester(const edm::ParameterSet& iConfig);
  ~MtdMuonIsoHarvester() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result);

  const std::string folder_;
  const bool optionalPlots_;

  // --- Histograms
  MonitorElement* mePtEffnoMTD_Sig_EB_;
  MonitorElement* meEtaEffnoMTD_Sig_EB_;
  MonitorElement* mePhiEffnoMTD_Sig_EB_;
  MonitorElement* mePtEffMTD_1_Sig_EB_;
  MonitorElement* meEtaEffMTD_1_Sig_EB_;
  MonitorElement* mePhiEffMTD_1_Sig_EB_;

  MonitorElement* mePtEffMTD_2_Sig_EB_;
  MonitorElement* meEtaEffMTD_2_Sig_EB_;
  MonitorElement* mePhiEffMTD_2_Sig_EB_;

  MonitorElement* mePtEffMTD_3_Sig_EB_;
  MonitorElement* meEtaEffMTD_3_Sig_EB_;
  MonitorElement* mePhiEffMTD_3_Sig_EB_;

  MonitorElement* mePtEffMTD_4_Sig_EB_;
  MonitorElement* meEtaEffMTD_4_Sig_EB_;
  MonitorElement* mePhiEffMTD_4_Sig_EB_;

  MonitorElement* mePtEffMTD_5_Sig_EB_;
  MonitorElement* meEtaEffMTD_5_Sig_EB_;
  MonitorElement* mePhiEffMTD_5_Sig_EB_;

  MonitorElement* mePtEffMTD_6_Sig_EB_;
  MonitorElement* meEtaEffMTD_6_Sig_EB_;
  MonitorElement* mePhiEffMTD_6_Sig_EB_;

  MonitorElement* mePtEffMTD_7_Sig_EB_;
  MonitorElement* meEtaEffMTD_7_Sig_EB_;
  MonitorElement* mePhiEffMTD_7_Sig_EB_;

  MonitorElement* mePtEffMTD_sim_1_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_2_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_3_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_4_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_5_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_6_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_7_Sig_EB_;

  MonitorElement* mePtEffgen_Sig_EB_;

  MonitorElement* mePtEffMTD_sim_4sigma_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_3sigma_Sig_EB_;
  MonitorElement* mePtEffMTD_sim_2sigma_Sig_EB_;

  MonitorElement* mePtEffMTD_4sigma_Sig_EB_;
  MonitorElement* meEtaEffMTD_4sigma_Sig_EB_;
  MonitorElement* mePhiEffMTD_4sigma_Sig_EB_;

  MonitorElement* mePtEffMTD_3sigma_Sig_EB_;
  MonitorElement* meEtaEffMTD_3sigma_Sig_EB_;
  MonitorElement* mePhiEffMTD_3sigma_Sig_EB_;

  MonitorElement* mePtEffMTD_2sigma_Sig_EB_;
  MonitorElement* meEtaEffMTD_2sigma_Sig_EB_;
  MonitorElement* mePhiEffMTD_2sigma_Sig_EB_;

  MonitorElement* mePtEffnoMTD_Sig_EE_;
  MonitorElement* meEtaEffnoMTD_Sig_EE_;
  MonitorElement* mePhiEffnoMTD_Sig_EE_;

  MonitorElement* mePtEffMTD_1_Sig_EE_;
  MonitorElement* meEtaEffMTD_1_Sig_EE_;
  MonitorElement* mePhiEffMTD_1_Sig_EE_;

  MonitorElement* mePtEffMTD_2_Sig_EE_;
  MonitorElement* meEtaEffMTD_2_Sig_EE_;
  MonitorElement* mePhiEffMTD_2_Sig_EE_;

  MonitorElement* mePtEffMTD_3_Sig_EE_;
  MonitorElement* meEtaEffMTD_3_Sig_EE_;
  MonitorElement* mePhiEffMTD_3_Sig_EE_;

  MonitorElement* mePtEffMTD_4_Sig_EE_;
  MonitorElement* meEtaEffMTD_4_Sig_EE_;
  MonitorElement* mePhiEffMTD_4_Sig_EE_;

  MonitorElement* mePtEffMTD_5_Sig_EE_;
  MonitorElement* meEtaEffMTD_5_Sig_EE_;
  MonitorElement* mePhiEffMTD_5_Sig_EE_;

  MonitorElement* mePtEffMTD_6_Sig_EE_;
  MonitorElement* meEtaEffMTD_6_Sig_EE_;
  MonitorElement* mePhiEffMTD_6_Sig_EE_;

  MonitorElement* mePtEffMTD_7_Sig_EE_;
  MonitorElement* meEtaEffMTD_7_Sig_EE_;
  MonitorElement* mePhiEffMTD_7_Sig_EE_;

  MonitorElement* mePtEffMTD_sim_1_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_2_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_3_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_4_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_5_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_6_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_7_Sig_EE_;

  MonitorElement* mePtEffgen_Sig_EE_;

  MonitorElement* mePtEffMTD_sim_4sigma_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_3sigma_Sig_EE_;
  MonitorElement* mePtEffMTD_sim_2sigma_Sig_EE_;

  MonitorElement* mePtEffMTD_4sigma_Sig_EE_;
  MonitorElement* meEtaEffMTD_4sigma_Sig_EE_;
  MonitorElement* mePhiEffMTD_4sigma_Sig_EE_;

  MonitorElement* mePtEffMTD_3sigma_Sig_EE_;
  MonitorElement* meEtaEffMTD_3sigma_Sig_EE_;
  MonitorElement* mePhiEffMTD_3sigma_Sig_EE_;

  MonitorElement* mePtEffMTD_2sigma_Sig_EE_;
  MonitorElement* meEtaEffMTD_2sigma_Sig_EE_;
  MonitorElement* mePhiEffMTD_2sigma_Sig_EE_;

  MonitorElement* mePtEffnoMTD_Bkg_EB_;
  MonitorElement* meEtaEffnoMTD_Bkg_EB_;
  MonitorElement* mePhiEffnoMTD_Bkg_EB_;
  MonitorElement* mePtEffMTD_1_Bkg_EB_;
  MonitorElement* meEtaEffMTD_1_Bkg_EB_;
  MonitorElement* mePhiEffMTD_1_Bkg_EB_;

  MonitorElement* mePtEffMTD_2_Bkg_EB_;
  MonitorElement* meEtaEffMTD_2_Bkg_EB_;
  MonitorElement* mePhiEffMTD_2_Bkg_EB_;

  MonitorElement* mePtEffMTD_3_Bkg_EB_;
  MonitorElement* meEtaEffMTD_3_Bkg_EB_;
  MonitorElement* mePhiEffMTD_3_Bkg_EB_;

  MonitorElement* mePtEffMTD_4_Bkg_EB_;
  MonitorElement* meEtaEffMTD_4_Bkg_EB_;
  MonitorElement* mePhiEffMTD_4_Bkg_EB_;

  MonitorElement* mePtEffMTD_5_Bkg_EB_;
  MonitorElement* meEtaEffMTD_5_Bkg_EB_;
  MonitorElement* mePhiEffMTD_5_Bkg_EB_;

  MonitorElement* mePtEffMTD_6_Bkg_EB_;
  MonitorElement* meEtaEffMTD_6_Bkg_EB_;
  MonitorElement* mePhiEffMTD_6_Bkg_EB_;

  MonitorElement* mePtEffMTD_7_Bkg_EB_;
  MonitorElement* meEtaEffMTD_7_Bkg_EB_;
  MonitorElement* mePhiEffMTD_7_Bkg_EB_;

  MonitorElement* mePtEffMTD_sim_1_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_2_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_3_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_4_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_5_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_6_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_7_Bkg_EB_;

  MonitorElement* mePtEffgen_Bkg_EB_;

  MonitorElement* mePtEffMTD_sim_4sigma_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_3sigma_Bkg_EB_;
  MonitorElement* mePtEffMTD_sim_2sigma_Bkg_EB_;

  MonitorElement* mePtEffMTD_4sigma_Bkg_EB_;
  MonitorElement* meEtaEffMTD_4sigma_Bkg_EB_;
  MonitorElement* mePhiEffMTD_4sigma_Bkg_EB_;

  MonitorElement* mePtEffMTD_3sigma_Bkg_EB_;
  MonitorElement* meEtaEffMTD_3sigma_Bkg_EB_;
  MonitorElement* mePhiEffMTD_3sigma_Bkg_EB_;

  MonitorElement* mePtEffMTD_2sigma_Bkg_EB_;
  MonitorElement* meEtaEffMTD_2sigma_Bkg_EB_;
  MonitorElement* mePhiEffMTD_2sigma_Bkg_EB_;

  MonitorElement* mePtEffnoMTD_Bkg_EE_;
  MonitorElement* meEtaEffnoMTD_Bkg_EE_;
  MonitorElement* mePhiEffnoMTD_Bkg_EE_;
  MonitorElement* mePtEffMTD_1_Bkg_EE_;
  MonitorElement* meEtaEffMTD_1_Bkg_EE_;
  MonitorElement* mePhiEffMTD_1_Bkg_EE_;

  MonitorElement* mePtEffMTD_2_Bkg_EE_;
  MonitorElement* meEtaEffMTD_2_Bkg_EE_;
  MonitorElement* mePhiEffMTD_2_Bkg_EE_;

  MonitorElement* mePtEffMTD_3_Bkg_EE_;
  MonitorElement* meEtaEffMTD_3_Bkg_EE_;
  MonitorElement* mePhiEffMTD_3_Bkg_EE_;

  MonitorElement* mePtEffMTD_4_Bkg_EE_;
  MonitorElement* meEtaEffMTD_4_Bkg_EE_;
  MonitorElement* mePhiEffMTD_4_Bkg_EE_;

  MonitorElement* mePtEffMTD_5_Bkg_EE_;
  MonitorElement* meEtaEffMTD_5_Bkg_EE_;
  MonitorElement* mePhiEffMTD_5_Bkg_EE_;

  MonitorElement* mePtEffMTD_6_Bkg_EE_;
  MonitorElement* meEtaEffMTD_6_Bkg_EE_;
  MonitorElement* mePhiEffMTD_6_Bkg_EE_;

  MonitorElement* mePtEffMTD_7_Bkg_EE_;
  MonitorElement* meEtaEffMTD_7_Bkg_EE_;
  MonitorElement* mePhiEffMTD_7_Bkg_EE_;

  MonitorElement* mePtEffMTD_sim_1_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_2_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_3_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_4_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_5_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_6_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_7_Bkg_EE_;

  MonitorElement* mePtEffgen_Bkg_EE_;

  MonitorElement* mePtEffMTD_sim_4sigma_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_3sigma_Bkg_EE_;
  MonitorElement* mePtEffMTD_sim_2sigma_Bkg_EE_;

  MonitorElement* mePtEffMTD_4sigma_Bkg_EE_;
  MonitorElement* meEtaEffMTD_4sigma_Bkg_EE_;
  MonitorElement* mePhiEffMTD_4sigma_Bkg_EE_;

  MonitorElement* mePtEffMTD_3sigma_Bkg_EE_;
  MonitorElement* meEtaEffMTD_3sigma_Bkg_EE_;
  MonitorElement* mePhiEffMTD_3sigma_Bkg_EE_;

  MonitorElement* mePtEffMTD_2sigma_Bkg_EE_;
  MonitorElement* meEtaEffMTD_2sigma_Bkg_EE_;
  MonitorElement* mePhiEffMTD_2sigma_Bkg_EE_;
};

// ------------ constructor and destructor --------------
MtdMuonIsoHarvester::MtdMuonIsoHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      optionalPlots_(iConfig.getParameter<bool>("option_plots")) {}

MtdMuonIsoHarvester::~MtdMuonIsoHarvester() {}

// auxiliary method to compute efficiency from the ratio of two 1D MonitorElement
void MtdMuonIsoHarvester::computeEfficiency1D(MonitorElement* num, MonitorElement* den, MonitorElement* result) {
  for (int ibin = 1; ibin <= den->getNbinsX(); ibin++) {
    double eff = num->getBinContent(ibin) / den->getBinContent(ibin);
    double bin_err = sqrt((num->getBinContent(ibin) * (den->getBinContent(ibin) - num->getBinContent(ibin))) /
                          pow(den->getBinContent(ibin), 3));
    if (den->getBinContent(ibin) == 0) {
      eff = 0;
      bin_err = 0;
    }
    result->setBinContent(ibin, eff);
    result->setBinError(ibin, bin_err);
  }
}

// ------------ endjob tasks ----------------------------
void MtdMuonIsoHarvester::dqmEndJob(DQMStore::IBooker& ibook, DQMStore::IGetter& igetter) {
  // --- Get the monitoring histograms
  // For promt (signal)
  MonitorElement* meMuon_pt_tot_Sig_EB_ = igetter.get(folder_ + "Muon_pT_tot_Sig_EB");
  MonitorElement* meMuon_pt_sim_tot_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_tot_Sig_EB");
  MonitorElement* meMuon_pt_noMTD_Sig_EB_ = igetter.get(folder_ + "Muon_pT_noMTD_Sig_EB");
  MonitorElement* meMuon_pt_gen_Sig_EB_ = igetter.get(folder_ + "Muon_pT_gen_Sig_EB");

  MonitorElement* meMuon_eta_tot_Sig_EB_ = igetter.get(folder_ + "Muon_eta_tot_Sig_EB");
  MonitorElement* meMuon_eta_noMTD_Sig_EB_ = igetter.get(folder_ + "Muon_eta_noMTD_Sig_EB");

  MonitorElement* meMuon_phi_tot_Sig_EB_ = igetter.get(folder_ + "Muon_phi_tot_Sig_EB");
  MonitorElement* meMuon_phi_noMTD_Sig_EB_ = igetter.get(folder_ + "Muon_phi_noMTD_Sig_EB");
  MonitorElement* meMuon_pt_MTD_1_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_1_Sig_EB");
  MonitorElement* meMuon_eta_MTD_1_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_1_Sig_EB");
  MonitorElement* meMuon_phi_MTD_1_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_1_Sig_EB");

  MonitorElement* meMuon_pt_MTD_2_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_2_Sig_EB");
  MonitorElement* meMuon_eta_MTD_2_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_2_Sig_EB");
  MonitorElement* meMuon_phi_MTD_2_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_2_Sig_EB");

  MonitorElement* meMuon_pt_MTD_3_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_3_Sig_EB");
  MonitorElement* meMuon_eta_MTD_3_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_3_Sig_EB");
  MonitorElement* meMuon_phi_MTD_3_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_3_Sig_EB");

  MonitorElement* meMuon_pt_MTD_4_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_4_Sig_EB");
  MonitorElement* meMuon_eta_MTD_4_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_4_Sig_EB");
  MonitorElement* meMuon_phi_MTD_4_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_4_Sig_EB");

  MonitorElement* meMuon_pt_MTD_5_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_5_Sig_EB");
  MonitorElement* meMuon_eta_MTD_5_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_5_Sig_EB");
  MonitorElement* meMuon_phi_MTD_5_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_5_Sig_EB");

  MonitorElement* meMuon_pt_MTD_6_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_6_Sig_EB");
  MonitorElement* meMuon_eta_MTD_6_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_6_Sig_EB");
  MonitorElement* meMuon_phi_MTD_6_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_6_Sig_EB");

  MonitorElement* meMuon_pt_MTD_7_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_7_Sig_EB");
  MonitorElement* meMuon_eta_MTD_7_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_7_Sig_EB");
  MonitorElement* meMuon_phi_MTD_7_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_7_Sig_EB");

  MonitorElement* meMuon_pt_sim_MTD_1_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_1_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_2_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_3_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_4_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_5_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_5_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_6_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_6_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_7_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_7_Sig_EB");

  MonitorElement* meMuon_pt_sim_MTD_2sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2sigma_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3sigma_Sig_EB");
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4sigma_Sig_EB");

  MonitorElement* meMuon_pt_MTD_2sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_2sigma_Sig_EB");
  MonitorElement* meMuon_eta_MTD_2sigma_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_2sigma_Sig_EB");
  MonitorElement* meMuon_phi_MTD_2sigma_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_2sigma_Sig_EB");

  MonitorElement* meMuon_pt_MTD_3sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_3sigma_Sig_EB");
  MonitorElement* meMuon_eta_MTD_3sigma_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_3sigma_Sig_EB");
  MonitorElement* meMuon_phi_MTD_3sigma_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_3sigma_Sig_EB");

  MonitorElement* meMuon_pt_MTD_4sigma_Sig_EB_ = igetter.get(folder_ + "Muon_pT_MTD_4sigma_Sig_EB");
  MonitorElement* meMuon_eta_MTD_4sigma_Sig_EB_ = igetter.get(folder_ + "Muon_eta_MTD_4sigma_Sig_EB");
  MonitorElement* meMuon_phi_MTD_4sigma_Sig_EB_ = igetter.get(folder_ + "Muon_phi_MTD_4sigma_Sig_EB");

  MonitorElement* meMuon_pt_tot_Sig_EE_ = igetter.get(folder_ + "Muon_pT_tot_Sig_EE");
  MonitorElement* meMuon_pt_sim_tot_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_tot_Sig_EE");
  MonitorElement* meMuon_pt_noMTD_Sig_EE_ = igetter.get(folder_ + "Muon_pT_noMTD_Sig_EE");
  MonitorElement* meMuon_pt_gen_Sig_EE_ = igetter.get(folder_ + "Muon_pT_gen_Sig_EE");

  MonitorElement* meMuon_eta_tot_Sig_EE_ = igetter.get(folder_ + "Muon_eta_tot_Sig_EE");
  MonitorElement* meMuon_eta_noMTD_Sig_EE_ = igetter.get(folder_ + "Muon_eta_noMTD_Sig_EE");

  MonitorElement* meMuon_phi_tot_Sig_EE_ = igetter.get(folder_ + "Muon_phi_tot_Sig_EE");
  MonitorElement* meMuon_phi_noMTD_Sig_EE_ = igetter.get(folder_ + "Muon_phi_noMTD_Sig_EE");
  MonitorElement* meMuon_pt_MTD_1_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_1_Sig_EE");
  MonitorElement* meMuon_eta_MTD_1_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_1_Sig_EE");
  MonitorElement* meMuon_phi_MTD_1_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_1_Sig_EE");

  MonitorElement* meMuon_pt_MTD_2_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_2_Sig_EE");
  MonitorElement* meMuon_eta_MTD_2_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_2_Sig_EE");
  MonitorElement* meMuon_phi_MTD_2_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_2_Sig_EE");

  MonitorElement* meMuon_pt_MTD_3_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_3_Sig_EE");
  MonitorElement* meMuon_eta_MTD_3_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_3_Sig_EE");
  MonitorElement* meMuon_phi_MTD_3_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_3_Sig_EE");

  MonitorElement* meMuon_pt_MTD_4_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_4_Sig_EE");
  MonitorElement* meMuon_eta_MTD_4_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_4_Sig_EE");
  MonitorElement* meMuon_phi_MTD_4_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_4_Sig_EE");

  MonitorElement* meMuon_pt_MTD_5_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_5_Sig_EE");
  MonitorElement* meMuon_eta_MTD_5_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_5_Sig_EE");
  MonitorElement* meMuon_phi_MTD_5_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_5_Sig_EE");

  MonitorElement* meMuon_pt_MTD_6_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_6_Sig_EE");
  MonitorElement* meMuon_eta_MTD_6_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_6_Sig_EE");
  MonitorElement* meMuon_phi_MTD_6_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_6_Sig_EE");

  MonitorElement* meMuon_pt_MTD_7_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_7_Sig_EE");
  MonitorElement* meMuon_eta_MTD_7_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_7_Sig_EE");
  MonitorElement* meMuon_phi_MTD_7_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_7_Sig_EE");

  MonitorElement* meMuon_pt_sim_MTD_1_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_1_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_2_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_3_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_4_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_5_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_5_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_6_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_6_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_7_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_7_Sig_EE");

  MonitorElement* meMuon_pt_sim_MTD_2sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2sigma_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3sigma_Sig_EE");
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4sigma_Sig_EE");

  MonitorElement* meMuon_pt_MTD_2sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_2sigma_Sig_EE");
  MonitorElement* meMuon_eta_MTD_2sigma_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_2sigma_Sig_EE");
  MonitorElement* meMuon_phi_MTD_2sigma_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_2sigma_Sig_EE");

  MonitorElement* meMuon_pt_MTD_3sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_3sigma_Sig_EE");
  MonitorElement* meMuon_eta_MTD_3sigma_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_3sigma_Sig_EE");
  MonitorElement* meMuon_phi_MTD_3sigma_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_3sigma_Sig_EE");

  MonitorElement* meMuon_pt_MTD_4sigma_Sig_EE_ = igetter.get(folder_ + "Muon_pT_MTD_4sigma_Sig_EE");
  MonitorElement* meMuon_eta_MTD_4sigma_Sig_EE_ = igetter.get(folder_ + "Muon_eta_MTD_4sigma_Sig_EE");
  MonitorElement* meMuon_phi_MTD_4sigma_Sig_EE_ = igetter.get(folder_ + "Muon_phi_MTD_4sigma_Sig_EE");

  // for non-promt (background)
  MonitorElement* meMuon_pt_tot_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_tot_Bkg_EB");
  MonitorElement* meMuon_pt_sim_tot_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_tot_Bkg_EB");
  MonitorElement* meMuon_pt_noMTD_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_noMTD_Bkg_EB");
  MonitorElement* meMuon_pt_gen_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_gen_Bkg_EB");

  MonitorElement* meMuon_eta_tot_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_tot_Bkg_EB");
  MonitorElement* meMuon_eta_noMTD_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_noMTD_Bkg_EB");

  MonitorElement* meMuon_phi_tot_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_tot_Bkg_EB");
  MonitorElement* meMuon_phi_noMTD_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_noMTD_Bkg_EB");
  MonitorElement* meMuon_pt_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_1_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_1_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_1_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_2_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_2_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_2_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_3_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_3_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_3_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_4_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_4_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_4_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_5_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_5_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_5_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_6_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_6_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_6_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_7_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_7_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_7_Bkg_EB");

  MonitorElement* meMuon_pt_sim_MTD_1_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_1_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_2_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_3_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_4_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_5_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_5_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_6_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_6_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_7_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_7_Bkg_EB");

  MonitorElement* meMuon_pt_sim_MTD_2sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2sigma_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3sigma_Bkg_EB");
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4sigma_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_2sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_2sigma_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_2sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_2sigma_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_2sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_2sigma_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_3sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_3sigma_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_3sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_3sigma_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_3sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_3sigma_Bkg_EB");

  MonitorElement* meMuon_pt_MTD_4sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_pT_MTD_4sigma_Bkg_EB");
  MonitorElement* meMuon_eta_MTD_4sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_eta_MTD_4sigma_Bkg_EB");
  MonitorElement* meMuon_phi_MTD_4sigma_Bkg_EB_ = igetter.get(folder_ + "Muon_phi_MTD_4sigma_Bkg_EB");

  MonitorElement* meMuon_pt_tot_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_tot_Bkg_EE");
  MonitorElement* meMuon_pt_sim_tot_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_tot_Bkg_EE");
  MonitorElement* meMuon_pt_noMTD_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_noMTD_Bkg_EE");
  MonitorElement* meMuon_pt_gen_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_gen_Bkg_EE");

  MonitorElement* meMuon_eta_tot_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_tot_Bkg_EE");
  MonitorElement* meMuon_eta_noMTD_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_noMTD_Bkg_EE");

  MonitorElement* meMuon_phi_tot_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_tot_Bkg_EE");
  MonitorElement* meMuon_phi_noMTD_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_noMTD_Bkg_EE");
  MonitorElement* meMuon_pt_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_1_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_1_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_1_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_2_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_2_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_2_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_3_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_3_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_3_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_4_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_4_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_4_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_5_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_5_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_5_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_6_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_6_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_6_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_7_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_7_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_7_Bkg_EE");

  MonitorElement* meMuon_pt_sim_MTD_1_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_1_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_2_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_3_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_4_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_5_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_5_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_6_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_6_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_7_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_7_Bkg_EE");

  MonitorElement* meMuon_pt_sim_MTD_2sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_2sigma_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_3sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_3sigma_Bkg_EE");
  MonitorElement* meMuon_pt_sim_MTD_4sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_sim_MTD_4sigma_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_2sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_2sigma_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_2sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_2sigma_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_2sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_2sigma_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_3sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_3sigma_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_3sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_3sigma_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_3sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_3sigma_Bkg_EE");

  MonitorElement* meMuon_pt_MTD_4sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_pT_MTD_4sigma_Bkg_EE");
  MonitorElement* meMuon_eta_MTD_4sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_eta_MTD_4sigma_Bkg_EE");
  MonitorElement* meMuon_phi_MTD_4sigma_Bkg_EE_ = igetter.get(folder_ + "Muon_phi_MTD_4sigma_Bkg_EE");
  if (optionalPlots_) {
    if (!meMuon_pt_tot_Sig_EB_ || !meMuon_pt_MTD_1_Sig_EB_ || !meMuon_pt_MTD_2_Sig_EB_ || !meMuon_pt_MTD_3_Sig_EB_ ||
        !meMuon_pt_MTD_4_Sig_EB_ || !meMuon_pt_MTD_5_Sig_EB_ || !meMuon_pt_MTD_6_Sig_EB_ || !meMuon_pt_MTD_7_Sig_EB_ ||
        !meMuon_pt_noMTD_Sig_EB_ || !meMuon_eta_tot_Sig_EB_ || !meMuon_eta_MTD_1_Sig_EB_ || !meMuon_eta_MTD_2_Sig_EB_ ||
        !meMuon_eta_MTD_3_Sig_EB_ || !meMuon_eta_MTD_4_Sig_EB_ || !meMuon_eta_MTD_5_Sig_EB_ || !meMuon_eta_MTD_6_Sig_EB_ ||
        !meMuon_eta_MTD_7_Sig_EB_ || !meMuon_eta_noMTD_Sig_EB_ || !meMuon_phi_tot_Sig_EB_ || !meMuon_phi_MTD_1_Sig_EB_ ||
        !meMuon_phi_MTD_2_Sig_EB_ || !meMuon_phi_MTD_3_Sig_EB_ || !meMuon_phi_MTD_4_Sig_EB_ || !meMuon_phi_MTD_5_Sig_EB_ ||
        !meMuon_phi_MTD_6_Sig_EB_ || !meMuon_phi_MTD_7_Sig_EB_ || !meMuon_phi_noMTD_Sig_EB_ || !meMuon_pt_tot_Sig_EE_ ||
        !meMuon_pt_MTD_1_Sig_EE_ || !meMuon_pt_MTD_2_Sig_EE_ || !meMuon_pt_MTD_3_Sig_EE_ || !meMuon_pt_MTD_4_Sig_EE_ ||
        !meMuon_pt_MTD_5_Sig_EE_ || !meMuon_pt_MTD_6_Sig_EE_ || !meMuon_pt_MTD_7_Sig_EE_ || !meMuon_pt_noMTD_Sig_EE_ ||
        !meMuon_eta_tot_Sig_EE_ || !meMuon_eta_MTD_1_Sig_EE_ || !meMuon_eta_MTD_2_Sig_EE_ || !meMuon_eta_MTD_3_Sig_EE_ ||
        !meMuon_eta_MTD_4_Sig_EE_ || !meMuon_eta_MTD_5_Sig_EE_ || !meMuon_eta_MTD_6_Sig_EE_ || !meMuon_eta_MTD_7_Sig_EE_ ||
        !meMuon_eta_noMTD_Sig_EE_ || !meMuon_phi_tot_Sig_EE_ || !meMuon_phi_MTD_1_Sig_EE_ || !meMuon_phi_MTD_2_Sig_EE_ ||
        !meMuon_phi_MTD_3_Sig_EE_ || !meMuon_phi_MTD_4_Sig_EE_ || !meMuon_phi_MTD_5_Sig_EE_ || !meMuon_phi_MTD_6_Sig_EE_ ||
        !meMuon_phi_MTD_7_Sig_EE_ || !meMuon_phi_noMTD_Sig_EE_ || !meMuon_pt_tot_Bkg_EB_ || !meMuon_pt_MTD_1_Bkg_EB_ ||
        !meMuon_pt_MTD_2_Bkg_EB_ || !meMuon_pt_MTD_3_Bkg_EB_ || !meMuon_pt_MTD_4_Bkg_EB_ || !meMuon_pt_MTD_5_Bkg_EB_ ||
        !meMuon_pt_MTD_6_Bkg_EB_ || !meMuon_pt_MTD_7_Bkg_EB_ || !meMuon_pt_noMTD_Bkg_EB_ || !meMuon_eta_tot_Bkg_EB_ ||
        !meMuon_eta_MTD_1_Bkg_EB_ || !meMuon_eta_MTD_2_Bkg_EB_ || !meMuon_eta_MTD_3_Bkg_EB_ || !meMuon_eta_MTD_4_Bkg_EB_ ||
        !meMuon_eta_MTD_5_Bkg_EB_ || !meMuon_eta_MTD_6_Bkg_EB_ || !meMuon_eta_MTD_7_Bkg_EB_ || !meMuon_eta_noMTD_Bkg_EB_ ||
        !meMuon_phi_tot_Bkg_EB_ || !meMuon_phi_MTD_1_Bkg_EB_ || !meMuon_phi_MTD_2_Bkg_EB_ || !meMuon_phi_MTD_3_Bkg_EB_ ||
        !meMuon_phi_MTD_4_Bkg_EB_ || !meMuon_phi_MTD_5_Bkg_EB_ || !meMuon_phi_MTD_6_Bkg_EB_ || !meMuon_phi_MTD_7_Bkg_EB_ ||
        !meMuon_phi_noMTD_Bkg_EB_ || !meMuon_pt_tot_Bkg_EE_ || !meMuon_pt_MTD_1_Bkg_EE_ || !meMuon_pt_MTD_2_Bkg_EE_ ||
        !meMuon_pt_MTD_3_Bkg_EE_ || !meMuon_pt_MTD_4_Bkg_EE_ || !meMuon_pt_MTD_5_Bkg_EE_ || !meMuon_pt_MTD_6_Bkg_EE_ ||
        !meMuon_pt_MTD_7_Bkg_EE_ || !meMuon_pt_noMTD_Bkg_EE_ || !meMuon_eta_tot_Bkg_EE_ || !meMuon_eta_MTD_1_Bkg_EE_ ||
        !meMuon_eta_MTD_2_Bkg_EE_ || !meMuon_eta_MTD_3_Bkg_EE_ || !meMuon_eta_MTD_4_Bkg_EE_ || !meMuon_eta_MTD_5_Bkg_EE_ ||
        !meMuon_eta_MTD_6_Bkg_EE_ || !meMuon_eta_MTD_7_Bkg_EE_ || !meMuon_eta_noMTD_Bkg_EE_ || !meMuon_phi_tot_Bkg_EE_ ||
        !meMuon_phi_MTD_1_Bkg_EE_ || !meMuon_phi_MTD_2_Bkg_EE_ || !meMuon_phi_MTD_3_Bkg_EE_ || !meMuon_phi_MTD_4_Bkg_EE_ ||
        !meMuon_phi_MTD_5_Bkg_EE_ || !meMuon_phi_MTD_6_Bkg_EE_ || !meMuon_phi_MTD_7_Bkg_EE_ || !meMuon_phi_noMTD_Bkg_EE_ ||
        !meMuon_pt_sim_MTD_1_Sig_EB_ || !meMuon_pt_sim_MTD_2_Sig_EB_ || !meMuon_pt_sim_MTD_3_Sig_EB_ ||
        !meMuon_pt_sim_MTD_4_Sig_EB_ || !meMuon_pt_sim_MTD_5_Sig_EB_ || !meMuon_pt_sim_MTD_6_Sig_EB_ ||
        !meMuon_pt_sim_MTD_7_Sig_EB_ || !meMuon_pt_sim_MTD_1_Sig_EE_ || !meMuon_pt_sim_MTD_2_Sig_EE_ ||
        !meMuon_pt_sim_MTD_3_Sig_EE_ || !meMuon_pt_sim_MTD_4_Sig_EE_ || !meMuon_pt_sim_MTD_5_Sig_EE_ ||
        !meMuon_pt_sim_MTD_6_Sig_EE_ || !meMuon_pt_sim_MTD_7_Sig_EE_ || !meMuon_pt_sim_MTD_1_Bkg_EB_ ||
        !meMuon_pt_sim_MTD_2_Bkg_EB_ || !meMuon_pt_sim_MTD_3_Bkg_EB_ || !meMuon_pt_sim_MTD_4_Bkg_EB_ ||
        !meMuon_pt_sim_MTD_5_Bkg_EB_ || !meMuon_pt_sim_MTD_6_Bkg_EB_ || !meMuon_pt_sim_MTD_7_Bkg_EB_ ||
        !meMuon_pt_sim_MTD_1_Bkg_EE_ || !meMuon_pt_sim_MTD_2_Bkg_EE_ || !meMuon_pt_sim_MTD_3_Bkg_EE_ ||
        !meMuon_pt_sim_MTD_4_Bkg_EE_ || !meMuon_pt_sim_MTD_5_Bkg_EE_ || !meMuon_pt_sim_MTD_6_Bkg_EE_ ||
        !meMuon_pt_sim_MTD_7_Bkg_EE_ || !meMuon_pt_MTD_4sigma_Sig_EB_ || !meMuon_pt_MTD_3sigma_Sig_EB_ ||
        !meMuon_pt_MTD_2sigma_Sig_EB_ || !meMuon_pt_MTD_4sigma_Sig_EE_ || !meMuon_pt_MTD_3sigma_Sig_EE_ ||
        !meMuon_pt_MTD_2sigma_Sig_EE_ || !meMuon_pt_MTD_4sigma_Bkg_EB_ || !meMuon_pt_MTD_3sigma_Bkg_EB_ ||
        !meMuon_pt_MTD_2sigma_Bkg_EB_ || !meMuon_pt_MTD_4sigma_Bkg_EE_ || !meMuon_pt_MTD_3sigma_Bkg_EE_ ||
        !meMuon_pt_MTD_2sigma_Bkg_EE_ || !meMuon_pt_sim_MTD_4sigma_Sig_EB_ || !meMuon_pt_sim_MTD_3sigma_Sig_EB_ ||
        !meMuon_pt_sim_MTD_2sigma_Sig_EB_ || !meMuon_pt_sim_MTD_4sigma_Sig_EE_ || !meMuon_pt_sim_MTD_3sigma_Sig_EE_ ||
        !meMuon_pt_sim_MTD_2sigma_Sig_EE_ || !meMuon_pt_sim_MTD_4sigma_Bkg_EB_ || !meMuon_pt_sim_MTD_3sigma_Bkg_EB_ ||
        !meMuon_pt_sim_MTD_2sigma_Bkg_EB_ || !meMuon_pt_sim_MTD_4sigma_Bkg_EE_ || !meMuon_pt_sim_MTD_3sigma_Bkg_EE_ ||
        !meMuon_pt_sim_MTD_2sigma_Bkg_EE_ || !meMuon_pt_gen_Sig_EB_ || !meMuon_pt_gen_Sig_EE_ || !meMuon_pt_gen_Bkg_EB_ ||
        !meMuon_pt_gen_Bkg_EE_ || !meMuon_eta_MTD_4sigma_Sig_EB_ || !meMuon_eta_MTD_3sigma_Sig_EB_ ||
        !meMuon_eta_MTD_2sigma_Sig_EB_ || !meMuon_phi_MTD_4sigma_Sig_EB_ || !meMuon_phi_MTD_3sigma_Sig_EB_ ||
        !meMuon_phi_MTD_2sigma_Sig_EB_ || !meMuon_eta_MTD_4sigma_Sig_EE_ || !meMuon_eta_MTD_3sigma_Sig_EE_ ||
        !meMuon_eta_MTD_2sigma_Sig_EE_ || !meMuon_phi_MTD_4sigma_Sig_EE_ || !meMuon_phi_MTD_3sigma_Sig_EE_ ||
        !meMuon_phi_MTD_2sigma_Sig_EE_ || !meMuon_eta_MTD_4sigma_Bkg_EB_ || !meMuon_eta_MTD_3sigma_Bkg_EB_ ||
        !meMuon_eta_MTD_2sigma_Bkg_EB_ || !meMuon_phi_MTD_4sigma_Bkg_EB_ || !meMuon_phi_MTD_3sigma_Bkg_EB_ ||
        !meMuon_phi_MTD_2sigma_Bkg_EB_ || !meMuon_eta_MTD_4sigma_Bkg_EE_ || !meMuon_eta_MTD_3sigma_Bkg_EE_ ||
        !meMuon_eta_MTD_2sigma_Bkg_EE_ || !meMuon_phi_MTD_4sigma_Bkg_EE_ || !meMuon_phi_MTD_3sigma_Bkg_EE_ ||
        !meMuon_phi_MTD_2sigma_Bkg_EE_) {
      edm::LogError("MtdMuonIsoHarvester") << "Monitoring histograms not found!" << std::endl;
      return;
    }
  } else {
    if (!meMuon_pt_tot_Sig_EB_ || !meMuon_pt_noMTD_Sig_EB_ || !meMuon_eta_tot_Sig_EB_ || !meMuon_eta_noMTD_Sig_EB_ ||
        !meMuon_phi_tot_Sig_EB_ || !meMuon_phi_noMTD_Sig_EB_ || !meMuon_pt_tot_Sig_EE_ || !meMuon_pt_noMTD_Sig_EE_ ||
        !meMuon_eta_tot_Sig_EE_ || !meMuon_eta_noMTD_Sig_EE_ || !meMuon_phi_tot_Sig_EE_ || !meMuon_phi_noMTD_Sig_EE_ ||
        !meMuon_pt_tot_Bkg_EB_ || !meMuon_pt_noMTD_Bkg_EB_ || !meMuon_eta_tot_Bkg_EB_ || !meMuon_eta_noMTD_Bkg_EB_ ||
        !meMuon_phi_tot_Bkg_EB_ || !meMuon_phi_noMTD_Bkg_EB_ || !meMuon_pt_tot_Bkg_EE_ || !meMuon_pt_noMTD_Bkg_EE_ ||
        !meMuon_eta_tot_Bkg_EE_ || !meMuon_eta_noMTD_Bkg_EE_ || !meMuon_phi_tot_Bkg_EE_ || !meMuon_phi_noMTD_Bkg_EE_ ||
        !meMuon_pt_MTD_4sigma_Sig_EB_ || !meMuon_pt_MTD_3sigma_Sig_EB_ || !meMuon_pt_MTD_2sigma_Sig_EB_ ||
        !meMuon_pt_MTD_4sigma_Sig_EE_ || !meMuon_pt_MTD_3sigma_Sig_EE_ || !meMuon_pt_MTD_2sigma_Sig_EE_ ||
        !meMuon_pt_MTD_4sigma_Bkg_EB_ || !meMuon_pt_MTD_3sigma_Bkg_EB_ || !meMuon_pt_MTD_2sigma_Bkg_EB_ ||
        !meMuon_pt_MTD_4sigma_Bkg_EE_ || !meMuon_pt_MTD_3sigma_Bkg_EE_ || !meMuon_pt_MTD_2sigma_Bkg_EE_ ||
        !meMuon_eta_MTD_4sigma_Sig_EB_ || !meMuon_eta_MTD_3sigma_Sig_EB_ || !meMuon_eta_MTD_2sigma_Sig_EB_ ||
        !meMuon_phi_MTD_4sigma_Sig_EB_ || !meMuon_phi_MTD_3sigma_Sig_EB_ || !meMuon_phi_MTD_2sigma_Sig_EB_ ||
        !meMuon_eta_MTD_4sigma_Sig_EE_ || !meMuon_eta_MTD_3sigma_Sig_EE_ || !meMuon_eta_MTD_2sigma_Sig_EE_ ||
        !meMuon_phi_MTD_4sigma_Sig_EE_ || !meMuon_phi_MTD_3sigma_Sig_EE_ || !meMuon_phi_MTD_2sigma_Sig_EE_ ||
        !meMuon_eta_MTD_4sigma_Bkg_EB_ || !meMuon_eta_MTD_3sigma_Bkg_EB_ || !meMuon_eta_MTD_2sigma_Bkg_EB_ ||
        !meMuon_phi_MTD_4sigma_Bkg_EB_ || !meMuon_phi_MTD_3sigma_Bkg_EB_ || !meMuon_phi_MTD_2sigma_Bkg_EB_ ||
        !meMuon_eta_MTD_4sigma_Bkg_EE_ || !meMuon_eta_MTD_3sigma_Bkg_EE_ || !meMuon_eta_MTD_2sigma_Bkg_EE_ ||
        !meMuon_phi_MTD_4sigma_Bkg_EE_ || !meMuon_phi_MTD_3sigma_Bkg_EE_ || !meMuon_phi_MTD_2sigma_Bkg_EE_) {
      edm::LogError("MtdMuonIsoHarvester") << "Monitoring histograms not found!" << std::endl;
      return;
    }
  }
  // --- Book  histograms
  ibook.cd(folder_);
  // Muon iso addition starts; MTD vs noMTD case
  /////////////////////////////////////////////////////// For promt (signal)
  if (optionalPlots_) {
    mePtEffMTD_1_Sig_EB_ = ibook.book1D("pTeffMTD_1_Sig_EB",
                                        " MTD isolation Efficiency 1 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_1_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_1_Sig_EB_);

    mePtEffMTD_2_Sig_EB_ = ibook.book1D("pTeffMTD_2_Sig_EB",
                                        " MTD isolation Efficiency 2 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_2_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_2_Sig_EB_);

    mePtEffMTD_3_Sig_EB_ = ibook.book1D("pTeffMTD_3_Sig_EB",
                                        " MTD isolation Efficiency 3 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_3_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_3_Sig_EB_);

    mePtEffMTD_4_Sig_EB_ = ibook.book1D("pTeffMTD_4_Sig_EB",
                                        " MTD isolation Efficiency 4 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_4_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_4_Sig_EB_);

    mePtEffMTD_5_Sig_EB_ = ibook.book1D("pTeffMTD_5_Sig_EB",
                                        " MTD isolation Efficiency 5 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_5_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_5_Sig_EB_);

    mePtEffMTD_6_Sig_EB_ = ibook.book1D("pTeffMTD_6_Sig_EB",
                                        " MTD isolation Efficiency 6 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_6_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_6_Sig_EB_);

    mePtEffMTD_7_Sig_EB_ = ibook.book1D("pTeffMTD_7_Sig_EB",
                                        " MTD isolation Efficiency 7 Signal Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_7_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_7_Sig_EB_);

    mePtEffMTD_sim_1_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_1_Sig_EB",
                                            " MTD isolation Efficiency 1 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_1_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_1_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_1_Sig_EB_);

    mePtEffMTD_sim_2_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_2_Sig_EB",
                                            " MTD isolation Efficiency 2 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_2_Sig_EB_);

    mePtEffMTD_sim_3_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_3_Sig_EB",
                                            " MTD isolation Efficiency 3 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_3_Sig_EB_);

    mePtEffMTD_sim_4_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_4_Sig_EB",
                                            " MTD isolation Efficiency 4 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_4_Sig_EB_);

    mePtEffMTD_sim_5_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_5_Sig_EB",
                                            " MTD isolation Efficiency 5 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_5_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_5_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_5_Sig_EB_);

    mePtEffMTD_sim_6_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_6_Sig_EB",
                                            " MTD isolation Efficiency 6 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_6_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_6_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_6_Sig_EB_);

    mePtEffMTD_sim_7_Sig_EB_ = ibook.book1D("mePtEffMTD_sim_7_Sig_EB",
                                            " MTD isolation Efficiency 7 SIM Signal barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_7_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_7_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_7_Sig_EB_);

    mePtEffMTD_sim_4sigma_Sig_EB_ =
        ibook.book1D("pTeffMTD_sim_4sigma_Sig_EB",
                     " MTD isolation Efficiency SIM - 4 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4sigma_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_4sigma_Sig_EB_);

    mePtEffMTD_sim_3sigma_Sig_EB_ =
        ibook.book1D("pTeffMTD_sim_3sigma_Sig_EB",
                     " MTD isolation Efficiency SIM - 3 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3sigma_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_3sigma_Sig_EB_);

    mePtEffMTD_sim_2sigma_Sig_EB_ =
        ibook.book1D("pTeffMTD_sim_2sigma_Sig_EB",
                     " MTD isolation Efficiency SIM - 2 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2sigma_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_sim_2sigma_Sig_EB_);

    meEtaEffMTD_1_Sig_EB_ = ibook.book1D("EtaEffMTD_1_Sig_EB",
                                         " MTD isolation Efficiency 1 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_1_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_1_Sig_EB_);

    meEtaEffMTD_2_Sig_EB_ = ibook.book1D("EtaEffMTD_2_Sig_EB",
                                         " MTD isolation Efficiency 2 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_2_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_2_Sig_EB_);

    meEtaEffMTD_3_Sig_EB_ = ibook.book1D("EtaEffMTD_3_Sig_EB",
                                         " MTD isolation Efficiency 3 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_3_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_3_Sig_EB_);

    meEtaEffMTD_4_Sig_EB_ = ibook.book1D("EtaEffMTD_4_Sig_EB",
                                         " MTD isolation Efficiency 4 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_4_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_4_Sig_EB_);

    meEtaEffMTD_5_Sig_EB_ = ibook.book1D("EtaEffMTD_5_Sig_EB",
                                         " MTD isolation Efficiency 5 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_5_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_5_Sig_EB_);

    meEtaEffMTD_6_Sig_EB_ = ibook.book1D("EtaEffMTD_6_Sig_EB",
                                         " MTD isolation Efficiency 6 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_6_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_6_Sig_EB_);

    meEtaEffMTD_7_Sig_EB_ = ibook.book1D("EtaEffMTD_7_Sig_EB",
                                         " MTD isolation Efficiency 7 Signal Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_7_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_7_Sig_EB_);

    mePhiEffMTD_1_Sig_EB_ = ibook.book1D("PhiEffMTD_1_Sig_EB",
                                         " MTD isolation Efficiency 1 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_1_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_1_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_1_Sig_EB_);

    mePhiEffMTD_2_Sig_EB_ = ibook.book1D("PhiEffMTD_2_Sig_EB",
                                         " MTD isolation Efficiency 2 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_2_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_2_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_2_Sig_EB_);

    mePhiEffMTD_3_Sig_EB_ = ibook.book1D("PhiEffMTD_3_Sig_EB",
                                         " MTD isolation Efficiency 3 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_3_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_3_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_3_Sig_EB_);

    mePhiEffMTD_4_Sig_EB_ = ibook.book1D("PhiEffMTD_4_Sig_EB",
                                         " MTD isolation Efficiency 4 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_4_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_4_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_4_Sig_EB_);

    mePhiEffMTD_5_Sig_EB_ = ibook.book1D("PhiEffMTD_5_Sig_EB",
                                         " MTD isolation Efficiency 5 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_5_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_5_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_5_Sig_EB_);

    mePhiEffMTD_6_Sig_EB_ = ibook.book1D("PhiEffMTD_6_Sig_EB",
                                         " MTD isolation Efficiency 6 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_6_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_6_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_6_Sig_EB_);

    mePhiEffMTD_7_Sig_EB_ = ibook.book1D("PhiEffMTD_7_Sig_EB",
                                         " MTD isolation Efficiency 7 Signal Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_7_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_7_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_7_Sig_EB_);
  }

  mePtEffMTD_4sigma_Sig_EB_ =
      ibook.book1D("pTeffMTD_4sigma_Sig_EB",
                   "MTD isolation Efficiency - 4 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EB_->getNbinsX(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_4sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_4sigma_Sig_EB_);

  mePtEffMTD_3sigma_Sig_EB_ =
      ibook.book1D("pTeffMTD_3sigma_Sig_EB",
                   " MTD isolation Efficiency - 3 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EB_->getNbinsX(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_3sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_3sigma_Sig_EB_);

  mePtEffMTD_2sigma_Sig_EB_ =
      ibook.book1D("pTeffMTD_2sigma_Sig_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EB_->getNbinsX(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_2sigma_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffMTD_2sigma_Sig_EB_);

  meEtaEffMTD_4sigma_Sig_EB_ =
      ibook.book1D("EtaEffMTD_4sigma_Sig_EB",
                   " MTD isolation Efficiency - 4 sigma compatibility - Signal Barrel VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EB_->getNbinsX(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_4sigma_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_4sigma_Sig_EB_);

  meEtaEffMTD_3sigma_Sig_EB_ =
      ibook.book1D("EtaEffMTD_3sigma_Sig_EB",
                   " MTD isolation Efficiency - 3 sigma compatibility - Signal Barrel VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EB_->getNbinsX(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_3sigma_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_3sigma_Sig_EB_);

  meEtaEffMTD_2sigma_Sig_EB_ =
      ibook.book1D("EtaEffMTD_2sigma_Sig_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Barrel VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EB_->getNbinsX(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_2sigma_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffMTD_2sigma_Sig_EB_);
  mePhiEffMTD_4sigma_Sig_EB_ =
      ibook.book1D("PhiEffMTD_4sigma_Sig_EB",
                   " MTD isolation Efficiency - 4 sigma compatibility - Signal Barrel VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EB_->getNbinsX(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_4sigma_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_4sigma_Sig_EB_);

  mePhiEffMTD_3sigma_Sig_EB_ =
      ibook.book1D("PhiEffMTD_3sigma_Sig_EB",
                   " MTD isolation Efficiency - 3  sigma compatibility - Signal Barrel VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EB_->getNbinsX(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_3sigma_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_3sigma_Sig_EB_);

  mePhiEffMTD_2sigma_Sig_EB_ =
      ibook.book1D("PhiEffMTD_2sigma_Sig_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Barrel VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EB_->getNbinsX(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2sigma_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_2sigma_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffMTD_2sigma_Sig_EB_);

  mePtEffnoMTD_Sig_EB_ = ibook.book1D("pTeffnoMTD_Sig_EB",
                                      " noMTD isolation Efficiency - Signal Barrel VS pT;p_{T};Efficiency",
                                      meMuon_pt_tot_Sig_EB_->getNbinsX(),
                                      meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_noMTD_Sig_EB_, meMuon_pt_tot_Sig_EB_, mePtEffnoMTD_Sig_EB_);

  meEtaEffnoMTD_Sig_EB_ = ibook.book1D("EtaEffnoMTD_Sig_EB",
                                       " noMTD isolation Efficiency - Signal Barrel VS Eta;#eta;Efficiency",
                                       meMuon_eta_tot_Sig_EB_->getNbinsX(),
                                       meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_eta_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_noMTD_Sig_EB_, meMuon_eta_tot_Sig_EB_, meEtaEffnoMTD_Sig_EB_);

  mePhiEffnoMTD_Sig_EB_ = ibook.book1D("PhiEffnoMTD_Sig_EB",
                                       " noMTD isolation Efficiency - Signal Barrel VS Phi;#phi;Efficiency",
                                       meMuon_phi_tot_Sig_EB_->getNbinsX(),
                                       meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_phi_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Sig_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_noMTD_Sig_EB_, meMuon_phi_tot_Sig_EB_, mePhiEffnoMTD_Sig_EB_);

  if (optionalPlots_) {
    mePtEffgen_Sig_EB_ = ibook.book1D("pTeffMTD_gen_Sig_EB",
                                      " MTD isolation Efficiency - genInfo - Signal Barrel VS pT;p_{T};Efficiency",
                                      meMuon_pt_sim_tot_Sig_EB_->getNbinsX(),
                                      meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_sim_tot_Sig_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffgen_Sig_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_gen_Sig_EB_, meMuon_pt_sim_tot_Sig_EB_, mePtEffgen_Sig_EB_);

    // Muon iso addition ends
    // For endcap now
    mePtEffMTD_1_Sig_EE_ = ibook.book1D("pTeffMTD_1_Sig_EE",
                                        " MTD isolation Efficiency - 1 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_1_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_1_Sig_EE_);

    mePtEffMTD_2_Sig_EE_ = ibook.book1D("pTeffMTD_2_Sig_EE",
                                        " MTD isolation Efficiency - 2 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_2_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_2_Sig_EE_);

    mePtEffMTD_3_Sig_EE_ = ibook.book1D("pTeffMTD_3_Sig_EE",
                                        " MTD isolation Efficiency - 3 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_3_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_3_Sig_EE_);

    mePtEffMTD_4_Sig_EE_ = ibook.book1D("pTeffMTD_4_Sig_EE",
                                        " MTD isolation Efficiency - 4 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_4_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_4_Sig_EE_);

    mePtEffMTD_5_Sig_EE_ = ibook.book1D("pTeffMTD_5_Sig_EE",
                                        " MTD isolation Efficiency - 5 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_5_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_5_Sig_EE_);

    mePtEffMTD_6_Sig_EE_ = ibook.book1D("pTeffMTD_6_Sig_EE",
                                        " MTD isolation Efficiency - 6 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_6_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_6_Sig_EE_);

    mePtEffMTD_7_Sig_EE_ = ibook.book1D("pTeffMTD_7_Sig_EE",
                                        " MTD isolation Efficiency - 7 Signal Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_7_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_7_Sig_EE_);

    mePtEffMTD_sim_1_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_1_Sig_EE_",
                                            " MTD isolation Efficiency - 1 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_1_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_1_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_1_Sig_EE_);

    mePtEffMTD_sim_2_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_2_Sig_EE_",
                                            " MTD isolation Efficiency - 2 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_2_Sig_EE_);

    mePtEffMTD_sim_3_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_3_Sig_EE_",
                                            " MTD isolation Efficiency - 3 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_3_Sig_EE_);

    mePtEffMTD_sim_4_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_4_Sig_EE_",
                                            " MTD isolation Efficiency - 4 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_4_Sig_EE_);

    mePtEffMTD_sim_5_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_5_Sig_EE_",
                                            " MTD isolation Efficiency - 5 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_5_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_5_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_5_Sig_EE_);

    mePtEffMTD_sim_6_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_6_Sig_EE_",
                                            " MTD isolation Efficiency - 6 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_6_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_6_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_6_Sig_EE_);

    mePtEffMTD_sim_7_Sig_EE_ = ibook.book1D("mePtEffMTD_sim_7_Sig_EE_",
                                            " MTD isolation Efficiency - 7 SIM Signal Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_7_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_7_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_7_Sig_EE_);

    mePtEffMTD_sim_4sigma_Sig_EE_ =
        ibook.book1D("pTeffMTD_sim_4sigma_Sig_EE",
                     " MTD isolation Efficiency SIM - 4 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4sigma_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_4sigma_Sig_EE_);

    mePtEffMTD_sim_3sigma_Sig_EE_ =
        ibook.book1D("pTeffMTD_sim_3sigma_Sig_EE",
                     " MTD isolation Efficiency SIM - 3 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3sigma_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_3sigma_Sig_EE_);

    mePtEffMTD_sim_2sigma_Sig_EE_ =
        ibook.book1D("pTeffMTD_sim_2sigma_Sig_EE",
                     " MTD isolation Efficiency SIM - 2 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2sigma_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_sim_2sigma_Sig_EE_);

    meEtaEffMTD_1_Sig_EE_ = ibook.book1D("EtaEffMTD_1_Sig_EE",
                                         " MTD isolation Efficiency 1 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_1_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_1_Sig_EE_);

    meEtaEffMTD_2_Sig_EE_ = ibook.book1D("EtaEffMTD_2_Sig_EE",
                                         " MTD isolation Efficiency 2 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_2_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_2_Sig_EE_);

    meEtaEffMTD_3_Sig_EE_ = ibook.book1D("EtaEffMTD_3_Sig_EE",
                                         " MTD isolation Efficiency 3 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_3_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_3_Sig_EE_);

    meEtaEffMTD_4_Sig_EE_ = ibook.book1D("EtaEffMTD_4_Sig_EE",
                                         " MTD isolation Efficiency 4 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_4_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_4_Sig_EE_);

    meEtaEffMTD_5_Sig_EE_ = ibook.book1D("EtaEffMTD_5_Sig_EE",
                                         " MTD isolation Efficiency 5 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_5_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_5_Sig_EE_);

    meEtaEffMTD_6_Sig_EE_ = ibook.book1D("EtaEffMTD_6_Sig_EE",
                                         " MTD isolation Efficiency 6 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_6_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_6_Sig_EE_);

    meEtaEffMTD_7_Sig_EE_ = ibook.book1D("EtaEffMTD_7_Sig_EE",
                                         " MTD isolation Efficiency 7 Signal Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_7_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_7_Sig_EE_);

    mePhiEffMTD_1_Sig_EE_ = ibook.book1D("PhiEffMTD_1_Sig_EE",
                                         " MTD isolation Efficiency - 1 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_1_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_1_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_1_Sig_EE_);

    mePhiEffMTD_2_Sig_EE_ = ibook.book1D("PhiEffMTD_2_Sig_EE",
                                         " MTD isolation Efficiency - 2 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_2_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_2_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_2_Sig_EE_);

    mePhiEffMTD_3_Sig_EE_ = ibook.book1D("PhiEffMTD_3_Sig_EE",
                                         " MTD isolation Efficiency - 3 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_3_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_3_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_3_Sig_EE_);

    mePhiEffMTD_4_Sig_EE_ = ibook.book1D("PhiEffMTD_4_Sig_EE",
                                         " MTD isolation Efficiency - 4 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_4_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_4_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_4_Sig_EE_);

    mePhiEffMTD_5_Sig_EE_ = ibook.book1D("PhiEffMTD_5_Sig_EE",
                                         " MTD isolation Efficiency - 5 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_5_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_5_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_5_Sig_EE_);

    mePhiEffMTD_6_Sig_EE_ = ibook.book1D("PhiEffMTD_6_Sig_EE",
                                         " MTD isolation Efficiency - 6 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_6_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_6_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_6_Sig_EE_);

    mePhiEffMTD_7_Sig_EE_ = ibook.book1D("PhiEffMTD_7_Sig_EE",
                                         " MTD isolation Efficiency - 7 Signal Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_7_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_7_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_7_Sig_EE_);

    mePtEffgen_Sig_EE_ = ibook.book1D("pTeffMTD_gen_Sig_EE",
                                      " MTD isolation Efficiency - genInfo - Signal Endcap VS pT;p_{T};Efficiency",
                                      meMuon_pt_sim_tot_Sig_EE_->getNbinsX(),
                                      meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_sim_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffgen_Sig_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_gen_Sig_EE_, meMuon_pt_sim_tot_Sig_EE_, mePtEffgen_Sig_EE_);
  }

  mePtEffMTD_4sigma_Sig_EE_ =
      ibook.book1D("pTeffMTD_4sigma_Sig_EE",
                   " MTD isolation Efficiency - 4 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EE_->getNbinsX(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_4sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_4sigma_Sig_EE_);

  mePtEffMTD_3sigma_Sig_EE_ =
      ibook.book1D("pTeffMTD_3sigma_Sig_EE",
                   " MTD isolation Efficiency - 3 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EE_->getNbinsX(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_3sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_3sigma_Sig_EE_);

  mePtEffMTD_2sigma_Sig_EE_ =
      ibook.book1D("pTeffMTD_2sigma_Sig_EE",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Sig_EE_->getNbinsX(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_2sigma_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffMTD_2sigma_Sig_EE_);

  meEtaEffMTD_4sigma_Sig_EE_ =
      ibook.book1D("EtaEffMTD_4sigma_Sig_EE",
                   " MTD isolation Efficiency - 4 sigma compatibility - Signal Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EE_->getNbinsX(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_4sigma_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_4sigma_Sig_EE_);

  meEtaEffMTD_3sigma_Sig_EE_ =
      ibook.book1D("EtaEffMTD_3sigma_Sig_EE",
                   " MTD isolation Efficiency - 3 sigma compatibility - Signal Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EE_->getNbinsX(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_3sigma_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_3sigma_Sig_EE_);

  meEtaEffMTD_2sigma_Sig_EE_ =
      ibook.book1D("EtaEffMTD_2sigma_Sig_EE",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Sig_EE_->getNbinsX(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_2sigma_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffMTD_2sigma_Sig_EE_);

  mePhiEffMTD_4sigma_Sig_EE_ =
      ibook.book1D("PhiEffMTD_4sigma_Sig_EE",
                   " MTD isolation Efficiency - 4 sigma compatibility - Signal Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EE_->getNbinsX(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_4sigma_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_4sigma_Sig_EE_);

  mePhiEffMTD_3sigma_Sig_EE_ =
      ibook.book1D("PhiEffMTD_3sigma_Sig_EE",
                   " MTD isolation Efficiency - 3 sigma compatibility - Signal Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EE_->getNbinsX(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_3sigma_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_3sigma_Sig_EE_);

  mePhiEffMTD_2sigma_Sig_EE_ =
      ibook.book1D("PhiEffMTD_2sigma_Sig_EE",
                   " MTD isolation Efficiency - 2 sigma compatibility - Signal Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Sig_EE_->getNbinsX(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2sigma_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_2sigma_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffMTD_2sigma_Sig_EE_);

  mePtEffnoMTD_Sig_EE_ = ibook.book1D("pTeffnoMTD_Sig_EE",
                                      " noMTD isolation Efficiency - Signal Endcap VS pT;p_{T};Efficiency",
                                      meMuon_pt_tot_Sig_EE_->getNbinsX(),
                                      meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_noMTD_Sig_EE_, meMuon_pt_tot_Sig_EE_, mePtEffnoMTD_Sig_EE_);

  meEtaEffnoMTD_Sig_EE_ = ibook.book1D("EtaEffnoMTD_Sig_EE",
                                       " noMTD isolation Efficiency - Signal Endcap VS Eta;#eta;Efficiency",
                                       meMuon_eta_tot_Sig_EE_->getNbinsX(),
                                       meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_eta_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_noMTD_Sig_EE_, meMuon_eta_tot_Sig_EE_, meEtaEffnoMTD_Sig_EE_);

  mePhiEffnoMTD_Sig_EE_ = ibook.book1D("PhiEffnoMTD_Sig_EE",
                                       " noMTD isolation Efficiency - Signal Endcap VS Phi;#phi;Efficiency",
                                       meMuon_phi_tot_Sig_EE_->getNbinsX(),
                                       meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_phi_tot_Sig_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Sig_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_noMTD_Sig_EE_, meMuon_phi_tot_Sig_EE_, mePhiEffnoMTD_Sig_EE_);

  /////////////////////////////////////////////////////// For non-promt (background)
  if (optionalPlots_) {
    mePtEffMTD_1_Bkg_EB_ = ibook.book1D("pTeffMTD_1_Bkg_EB",
                                        " MTD isolation Efficiency - 1 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_1_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_1_Bkg_EB_);

    mePtEffMTD_2_Bkg_EB_ = ibook.book1D("pTeffMTD_2_Bkg_EB",
                                        " MTD isolation Efficiency - 2 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_2_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_2_Bkg_EB_);

    mePtEffMTD_3_Bkg_EB_ = ibook.book1D("pTeffMTD_3_Bkg_EB",
                                        " MTD isolation Efficiency - 3 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_3_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_3_Bkg_EB_);

    mePtEffMTD_4_Bkg_EB_ = ibook.book1D("pTeffMTD_4_Bkg_EB",
                                        " MTD isolation Efficiency - 4 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_4_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_4_Bkg_EB_);

    mePtEffMTD_5_Bkg_EB_ = ibook.book1D("pTeffMTD_5_Bkg_EB",
                                        " MTD isolation Efficiency - 5 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_5_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_5_Bkg_EB_);

    mePtEffMTD_6_Bkg_EB_ = ibook.book1D("pTeffMTD_6_Bkg_EB",
                                        " MTD isolation Efficiency - 6 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_6_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_6_Bkg_EB_);

    mePtEffMTD_7_Bkg_EB_ = ibook.book1D("pTeffMTD_7_Bkg_EB",
                                        " MTD isolation Efficiency - 7 Bkg Barrel VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_7_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_7_Bkg_EB_);

    mePtEffMTD_sim_1_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_1_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 1 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_1_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_1_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_1_Bkg_EB_);

    mePtEffMTD_sim_2_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_2_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 2 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_2_Bkg_EB_);

    mePtEffMTD_sim_3_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_3_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 3 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_3_Bkg_EB_);

    mePtEffMTD_sim_4_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_4_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 4 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_4_Bkg_EB_);

    mePtEffMTD_sim_5_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_5_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 5 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_5_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_5_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_5_Bkg_EB_);

    mePtEffMTD_sim_6_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_6_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 6 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_6_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_6_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_6_Bkg_EB_);

    mePtEffMTD_sim_7_Bkg_EB_ = ibook.book1D("mePtEffMTD_sim_7_Bkg_EB_",
                                            " MTD isolation Efficiency SIM - 7 Bkg Barrel VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_7_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_7_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_7_Bkg_EB_);

    mePtEffMTD_sim_4sigma_Bkg_EB_ =
        ibook.book1D("pTeffMTD_sim_4sigma_Bkg_EB",
                     " MTD isolation Efficiency SIM - 4 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_4sigma_Bkg_EB_);

    mePtEffMTD_sim_3sigma_Bkg_EB_ =
        ibook.book1D("pTeffMTD_sim_3sigma_Bkg_EB",
                     " MTD isolation Efficiency SIM - 3 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_3sigma_Bkg_EB_);

    mePtEffMTD_sim_2sigma_Bkg_EB_ =
        ibook.book1D("pTeffMTD_sim_2sigma_Bkg_EB",
                     " MTD isolation Efficiency SIM - 2 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_sim_2sigma_Bkg_EB_);

    meEtaEffMTD_1_Bkg_EB_ = ibook.book1D("EtaEffMTD_1_Bkg_EB",
                                         " MTD isolation Efficiency - 1 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_1_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_1_Bkg_EB_);

    meEtaEffMTD_2_Bkg_EB_ = ibook.book1D("EtaEffMTD_2_Bkg_EB",
                                         " MTD isolation Efficiency - 2 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_2_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_2_Bkg_EB_);

    meEtaEffMTD_3_Bkg_EB_ = ibook.book1D("EtaEffMTD_3_Bkg_EB",
                                         " MTD isolation Efficiency - 3 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_3_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_3_Bkg_EB_);

    meEtaEffMTD_4_Bkg_EB_ = ibook.book1D("EtaEffMTD_4_Bkg_EB",
                                         " MTD isolation Efficiency - 4 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_4_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_4_Bkg_EB_);

    meEtaEffMTD_5_Bkg_EB_ = ibook.book1D("EtaEffMTD_5_Bkg_EB",
                                         " MTD isolation Efficiency - 5 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_5_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_5_Bkg_EB_);

    meEtaEffMTD_6_Bkg_EB_ = ibook.book1D("EtaEffMTD_6_Bkg_EB",
                                         " MTD isolation Efficiency - 6 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_6_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_6_Bkg_EB_);

    meEtaEffMTD_7_Bkg_EB_ = ibook.book1D("EtaEffMTD_7_Bkg_EB",
                                         " MTD isolation Efficiency - 7 Bkg Barrel VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_7_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_7_Bkg_EB_);

    mePhiEffMTD_1_Bkg_EB_ = ibook.book1D("PhiEffMTD_1_Bkg_EB",
                                         " MTD isolation Efficiency - 1 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_1_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_1_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_1_Bkg_EB_);

    mePhiEffMTD_2_Bkg_EB_ = ibook.book1D("PhiEffMTD_2_Bkg_EB",
                                         " MTD isolation Efficiency - 2 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_2_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_2_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_2_Bkg_EB_);

    mePhiEffMTD_3_Bkg_EB_ = ibook.book1D("PhiEffMTD_3_Bkg_EB",
                                         " MTD isolation Efficiency - 3 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_3_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_3_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_3_Bkg_EB_);

    mePhiEffMTD_4_Bkg_EB_ = ibook.book1D("PhiEffMTD_4_Bkg_EB",
                                         " MTD isolation Efficiency - 4 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_4_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_4_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_4_Bkg_EB_);

    mePhiEffMTD_5_Bkg_EB_ = ibook.book1D("PhiEffMTD_5_Bkg_EB",
                                         " MTD isolation Efficiency - 5 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_5_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_5_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_5_Bkg_EB_);

    mePhiEffMTD_6_Bkg_EB_ = ibook.book1D("PhiEffMTD_6_Bkg_EB",
                                         " MTD isolation Efficiency - 6 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_6_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_6_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_6_Bkg_EB_);

    mePhiEffMTD_7_Bkg_EB_ = ibook.book1D("PhiEffMTD_7_Bkg_EB",
                                         " MTD isolation Efficiency - 7 Bkg Barrel VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_7_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_7_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_7_Bkg_EB_);
  }
  meEtaEffMTD_4sigma_Bkg_EB_ =
      ibook.book1D("EtaEffMTD_4sigma_Bkg_EB",
                   " MTD isolation Efficiency - 4 sigma compatibility - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_4sigma_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_4sigma_Bkg_EB_);

  meEtaEffMTD_3sigma_Bkg_EB_ =
      ibook.book1D("EtaEffMTD_3sigma_Bkg_EB",
                   " MTD isolation Efficiency - 3 sigma compatibility - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_3sigma_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_3sigma_Bkg_EB_);

  meEtaEffMTD_2sigma_Bkg_EB_ =
      ibook.book1D("EtaEffMTD_2sigma_Bkg_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_2sigma_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffMTD_2sigma_Bkg_EB_);

  mePhiEffMTD_4sigma_Bkg_EB_ =
      ibook.book1D("PhiEffMTD_4sigma_Bkg_EB",
                   " MTD isolation Efficiency - 4 sigma compatibility - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_4sigma_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_4sigma_Bkg_EB_);

  mePhiEffMTD_3sigma_Bkg_EB_ =
      ibook.book1D("PhiEffMTD_3sigma_Bkg_EB",
                   " MTD isolation Efficiency - 3 sigma compatibility - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_3sigma_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_3sigma_Bkg_EB_);

  mePhiEffMTD_2sigma_Bkg_EB_ =
      ibook.book1D("PhiEffMTD_2sigma_Bkg_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_2sigma_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffMTD_2sigma_Bkg_EB_);

  mePtEffMTD_4sigma_Bkg_EB_ =
      ibook.book1D("pTeffMTD_4sigma_Bkg_EB",
                   " MTD isolation Efficiency - 4 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_4sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_4sigma_Bkg_EB_);

  mePtEffMTD_3sigma_Bkg_EB_ =
      ibook.book1D("pTeffMTD_3sigma_Bkg_EB",
                   " MTD isolation Efficiency - 3 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_3sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_3sigma_Bkg_EB_);

  mePtEffMTD_2sigma_Bkg_EB_ =
      ibook.book1D("pTeffMTD_2sigma_Bkg_EB",
                   " MTD isolation Efficiency - 2 sigma compatibility - Bkg Barrel VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2sigma_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_2sigma_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffMTD_2sigma_Bkg_EB_);
  mePtEffnoMTD_Bkg_EB_ = ibook.book1D("pTeffnoMTD_Bkg_EB",
                                      " noMTD isolation Efficiency - Bkg Barrel VS pT;p_{T};Efficiency",
                                      meMuon_pt_tot_Bkg_EB_->getNbinsX(),
                                      meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_noMTD_Bkg_EB_, meMuon_pt_tot_Bkg_EB_, mePtEffnoMTD_Bkg_EB_);

  meEtaEffnoMTD_Bkg_EB_ = ibook.book1D("EtaEffnoMTD_Bkg_EB",
                                       " noMTD isolation Efficiency - Bkg Barrel VS pT;p_{T};Efficiency",
                                       meMuon_eta_tot_Bkg_EB_->getNbinsX(),
                                       meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_eta_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_noMTD_Bkg_EB_, meMuon_eta_tot_Bkg_EB_, meEtaEffnoMTD_Bkg_EB_);

  mePhiEffnoMTD_Bkg_EB_ = ibook.book1D("PhiEffnoMTD_Bkg_EB",
                                       " noMTD isolation Efficiency - Bkg Barrel VS pT;p_{T};Efficiency",
                                       meMuon_phi_tot_Bkg_EB_->getNbinsX(),
                                       meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_phi_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Bkg_EB_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_noMTD_Bkg_EB_, meMuon_phi_tot_Bkg_EB_, mePhiEffnoMTD_Bkg_EB_);

  if (optionalPlots_) {
    mePtEffgen_Bkg_EB_ = ibook.book1D("pTeffMTD_gen_Bkg_EB",
                                      " MTD isolation Efficiency - genInfo - Bkg Barrel VS pT;p_{T};Efficiency",
                                      meMuon_pt_sim_tot_Bkg_EB_->getNbinsX(),
                                      meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_sim_tot_Bkg_EB_->getTH1()->GetXaxis()->GetXmax());
    mePtEffgen_Bkg_EB_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_gen_Bkg_EB_, meMuon_pt_sim_tot_Bkg_EB_, mePtEffgen_Bkg_EB_);

    // Muon iso addition ends
    // For endcap now
    mePtEffMTD_1_Bkg_EE_ = ibook.book1D("pTeffMTD_1_Bkg_EE",
                                        " MTD isolation Efficiency - 1 BkgL Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_1_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_1_Bkg_EE_);

    mePtEffMTD_2_Bkg_EE_ = ibook.book1D("pTeffMTD_2_Bkg_EE",
                                        " MTD isolation Efficiency - 2 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_2_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_2_Bkg_EE_);

    mePtEffMTD_3_Bkg_EE_ = ibook.book1D("pTeffMTD_3_Bkg_EE",
                                        " MTD isolation Efficiency - 3 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_3_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_3_Bkg_EE_);

    mePtEffMTD_4_Bkg_EE_ = ibook.book1D("pTeffMTD_4_Bkg_EE",
                                        " MTD isolation Efficiency - 4 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_4_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_4_Bkg_EE_);

    mePtEffMTD_5_Bkg_EE_ = ibook.book1D("pTeffMTD_5_Bkg_EE",
                                        " MTD isolation Efficiency - 5 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_5_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_5_Bkg_EE_);

    mePtEffMTD_6_Bkg_EE_ = ibook.book1D("pTeffMTD_6_Bkg_EE",
                                        " MTD isolation Efficiency - 6 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_6_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_6_Bkg_EE_);

    mePtEffMTD_7_Bkg_EE_ = ibook.book1D("pTeffMTD_7_Bkg_EE",
                                        " MTD isolation Efficiency - 7 Bkg Endcap VS pT;p_{T};Efficiency",
                                        meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                        meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_MTD_7_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_7_Bkg_EE_);

    mePtEffMTD_sim_1_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_1_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 1 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_1_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_1_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_1_Bkg_EE_);

    mePtEffMTD_sim_2_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_2_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 2 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_2_Bkg_EE_);

    mePtEffMTD_sim_3_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_3_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 3 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_3_Bkg_EE_);

    mePtEffMTD_sim_4_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_4_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 4 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_4_Bkg_EE_);

    mePtEffMTD_sim_5_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_5_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 5 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_5_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_5_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_5_Bkg_EE_);

    mePtEffMTD_sim_6_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_6_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 6 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_6_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_6_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_6_Bkg_EE_);

    mePtEffMTD_sim_7_Bkg_EE_ = ibook.book1D("mePtEffMTD_sim_7_Bkg_EE_",
                                            " MTD isolation Efficiency SIM - 7 Bkg Endcap VS pT;p_{T};Efficiency",
                                            meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                            meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_7_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_7_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_7_Bkg_EE_);
  }
  mePtEffMTD_4sigma_Bkg_EE_ =
      ibook.book1D("pTeffMTD_4sigma_Bkg_EE",
                   " MTD isolation Efficiency - 4 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_4sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_4sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_4sigma_Bkg_EE_);

  mePtEffMTD_3sigma_Bkg_EE_ =
      ibook.book1D("pTeffMTD_3sigma_Bkg_EE",
                   " MTD isolation Efficiency - 3 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_3sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_3sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_3sigma_Bkg_EE_);

  mePtEffMTD_2sigma_Bkg_EE_ =
      ibook.book1D("pTeffMTD_2sigma_Bkg_EE",
                   " MTD isolation Efficiency - 2 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                   meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffMTD_2sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_MTD_2sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_2sigma_Bkg_EE_);

  if (optionalPlots_) {
    mePtEffMTD_sim_4sigma_Bkg_EE_ =
        ibook.book1D("pTeffMTD_sim_4sigma_Bkg_EE",
                     " MTD isolation Efficiency SIM - 4 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_4sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_4sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_4sigma_Bkg_EE_);

    mePtEffMTD_sim_3sigma_Bkg_EE_ =
        ibook.book1D("pTeffMTD_sim_3sigma_Bkg_EE",
                     " MTD isolation Efficiency SIM - 3 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_3sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_3sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_3sigma_Bkg_EE_);

    mePtEffMTD_sim_2sigma_Bkg_EE_ =
        ibook.book1D("pTeffMTD_sim_2sigma_Bkg_EE",
                     " MTD isolation Efficiency SIM - 2 sigma compatibility - Bkg Endcap VS pT;p_{T};Efficiency",
                     meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                     meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffMTD_sim_2sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_sim_MTD_2sigma_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffMTD_sim_2sigma_Bkg_EE_);

    meEtaEffMTD_1_Bkg_EE_ = ibook.book1D("EtaEffMTD_1_Bkg_EE",
                                         " MTD isolation Efficiency - 1 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_1_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_1_Bkg_EE_);

    meEtaEffMTD_2_Bkg_EE_ = ibook.book1D("EtaEffMTD_2_Bkg_EE",
                                         " MTD isolation Efficiency - 2 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_2_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_2_Bkg_EE_);

    meEtaEffMTD_3_Bkg_EE_ = ibook.book1D("EtaEffMTD_3_Bkg_EE",
                                         " MTD isolation Efficiency - 3 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_3_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_3_Bkg_EE_);

    meEtaEffMTD_4_Bkg_EE_ = ibook.book1D("EtaEffMTD_4_Bkg_EE",
                                         " MTD isolation Efficiency - 4 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_4_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_4_Bkg_EE_);

    meEtaEffMTD_5_Bkg_EE_ = ibook.book1D("EtaEffMTD_5_Bkg_EE",
                                         " MTD isolation Efficiency - 5 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_5_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_5_Bkg_EE_);

    meEtaEffMTD_6_Bkg_EE_ = ibook.book1D("EtaEffMTD_6_Bkg_EE",
                                         " MTD isolation Efficiency - 6 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_6_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_6_Bkg_EE_);

    meEtaEffMTD_7_Bkg_EE_ = ibook.book1D("EtaEffMTD_7_Bkg_EE",
                                         " MTD isolation Efficiency - 7 Bkg Endcap VS Eta;#eta;Efficiency",
                                         meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    meEtaEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_eta_MTD_7_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_7_Bkg_EE_);

    mePhiEffMTD_1_Bkg_EE_ = ibook.book1D("PhiEffMTD_1_Bkg_EE",
                                         " MTD isolation Efficiency - 1 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_1_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_1_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_1_Bkg_EE_);

    mePhiEffMTD_2_Bkg_EE_ = ibook.book1D("PhiEffMTD_2_Bkg_EE",
                                         " MTD isolation Efficiency - 2 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_2_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_2_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_2_Bkg_EE_);

    mePhiEffMTD_3_Bkg_EE_ = ibook.book1D("PhiEffMTD_3_Bkg_EE",
                                         " MTD isolation Efficiency - 3 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_3_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_3_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_3_Bkg_EE_);

    mePhiEffMTD_4_Bkg_EE_ = ibook.book1D("PhiEffMTD_4_Bkg_EE",
                                         " MTD isolation Efficiency - 4 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_4_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_4_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_4_Bkg_EE_);

    mePhiEffMTD_5_Bkg_EE_ = ibook.book1D("PhiEffMTD_5_Bkg_EE",
                                         " MTD isolation Efficiency - 5 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_5_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_5_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_5_Bkg_EE_);

    mePhiEffMTD_6_Bkg_EE_ = ibook.book1D("PhiEffMTD_6_Bkg_EE",
                                         " MTD isolation Efficiency - 6 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_6_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_6_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_6_Bkg_EE_);

    mePhiEffMTD_7_Bkg_EE_ = ibook.book1D("PhiEffMTD_7_Bkg_EE",
                                         " MTD isolation Efficiency - 7 Bkg Endcap VS Phi;#phi;Efficiency",
                                         meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                         meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePhiEffMTD_7_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_phi_MTD_7_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_7_Bkg_EE_);
  }

  meEtaEffMTD_4sigma_Bkg_EE_ =
      ibook.book1D("EtaEffMTD_4sigma_Bkg_EE",
                   " MTD isolation Efficiency - 4 sigma compatibiliy - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_4sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_4sigma_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_4sigma_Bkg_EE_);

  meEtaEffMTD_3sigma_Bkg_EE_ =
      ibook.book1D("EtaEffMTD_3sigma_Bkg_EE",
                   " MTD isolation Efficiency - 3 sigma compatibiliy - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_3sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_3sigma_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_3sigma_Bkg_EE_);

  meEtaEffMTD_2sigma_Bkg_EE_ =
      ibook.book1D("EtaEffMTD_2sigma_Bkg_EE",
                   " MTD isolation Efficiency - 2 sigma compatibiliy - Bkg Endcap VS Eta;#eta;Efficiency",
                   meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffMTD_2sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_MTD_2sigma_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffMTD_2sigma_Bkg_EE_);

  mePhiEffMTD_4sigma_Bkg_EE_ =
      ibook.book1D("PhiEffMTD_4sigma_Bkg_EE",
                   " MTD isolation Efficiency - 4 sigma compatibiliy - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_4sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_4sigma_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_4sigma_Bkg_EE_);

  mePhiEffMTD_3sigma_Bkg_EE_ =
      ibook.book1D("PhiEffMTD_3sigma_Bkg_EE",
                   " MTD isolation Efficiency - 3 sigma compatibiliy - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_3sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_3sigma_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_3sigma_Bkg_EE_);

  mePhiEffMTD_2sigma_Bkg_EE_ =
      ibook.book1D("PhiEffMTD_2sigma_Bkg_EE",
                   " MTD isolation Efficiency - 2 sigma compatibiliy - Bkg Endcap VS Phi;#phi;Efficiency",
                   meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                   meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffMTD_2sigma_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_MTD_2sigma_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffMTD_2sigma_Bkg_EE_);

  mePtEffnoMTD_Bkg_EE_ = ibook.book1D("pTeffnoMTD_Bkg_EE",
                                      " noMTD isolation Efficiency - Bkg Endcap VS pT;p_{T};Efficiency",
                                      meMuon_pt_tot_Bkg_EE_->getNbinsX(),
                                      meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePtEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_pt_noMTD_Bkg_EE_, meMuon_pt_tot_Bkg_EE_, mePtEffnoMTD_Bkg_EE_);

  meEtaEffnoMTD_Bkg_EE_ = ibook.book1D("EtaEffnoMTD_Bkg_EE",
                                       " noMTD isolation Efficiency - Bkg Endcap VS Eta;#eta;Efficiency",
                                       meMuon_eta_tot_Bkg_EE_->getNbinsX(),
                                       meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_eta_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  meEtaEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_eta_noMTD_Bkg_EE_, meMuon_eta_tot_Bkg_EE_, meEtaEffnoMTD_Bkg_EE_);

  mePhiEffnoMTD_Bkg_EE_ = ibook.book1D("PhiEffnoMTD_Bkg_EE",
                                       " noMTD isolation Efficiency - Bkg Endcap VS Phi;#phi;Efficiency",
                                       meMuon_phi_tot_Bkg_EE_->getNbinsX(),
                                       meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                       meMuon_phi_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
  mePhiEffnoMTD_Bkg_EE_->getTH1()->SetMinimum(0.);
  computeEfficiency1D(meMuon_phi_noMTD_Bkg_EE_, meMuon_phi_tot_Bkg_EE_, mePhiEffnoMTD_Bkg_EE_);

  if (optionalPlots_) {
    mePtEffgen_Bkg_EE_ = ibook.book1D("pTeffMTD_gen_Bkg_EE",
                                      " MTD isolation Efficiency - genInfo - Bkg Endcap VS pT;p_{T};Efficiency",
                                      meMuon_pt_sim_tot_Bkg_EE_->getNbinsX(),
                                      meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmin(),
                                      meMuon_pt_sim_tot_Bkg_EE_->getTH1()->GetXaxis()->GetXmax());
    mePtEffgen_Bkg_EE_->getTH1()->SetMinimum(0.);
    computeEfficiency1D(meMuon_pt_gen_Bkg_EE_, meMuon_pt_sim_tot_Bkg_EE_, mePtEffgen_Bkg_EE_);
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ----------
void MtdMuonIsoHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/MuonIso/");
  //desc.add<bool>("option_plots", false);
  desc.add<bool>("option_plots", true);

  descriptions.add("MtdMuonIsoPostProcessor", desc);
}

DEFINE_FWK_MODULE(MtdMuonIsoHarvester);
