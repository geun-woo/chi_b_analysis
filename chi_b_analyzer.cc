// C++
#include <iostream>
#include <cmath>
// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"



using namespace std;


double getDPHI( double phi1, double phi2);

void chi_b_analyzer() {
    TFile *fIn = new TFile("/home/gekim/data/oniaTree_Pbp+pPb_merged.root", "read");
    // TFile *fOut = new TFile("output.root", "recreate");

    TTree *oldTree = (TTree *)fIn->Get("hionia/myTree");
    // TTree *newTree = new TTree("newTree", "");
    
    // TClonesArrays
    TClonesArray* Reco_conv_4mom = nullptr;
    TClonesArray* Reco_conv_vertex = nullptr;
    TClonesArray* Reco_QQ_4mom = nullptr;
    TClonesArray* Reco_QQ_vertex = nullptr;
    TClonesArray* Reco_QQ_mupl_4mom = nullptr;
    TClonesArray* Reco_QQ_mumi_4mom = nullptr;
    
    // Conversion
    int Reco_conv_size;
    double Reco_conv_pt[500];
    double Reco_conv_eta[500];
    double Reco_conv_phi[500];
    // Conversion track
    double Reco_conv_chi2_probability[500];
    double Reco_conv_dz_innerPosition[500];
    double Reco_conv_distOfMinimumApproach[500];
    double Reco_conv_zOfPrimaryVertexFromTracks[500];
    // Conversion track 1
    int Reco_conv_trk1_charge[500];
    int Reco_conv_trk1_nhit[500];
    double Reco_conv_trk1_normalizedChi2[500];
    double Reco_conv_trk1_theta[500];
    double Reco_conv_trk1_phi[500];
    double Reco_conv_trk1_d0[500];
    // Conversion track 2
    int Reco_conv_trk2_charge[500];
    int Reco_conv_trk2_nhit[500];
    double Reco_conv_trk2_normalizedChi2[500];
    double Reco_conv_trk2_theta[500];
    double Reco_conv_trk2_phi[500];
    double Reco_conv_trk2_d0[500];
    
    // Dimuon
    int Reco_QQ_size;
    float Reco_QQ_VtxProb[500];
    // Anti-muon
    int Reco_QQ_mupl_nTrkWMea[500];
    int Reco_QQ_mupl_nPixWMea[500];
    bool Reco_QQ_mupl_TMOneStaTight[500];
    float Reco_QQ_mupl_normChi2_inner[500];
    float Reco_QQ_mupl_dxy[500];
    float Reco_QQ_mupl_dz[500];
    // Muon
    int Reco_QQ_mumi_nPixWMea[500];
    int Reco_QQ_mumi_nTrkWMea[500];
    bool Reco_QQ_mumi_TMOneStaTight[500];
    float Reco_QQ_mumi_normChi2_inner[500];
    float Reco_QQ_mumi_dxy[500];
    float Reco_QQ_mumi_dz[500];


    // Conversion
    oldTree->SetBranchAddress("Reco_conv_size", &Reco_conv_size);
    oldTree->SetBranchAddress("Reco_conv_4mom", &Reco_conv_4mom);
    oldTree->SetBranchAddress("Reco_conv_vertex", &Reco_conv_vertex);
    //
    oldTree->SetBranchAddress("Reco_conv_pt", Reco_conv_pt);
    oldTree->SetBranchAddress("Reco_conv_eta", Reco_conv_eta);
    oldTree->SetBranchAddress("Reco_conv_phi", Reco_conv_phi);
    //
    oldTree->SetBranchAddress("Reco_conv_chi2_probability", Reco_conv_chi2_probability);
    oldTree->SetBranchAddress("Reco_conv_dz_innerPosition", Reco_conv_dz_innerPosition);
    oldTree->SetBranchAddress("Reco_conv_distOfMinimumApproach", Reco_conv_distOfMinimumApproach);
    oldTree->SetBranchAddress("Reco_conv_zOfPrimaryVertexFromTracks", Reco_conv_zOfPrimaryVertexFromTracks);
    // Conversion track 1
    oldTree->SetBranchAddress("Reco_conv_trk1_charge", Reco_conv_trk1_charge);
    oldTree->SetBranchAddress("Reco_conv_trk1_nhit", Reco_conv_trk1_nhit);
    oldTree->SetBranchAddress("Reco_conv_trk1_normalizedChi2", Reco_conv_trk1_normalizedChi2);
    oldTree->SetBranchAddress("Reco_conv_trk1_theta", Reco_conv_trk1_theta);
    oldTree->SetBranchAddress("Reco_conv_trk1_phi", Reco_conv_trk1_phi);
    oldTree->SetBranchAddress("Reco_conv_trk1_d0", Reco_conv_trk1_d0);
    // Conversion track 2
    oldTree->SetBranchAddress("Reco_conv_trk2_charge", Reco_conv_trk2_charge);
    oldTree->SetBranchAddress("Reco_conv_trk2_nhit", Reco_conv_trk2_nhit);
    oldTree->SetBranchAddress("Reco_conv_trk2_normalizedChi2", Reco_conv_trk2_normalizedChi2);
    oldTree->SetBranchAddress("Reco_conv_trk2_theta", Reco_conv_trk2_theta);
    oldTree->SetBranchAddress("Reco_conv_trk2_phi", Reco_conv_trk2_phi);
    oldTree->SetBranchAddress("Reco_conv_trk2_d0", Reco_conv_trk2_d0);
    
    // Dimuon
    oldTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
    oldTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
    oldTree->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vertex);
    oldTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
    oldTree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom);
    oldTree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom);
    // Anti-muon
    oldTree->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner);
    oldTree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea);
    oldTree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea);
    oldTree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy);
    oldTree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz);
    oldTree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight);
    // Muon
    oldTree->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner);
    oldTree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea);
    oldTree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea);
    oldTree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy);
    oldTree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz);
    oldTree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight);
    
    TH1D *h_pt_org = new TH1D("h_pt_org", "Reco_conv_Pt;Pt;Events", 10, 0, 20);
    TH1D *h_pt = new TH1D("h_pt", ";Pt;Events", 10, 0, 20);
    TH1D *h_eta_org = new TH1D("h_eta_org", "Reco_conv_Eta;Eta;Events", 100, -3., +3.);
    TH1D *h_eta = new TH1D("h_eta", ";Eta;Events", 100, -3., +3.);
    TH1D *h_phi_org = new TH1D("h_phi_org", "Reco_conv_phi;Phi;Events", 20, -3.15, +3.15);
    TH1D *h_phi = new TH1D("h_phi", ";Phi;Events", 20, -3.15, +3.15);
    TH1F *h_chi_b_mass = new TH1F("h_chi_b_mass", "Chi_b mass;GeV/c^2;Candidates / (25 MeV/c^2)", 100, 8.5, 11);
    // TH1F *h_dimuon_mass = new TH1F("h_dimuon_b_mass", "Chi_b mass;GeV/c^2;Candidates / (25 MeV/c^2)", 50, 8.5, 11);

    TLorentzVector *conv_4mom = new TLorentzVector;
    TLorentzVector *mupl_4mom = new TLorentzVector;
    TLorentzVector *mumi_4mom = new TLorentzVector;
    TLorentzVector *dimu_4mom = new TLorentzVector;
    TLorentzVector *chi_b_4mom = new TLorentzVector;

    TVector3 *conv_vtx = new TVector3;
    TVector3 *dimu_vtx = new TVector3;

    int nEvents = oldTree->GetEntries();
    for(int iEvent = 0; iEvent < nEvents; ++iEvent) {
        cout << "Event " << iEvent + 1 << "/" << nEvents << endl;
        oldTree->GetEntry(iEvent);
        for(int iConv = 0; iConv < Reco_conv_size; ++iConv) {
            if( Reco_conv_trk1_charge[iConv] + Reco_conv_trk2_charge[iConv] != 0 ) continue;
            // Electron track hits >= 4, 3
            if( ( Reco_conv_trk1_nhit[iConv] <= 3 ) && ( Reco_conv_trk2_nhit[iConv] <= 3 ) ) continue;
            // Electron track fit chi squared / ndof < 10
            if( ( Reco_conv_trk1_normalizedChi2[iConv] >= 10 ) || ( Reco_conv_trk2_normalizedChi2[iConv] >= 10 ) ) continue;
            // Distance of approach : -0.25 cm < d_m < 1 cm
            if( ( Reco_conv_distOfMinimumApproach[iConv] <= -.25 ) || ( Reco_conv_distOfMinimumApproach[iConv] >= 1. ) ) continue;
            // Signed impact parameter : q * d_0 > 0
            if( ( Reco_conv_trk1_d0[iConv] * Reco_conv_trk1_charge[iConv] <= 0 ) || ( Reco_conv_trk2_d0[iConv] * Reco_conv_trk2_charge[iConv] <= 0 ) ) continue;
            // z distance of inner hits < 5 cm
            if( fabs( Reco_conv_dz_innerPosition[iConv] ) > 5) continue;
            // positron - electron vertex fit probablilty > 5 * 10^(-4)
            if( Reco_conv_chi2_probability[iConv] < 0.0005 ) continue;
            // Radius of conversion > 1.5 cm
            conv_vtx = (TVector3*)Reco_conv_vertex->At(iConv);
            if( TMath::Sqrt( conv_vtx->x()*conv_vtx->x() + conv_vtx->y()*conv_vtx->y() ) <= 1.5) continue;
            // delta phi (positron-electron) < 0.2
            if( fabs( getDPHI( Reco_conv_trk1_phi[iConv], Reco_conv_trk2_phi[iConv] ) ) >= 0.2 ) continue;
            // delta cot(theta) (positron-electron) < 0.1
            if( fabs( 1./tan( Reco_conv_trk1_theta[iConv] ) - 1./tan( Reco_conv_trk2_theta[iConv] ) ) >= 0.1 ) continue;
            // | Eta(gamma) | < 1.0
            if( fabs( Reco_conv_eta[iConv] ) >= 1. ) continue;
            //
            for(int iMuMu = 0; iMuMu < Reco_QQ_size; ++iMuMu) {
                // track fit chi squared / ndof < 1.8
                if( Reco_QQ_mupl_normChi2_inner[iMuMu] >= 1.8) continue;
                if( Reco_QQ_mumi_normChi2_inner[iMuMu] >= 1.8) continue;
                // hits in Pixel > 1
                if( Reco_QQ_mupl_nPixWMea[iMuMu] <= 1) continue;
                if( Reco_QQ_mumi_nPixWMea[iMuMu] <= 1) continue;
                // hits in Tracker > 5
                if( Reco_QQ_mupl_nTrkWMea[iMuMu] <= 5) continue;
                if( Reco_QQ_mumi_nTrkWMea[iMuMu] <= 5) continue;
                // Fiducial cylinder : 3 cm (r) x 30 cm (z)
                if( Reco_QQ_mupl_dxy[iMuMu] >= 3) continue;
                if( Reco_QQ_mumi_dxy[iMuMu] >= 3) continue;
                if( Reco_QQ_mumi_dz[iMuMu] >= 30) continue;
                if( Reco_QQ_mumi_dz[iMuMu] >= 30) continue;
                // mu-mu vertex fit probalbility > 0.01
                if( Reco_QQ_VtxProb[iMuMu] <= 0.01 ) continue;
                // Muon ID = "TMuonOneStationTight"
                if( !(Reco_QQ_mupl_TMOneStaTight[iMuMu]) ) continue;
                if( !(Reco_QQ_mupl_TMOneStaTight[iMuMu]) ) continue;
                // Pt(mu) > 3.5 GeV
                mupl_4mom = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iMuMu);
                mumi_4mom = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iMuMu);
                if( mupl_4mom->Pt() <= 3.5 ) continue;
                if( mumi_4mom->Pt() <= 3.5 ) continue;
                // | Eta(mu) | < 1.9
                if( fabs( mupl_4mom->Eta() ) >= 1.9 ) continue;
                if( fabs( mumi_4mom->Eta() ) >= 1.9 ) continue;
                conv_4mom = (TLorentzVector*)Reco_conv_4mom->At(iConv);
                dimu_4mom = (TLorentzVector*)Reco_QQ_4mom->At(iMuMu);
                // mass_mu-mu : 8.5 - 11 GeV
                if( dimu_4mom->M() < 8.5 || dimu_4mom->M() > 11) continue;
                // | y(mu-mu) | < 1.5
                if( fabs( dimu_4mom->Rapidity() ) >= 1.5 ) continue;
                // mass(mu-mu-gamma) - mass(mu-mu) < 2 GeV
                *chi_b_4mom = *conv_4mom + *dimu_4mom;
                if( chi_b_4mom->M() - dimu_4mom->M() >= 2) continue;
                // Fill chi_b
                h_chi_b_mass->Fill( chi_b_4mom->M() );
            }
        }
    }
    // h_pt_org->SetStats(0);
    // h_pt->SetStats(0);
    // h_pt->SetLineColor(kRed);
    // tree->Draw("Reco_conv_pt>>h_pt_org");
    // h_pt->Draw("same");

    // h_eta_org->SetStats(0);
    // h_eta->SetStats(0);
    // h_eta->SetLineColor(kRed);
    // tree->Draw("Reco_conv_eta>>h_eta_org");
    // h_eta->Draw("same");

    // h_phi_org->SetStats(0);
    // h_phi->SetStats(0);
    // h_phi->SetLineColor(kRed);
    // tree->Draw("Reco_conv_phi>>h_phi_org");
    // h_phi->Draw("same");

}

double getDPHI( double phi1, double phi2) {
    double dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ){
        dphi = dphi - 2. * TMath::Pi();
    }
    if ( dphi <= -TMath::Pi() ){
        dphi = dphi + 2. * TMath::Pi();
    }
    if ( fabs(dphi) > TMath::Pi() ) {
        cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
    }

    return dphi;
}