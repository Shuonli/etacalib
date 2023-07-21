#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream>	//std::ostringstsream
#include <fstream>	//std::ifstream
#include <iostream> //std::cout, std::endl
// #include <math.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <vector>
#include <memory>
#include <TF1.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include "RollingBuffer.h"

RollingBuffer::RollingBuffer(int maxLength) : maxLength(maxLength) {}

void RollingBuffer::AddFront(const std::vector<TLorentzVector> &vec)
{
	buffer.insert(buffer.begin(), vec);
	if ((int)buffer.size() > maxLength)
	{
		PopBack();
	}
}

int RollingBuffer::GetSize() const
{
	return buffer.size();
}

void RollingBuffer::Flush()
{
	buffer.clear();
}

std::vector<TLorentzVector> RollingBuffer::GetElement(int index) const
{
	if (index >= 0 && index < (int)buffer.size())
	{
		return buffer[index];
	}
	return std::vector<TLorentzVector>(); // Return an empty vector if index is out of bounds
}

void RollingBuffer::PopBack()
{
	buffer.pop_back();
}

float DeltaR(TLorentzVector photon1, TLorentzVector photon2)
{
	float deta = photon1.PseudoRapidity() - photon2.PseudoRapidity();
	float dphi = abs(photon1.Phi() - photon2.Phi());

	if (dphi > M_PI)
		dphi = 2 * M_PI - dphi;
	float dr = sqrt(deta * deta + dphi * dphi);
	return dr;
}

bool CutPhoton(TLorentzVector photon1, TLorentzVector photon2)
{
	// lead photon energy cut
	float leadingphotoncut = 2;
	float asymcut = 0.5;
	// drcut
	// float dr = DeltaR(photon1, photon2);
	// if (dr < 0.1) return false;
	// if (dr > 0.8) return false;

	if (photon1.Energy() < leadingphotoncut && photon2.Energy() < leadingphotoncut)
		return false;
	// asym cut
	float asym = abs(photon1.Energy() - photon2.Energy()) / (photon1.Energy() + photon2.Energy());
	if (asym > asymcut)
		return false;
	// cut on the pt of the sum of the two photons
	float ptsum = (photon1 + photon2).Pt();
	if (ptsum < 4)
		return false;
	if (ptsum > 100)
		return false;
	return true;
}

int findbin(float value, float *binarray, int nbin)
{
	for (int i = 0; i < nbin; i++)
	{
		if (value >= binarray[i] && value < binarray[i + 1])
			return i;
	}
	return -1;
}

void MixEvent(RollingBuffer *rb, TH1F *hm[5], TRandom3 randomGenerator, float *binarray)
{

	if (rb->GetSize() < 2)
	{
		std::cout << "There must be at least 2 elements in the buffer." << std::endl;
		return;
	}
	std::vector<TLorentzVector> headerevent = rb->GetElement(0);
	for (int i = 1; i < rb->GetSize(); i++)
	{
		std::vector<TLorentzVector> poolevent = rb->GetElement(i);
		for (TLorentzVector headervalue : headerevent)
		{
			if ((int)poolevent.size() < 1)
				continue;

			int randomIndex = randomGenerator.Integer(poolevent.size());
			// for (TLorentzVector poolvalue : poolevent) {
			TLorentzVector poolvalue = poolevent.at(randomIndex);
			if (!CutPhoton(headervalue, poolvalue))
				continue;
			TLorentzVector mixedeta = headervalue + poolvalue;

			int ibin = findbin(mixedeta.Pt(), binarray, 5);
			if (ibin >= 0)
				hm[ibin]->Fill(mixedeta.M());
			//}
		}
	}
}

void simpleetamassreco()
{
	TRandom3 randomGenerator;
	randomGenerator.SetSeed(std::time(0));

	const int ncentbin = 10;
	const int nvtxbin = 6;
	const int depth = 50;
	float ptbin[6] = {4, 5, 6, 7, 8, 100};
	RollingBuffer *rblist[ncentbin][nvtxbin];

	for (int i = 0; i < ncentbin; i++)
	{
		for (int j = 0; j < nvtxbin; j++)
		{
			rblist[i][j] = new RollingBuffer(depth);
		}
	}
	TFile *out = new TFile("etamassrecoout.root", "RECREATE");
	TH1F *hclusterE = new TH1F("hclusterE", "", 100, 0, 30);
	TH1F *hclusterN = new TH1F("hclusterN", "", 20, 0, 20);
	TH1F *hetapt = new TH1F("hetapt", "", 20, 0, 20);
	TH2F *hclusterEN = new TH2F("hclusterEN", "", 30, 0, 40, 100, 0, 40);
	TH1F *hetamass = new TH1F("hetamass", "", 200, 0, 2);
	// make histograms for eta at different pt ranges
	TH1F *hetamassrealpt[5];
	hetamassrealpt[0] = new TH1F("hetamasspt45", "", 200, 0, 2);
	hetamassrealpt[1] = new TH1F("hetamasspt56", "", 200, 0, 2);
	hetamassrealpt[2] = new TH1F("hetamasspt67", "", 200, 0, 2);
	hetamassrealpt[3] = new TH1F("hetamasspt78", "", 200, 0, 2);
	hetamassrealpt[4] = new TH1F("hetamasspt8u", "", 200, 0, 2);
	TH1F *hetamassmixpt[5];
	hetamassmixpt[0] = new TH1F("hetamassptmix45", "", 200, 0, 2);
	hetamassmixpt[1] = new TH1F("hetamassptmix56", "", 200, 0, 2);
	hetamassmixpt[2] = new TH1F("hetamassptmix67", "", 200, 0, 2);
	hetamassmixpt[3] = new TH1F("hetamassptmix78", "", 200, 0, 2);
	hetamassmixpt[4] = new TH1F("hetamassptmix8u", "", 200, 0, 2);

	TH1F *hphotonetamass = new TH1F("hphotonetamass", "", 100, 0, 2);
	TH1F *hetamassmix = new TH1F("hetamassmix", "", 200, 0, 2);
	TH2F *hetamasspt = new TH2F("hetamasspt", "", 100, 0, 40, 100, 0, 2);
	TH2F *hdrmass = new TH2F("hdrmass", "", 100, 0, 6, 100, 0, 2);
	TH2F *hdrtruth = new TH2F("hdrtruth", "", 100, 0, 3, 100, 0, 3);
	TH1F *hngoodcluster = new TH1F("hngoodcluster", "", 100, 0, 100);
	TH2F *hmassncluster = new TH2F("hmassncluster", "", 10, 0, 10, 100, 0, 2);
	TH2F *hmassntower = new TH2F("hmassntower", "", 30, 0, 30, 100, 0, 2);
	TH1F *hdr = new TH1F("hdr", "", 100, 0, 6);
	TBranch *b_eta_center = nullptr;
	TBranch *b_phi_center = nullptr;
	TBranch *b_tower_energy = nullptr;
	TBranch *b_cluster_pt = nullptr;
	TBranch *b_cluster_eta = nullptr;
	TBranch *b_cluster_phi = nullptr;
	TBranch *b_cluster_e = nullptr;
	TBranch *b_cluster_chi2 = nullptr;
	TBranch *b_cluster_prob = nullptr;
	TBranch *b_cluster_nTowers = nullptr;
	TBranch *b_cent = nullptr;
	TBranch *b_vtxz = nullptr;
	TBranch *b_asym = nullptr;
	TBranch *b_deltaR = nullptr;
	TBranch *b_lead_E = nullptr;
	TBranch *b_sublead_E = nullptr;
	TBranch *b_lead_phi = nullptr;
	TBranch *b_lead_eta = nullptr;
	TBranch *b_sublead_phi = nullptr;
	TBranch *b_sublead_eta = nullptr;
	TBranch *b_pi0_E = nullptr;
	TBranch *b_pi0_eta = nullptr;
	TBranch *b_pi0_phi = nullptr;
	TBranch *b_pi0_pt = nullptr;
	TBranch *b_pi0_M = nullptr;
	TBranch *b_tr_p = nullptr;
	TBranch *b_tr_pt = nullptr;
	TBranch *b_tr_eta = nullptr;
	TBranch *b_tr_phi = nullptr;
	TBranch *b_tr_charge = nullptr;
	TBranch *b_tr_chisq = nullptr;
	TBranch *b_tr_ndf = nullptr;
	TBranch *b_tr_cemc_eta = nullptr;
	TBranch *b_tr_cemc_phi = nullptr;

	std::vector<float> *m_eta_center = nullptr;
	std::vector<float> *m_phi_center = nullptr;
	std::vector<float> *m_tower_energy = nullptr;
	std::vector<float> *m_cluster_pt = nullptr;
	std::vector<float> *m_cluster_eta = nullptr;
	std::vector<float> *m_cluster_phi = nullptr;
	std::vector<float> *m_cluster_e = nullptr;
	std::vector<float> *m_cluster_chi2 = nullptr;
	std::vector<float> *m_cluster_prob = nullptr;
	std::vector<float> *m_cluster_nTowers = nullptr;
	std::vector<float> *m_cent = nullptr;
	std::vector<float> *m_vtxz = nullptr;

	std::vector<float> *m_asym = nullptr;
	std::vector<float> *m_deltaR = nullptr;
	std::vector<float> *m_lead_E = nullptr;
	std::vector<float> *m_sublead_E = nullptr;
	std::vector<float> *m_lead_phi = nullptr;
	std::vector<float> *m_lead_eta = nullptr;
	std::vector<float> *m_sublead_phi = nullptr;
	std::vector<float> *m_sublead_eta = nullptr;

	std::vector<float> *m_pi0_E = nullptr;
	std::vector<float> *m_pi0_eta = nullptr;
	std::vector<float> *m_pi0_phi = nullptr;
	std::vector<float> *m_pi0_pt = nullptr;
	std::vector<float> *m_pi0_M = nullptr;

	std::vector<float> *m_tr_p = nullptr;
	std::vector<float> *m_tr_pt = nullptr;
	std::vector<float> *m_tr_eta = nullptr;
	std::vector<float> *m_tr_phi = nullptr;
	std::vector<float> *m_tr_charge = nullptr;
	std::vector<float> *m_tr_chisq = nullptr;
	std::vector<float> *m_tr_ndf = nullptr;
	std::vector<float> *m_tr_cemc_eta = nullptr;
	std::vector<float> *m_tr_cemc_phi = nullptr;

	TChain *t0 = new TChain("TreeTruthEta");
	// t0->Add("/sphenix/user/shuhangli/etacalib/analysis/pi0ClusterAna/rootFiles/eta_simple_0/all.root");
	t0->Add("rootfile/MBhijing/all*.root");
	TChain *t1 = new TChain("TreeTruthPhoton");
	// t1->Add("/sphenix/user/shuhangli/etacalib/analysis/pi0ClusterAna/rootFiles/eta_simple_0/all.root");
	t1->Add("rootfile/MBhijing/all*.root");
	TChain *t2 = new TChain("TreeClusterTower");
	// t2->Add("/sphenix/user/shuhangli/etacalib/analysis/pi0ClusterAna/rootFiles/eta_simple_0/all.root");
	t2->Add("rootfile/MBhijing/all*.root");
	TChain *t3 = new TChain("TreeTracks");
	t3->Add("rootfile/MBhijing/all*.root");

	t2->SetBranchAddress("towerEtaCEnter", &m_eta_center, &b_eta_center);
	t2->SetBranchAddress("towerPhiCenter", &m_phi_center, &b_phi_center);
	t2->SetBranchAddress("towerEnergy", &m_tower_energy, &b_tower_energy);
	t2->SetBranchAddress("clusterpt", &m_cluster_pt, &b_cluster_pt);
	t2->SetBranchAddress("clusterEta", &m_cluster_eta, &b_cluster_eta);
	t2->SetBranchAddress("clusterPhi", &m_cluster_phi, &b_cluster_phi);
	t2->SetBranchAddress("clusterE", &m_cluster_e, &b_cluster_e);
	t2->SetBranchAddress("clusterChi2", &m_cluster_chi2, &b_cluster_chi2);
	t2->SetBranchAddress("clusterProb", &m_cluster_prob, &b_cluster_prob);
	// t2->SetBranchAddress("clusterNTowers", &m_cluster_nTowers, &b_cluster_nTowers);
	t2->SetBranchAddress("cent", &m_cent, &b_cent);
	t2->SetBranchAddress("vtxz", &m_vtxz, &b_vtxz);

	t1->SetBranchAddress("pairAsym", &m_asym, &b_asym);
	t1->SetBranchAddress("pairDeltaR", &m_deltaR, &b_deltaR);
	t1->SetBranchAddress("leadPhotonPhi", &m_lead_E, &b_lead_E);
	t1->SetBranchAddress("leadPhotonEta", &m_sublead_E, &b_sublead_E);
	t1->SetBranchAddress("subleadPhotonPhi", &m_lead_phi, &b_lead_phi);
	t1->SetBranchAddress("subleadPhotonEta", &m_lead_eta, &b_lead_eta);
	t1->SetBranchAddress("leadPhotonE", &m_sublead_phi, &b_sublead_phi);
	t1->SetBranchAddress("subLeadPhotonE", &m_sublead_eta, &b_sublead_eta);

	t0->SetBranchAddress("EtaE", &m_pi0_E, &b_pi0_E);
	t0->SetBranchAddress("Eta_eta", &m_pi0_eta, &b_pi0_eta);
	t0->SetBranchAddress("Eta_phi", &m_pi0_phi, &b_pi0_phi);
	t0->SetBranchAddress("Eta_pt", &m_pi0_pt, &b_pi0_pt);
	t0->SetBranchAddress("Eta_M", &m_pi0_M, &b_pi0_M);

	t3->SetBranchAddress("tr_p", &m_tr_p, &b_tr_p);
	t3->SetBranchAddress("tr_pt", &m_tr_pt, &b_tr_pt);
	t3->SetBranchAddress("tr_eta", &m_tr_eta, &b_tr_eta);
	t3->SetBranchAddress("tr_phi", &m_tr_phi, &b_tr_phi);
	t3->SetBranchAddress("tr_charge", &m_tr_charge, &b_tr_charge);
	t3->SetBranchAddress("tr_chisq", &m_tr_chisq, &b_tr_chisq);
	t3->SetBranchAddress("tr_ndf", &m_tr_ndf, &b_tr_ndf);
	t3->SetBranchAddress("tr_cemc_eta", &m_tr_cemc_eta, &b_tr_cemc_eta);
	t3->SetBranchAddress("tr_cemc_phi", &m_tr_cemc_phi, &b_tr_cemc_phi);

	int ne = t2->GetEntries();
	std::cout << "total events: " << ne << std::endl;
	if (!(t0->GetEntries() == t1->GetEntries()) && (t0->GetEntries() == t2->GetEntries()))
		std::cout << "event numbers do not match" << std::endl;

	for (int e = 0; e < ne; e++)
	{
		// for(int e=0;e<20000;e++){
		if (e % 10000 == 0)
			std::cout << "number " << e << "/" << ne << std::endl;
		//	std::cout << e << std::endl;

		t2->GetEntry(e);
		int centbin = (int)m_cent->at(0) / (100 / (ncentbin - 1));
		int vtxbin = (int)(m_vtxz->at(0) + 10) / (20 / (nvtxbin - 1));
		if (vtxbin < 0)
			continue;
		// std::cout<< m_vtxz->at(0)<<" vtxbin: "<<vtxbin<<std::endl;
		// std::cout<< m_cent->at(0)<<std::endl;
		// if (m_cent->at(0) < 60.)
		// continue;
		// if (m_cent->at(0) > 1000.) continue;
		t0->GetEntry(e);
		// t1->GetEntry(e);
		t3->GetEntry(e);
		TLorentzVector photon4;

		std::vector<TLorentzVector> goodphoton;
		for (int i = 0; i < (int)m_pi0_E->size(); i++)
		{

			if (m_pi0_E->at(i) < 1)
				continue;

			photon4.SetPtEtaPhiE(m_pi0_pt->at(i), m_pi0_eta->at(i), m_pi0_phi->at(i), m_pi0_E->at(i));
			goodphoton.push_back(photon4);

			// std::cout<<photon4.M()<<" "<<m_pi0_eta->at(i)<<" "<<photon4.Phi()<<" "<<photon4.E()<<std::endl;
		}

		// save track info to veto on the cluster

		std::vector<float> trackemphi;
		std::vector<float> trackemeta;
		for (int i = 0; i < (int)m_tr_p->size(); i++)
		{
			if (m_tr_p->at(i) < 1)
				continue;
			if (abs(m_tr_cemc_eta->at(i)) > 1.1)
				continue;
			if ((m_tr_chisq->at(i) / m_tr_ndf->at(i)) > 10)
				continue;
			trackemeta.push_back(m_tr_cemc_eta->at(i));
			trackemphi.push_back(m_tr_cemc_phi->at(i));
		}

		TLorentzVector cluster4;

		std::vector<TLorentzVector> goodcluster;
		// std::vector<float> clusterntower;
		// selected cluster that matched with primary photon here just to make sure
		for (int i = 0; i < (int)m_cluster_e->size(); i++)
		{
			// std::cout<<i<<"/"<<(int)m_cluster_e->size();
			float Cpt = m_cluster_pt->at(i);
			if (Cpt < 3)
				continue;
			// if (Cpt > 8)
			//	continue;
			float chi2 = m_cluster_chi2->at(i);
			if (chi2 > 1)
				continue;
			float Ceta = m_cluster_eta->at(i);
			if (abs(Ceta) > 0.9)
				continue;
			// if (abs(m_cluster_eta->at(i)) > 0.3) continue;
			float Cphi = m_cluster_phi->at(i);
			// veto cut on tracks

			bool closetrack = 0;
			for (int j = 0; j < (int)trackemphi.size(); j++)
			{
				float deta = abs(trackemeta.at(j) - Ceta);
				float dphi = abs(trackemphi.at(j) - Cphi);

				if (dphi > M_PI)
					dphi = 2 * M_PI - dphi;
				float dr = sqrt(deta * deta + dphi * dphi);

				if (dr < 0.2)
				{
					closetrack = 1;
					break;
				}
			}
			if (closetrack)
				continue;

			cluster4.SetPtEtaPhiE(m_cluster_pt->at(i), m_cluster_eta->at(i), m_cluster_phi->at(i), m_cluster_e->at(i));
			/*
			bool photonmatch = 0;
			for (int j = 0; j < (int)goodphoton.size(); j++) {
				TLorentzVector photon1 = goodphoton.at(j);
				if (DeltaR(photon1, cluster4) < 0.1) {
					photonmatch = true;
					break;
				}
			}
			if (!photonmatch) continue;
			*/
			goodcluster.push_back(cluster4);
			// clusterntower.push_back(m_cluster_nTowers->at(i));
		}
		// std::cout<<(int)goodcluster.size()<<std::endl;
		hngoodcluster->Fill((int)goodcluster.size());

		rblist[centbin][vtxbin]->AddFront(goodcluster);

		// if((int)goodcluster.size()<=2) continue;

		// rblist[centbin][vtxbin]->AddFront(goodphoton);
		MixEvent(rblist[centbin][vtxbin], hetamassmixpt, randomGenerator, ptbin);

		for (int i = 0; i < (int)goodphoton.size(); i++)
		{
			hclusterE->Fill(goodphoton.at(i).Energy());

			for (int j = 0; j < (int)goodphoton.size(); j++)
			{
				if (j <= i)
					continue;
				TLorentzVector photon1 = goodphoton.at(i);
				TLorentzVector photon2 = goodphoton.at(j);
				if (!CutPhoton(photon1, photon2))
					continue;
				TLorentzVector eta = photon1 + photon2;
				if (eta.Pt() < 3)
					continue;
				hphotonetamass->Fill(eta.M());
			}
		}

		for (int i = 0; i < (int)goodcluster.size(); i++)
		{
			// std::cout<<i<<std::endl;

			hclusterE->Fill(goodcluster.at(i).Energy());
			for (int j = 0; j < (int)goodcluster.size(); j++)
			{
				if (j <= i)
					continue;

				TLorentzVector photon1 = goodcluster.at(i);
				TLorentzVector photon2 = goodcluster.at(j);
				if (!CutPhoton(photon1, photon2))
					continue;
				TLorentzVector eta = photon1 + photon2;

				float dr = DeltaR(photon1, photon2);
				// if(dr>0.8)continue;
				// if(dr<0.025)continue;
				// if((photon1.PseudoRapidity() / photon2.PseudoRapidity())<0) continue;

				// hdrtruth->Fill(m_deltaR->at(0), dr);
				// hdrmass->Fill(dr, eta.M());
				int ibin = findbin(eta.Pt(), ptbin, 5);
				if (ibin >= 0)
					hetamassrealpt[ibin]->Fill(eta.M());
				hetamasspt->Fill(eta.Pt(), eta.M());
				hmassncluster->Fill((int)goodcluster.size(), eta.M());
				// hmassntower->Fill(clusterntower.at(i)+clusterntower.at(j),eta.M());
				hetapt->Fill(eta.Pt());
				hdr->Fill(dr);
			}
		}
	}
	out->Write();
	out->Close();
}
