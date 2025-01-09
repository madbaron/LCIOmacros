from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree, TVector3
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_kaons.root',
                  type=str, default='ntup_kaons.root')
(options, args) = parser.parse_args()

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000., 2500., 5000.))

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_matchedEMpfo_E = TH1D('matchedEMpfo_E', 'matchedEMpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_theta = TH1D('matchedEMpfo_theta', 'matchedEMpfo_theta',
                            len(arrBins_theta)-1, arrBins_theta)
h_matchedHADpfo_E = TH1D('matchedHADpfo_E', 'matchedHADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedHADpfo_theta = TH1D('matchedHADpfo_theta', 'matchedHADpfo_theta',
                             len(arrBins_theta)-1, arrBins_theta)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)

h_deltaEM_E = TH1D('deltaEM_E', 'deltaEM_E', 250, -1000, 1000)
h_deltaHAD_E = TH1D('deltaHAD_E', 'deltaHAD_E', 250, -1000, 1000)

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta, 
               h_matched_theta, h_matched_E,
               h_matchedEMpfo_E, h_matchedEMpfo_theta,
               h_matchedHADpfo_E, h_matchedHADpfo_theta,
               h_deltaEM_E, h_deltaHAD_E,
               h_Npfo, h_pfo_type]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
kaon_tree = TTree("kaon_tree", "kaon_tree")
E = array('d', [0])
theta = array('d', [0])
pid = array('i', [0])
E_jet = array('d', [0])
theta_jet = array('d', [0])
E_truth = array('d', [0])
theta_truth = array('d', [0])
sumECAL = array('d', [0])
sumHCAL = array('d', [0])
kaon_tree.Branch("E",  E,  'var/D')
kaon_tree.Branch("theta", theta, 'var/D')
kaon_tree.Branch("pid", pid, 'var/I')
kaon_tree.Branch("E_jet",  E_jet,  'var/D')
kaon_tree.Branch("theta_jet", theta_jet, 'var/D')
kaon_tree.Branch("E_truth",  E_truth,  'var/D')
kaon_tree.Branch("theta_truth", theta_truth, 'var/D')
kaon_tree.Branch("sumECAL", sumECAL, 'var/D')
kaon_tree.Branch("sumHCAL", sumHCAL, 'var/D')

to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # loop over all events in the file
    for ievt, event in enumerate(reader):

        if ievt % 1 == 0:
            print(" ")
            print("Processing event " + str(ievt))

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        h_truth_E.Fill(mcpCollection[0].getEnergy())
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        E_truth[0] = mcpCollection[0].getEnergy()
        theta_truth[0] = tlv.Theta()

        print("True energy", E_truth[0])

        '''
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus() == 1 and len(mcp.getParents()) == 0:
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                h_truth_E.Fill(mcp.getEnergy())
                h_truth_theta.Fill(tlv.Theta())
        '''
        
        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')
        h_Npfo.Fill(len(pfoCollection))

        # Match true pfo with closest reco PFO in deltaR
        matchedEM_E = -1.
        matchedEM_theta = -1.
        matchedHAD_E = -1.
        matchedHAD_theta = -1.

        minDREM = 999999.
        minDRHAD = 999999.

        E[0] = -1
        theta[0] = -1

        which_pfo = []

        for pfo in pfoCollection:
            h_pfo_type.Fill(abs(pfo.getType()))
            dp3 = pfo.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)
            which_pfo.append(pfo.getType())

            if pfo.getEnergy()/mcpCollection[0].getEnergy() > 0.: #change threshold here if needed
                if (dR < minDREM) and abs(pfo.getType()) == 22:
                    minDREM = dR
                    matchedEM_E = pfo.getEnergy()
                    matchedEM_theta = tlv_pfo.Theta()
                if (dR < minDRHAD) and abs(pfo.getType()) != 22:
                    minDRHAD = dR
                    matchedHAD_E = pfo.getEnergy()
                    matchedHAD_theta = tlv_pfo.Theta()
                    E[0] = pfo.getEnergy()
                    theta[0] = tlv_pfo.Theta()

        pid[0] = 22
        if minDRHAD < minDREM:
            pid[0] = 2212

        print("Matched EM, HAD:", matchedEM_E, matchedHAD_E)
        print("Which PFO:", which_pfo)

        if matchedEM_E > 0:
            h_matchedEMpfo_E.Fill(matchedEM_E)
            h_matchedEMpfo_theta.Fill(matchedEM_theta)
            h_deltaEM_E.Fill(matchedEM_E-mcpCollection[0].getEnergy())
        if matchedHAD_E > 0:
            h_matchedHADpfo_E.Fill(matchedHAD_E)
            h_matchedHADpfo_theta.Fill(matchedHAD_theta)
            h_deltaHAD_E.Fill(matchedHAD_E-mcpCollection[0].getEnergy())

        '''
        clusterCollection = event.getCollection('PandoraClusters')
        clus_sumE = 0
        minDRclus = 99999.
        for cluster in clusterCollection:

            px = cluster.getEnergy()*sin(cluster.getITheta())*cos(cluster.getIPhi())
            py = cluster.getEnergy()*sin(cluster.getITheta())*sin(cluster.getIPhi())
            pz = cluster.getEnergy()*cos(cluster.getITheta())

            tlv_clus = TLorentzVector()
            tlv_clus.SetPxPyPzE(px, py, pz, cluster.getEnergy())

            dR = tlv_clus.DeltaR(tlv)

            if dR < minDRclus:
                minDRclus = dR
            clus_sumE = clus_sumE + cluster.getEnergy()
        '''

        jetCollection = event.getCollection('JetOut')
        minDRjet = 999999.
        E_jet[0] = -1
        theta_jet[0] = -1

        for jet in jetCollection:
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)

            if dR < minDRjet:
                minDRjet = dR
                E_jet[0] =jet.getEnergy()
                theta_jet[0] = tlv_pfo.Theta()


        # Fill the hit-level histos and aggregated energy
        ECAL_sumE = 0.
        #ecal_coll = ['EcalBarrelCollectionSel','EcalEndcapCollectionSel']
        ecal_coll = ['EcalBarrelCollectionRec','EcalEndcapCollectionRec']

        for icoll, coll in enumerate(ecal_coll):

            try:
                ECALhitCollection = event.getCollection(coll)
                for hit in ECALhitCollection:
                    ECAL_sumE = ECAL_sumE + hit.getEnergy()
            except:
                print("No", coll, "found")

        sumECAL[0] = ECAL_sumE

        HCAL_sumE = 0.
        #hcal_coll = ['HcalBarrelCollectionSel','HcalEndcapCollectionSel']
        hcal_coll = ['HcalBarrelCollectionRec','HcalEndcapCollectionRec']

        for icoll, coll in enumerate(hcal_coll):

            try:
                HCALhitCollection = event.getCollection(coll)
                for hit in HCALhitCollection:
                    HCAL_sumE = HCAL_sumE + hit.getEnergy()
            except:
                print("No", coll, "found")

        sumHCAL[0] = HCAL_sumE

        #print("Sum clusters:", clus_sumE)
        print("SumCAL ECAL, HCAL:",ECAL_sumE, HCAL_sumE)

        kaon_tree.Fill()

    reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
kaon_tree.Write()
output_file.Close()
