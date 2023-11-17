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
parser.add_option('-o', '--outFile', help='--outFile histos_photons.root',
                  type=str, default='histos_photons.root')
(options, args) = parser.parse_args()

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000., 1500., 2500., 5000.))

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta, h_matched_theta, h_matched_E,
               h_Npfo, h_pfo_type
               ]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
photon_tree = TTree("photon_tree", "photon_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("E_truth",  E_truth,  'var/D')
photon_tree.Branch("phi_truth", phi_truth, 'var/D')
photon_tree.Branch("theta_truth", theta_truth, 'var/D')

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
        pfoCollection = event.getCollection('PandoraPFOs')

        h_Npfo.Fill(len(pfoCollection))

        for mcp in mcpCollection:

            E_truth[0] = -1.
            phi_truth[0] = -99.
            theta_truth[0] = -1.
            E[0] =  -1.
            phi[0] = -99.
            theta[0] = -1.

            if (mcp.getGeneratorStatus() == 1) and (mcp.getPDG()==22):
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                h_truth_E.Fill(mcp.getEnergy())
                h_truth_theta.Fill(tlv.Theta())
                E_truth[0] = mcp.getEnergy()
                phi_truth[0] = tlv.Phi()
                theta_truth[0] = tlv.Theta()

                # Match true pfo with closest reco PFO in deltaR
                matched_tlv = TLorentzVector()
                minDR = 999999.

                for pfo in pfoCollection:
                    h_pfo_type.Fill(abs(pfo.getType()))

                    if (pfo.getType() == 22) or (pfo.getType() == 2112):
                        dp3 = pfo.getMomentum()
                        tlv_pfo = TLorentzVector()
                        tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                        dR = tlv_pfo.DeltaR(tlv)

                        if dR < minDR:
                            minDR = dR
                            matched_tlv = tlv_pfo

                #look for hadronic punch-through
                for pfo in pfoCollection:
                    if abs(pfo.getType()) == 2112:
                        dp3 = pfo.getMomentum()
                        tlv_pfo = TLorentzVector()
                        tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                        dR = tlv_pfo.DeltaR(matched_tlv)

                        if dR < 0.4:
                            matched_tlv = matched_tlv+tlv_pfo
  
                if minDR<0.4:
                    E[0] = matched_tlv.E()
                    phi[0] = matched_tlv.Phi()
                    theta[0] = matched_tlv.Theta()
                    h_matched_E.Fill(mcp.getEnergy())
                    h_matched_theta.Fill(tlv.Theta())

                photon_tree.Fill()

    reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()
output_file.Close()
