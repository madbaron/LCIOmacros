from pyLCIO import IOIMPL
from ROOT import TH1D, TFile, TMath, TTree
from ROOT.Math import PtEtaPhiEVector, VectorUtil
from optparse import OptionParser
from array import array
import os

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_photons.root',
                  type=str, default='ntup_photons.root')
(options, args) = parser.parse_args()

#Global variables
matching_threshold = 0.01
hadfrac_cleaning = 0.1

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 75., 100., 150., 200., 250., 300., 400., 500., 750., 1000.))

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)

h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)

h_alt_E = TH1D('alt_E', 'alt_E', len(arrBins_E)-1, arrBins_E)
h_alt_theta = TH1D('alt_theta', 'alt_theta', len(arrBins_theta)-1, arrBins_theta)

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta,
               h_matched_E, h_matched_theta,
               h_alt_E, h_alt_theta,
               h_Npfo, h_pfo_type
               ]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
photon_tree = TTree("photon_tree", "photon_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
hfrac = array('d', [0])
E_alt = array('d', [0])
phi_alt = array('d', [0])
theta_alt = array('d', [0])
hfrac_alt = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("hfrac", hfrac, 'var/D')
photon_tree.Branch("E_alt",  E_alt,  'var/D')
photon_tree.Branch("phi_alt", phi_alt, 'var/D')
photon_tree.Branch("theta_alt", theta_alt, 'var/D')
photon_tree.Branch("hfrac_alt", hfrac_alt, 'var/D')
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
        h_truth_E.Fill(mcpCollection[0].getEnergy())
        dp3 = mcpCollection[0].getMomentum()
        tlv = PtEtaPhiEVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        E_truth[0] = mcpCollection[0].getEnergy()
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()
        
        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')

        # Alternative photon candidate first from reco only info
        allEM_E = 0.
        allHAD_E = 0.
        max_E = -1.
        max_tlv = PtEtaPhiEVector()

        for pfo in pfoCollection:
            if pfo.getEnergy() > max_E:
                if abs(pfo.getType()) == 22:
                    max_E = pfo.getEnergy()
                    dp3 = pfo.getMomentum()
                    max_tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

        E_alt[0] = -1
        phi_alt[0] = -4
        theta_alt[0] = -1
        hfrac_alt[0] = 2
        tlv_sum = PtEtaPhiEVector()
        n_pfos_added = 0
        for pfo in pfoCollection:

            dp3 = pfo.getMomentum()
            tlv_pfo = PtEtaPhiEVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            dR = VectorUtil.DeltaR(tlv_pfo,max_tlv)

            if dR < matching_threshold:
                tlv_sum = tlv_sum + tlv_pfo

                n_pfos_added = n_pfos_added + 1

                if abs(pfo.getType()) == 22:
                    allEM_E = allEM_E + pfo.getEnergy()
                elif abs(pfo.getType()) == 2112:
                    allHAD_E = allHAD_E + pfo.getEnergy()
        
        h_Npfo.Fill(n_pfos_added)

        E_alt[0] = tlv_sum.E()
        phi_alt[0] = tlv_sum.Phi()
        theta_alt[0] = tlv_sum.Theta()
        if allEM_E+allHAD_E>0:
            hfrac_alt[0] = allHAD_E/(allEM_E+allHAD_E)

        if hfrac_alt[0] < hadfrac_cleaning:
            h_alt_E.Fill(E_truth[0]) #for efficiency calculation
            h_alt_theta.Fill(theta_truth[0])

        # Build a photon candidate summing up all PFOs within matching_threshold of truth
        allEM_E = 0.
        allHAD_E = 0.

        E[0] = -1
        phi[0] = -4
        theta[0] = -1
        hfrac[0] = 2
        tlv_sum = PtEtaPhiEVector()
        for pfo in pfoCollection:

            dp3 = pfo.getMomentum()
            tlv_pfo = PtEtaPhiEVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            dR = VectorUtil.DeltaR(tlv_pfo,tlv)

            if dR < matching_threshold:
                h_pfo_type.Fill(abs(pfo.getType()))
                tlv_sum = tlv_sum + tlv_pfo
            
                if abs(pfo.getType()) == 22:
                    allEM_E = allEM_E + pfo.getEnergy()
                elif abs(pfo.getType()) == 2112:
                    allHAD_E = allHAD_E + pfo.getEnergy()

        if allEM_E+allHAD_E>0:
            hfrac[0] = allHAD_E/(allEM_E+allHAD_E)
        E[0] = tlv_sum.E()
        phi[0] = tlv_sum.Phi()
        theta[0] = tlv_sum.Theta()

        dR_match = VectorUtil.DeltaR(tlv_sum,tlv)

        if hfrac[0] < hadfrac_cleaning:
            if dR_match < 0.1:
                h_matched_E.Fill(E_truth[0]) #for efficiency calculation
                h_matched_theta.Fill(theta_truth[0])

        photon_tree.Fill()

    reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()
output_file.Close()
