from pyLCIO import IOIMPL
from ROOT import TH1D, TFile, TLorentzVector, TProfile2D, TMath
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_jets.root',
                  type=str, default='ntup_jets.root')
parser.add_option('-p', '--photonCalibFile', help='--photonCalibFile ResponseMap_reco_photons_v5_noBIB.root',
                  type=str, default='ResponseMap_reco_photons_v5_noBIB.root')
parser.add_option('-n', '--neutronCalibFile', help='--neutronCalibFile ResponseMap_reco_neutrons_v5_noBIB.root',
                  type=str, default='ResponseMap_reco_neutrons_v5_noBIB.root')
(options, args) = parser.parse_args()

arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 35.,
                         40., 45., 50., 75., 100., 250., 500.))

# declare histograms
h_mjj = TH1D('mjj', 'mjj', 150, 0, 150)
h_truth_mjj = TH1D('truth_mjj', 'truth_mjj', 150, 0, 150)

h_correction_visible = TProfile2D('h_vis', 'h_vis',
                          20, 0, TMath.Pi(), len(arrBins_E)-1, arrBins_E,
                          's')
h_correction_equalbins = TProfile2D('h_equalbins', 'h_equalbins',
                          20, 0, TMath.Pi(), 20, 0., 200., 0., 2.,
                          's')

# Histo list for writing to outputs
histos_list = [h_mjj, h_truth_mjj, h_correction_visible, h_correction_equalbins]

for histo in histos_list:
    histo.SetDirectory(0)

calibFile_photons = TFile(options.photonCalibFile, "READ")
calibMap_photons = calibFile_photons.Get("calib_2d")
calibFile_neutrons = TFile(options.neutronCalibFile, "READ")
calibMap_neutrons = calibFile_neutrons.Get("calib_2d")

def get_calibrated_tlv(particle):
    tlv = TLorentzVector()
    E = particle.getEnergy()
    dp3 = particle.getMomentum()
    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], E)

    theta = tlv.Theta()
    E_cap = E
    correction = 1.0

    if particle.getType() == 22:
        if E_cap > 1000.:
            E_cap = 999.
        correction = calibMap_photons.GetBinContent(calibMap_photons.FindBin(theta, E_cap))
    elif particle.getType() == 2112:
        if E_cap > 250.:
            E_cap = 249.        
        correction = calibMap_neutrons.GetBinContent(calibMap_neutrons.FindBin(theta, E_cap))
    else:
        correction = 1.0

    if correction < 0.01:
        correction = 1.0

    E = E*correction
    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], E)

    return tlv


# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))

    # get the truth jets
    truthJetCollection = event.getCollection('TruthValenciaJetOut')

    if len(truthJetCollection) != 2:
        print("Skipping event " + str(ievt) + " because it does not have exactly two truth jets")
        continue

    tlv_truthJet1 = TLorentzVector()
    dp3 = truthJetCollection[0].getMomentum()
    tlv_truthJet1.SetPxPyPzE(dp3[0], dp3[1], dp3[2], truthJetCollection[0].getEnergy())

    tlv_truthJet2 = TLorentzVector()
    dp3 = truthJetCollection[1].getMomentum()
    tlv_truthJet2.SetPxPyPzE(dp3[0], dp3[1], dp3[2], truthJetCollection[1].getEnergy())

    h_truth_mjj.Fill((tlv_truthJet1 + tlv_truthJet2).M())

    ij1 = -1
    ij2 = -1

    jetCollection = event.getCollection('ValenciaJetOut')
    if len(jetCollection) < 2:
        print("Skipping event " + str(ievt) + " because it does not have at least two jets")
        continue

    deltaR_j = 99
    for ijet, jet in enumerate(jetCollection):
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet1)

        if deltaR < deltaR_j:
            ij1 = ijet
            deltaR_j = deltaR

    deltaR_j = 99
    for ijet, jet in enumerate(jetCollection):
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet2)

        if deltaR < deltaR_j:
            ij2 = ijet
            deltaR_j = deltaR
    
    tlv_j1 = TLorentzVector()
    tlv_j2 = TLorentzVector()

    jet1 = jetCollection[ij1]
    for constituent in jet1.getParticles():
        tlv_j1 += get_calibrated_tlv(constituent)

    jet2 = jetCollection[ij2]
    for constituent in jet2.getParticles():
        tlv_j2 += get_calibrated_tlv(constituent)

    h_correction_visible.Fill(tlv_truthJet1.Theta(), tlv_truthJet1.E(), tlv_truthJet1.E() / tlv_j1.E())
    h_correction_visible.Fill(tlv_truthJet2.Theta(), tlv_truthJet2.E(), tlv_truthJet2.E() / tlv_j2.E())

    h_correction_equalbins.Fill(tlv_truthJet1.Theta(), tlv_truthJet1.E(), tlv_truthJet1.E() / tlv_j1.E())
    h_correction_equalbins.Fill(tlv_truthJet2.Theta(), tlv_truthJet2.E(), tlv_truthJet2.E() / tlv_j2.E())

    h_mjj.Fill((tlv_j1 + tlv_j2).M())

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
