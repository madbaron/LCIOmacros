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
(options, args) = parser.parse_args()

arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 35.,
                         40., 45., 50., 75., 100., 250., 500.))

# declare histograms
h_mjj = TH1D('mjj', 'mjj', 100, 50, 150)
h_truth_mjj = TH1D('truth_mjj', 'truth_mjj', 100, 50, 150)

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

    tlv_j1 = TLorentzVector()
    tlv_j2 = TLorentzVector()

    jetCollection = event.getCollection('ValenciaJetOut')
    if len(jetCollection) < 2:
        print("Skipping event " + str(ievt) + " because it does not have at least two jets")
        continue

    deltaR_j = 99
    for jet in jetCollection:
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet1)

        if deltaR < deltaR_j:
            tlv_j1 = tlv
            deltaR_j = deltaR

    deltaR_j = 99
    for jet in jetCollection:
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet2)

        if deltaR < deltaR_j:
            tlv_j2 = tlv
            deltaR_j = deltaR
    
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
