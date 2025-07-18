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
h_mjj = TH1D('mjj', 'mjj', 100, 0, 1000)
h_truth_mjj = TH1D('truth_mjj', 'truth_mjj', 100, 0, 1000)

h_correction_visible = TProfile2D('h_vis', 'h_vis',
                          20, 0, TMath.Pi(), len(arrBins_E)-1, arrBins_E,
                          's')
h_correction_total = TProfile2D('h_total', 'h_total',
                          20, 0, TMath.Pi(), len(arrBins_E)-1, arrBins_E,
                          's')
h_correction_equalbins = TProfile2D('h_equalbins', 'h_equalbins',
                          20, 0, TMath.Pi(), 20, 0., 200., 0., 2.,
                          's')


# Histo list for writing to outputs
histos_list = [h_mjj, h_truth_mjj, h_correction_visible, h_correction_total, h_correction_equalbins]

for histo in histos_list:
    histo.SetDirectory(0)

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))

    '''
    jetCollection = event.getCollection('ValenciaJetOut')

    dp3 = jetCollection[0].getMomentum()
    tlv = TLorentzVector()
    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jetCollection[0].getEnergy())

    print(str(tlv.Perp()) + " " + str(len(jetCollection[0].getParticles())))
    
    dp32 = jetCollection[1].getMomentum()
    tlv2 = TLorentzVector()
    tlv2.SetPxPyPzE(dp32[0], dp32[1], dp32[2], jetCollection[1].getEnergy())

    print(str(tlv2.Perp()) + " " + str(len(jetCollection[1].getParticles())))

    print("Mass", (tlv+tlv2).M())

    h_mjj.Fill((tlv+tlv2).M())
    '''

    #Build the truth jets
    MCparticleCollection = event.getCollection('MCParticle')

    out_quark1 = MCparticleCollection[0]
    out_quark2 = MCparticleCollection[0]

    #find first particle with 2 daughters
    for mcparticle in MCparticleCollection:
        if len(mcparticle.getDaughters()) == 2:
            out_quark1 = mcparticle.getDaughters()[0]
            out_quark2 = mcparticle.getDaughters()[1]  
            break

    tlv_d1 = TLorentzVector()
    tlv_d2 = TLorentzVector()
    tlv_d1.SetPxPyPzE(out_quark1.getMomentum()[0],
                     out_quark1.getMomentum()[1],
                     out_quark1.getMomentum()[2],
                     out_quark1.getEnergy())
    tlv_d2.SetPxPyPzE(out_quark2.getMomentum()[0],
                     out_quark2.getMomentum()[1],
                     out_quark2.getMomentum()[2],
                     out_quark2.getEnergy())
    
    tlv_truthJet1 = TLorentzVector()
    tlv_truthJet2 = TLorentzVector()

    nu_pdg = [12, 14, 16]  # neutrino PDGs
    for mcparticle in MCparticleCollection:

        if mcparticle.getGeneratorStatus() == 1:
            if fabs(mcparticle.getPDG()) in nu_pdg:
                continue

            tlv = TLorentzVector()
            tlv.SetPxPyPzE(mcparticle.getMomentum()[0],
                          mcparticle.getMomentum()[1],
                          mcparticle.getMomentum()[2],
                          mcparticle.getEnergy())

            if tlv.DeltaR(tlv_d1) < 1.2:
                tlv_truthJet1 += tlv
            if tlv.DeltaR(tlv_d2) < 1.2:
                tlv_truthJet2 += tlv

    #print("Visible mjj", (tlv_truthJet1 + tlv_truthJet2).M())
    h_truth_mjj.Fill((tlv_truthJet1 + tlv_truthJet2).M())

    tlv_j1 = TLorentzVector()
    tlv_j2 = TLorentzVector()

    jetCollection = event.getCollection('ValenciaJetOut')
    
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
    
    h_correction_total.Fill(tlv_d1.Theta(), tlv_d1.E(), tlv_d1.E() / tlv_j1.E())
    h_correction_total.Fill(tlv_d2.Theta(), tlv_d2.E(), tlv_d2.E() / tlv_j2.E())

    h_correction_equalbins.Fill(tlv_truthJet1.Theta(), tlv_truthJet1.E(), tlv_truthJet1.E() / tlv_j1.E())
    h_correction_equalbins.Fill(tlv_truthJet2.Theta(), tlv_truthJet2.E(), tlv_truthJet2.E() / tlv_j2.E())

    h_mjj.Fill((tlv_j1 + tlv_j2).M())

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
