from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
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

# declare histograms
h_mjj = TH1D('mjj', 'mjj', 100, 50, 150)

# Histo list for writing to outputs
histos_list = [h_mjj]

for histo in histos_list:
    histo.SetDirectory(0)

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))

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

    MCparticleCollection = event.getCollection('MCParticle')
    tlv_d1 = TLorentzVector()
    tlv_d2 = TLorentzVector()
    for mcparticle in MCparticleCollection:
        if mcparticle.getPDG() == 23:
            print("Z", mcparticle.getEnergy())
            daughters = mcparticle.getDaughters()
            tlv_d1.SetPxPyPzE(daughters[0].getMomentum()[0], daughters[0].getMomentum()[1], daughters[0].getMomentum()[2], daughters[0].getEnergy())
            tlv_d2.SetPxPyPzE(daughters[1].getMomentum()[0], daughters[1].getMomentum()[1], daughters[1].getMomentum()[2], daughters[1].getEnergy())
            break

    print("D1", tlv_d1.Perp())
    print("D2", tlv_d2.Perp())

    tlv_j1 = TLorentzVector()
    tlv_j2 = TLorentzVector()

    jetCollection = event.getCollection('JetOut')
    
    deltaR_j = 99
    for jet in jetCollection:
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        has_charged_pfo = False
        for constituent in jet.getParticles():
            if constituent.getCharge() != 0:
                has_charged_pfo = True
                break

        if not has_charged_pfo:
            continue

        deltaR = tlv.DeltaR(tlv_d1)

        if deltaR < deltaR_j:
            tlv_j1 = tlv
            deltaR_j = deltaR

    deltaR_j = 99
    for jet in jetCollection:
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        has_charged_pfo = False
        for constituent in jet.getParticles():
            if constituent.getCharge() != 0:
                has_charged_pfo = True
                break

        if not has_charged_pfo:
            continue

        deltaR = tlv.DeltaR(tlv_d2)

        if deltaR < deltaR_j:
            tlv_j2 = tlv
            deltaR_j = deltaR

    print(tlv_j1.Perp(), tlv_j2.Perp(), (tlv_j1+tlv_j2).M())

        #for constituent in jet.getParticles():
        #    ids = constituent.getParticleIDs()
        #    print(constituent.getEnergy(), constituent.getCharge(),
        #          len(constituent.getTracks()))

            # if len(ids) > 0:
            #    print(str(ids[0]))
    

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
