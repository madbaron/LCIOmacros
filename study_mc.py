from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TMath
from math import *
from optparse import OptionParser
import os
import fnmatch

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_mc.root',
                  type=str, default='histos_mc.root')
(options, args) = parser.parse_args()

#########################

# declare histograms
h_photon_E = TH1D('photon_E', 'photon_E', 100, 0, 1)  # GeV
h_photon_theta = TH1D('photon_theta', 'photon_theta', 100, 0, TMath.Pi())
h_photon_z = TH1D('photon_z', 'photon_z', 100, -100, 100)

h_electron_E = TH1D('electron_E', 'electron_E', 100, 0, 1)  # GeV
h_electron_theta = TH1D('electron_theta', 'electron_theta', 100, 0, TMath.Pi())
h_electron_z = TH1D('electron_z', 'electron_z', 100, -100, 100)

h_hadron_E = TH1D('hadron_E', 'hadron_E', 100, 0, 10)  # GeV
h_hadron_theta = TH1D('hadron_theta', 'hadron_theta', 100, 0, TMath.Pi())
h_hadron_z = TH1D('hadron_z', 'hadron_z', 100, -100, 100)

h_muon_E = TH1D('muon_E', 'muon_E', 100, 0, 1)  # GeV
h_muon_theta = TH1D('muon_theta', 'muon_theta', 100, 0, TMath.Pi())
h_muon_z = TH1D('muon_z', 'muon_z', 100, -100, 100)

h_entry_point = TH2D('entry_point', 'entry_point',
                     400, -200., 200., 50, 0., 50.)
h_entry_point_z = TH1D('entry_point_z', 'entry_point_z', 400, -200., 200.)
h_entry_point_R = TH1D('entry_point_R', 'entry_point_R', 50, 0., 50.)

h_species = TH1D('species', 'species', 50, 0., 50.)
h_species_time = TH1D('species_time', 'species_time', 50, 0., 50.)

# Histo list for writing to outputs
histos_list = [h_photon_E, h_photon_theta, h_photon_z,
               h_electron_E, h_electron_theta, h_electron_z,
               h_hadron_E, h_hadron_theta, h_hadron_z,
               h_muon_E, h_muon_theta, h_muon_z,
               h_entry_point, h_entry_point_z, h_entry_point_R,
               h_species, h_species_time
               ]

for histo in histos_list:
    histo.SetDirectory(0)

#########################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    print(" ")
    print("New event, #" + str(ievt))

    # get mc particle collection and loop over it
    mcpCollection = event.getCollection('MCParticle')
    for part in mcpCollection:

        mom = part.getMomentum()  # GeV
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(mom[0], mom[1], mom[2], part.getEnergy())

        pos = part.getVertex()  # mm

        h_entry_point.Fill(pos[2], sqrt(pos[0]*pos[0]+pos[1]*pos[1]))
        h_entry_point_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))
        h_entry_point_z.Fill(pos[2])

        if part.getEnergy() > 0.0001:

            pdgId = fabs(part.getPDG())
            if pdgId > 22:
                if part.getCharge() > 0:
                    pdgId = 0

            h_species.Fill(pdgId)
            time = part.getTime()
            if time > -0.5:
                if time < 15.:
                    h_species_time.Fill(pdgId)

            if fabs(part.getPDG()) == 22:
                h_photon_E.Fill(part.getEnergy())
                h_photon_theta.Fill(tlv.Theta())
                h_photon_z.Fill(pos[2])
            elif fabs(part.getPDG()) == 11:
                h_electron_E.Fill(part.getEnergy())
                h_electron_theta.Fill(tlv.Theta())
                h_electron_z.Fill(pos[2])
            elif fabs(part.getPDG()) == 13:
                h_muon_E.Fill(part.getEnergy())
                h_muon_theta.Fill(tlv.Theta())
                h_muon_z.Fill(pos[2])
            else:
                if part.getCharge() > 0:
                    h_hadron_E.Fill(part.getEnergy())
                    h_hadron_theta.Fill(tlv.Theta())
                    h_hadron_z.Fill(pos[2])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
