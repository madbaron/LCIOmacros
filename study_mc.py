from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath
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

Bfield = 5.00  # T

#########################

# declare histograms
h_photon_pt = TH1D('photon_pt', 'photon_pt', 100, 0, 10)  # GeV
h_photon_theta = TH1D('photon_theta', 'photon_theta', 100, 0, TMath.Pi())
h_photon_z = TH1D('photon_z', 'photon_z', 100, -100, 100)

h_electron_pt = TH1D('electron_pt', 'electron_pt', 100, 0, 1)  # GeV
h_electron_theta = TH1D('electron_theta', 'electron_theta', 100, 0, TMath.Pi())
h_electron_z = TH1D('electron_z', 'electron_z', 100, -100, 100)

h_positron_pt = TH1D('positron_pt', 'positron_pt', 100, 0, 1)  # GeV
h_positron_theta = TH1D('positron_theta', 'positron_theta', 100, 0, TMath.Pi())
h_positron_z = TH1D('positron_z', 'positron_z', 100, -100, 100)

h_hadron_pt = TH1D('hadron_pt', 'hadron_pt', 100, 0, 1)  # GeV
h_hadron_theta = TH1D('hadron_theta', 'hadron_theta', 100, 0, TMath.Pi())
h_hadron_z = TH1D('hadron_z', 'hadron_z', 100, -100, 100)

h_entry_point = TH2D('h_entry_point', 'h_entry_point',
                     400, -200., 200., 50, 0., 50.)

# Histo list for writing to outputs
histos_list = [h_photon_pt, h_photon_theta, h_photon_z,
               h_electron_pt, h_electron_theta, h_electron_z,
               h_positron_pt, h_positron_theta, h_positron_z,
               h_hadron_pt, h_hadron_theta, h_hadron_z,
               h_entry_point
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

        if fabs(part.getPDG()) == 22:
            h_photon_pt.Fill(tlv.Perp())
            h_photon_theta.Fill(tlv.Theta())
            h_photon_z.Fill(pos[2])
        elif part.getPDG() == 11:
            h_electron_pt.Fill(tlv.Perp())
            h_electron_theta.Fill(tlv.Theta())
            h_electron_z.Fill(pos[2])
        elif part.getPDG() == -11:
            h_positron_pt.Fill(tlv.Perp())
            h_positron_theta.Fill(tlv.Theta())
            h_positron_z.Fill(pos[2])
        else:
            if part.getCharge() > 0:
                h_hadron_pt.Fill(tlv.Perp())
                h_hadron_theta.Fill(tlv.Theta())
                h_hadron_z.Fill(pos[2])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
