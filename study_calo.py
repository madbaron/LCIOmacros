from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TColor, TH1D, TH2D, TFile, TCanvas, gROOT, gStyle, TLegend, TVector3, TMath
from math import *
from optparse import OptionParser
import os
import fnmatch

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_hits.root',
                  type=str, default='histos_hits.root')
(options, args) = parser.parse_args()

hit_E = TH1D('hitEnergy', 'hitEnergy', 100, 0., 0.001)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    # setting decoder
    muonHitsCollection = event.getCollection('MUON')

    for hit in muonHitsCollection:
        hit_E.Fill(hit.getEnergy())
        # print(hit.getEnergy())

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
hit_E.Write()
output_file.Close()
