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
parser.add_option('-o', '--outFile', help='--outFile histos_BIB.root',
                  type=str, default='histos_BIB.root')
(options, args) = parser.parse_args()


# h_entry_point = TH2D('h_entry_point', 'h_entry_point',
#                     400, -200., 200., 50, 0., 50.)

h_nHits = TH1D('n_simHits', 'n_simHits', 100, 0, 100)  # GeV

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

# Initialize an empty dictionary with counters
my_dict = {}

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    VBTrackerHitsCollection = event.getCollection('VertexBarrelCollection')
    for hit in VBTrackerHitsCollection:

        part = hit.getMCParticle()

        # Check if the element is already in the dictionary
        if part not in my_dict:
            # If not, add it to the dictionary with a counter initialized to 1
            my_dict[part] = {
                "message": f"Element {part.getPDG()} found!", "counter": 1}
        else:
            # If it's already in the dictionary, increment the counter
            my_dict[part]["counter"] += 1

    for key, value in my_dict.items():
        # print(f"Counter for element {key}: {value['counter']}")
        h_nHits.Fill(value['counter'])


reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_nHits.Write()
output_file.Close()
