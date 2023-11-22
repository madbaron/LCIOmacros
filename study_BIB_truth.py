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

h_nHits_VXB = TH1D('n_simHits_VXB', 'n_simHits_VXB', 100, 0, 100)
h_nHits_VXE = TH1D('n_simHits_VXE', 'n_simHits_VXE', 100, 0, 100)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

# Initialize an empty dictionary with counters
my_vb_dict = {}
my_ve_dict = {}

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    VBTrackerHitsCollection = event.getCollection('VertexBarrelCollection')
    encoding = VBTrackerHitsCollection.getParameters(
    ).getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    for hit in VBTrackerHitsCollection:

        cellID = int(hit.getCellID0())
        decoder.setValue(cellID)
        layer = decoder['layer'].value()

        if layer == 0:
            part = hit.getMCParticle()
            # Check if the element is already in the dictionary
            if part not in my_vb_dict:
                # If not, add it to the dictionary with a counter initialized to 1
                my_vb_dict[part] = {
                    "PDG": part.getPDG(), "counter": 1}
            else:
                # If it's already in the dictionary, increment the counter
                my_vb_dict[part]["counter"] += 1

    for key, value in my_vb_dict.items():
        # print(f"Counter for element {key}: {value['counter']}")
        h_nHits_VXB.Fill(value['counter'])

    VETrackerHitsCollection = event.getCollection('VertexEndcapCollection')
    encoding = VETrackerHitsCollection.getParameters(
    ).getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    for hit in VETrackerHitsCollection:

        cellID = int(hit.getCellID0())
        decoder.setValue(cellID)
        layer = decoder['layer'].value()

        if layer == 0:
            part = hit.getMCParticle()

            # Check if the element is already in the dictionary
            if part not in my_ve_dict:
                # If not, add it to the dictionary with a counter initialized to 1
                my_ve_dict[part] = {
                    "message": f"Element {part.getPDG()} found!", "counter": 1}
            else:
                # If it's already in the dictionary, increment the counter
                my_ve_dict[part]["counter"] += 1

    for key, value in my_ve_dict.items():
        # print(f"Counter for element {key}: {value['counter']}")
        h_nHits_VXE.Fill(value['counter'])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_nHits_VXB.Write()
h_nHits_VXE.Write()
output_file.Close()
