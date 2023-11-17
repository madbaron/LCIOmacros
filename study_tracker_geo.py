from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TColor, TH1D, TH2D, TGraph, TFile, TCanvas, gROOT, gStyle, TLegend, TVector3, TMath
from math import *
from optparse import OptionParser
import os
from array import array

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_geo.root',
                  type=str, default='histos_geo.root')
(options, args) = parser.parse_args()

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

x_coord = array('d', [0])
y_coord = array('d', [0])
r_coord = array('d', [0])
z_coord = array('d', [0])

for ievt, event in enumerate(reader):

    hit_collections = []
    IBTrackerHitsCollection = event.getCollection('IBTrackerHits')
    hit_collections.append(IBTrackerHitsCollection)
    VBTrackerHitsCollection = event.getCollection('VBTrackerHits')
    hit_collections.append(VBTrackerHitsCollection)

    # Encoding equal for all detectors
    encoding = VBTrackerHitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    for coll in hit_collections:

        print("Nhits ",len(coll))

        for hit in coll:

            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            detector = decoder["system"].value()

            position = hit.getPosition()
            x_coord.append(position[0])
            y_coord.append(position[1])
            z_coord.append(position[2])
            r_coord.append(sqrt(position[0]*position[0]+position[1]*position[1]))

reader.close()

print("Npoints", len(x_coord))

my_vertex_xy = TGraph(len(x_coord),x_coord,y_coord)
my_vertex_xy.SetMarkerSize(0.03)
my_vertex_xy.SetName("VXD_xy")

my_vertex_rz = TGraph(len(r_coord),z_coord,r_coord)
my_vertex_rz.SetMarkerSize(0.03)
my_vertex_rz.SetName("VXD_rz")

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
my_vertex_xy.Write()
my_vertex_rz.Write()
output_file.Close()
