from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TMath, TFile, TLorentzVector, TEfficiency, TH2D
from math import *
from optparse import OptionParser
from itertools import combinations

is_ten_TeV = 3


def getBin(detector, side, layer):
    """
    docstring
    """
    weight = 1.

    '''
    layer_area = [270.40, 270.40, 448.50, 448.50, 655.20, 655.20, 904.80, 904.80,  # VXD barrel
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48,  # VXD endcaps
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48,
                  8117.85, 22034.16, 51678.81,  # IT barrel
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56,  # IT endcaps
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56,
                  140032.91, 194828.39, 249623.88,  # OT barrel
                  69545.45, 69545.45, 69545.45, 69545.45,  # OT endcaps
                  69545.45, 69545.45, 69545.45, 69545.45]
    '''
    layer_area = [270.40, 270.40, 448.50, 655.20, 904.80,  # VXD barrel
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48,  # VXD endcaps
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48,
                  8117.85, 22034.16, 51678.81,  # IT barrel
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56,  # IT endcaps
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56,
                  140032.91, 194828.39, 249623.88,  # OT barrel
                  69545.45, 69545.45, 69545.45, 69545.45,  # OT endcaps
                  69545.45, 69545.45, 69545.45, 69545.45]

    bin_n = 0

    if detector == 1:
        if layer < 3:
            bin_n = layer
        elif layer == 4:
            bin_n = 3
        elif layer == 6:
            bin_n = 4
        elif layer == 8:
            bin_n = 5
    elif detector == 2:
        if side > 0:
            bin_n = layer+8-is_ten_TeV
        else:
            bin_n = layer+16-is_ten_TeV
    elif detector == 3:
        bin_n = layer+24-is_ten_TeV
    elif detector == 4:
        if side > 0:
            bin_n = layer+27-is_ten_TeV
        else:
            bin_n = layer+34-is_ten_TeV
    elif detector == 5:
        bin_n = layer+41-is_ten_TeV
    elif detector == 6:
        if side > 0:
            bin_n = layer+44-is_ten_TeV
        else:
            bin_n = layer+48-is_ten_TeV

    weight = 1./layer_area[bin_n]

    return bin_n, weight


#########################
# parameters
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_occupancy.root',
                  type=str, default='histos_occupancy.root')
(options, args) = parser.parse_args()
#########################

h_nhits_nowei = TH1D('h_nhits_nowei', 'h_nhits_nowei',
                     52-is_ten_TeV, 0., 52-is_ten_TeV)
h_nhits = TH1D('h_nhits', 'h_nhits', 52-is_ten_TeV, 0., 52-is_ten_TeV)
h_ntimehits = TH1D('h_ntimehits', 'h_ntimehits',
                   52-is_ten_TeV, 0., 52-is_ten_TeV)

h_nhits_endcap_2D = TH2D('h_nhits_endcap_2D', 'h_nhits_endcap_2D', 300, -150., 150., 300, -150., 150.)
h_nhits_endcap_R = TH1D('h_nhits_endcap_R', 'h_nhits_endcap_R', 150, 0., 150.)
h_nhits_endcap_R_time = TH1D('h_nhits_endcap_R_time', 'h_nhits_endcap_R_time', 150, 0., 150.)

h_nhits_barrel_z = TH1D('h_nhits_barrel_z', 'h_nhits_barrel_z', 140, -70., 70.)
h_nhits_barrel_z_time = TH1D('h_nhits_barrel_z_time', 'h_nhits_barrel_z_time', 140, -70., 70.)

#get TH1D histogram bin width
z_binwidth = h_nhits_barrel_z.GetBinLowEdge(2) - h_nhits_barrel_z.GetBinLowEdge(1)
R_binwidth = h_nhits_endcap_R.GetBinLowEdge(2) - h_nhits_endcap_R.GetBinLowEdge(1)
# get TH1D range
z_min = h_nhits_barrel_z.GetBinLowEdge(1)
z_max = h_nhits_barrel_z.GetBinLowEdge(h_nhits_barrel_z.GetNbinsX()+1)
z_range = z_max - z_min
R_min = h_nhits_endcap_R.GetBinLowEdge(1)
R_max = h_nhits_endcap_R.GetBinLowEdge(h_nhits_endcap_R.GetNbinsX()+1)
R_range = R_max - R_min

#########################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

allCollections = [
    "OuterTrackerEndcapCollection",
    "OuterTrackerBarrelCollection",
    "InnerTrackerEndcapCollection",
    "InnerTrackerBarrelCollection",
    "VertexEndcapCollection",
    "VertexBarrelCollection"]

timeCollections = [
    "OETrackerHits",
    "OBTrackerHits",
    "IETrackerHits",
    "IBTrackerHits",
    "VETrackerHits",
    "VBTrackerHits"]

totEv = 0

# loop over all events in the file
for ievt, event in enumerate(reader):

    totEv = totEv+1

    if ievt % 10 == 0:
        print("Processing "+str(ievt))

    # do all first
    for coll in allCollections:

        print("Collection: ", coll)

        # setting decoder
        hitsCollection = event.getCollection(coll)
        encoding = hitsCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        # looping over hit collections
        for ihit, hit in enumerate(hitsCollection):
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            detector = decoder["system"].value()
            side = decoder["side"].value()

            bin, wei = getBin(detector, side, layer)
            h_nhits.Fill(bin, wei)
            h_nhits_nowei.Fill(bin)

            if detector == 2:
                if layer == 3:
                    pos = hit.getPosition()
                    h_nhits_endcap_2D.Fill(pos[0], pos[1])
                    h_nhits_endcap_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]), R_range*wei/R_binwidth)

            if detector == 1:
                if layer == 0:
                    pos = hit.getPosition()
                    h_nhits_barrel_z.Fill(pos[2], z_range*wei/z_binwidth)

    # now time
    for coll in timeCollections:

        # setting decoder
        hitsCollection = event.getCollection(coll)
        encoding = hitsCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        # looping over hit collections
        for ihit, hit in enumerate(hitsCollection):
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            detector = decoder["system"].value()
            side = decoder["side"].value()

            bin, wei = getBin(detector, side, layer)
            h_ntimehits.Fill(bin, wei)

            if detector == 2:
                if layer == 3:
                    pos = hit.getPosition()
                    h_nhits_endcap_R_time.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]), R_range*wei/R_binwidth)

            if detector == 1:
                if layer == 0:
                    pos = hit.getPosition()
                    h_nhits_barrel_z_time.Fill(pos[2], z_range*wei/z_binwidth)

h_nhits.Scale(1./totEv)
h_ntimehits.Scale(1./totEv)

##################
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_nhits.Write()
h_ntimehits.Write()
h_nhits_nowei.Write()
h_nhits_endcap_2D.Write()
h_nhits_endcap_R.Write()
h_nhits_endcap_R_time.Write()
h_nhits_barrel_z.Write()
h_nhits_barrel_z_time.Write()
output_file.Close()
