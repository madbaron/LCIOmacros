from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TColor, TH1D, TH2D, TFile, TCanvas, gROOT, gStyle, TLegend, TVector3, TMath, TRandom3, TGraph2D
from math import *
from optparse import OptionParser
from array import array

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_BIB.root',
                  type=str, default='histos_BIB.root')
(options, args) = parser.parse_args()

speedoflight = 299792458/1000000  # mm/ns

# h_entry_point = TH2D('h_entry_point', 'h_entry_point',
#                     400, -200., 200., 50, 0., 50.)

h_nHits_VXB = TH1D('n_simHits_VXB', 'n_simHits_VXB', 100, 0, 100)
h_nHits_VXE = TH1D('n_simHits_VXE', 'n_simHits_VXE', 100, 0, 100)

h_nHits_VXB_time = TH1D('n_simHits_VXB_time',
                        'n_simHits_VXB_time', 100, 0, 100)
h_nHits_VXE_time = TH1D('n_simHits_VXE_time',
                        'n_simHits_VXE_time', 100, 0, 100)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

for ievt, event in enumerate(reader):

    rndm = TRandom3(ievt)

    # Initialize an empty dictionary with counters
    my_vb_dict = {}
    my_vb_dict_time = {}
    my_ve_dict = {}
    my_ve_dict_time = {}

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
                hit_list = [hit]
                my_vb_dict[part] = {
                    "PDG": part.getPDG(), "counter": 1, "hits": hit_list}
            else:
                # If it's already in the dictionary, increment the counter
                my_vb_dict[part]["counter"] += 1
                my_vb_dict[part]["hits"].append(hit)

            pos = hit.getPosition()  # mm
            d = sqrt(pos[0]*pos[0] + pos[1]* pos[1] + pos[2]*pos[2])
            tof = d/speedoflight
            corrected_time = hit.getTime()*(1.+rndm.Gaus(0., 0.03)) - tof

            if (corrected_time > -0.09) and (corrected_time < 0.15):
                # Check if the element is already in the dictionary
                if part not in my_vb_dict_time:
                    # If not, add it to the dictionary with a counter initialized to 1
                    my_vb_dict_time[part] = {
                        "message": f"Element {part.getPDG()} found!", "counter": 1}
                else:
                    # If it's already in the dictionary, increment the counter
                    my_vb_dict_time[part]["counter"] += 1


    max_hits = 0
    the_hits = []
    for key, value in my_vb_dict.items():
        h_nHits_VXB.Fill(value['counter'])
        if value['counter'] > max_hits:
            max_hits = value['counter']
            the_hits = value['hits']

    print("Max hits", len(the_hits))
    x = array('d')
    y = array('d')
    z = array('d')
    for hit in the_hits:
        x.append(hit.getPosition()[0])
        y.append(hit.getPosition()[1])
        z.append(hit.getPosition()[2])

    print(z)        
    print(x)        
    print(y)        
    my_trajectory = TGraph2D(len(the_hits), z, x, y)

    for key, value in my_vb_dict_time.items():
        h_nHits_VXB_time.Fill(value['counter'])

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

            pos = hit.getPosition()  # mm
            d = sqrt(pos[0]*pos[0] + pos[1]* pos[1] + pos[2]*pos[2])
            tof = d/speedoflight
            corrected_time = hit.getTime()*(1.+rndm.Gaus(0., 0.03)) - tof

            if (corrected_time > -0.09) and (corrected_time < 0.15):
                # Check if the element is already in the dictionary
                if part not in my_ve_dict_time:
                    # If not, add it to the dictionary with a counter initialized to 1
                    my_ve_dict_time[part] = {
                        "message": f"Element {part.getPDG()} found!", "counter": 1}
                else:
                    # If it's already in the dictionary, increment the counter
                    my_ve_dict_time[part]["counter"] += 1

    for key, value in my_ve_dict.items():
        h_nHits_VXE.Fill(value['counter'])

    for key, value in my_ve_dict_time.items():
        h_nHits_VXE_time.Fill(value['counter'])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_nHits_VXB.Write()
h_nHits_VXE.Write()
h_nHits_VXB_time.Write()
h_nHits_VXE_time.Write()
my_trajectory.Write()
output_file.Close()
