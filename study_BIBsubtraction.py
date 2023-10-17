from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree, TVector3
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50.,
                  100., 250., 500., 1000., 2500., 5000.))

# Low-level digitised hit distributions
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)

# Histo list for writing to outputs
histos_list = [
    h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
    h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
    h_ECAL_sumE, h_HCAL_sumE,
    h_EMfrac,
    h_ECAL_hit_layer, h_HCAL_hit_layer,
]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
hit_tree = TTree("hit_tree", "hit_tree")
E_hit = array('d', [0])
phi_hit = array('d', [0])
theta_hit = array('d', [0])
time_hit = array('d', [0])
layer_hit = array('i', [0])
isECAL_hit = array('i', [0])
hit_tree.Branch("E",  E_hit,  'var/D')
hit_tree.Branch("phi", phi_hit, 'var/D')
hit_tree.Branch("theta", theta_hit, 'var/D')
hit_tree.Branch("time", time_hit, 'var/D')
hit_tree.Branch("layer", layer_hit, 'var/I')
hit_tree.Branch("isECAL", isECAL_hit, 'var/I')

simhit_tree = TTree("simhit_tree", "hit_tree")
E_simhit = array('d', [0])
phi_simhit = array('d', [0])
theta_simhit = array('d', [0])
layer_simhit = array('i', [0])
isECAL_simhit = array('i', [0])
simhit_tree.Branch("E",  E_simhit,  'var/D')
simhit_tree.Branch("phi", phi_simhit, 'var/D')
simhit_tree.Branch("theta", theta_simhit, 'var/D')
simhit_tree.Branch("layer", layer_simhit, 'var/I')
simhit_tree.Branch("isECAL", isECAL_simhit, 'var/I')

to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # loop over all events in the file
    for ievt, event in enumerate(reader):

        if ievt % 1 == 0:
            print(" ")
            print("Processing event " + str(ievt))

        # Fill the hit-level histos and aggregated energy
        ECAL_sumE = 0.
        ecal_coll = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']
        for coll in ecal_coll:

            try:
                ECALhitCollection = event.getCollection(coll)
            except:
                print("No ", coll, " found!")
                continue

            encoding = ECALhitCollection.getParameters(
            ).getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)

            for hit in ECALhitCollection:
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder["layer"].value()

                h_ECAL_hit_time.Fill(hit.getTime())
                h_ECAL_hit_E.Fill(hit.getEnergy())
                h_ECAL_hit_layer.Fill(layer, hit.getEnergy())

                ECAL_sumE = ECAL_sumE + hit.getEnergy()
                pos = hit.getPosition()
                h_ECAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                posvec = TVector3(pos[0], pos[1], pos[2])

                E_hit[0] = hit.getEnergy()
                phi_hit[0] = posvec.Phi()
                theta_hit[0] = posvec.Theta()
                time_hit[0] = hit.getTime()
                layer_hit[0] = layer
                isECAL_hit[0] = 1
                hit_tree.Fill()

            h_ECAL_sumE.Fill(ECAL_sumE)

        HCAL_sumE = 0.
        hcal_coll = ['HcalBarrelsCollectionRec', 'HcalEndcapsCollectionRec']
        for coll in hcal_coll:

            try:
                HCALhitCollection = event.getCollection(coll)
            except:
                print("No ", coll, " found!")
                continue

            encoding = HCALhitCollection.getParameters(
            ).getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)

            for hit in HCALhitCollection:
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder["layer"].value()

                h_HCAL_hit_time.Fill(hit.getTime())
                h_HCAL_hit_E.Fill(hit.getEnergy())
                h_HCAL_hit_layer.Fill(layer, hit.getEnergy())
                HCAL_sumE = HCAL_sumE + hit.getEnergy()
                pos = hit.getPosition()
                h_HCAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                posvec = TVector3(pos[0], pos[1], pos[2])

                E_hit[0] = hit.getEnergy()
                phi_hit[0] = posvec.Phi()
                theta_hit[0] = posvec.Theta()
                time_hit[0] = hit.getTime()
                layer_hit[0] = layer
                isECAL_hit[0] = 0
                hit_tree.Fill()

            h_HCAL_sumE.Fill(HCAL_sumE)

        print(ECAL_sumE, HCAL_sumE)

        h_sumE.Fill(ECAL_sumE+HCAL_sumE)

        if ECAL_sumE+HCAL_sumE > 0:
            h_EMfrac.Fill(ECAL_sumE/(ECAL_sumE+HCAL_sumE))
        else:
            h_EMfrac.Fill(0)

        # SimHit
        ecal_coll = ['ECalBarrelCollection', 'ECalEndcapCollection']
        for coll in ecal_coll:

            try:
                simECALhitCollection = event.getCollection(coll)
            except:
                print("No ", coll, " found!")
                continue

            encoding = simECALhitCollection.getParameters(
            ).getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)

            for hit in simECALhitCollection:
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder["layer"].value()

                pos = hit.getPosition()
                posvec = TVector3(pos[0], pos[1], pos[2])

                E_simhit[0] = hit.getEnergy()
                phi_simhit[0] = posvec.Phi()
                theta_simhit[0] = posvec.Theta()
                layer_simhit[0] = layer
                isECAL_simhit[0] = 1
                simhit_tree.Fill()

        hcal_coll = ['HCalBarrelCollection', 'HCalEndcapCollection']
        for coll in hcal_coll:

            try:
                simHCALhitCollection = event.getCollection(coll)
            except:
                print("No ", coll, " found!")
                continue

            encoding = simHCALhitCollection.getParameters(
            ).getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)

            for hit in simHCALhitCollection:
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder["layer"].value()

                pos = hit.getPosition()
                posvec = TVector3(pos[0], pos[1], pos[2])

                E_simhit[0] = hit.getEnergy()
                phi_simhit[0] = posvec.Phi()
                theta_simhit[0] = posvec.Theta()
                layer_simhit[0] = layer
                isECAL_simhit[0] = 0
                simhit_tree.Fill()

    reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_BIBsub.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
hit_tree.Write()
simhit_tree.Write()
output_file.Close()
