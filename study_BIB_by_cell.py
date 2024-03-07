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
parser.add_option('-o', '--outFile', help='--outFile histos_calo.root',
                  type=str, default='histos_calo.root')
(options, args) = parser.parse_args()

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    # ECAL barrel
    ECALhitCollection = event.getCollection("ECalBarrelCollection")
    
    encoding = ECALhitCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    for ihit, hit in enumerate(ECALhitCollection):
        
        n_contributions =  hit.getNMCContributions()

        if n_contributions>1:

            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            print(hit.getEnergy(), layer, n_contributions)

            for icont in range(n_contributions):
                print(hit.getPDGCont(icont), hit.getEnergyCont(icont))

            print("")

            break



