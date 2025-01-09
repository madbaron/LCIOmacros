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
                  type=str, default='histos_calo_cell.root')
(options, args) = parser.parse_args()

h_ECAL_L0_photon = TH1D('ECAL_L0_photon', 'ECAL_L0_photon', 100, 0, 0.1)  # MeV
h_ECAL_L0_neutron = TH1D('ECAL_L0_neutron', 'ECAL_L0_neutron', 100, 0, 0.1)  # MeV
h_ECAL_L0_other = TH1D('ECAL_L0_other', 'ECAL_L0_other', 100, 0, 0.1)  # MeV

h_ECAL_L15_photon = TH1D('ECAL_L15_photon', 'ECAL_L15_photon', 100, 0, 0.1)  # MeV
h_ECAL_L15_neutron = TH1D('ECAL_L15_neutron', 'ECAL_L15_neutron', 100, 0, 0.1)  # MeV
h_ECAL_L15_other = TH1D('ECAL_L15_other', 'ECAL_L15_other', 100, 0, 0.1)  # MeV

h_ECAL_L40_photon = TH1D('ECAL_L40_photon', 'ECAL_L40_photon', 100, 0, 0.1)  # MeV
h_ECAL_L40_neutron = TH1D('ECAL_L40_neutron', 'ECAL_L40_neutron', 100, 0, 0.1)  # MeV
h_ECAL_L40_other = TH1D('ECAL_L40_other', 'ECAL_L40_other', 100, 0, 0.1)  # MeV

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    # ECAL barrel
    ECALhitCollection = event.getCollection("EcalBarrelCollectionRec")
    relECALhitCollection = event.getCollection("EcalBarrelRelationsSimRec")
    relation = UTIL.LCRelationNavigator(relECALhitCollection)
    
    encoding = ECALhitCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    for ihit, hit in enumerate(ECALhitCollection):

        cellID = int(hit.getCellID0())
        decoder.setValue(cellID)
        layer = decoder["layer"].value()

        if layer==0 or layer==15 or layer==40:

            pos = hit.getPosition()
            posvec = TVector3(pos[0], pos[1], pos[2])
            folded = posvec.Theta()
            if (folded > TMath.Pi()/2):
                folded = TMath.Pi()-posvec.Theta()

            if fabs(TMath.Pi()/2 - folded)*180./TMath.Pi() > 5:
                continue

            rel_simhits = relation.getRelatedToObjects(hit)
            for simhit in rel_simhits:

                n_contributions =  simhit.getNMCContributions()

                photon_frac = 0.
                neutron_frac = 0.
                other_frac = 0.

                for icont in range(n_contributions):

                    pdg_shower = simhit.getParticleCont(icont).getPDG()

                    if pdg_shower == 2112:
                        neutron_frac = neutron_frac + simhit.getEnergyCont(icont)
                    elif pdg_shower == 22:
                        photon_frac = photon_frac + simhit.getEnergyCont(icont)
                    else:
                        other_frac = other_frac + simhit.getEnergyCont(icont)

                sumsimE = photon_frac + neutron_frac + other_frac
                photon_frac = photon_frac/sumsimE
                neutron_frac = neutron_frac/sumsimE
                other_frac = other_frac/sumsimE

                if layer == 0:
                    h_ECAL_L0_photon.Fill(hit.getEnergy(), photon_frac)
                    h_ECAL_L0_neutron.Fill(hit.getEnergy(), neutron_frac)
                    h_ECAL_L0_other.Fill(hit.getEnergy(), other_frac)
                elif layer==15:
                    h_ECAL_L15_photon.Fill(hit.getEnergy(), photon_frac)
                    h_ECAL_L15_neutron.Fill(hit.getEnergy(), neutron_frac)
                    h_ECAL_L15_other.Fill(hit.getEnergy(), other_frac)
                else:
                    h_ECAL_L40_photon.Fill(hit.getEnergy(), photon_frac)
                    h_ECAL_L40_neutron.Fill(hit.getEnergy(), neutron_frac)
                    h_ECAL_L40_other.Fill(hit.getEnergy(), other_frac)

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_ECAL_L0_photon.Write()
h_ECAL_L0_neutron.Write()
h_ECAL_L0_other.Write()
h_ECAL_L15_photon.Write()
h_ECAL_L15_neutron.Write()
h_ECAL_L15_other.Write()
h_ECAL_L40_photon.Write()
h_ECAL_L40_neutron.Write()
h_ECAL_L40_other.Write()
output_file.Close()

