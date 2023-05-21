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

ecal_hit_E = TH1D('ECAL_hitEnergy', 'ECAL_hitEnergy', 100, 0., 0.1)
hcal_hit_E = TH1D('HCAL_hitEnergy', 'HCAL_hitEnergy', 100, 0., 0.5)
ecal_hit_r_barrel = TH1D('ECAL_barrel', 'ECAL_barrel', 100, 1590, 1850)
ecal_hit_rz = TH2D('ECAL_rz', 'ECAL_rz', 1000, -2800, 2800, 1000, 0, 2000)
ecal_hit_depth_barrel = TH1D(
    'ECAL_barrel_depth', 'ECAL_barrel_depth', 100, 0, 300)
ecal_hit_r_endcap = TH1D('ECAL_endcap', 'ECAL_endcap', 100, 310, 1700)
hcal_hit_r_barrel = TH1D('HCAL_barrel', 'HCAL_barrel', 100, 1852, 3852)
hcal_hit_r_endcap = TH1D('HCAL_endcap', 'HCAL_endcap', 100, 307, 3246)
hcal_hit_depth_barrel = TH1D(
    'HCAL_barrel_depth', 'HCAL_barrel_depth', 100, 0, 2000)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

doTruth = True
r_min_ecal = 1590
r_min_hcal = 1852

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    # ECAL barrel
    hitsCollection = event.getCollection("ECALBarrel")

    if doTruth:
        relationCollection = event.getCollection('RelationCaloHit')
        relation = UTIL.LCRelationNavigator(relationCollection)
        for ihit, hit in enumerate(hitsCollection):
            sim_vec = relation.getRelatedToObjects(hit)
            print(len(sim_vec))
            '''
            for simhit in sim_vec:
                mystring = ""
                for ipart in range(0, simhit.getNMCParticles()):
                    mcpart = simhit.getPDGCont(ipart)
                    if mcpart != 0:
                        mystring = mystring + " " + str(mcpart)
                print(mystring)
            '''

    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        ecal_hit_E.Fill(hit.getEnergy())
        # [1/cm^2], I hope
        density_weight = hit.getEnergy()*100./(2.*pi*r*2210.*2.)
        ecal_hit_depth_barrel.Fill(r-r_min_ecal, density_weight)
        ecal_hit_r_barrel.Fill(r, density_weight)
        ecal_hit_rz.Fill(hit.getPosition()[2], r)

    # ECAL endcap
    hitsCollection = event.getCollection("ECALEndcap")
    for ihit, hit in enumerate(hitsCollection):
        ecal_hit_E.Fill(hit.getEnergy())
        ecal_hit_r_endcap.Fill(sqrt(hit.getPosition()[
            0]*hit.getPosition()[0] + hit.getPosition()[1]*hit.getPosition()[1]))

    # HCAL barrel
    hitsCollection = event.getCollection("HCALBarrel")
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        hcal_hit_E.Fill(hit.getEnergy())
        # [1/cm^2], I hope
        density_weight = hit.getEnergy()*100./(2.*pi*r*2569.*2.)
        hcal_hit_r_barrel.Fill(r, density_weight)
        hcal_hit_depth_barrel.Fill(r-r_min_hcal, density_weight)

    # HCAL endcap
    hitsCollection = event.getCollection("HCALEndcap")
    for ihit, hit in enumerate(hitsCollection):
        hcal_hit_E.Fill(hit.getEnergy())
        hcal_hit_r_endcap.Fill(sqrt(hit.getPosition()[
            0]*hit.getPosition()[0] + hit.getPosition()[1]*hit.getPosition()[1]))

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
ecal_hit_E.Write()
hcal_hit_E.Write()
ecal_hit_r_barrel.Write()
ecal_hit_depth_barrel.Write()
ecal_hit_rz.Write()
ecal_hit_r_endcap.Write()
hcal_hit_r_barrel.Write()
hcal_hit_r_endcap.Write()
hcal_hit_depth_barrel.Write()
output_file.Close()
