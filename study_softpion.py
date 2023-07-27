from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TMath, TTree
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

# declare histograms
arrBins_pT = array('d', (0., 0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3.))
arrBins_theta = array('d', (0, 10.*TMath.Pi()/180.,20.*TMath.Pi()/180.,30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160.*TMath.Pi()/180., 170.*TMath.Pi()/180.,TMath.Pi()))

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)

h_track_d0 = TH1D('track_d0', 'track_d0', 100, -5., 5.)
h_track_z0 = TH1D('track_z0', 'track_z0', 100, -20., 20.)
h_track_Nhits = TH1D('track_Nhits', 'track_Nhits', 30, 0, 30)
h_track_pT = TH1D('track_pT', 'track_pT', len(arrBins_pT)-1, arrBins_pT)
h_track_theta = TH1D('track_theta', 'track_theta',
                     len(arrBins_theta)-1, arrBins_theta)

h_bg_d0 = TH1D('bg_d0', 'bg_d0', 100, -5., 5.)
h_bg_z0 = TH1D('bg_z0', 'bg_z0', 100, -20., 20.)
h_bg_Nhits = TH1D('bg_Nhits', 'bg_Nhits', 30, 0, 30)
h_bg_Nholes = TH1D('bg_Nholes', 'bg_Nholes', 30, 0, 30)
h_bg_pT = TH1D('bg_pT', 'bg_pT', len(arrBins_pT)-1, arrBins_pT)
h_bg_theta = TH1D('bg_theta', 'bg_theta',
                     len(arrBins_theta)-1, arrBins_theta)

h_truth_theta_pT = TH2D(
    'h_truth_theta_pT', 'h_truth_theta_pT', len(arrBins_theta)-1, arrBins_theta, len(arrBins_pT)-1, arrBins_pT)
h_track_theta_pT = TH2D(
    'h_track_theta_pT', 'h_track_theta_pT', len(arrBins_theta)-1, arrBins_theta, len(arrBins_pT)-1, arrBins_pT)

h_truth_theta_pT_highR = TH2D(
    'h_truth_theta_pT_highR', 'h_truth_theta_pT_highR', len(arrBins_theta)-1, arrBins_theta, len(arrBins_pT)-1, arrBins_pT)
h_track_theta_pT_highR = TH2D(
    'h_track_theta_pT_highR', 'h_track_theta_pT_highR', len(arrBins_theta)-1, arrBins_theta, len(arrBins_pT)-1, arrBins_pT)


histos_list = [h_truth_pT, h_truth_theta, h_track_d0, h_track_z0, h_track_Nhits, h_track_pT, h_track_theta, h_bg_d0, h_bg_z0, h_bg_Nhits, h_bg_pT, h_bg_theta, h_bg_Nholes, h_truth_theta_pT, h_track_theta_pT, h_truth_theta_pT_highR, h_track_theta_pT_highR]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
BIB_tree = TTree("BIB_tree", "BIB_tree")
pt_BIB = array('d', [0])
phi_BIB = array('d', [0])
theta_BIB = array('d', [0])
d0_BIB = array('d', [0])
z0_BIB = array('d', [0])
Nhits_BIB = array('i', [0])
Nholes_BIB = array('i', [0])
BIB_tree.Branch("pT",  pt_BIB,  'var/D')
BIB_tree.Branch("phi", phi_BIB, 'var/D')
BIB_tree.Branch("theta", theta_BIB, 'var/D')
BIB_tree.Branch("d0", d0_BIB, 'var/D')
BIB_tree.Branch("z0", z0_BIB, 'var/D')
BIB_tree.Branch("Nhits", Nhits_BIB, 'var/I')
BIB_tree.Branch("Nholes", Nholes_BIB, 'var/I')

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

Bfield = 3.57  # T

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print(" ")
        print("Processing event " + str(ievt))

    #print("MCtruth pions (pT, theta)")
    mcpCollection = event.getCollection('MCParticle')
    relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
    relation = UTIL.LCRelationNavigator(relationCollection)

    for mcp in mcpCollection:
        if fabs(mcp.getPDG())==211:
            dp3 = mcp.getMomentum()
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

            #print(tlv.Perp(), tlv.Theta())
            h_truth_pT.Fill(tlv.Perp())
            h_truth_theta.Fill(tlv.Theta())
            h_truth_theta_pT.Fill(tlv.Theta(),tlv.Perp())

            vpos = mcp.getVertex()
            r = sqrt(vpos[0]*vpos[0]+vpos[1]*vpos[1])
            if r > 3.1:
                h_truth_theta_pT_highR.Fill(tlv.Theta(),tlv.Perp())

            tracks = relation.getRelatedToObjects(mcp)
            for track in tracks:
                h_track_d0.Fill(track.getD0())
                h_track_z0.Fill(track.getZ0())
                h_track_Nhits.Fill(len(track.getTrackerHits()))
                h_track_pT.Fill(tlv.Perp())
                h_track_theta.Fill(tlv.Theta())
                h_track_theta_pT.Fill(tlv.Theta(),tlv.Perp())
                if r > 3.1:
                    h_track_theta_pT_highR.Fill(tlv.Theta(),tlv.Perp())


    #print("PFOs (PDG, pT, theta)")
    #pfoCollection = event.getCollection('PandoraPFOs')
    #for pfo in pfoCollection:
    #    dp3 = pfo.getMomentum()
    #    tlv = TLorentzVector()
    #    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getMass())
    #
    #    print(pfo.getType(), tlv.Perp(), tlv.Theta())

    #print("Tracks (pT, theta, d0, z0, Nhits)")
    tracks = event.getCollection('SiTracks_Refitted')
    for itrack, track in enumerate(tracks):
        try:
            mcp = relation.getRelatedFromObjects(track)[0]
            #print(0.3 * Bfield / fabs(track.getOmega() * 1000.), TMath.Pi()/2-atan(track.getTanLambda()), track.getD0(), track.getZ0(), len(track.getTrackerHits()))
            #print(" matched to ", mcp.getPDG())
        except:
            h_bg_d0.Fill(track.getD0())
            h_bg_z0.Fill(track.getZ0())
            h_bg_Nhits.Fill(len(track.getTrackerHits()))
            h_bg_Nholes.Fill(track.getNholes())
            h_bg_pT.Fill(0.3 * Bfield / fabs(track.getOmega() * 1000.))
            h_bg_theta.Fill(TMath.Pi()/2-atan(track.getTanLambda()))

            pt_BIB[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            phi_BIB[0] = track.getPhi()
            theta_BIB[0] = TMath.Pi()/2-atan(track.getTanLambda())
            d0_BIB[0] = track.getD0()
            z0_BIB[0] = track.getZ0()
            Nhits_BIB[0] = len(track.getTrackerHits())
            Nholes_BIB[0] = track.getNholes()
            BIB_tree.Fill()

reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_softpion.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
BIB_tree.Write()
output_file.Close()
