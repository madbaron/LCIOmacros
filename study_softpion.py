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
parser.add_option('-o', '--outFile', help='--outFile histos_softpion.root',
                  type=str, default='histos_softpion.root')
parser.add_option('-t', '--doTree', dest='doTree', help='Fill BIB tree', action='store_true', default=False)
parser.add_option('-g', '--doHoles', dest='doHoles', help='Check for holes', action='store_true', default=False)
(options, args) = parser.parse_args()

#Bfield = 5  # T
#if "3TeV" in options.inFile:
#    Bfield = 3.57  # T
Bfield = 3.57  # T

# declare histograms
#arrBins_pT = array('d', (0., 0.2, 0.5, 0.75, 1., 1.5, 2., 3., 4., 5.))
arrBins_pT = array('d', (0., 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1., 2., 3., 4., 5.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180.))

h_Cutflow = TH1D("CutFlowHist", "CutFlowHist", 30, 0, 30)

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)
h_truth_d0 = TH1D('truth_d0', 'truth_d0', 100, -5., 5.)
h_truth_z0 = TH1D('truth_z0', 'truth_z0', 100, -20., 20.)

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


histos_list = [h_truth_pT, h_truth_theta, h_truth_d0, h_truth_z0, h_track_d0, h_track_z0, h_track_Nhits, h_track_pT, h_track_theta, h_bg_d0, h_bg_z0, h_bg_Nhits, h_bg_pT, h_bg_theta, h_bg_Nholes, h_truth_theta_pT, h_track_theta_pT, h_truth_theta_pT_highR, h_track_theta_pT_highR]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
SIG_tree = TTree("SIG_tree", "SIG_tree")
pt_SIG = array('d', [0])
phi_SIG = array('d', [0])
theta_SIG = array('d', [0])
d0_SIG = array('d', [0])
z0_SIG = array('d', [0])
Nhits_SIG = array('i', [0])
Nholes_SIG = array('i', [0])
charge_SIG = array('i', [0])
chi2_SIG = array('d', [0])
ndf_SIG = array('i', [0])
SIG_tree.Branch("pT",  pt_SIG,  'var/D')
SIG_tree.Branch("phi", phi_SIG, 'var/D')
SIG_tree.Branch("theta", theta_SIG, 'var/D')
SIG_tree.Branch("d0", d0_SIG, 'var/D')
SIG_tree.Branch("z0", z0_SIG, 'var/D')
SIG_tree.Branch("chi2", chi2_SIG, 'var/D')
SIG_tree.Branch("ndf", ndf_SIG, 'var/I')
SIG_tree.Branch("Nhits", Nhits_SIG, 'var/I')
SIG_tree.Branch("Nholes", Nholes_SIG, 'var/I')
SIG_tree.Branch("charge", charge_SIG, 'var/I')

BIB_tree = TTree("BIB_tree", "BIB_tree")
pt_BIB = array('d', [0])
phi_BIB = array('d', [0])
theta_BIB = array('d', [0])
d0_BIB = array('d', [0])
z0_BIB = array('d', [0])
Nhits_BIB = array('i', [0])
Nholes_BIB = array('i', [0])
charge_BIB = array('i', [0])
chi2_BIB = array('d', [0])
ndf_BIB = array('i', [0])
BIB_tree.Branch("pT",  pt_BIB,  'var/D')
BIB_tree.Branch("phi", phi_BIB, 'var/D')
BIB_tree.Branch("theta", theta_BIB, 'var/D')
BIB_tree.Branch("d0", d0_BIB, 'var/D')
BIB_tree.Branch("z0", z0_BIB, 'var/D')
BIB_tree.Branch("chi2", chi2_BIB, 'var/D')
BIB_tree.Branch("ndf", ndf_BIB, 'var/I')
BIB_tree.Branch("Nhits", Nhits_BIB, 'var/I')
BIB_tree.Branch("Nholes", Nholes_BIB, 'var/I')
BIB_tree.Branch("charge", charge_BIB, 'var/I')

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

print("Bfield ", Bfield)
if options.doTree:
    print("Filling bg track tree")

# loop over all events in the file
for ievt, event in enumerate(reader):

    decoder = UTIL.BitField64('system:5,side:-2,layer:6,module:11,sensor:8')

    #print("MCtruth pions (pT, theta)")
    mcpCollection = event.getCollection('MCParticle')
    relationCollection = event.getCollection('MCParticle_SiTracks')
    relation = UTIL.LCRelationNavigator(relationCollection)
    tracks = event.getCollection('SiTracks')

    if ievt % 100 == 0:
        print(" ")
        print("Processing event ", ievt)
        print("MCP:", len(mcpCollection))
        print("Tracks:", len(tracks))
        print("Now loop on MCP")
        
    for mcp in mcpCollection:
        if fabs(mcp.getPDG())==211 or fabs(mcp.getPDG())==13:
            dp3 = mcp.getMomentum()
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

            folded = tlv.Theta()
            if (tlv.Theta() > TMath.Pi()/2):
                folded = TMath.Pi()-tlv.Theta()
            #print(tlv.Perp(), tlv.Theta())
            h_truth_pT.Fill(tlv.Perp())
            h_truth_theta.Fill(folded)
            h_truth_theta_pT.Fill(folded,tlv.Perp())

            vpos = mcp.getVertex()
            r = sqrt(vpos[0]*vpos[0]+vpos[1]*vpos[1])

            h_truth_d0.Fill(r)
            h_truth_z0.Fill(vpos[2])

            if r > 31.:
                h_truth_theta_pT_highR.Fill(folded,tlv.Perp())

            rel_tracks = relation.getRelatedToObjects(mcp)
            for track in rel_tracks:

                hits = track.getTrackerHits()

                Nhits = len(hits)
                h_track_Nhits.Fill(Nhits)

                if Nhits>7:
                    h_track_d0.Fill(r)
                    h_track_z0.Fill(vpos[2])
                    h_track_pT.Fill(tlv.Perp())
                    h_track_theta.Fill(folded)
                    h_track_theta_pT.Fill(folded,tlv.Perp())
                    if r > 31.:
                        h_track_theta_pT_highR.Fill(folded,tlv.Perp())

                if folded > TMath.Pi()/3:            
                    pt_SIG[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                    phi_SIG[0] = track.getPhi()
                    theta_SIG[0] = TMath.Pi()/2-atan(track.getTanLambda())
                    d0_SIG[0] = track.getD0()
                    z0_SIG[0] = track.getZ0()
                    Nhits_SIG[0] = Nhits
                    charge_SIG[0] = int(copysign(1, track.getOmega()))
                    Nholes_SIG[0] = 0
                    chi2_SIG[0] = track.getChi2()
                    ndf_SIG[0] = track.getNdf()

                    if options.doHoles:

                        the_hits = []
                        for ihit, hit in enumerate(hits):
                            cellID = int(hit.getCellID0())
                            decoder.setValue(cellID)
                            layer = decoder['layer'].value()
                            detector = decoder["system"].value()

                            this_hit = [detector, layer]
                            the_hits.append(this_hit)

                        previous_layer = 0
                        previous_detector = 1
                        for hit in the_hits:

                            new_holes = 0

                            if hit[0] > previous_detector:
                                #print("det jump", hit[0], "det,layer ", previous_detector, previous_layer)
                                extra = 0
                                if hit[1] > 0:
                                    extra = hit[1]

                                if hit[0] == 3:
                                    diff = 7 - previous_layer + extra
                                elif hit[0] == 5:
                                    if previous_detector==3:
                                        diff = 2 - previous_layer + extra
                                    else:
                                        diff = 7 - previous_layer + 3 + extra
                                #print("diff", diff)
                                new_holes = diff
                                previous_detector = hit[0]
                                previous_layer = hit[1]
                            elif hit[0] < previous_detector:
                                #going backwards - everything goes here
                                pass
                            else:
                                diff = hit[1] - previous_layer
                                previous_layer = hit[1]
                                previous_detector = hit[0]
                                if diff > 1:
                                    new_holes = diff-1

                            if new_holes>0:
                                Nholes_SIG[0] = Nholes_SIG[0]+new_holes

                    #if Nholes_SIG[0]>0:
                    #    print("MEGAHOLES", Nholes_SIG[0])
                    #    print(the_hits)

                    SIG_tree.Fill()   

    #print("PFOs (PDG, pT, theta)")
    #pfoCollection = event.getCollection('PandoraPFOs')
    #for pfo in pfoCollection:
    #    dp3 = pfo.getMomentum()
    #    tlv = TLorentzVector()
    #    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getMass())
    #
    #    print(pfo.getType(), tlv.Perp(), tlv.Theta())

    if ievt % 100 == 0:
        print("Now loop on tracks")

    n_ST_BIB = 0
    lead_ST = 0
    sublead_ST = 0
    lead_charge = 0
    sublead_charge = 0

    #print("Tracks (pT, theta, d0, z0, Nhits)")
    for itrack, track in enumerate(tracks):
        #print(0.3 * Bfield / fabs(track.getOmega() * 1000.), TMath.Pi()/2-atan(track.getTanLambda()), track.getD0(), track.getZ0(), len(track.getTrackerHits()))

        pdgId = 0
        try:
            mcp = relation.getRelatedFromObjects(track)[0]
            pdgId = mcp.getPDG()
        except:
            pdgId = 0
            
        if fabs(pdgId)==13:
            pass
        else:
            hits = track.getTrackerHits()

            Nhits = len(hits)
            h_bg_Nhits.Fill(Nhits)

            pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            folded = TMath.Pi()/2-atan(track.getTanLambda())
            if (folded > TMath.Pi()/2):
                folded = TMath.Pi()-folded
            
            if Nhits>7:
                if folded > TMath.Pi()/3:
                    if pT>lead_ST:
                        sublead_ST = lead_ST
                        sublead_charge = lead_charge
                        lead_ST = pT
                        lead_charge = int(copysign(1, track.getOmega()))
                    elif pT>sublead_ST:
                        sublead_ST = pT
                        sublead_charge = int(copysign(1, track.getOmega()))
                    else:
                        pass

                    n_ST_BIB = n_ST_BIB+1

                if pT < 1.:
                    h_bg_d0.Fill(track.getD0())
                    h_bg_z0.Fill(track.getZ0())
                    h_bg_Nholes.Fill(track.getNholes())
                    h_bg_theta.Fill(TMath.Pi()/2-atan(track.getTanLambda()))

            if options.doTree:

                if folded > TMath.Pi()/3:            
                    pt_BIB[0] = pT
                    phi_BIB[0] = track.getPhi()
                    theta_BIB[0] = TMath.Pi()/2-atan(track.getTanLambda())
                    d0_BIB[0] = track.getD0()
                    z0_BIB[0] = track.getZ0()
                    Nhits_BIB[0] = Nhits
                    charge_BIB[0] = int(copysign(1, track.getOmega()))
                    chi2_BIB[0] = track.getChi2()
                    ndf_BIB[0] = track.getNdf()

                    Nholes_BIB[0] = 0

                    if options.doHoles:
                        the_hits = []
                        for ihit, hit in enumerate(hits):
                            cellID = int(hit.getCellID0())
                            decoder.setValue(cellID)
                            layer = decoder['layer'].value()
                            detector = decoder["system"].value()

                            this_hit = [detector, layer]
                            the_hits.append(this_hit)

                        previous_layer = 0
                        previous_detector = 1
                        for hit in the_hits:
                        
                            new_holes = 0

                            if hit[0] > previous_detector:
                                extra = 0
                                if hit[1] > 0:
                                    extra = hit[1]

                                if hit[0] == 3:
                                    diff = 7 - previous_layer + extra
                                elif hit[0] == 5:
                                    if previous_detector==3:
                                        diff = 2 - previous_layer + extra
                                    else:
                                        diff = 7 - previous_layer + 3 + extra
                                new_holes = diff
                                previous_detector = hit[0]
                                previous_layer = hit[1]
                            elif hit[0] < previous_detector:
                                pass
                            else:
                                diff = hit[1] - previous_layer
                                previous_layer = hit[1]
                                previous_detector = hit[0]
                                if diff > 1:
                                    new_holes = diff-1

                            if new_holes>0:
                                Nholes_BIB[0] = Nholes_BIB[0]+new_holes
    
                    BIB_tree.Fill()

    if n_ST_BIB == 2 and lead_charge*sublead_charge<0:
        h_bg_pT.Fill(lead_ST)   

    if ievt % 100 == 0:
        print("")
    h_Cutflow.Fill(0)

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_Cutflow.Write()
for histo in histos_list:
    histo.Write()
if options.doTree:
    BIB_tree.Write()
SIG_tree.Write()
output_file.Close()
