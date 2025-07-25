from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TVector3, TTree, TMath
from math import *
from optparse import OptionParser

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_bg_tracks.root',
                  type=str, default='ntup_bg_tracks.root')
(options, args) = parser.parse_args()


tree = TTree("tracks_tree", "tracks_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
d0 = array('d', [0])
z0 = array('d', [0])
chi2 = array('d', [0])
ndf = array('i', [0])
nhits = array('i', [0])
nholes = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt,  'var/D')
tree.Branch("phi", phi, 'var/D')
tree.Branch("theta", theta, 'var/D')
tree.Branch("d0", d0, 'var/D')
tree.Branch("z0", z0, 'var/D')
tree.Branch("chi2", chi2, 'var/D')
tree.Branch("ndf", ndf, 'var/I')
tree.Branch("nhits", nhits, 'var/I')
tree.Branch("nholes", nholes, 'var/I')

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    # getting tracks and relative hits
    tracks = event.getCollection("SiTracks")

    for itrack, track in enumerate(tracks):

        if itrack % 100 == 0:
            print("Processing track " + str(itrack))

        pt[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi[0] = track.getPhi()
        theta[0] = TMath.Pi()/2-atan(track.getTanLambda())
        d0[0] = track.getD0()
        z0[0] = track.getZ0()
        chi2[0] = track.getChi2()
        ndf[0] = track.getNdf()
        nhits[0] = len(track.getTrackerHits())
        nholes[0] = int(track.getdEdxError())  # BADHACK

        tree.Fill()

reader.close()

# write tree
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
output_file.Close()
