from pyLCIO import IOIMPL
from ROOT import TH1D, TFile, TLorentzVector, TProfile2D, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_jets.root',
                  type=str, default='ntup_jets.root')
parser.add_option('-p', '--photonCalibFile', help='--photonCalibFile ResponseMap_reco_photons_v5_noBIB.root',
                  type=str, default='ResponseMap_reco_photons_v5_noBIB.root')
parser.add_option('-n', '--neutronCalibFile', help='--neutronCalibFile ResponseMap_reco_neutrons_v5_noBIB.root',
                  type=str, default='ResponseMap_reco_neutrons_v5_noBIB.root')
parser.add_option('-j', '--jetCalibFile', help='--jetCalibFile ResponseMap_reco_jets_v5_noBIB.root',
                  type=str, default='ResponseMap_reco_jets_v5_noBIB.root')
(options, args) = parser.parse_args()

arrBins_pT = array('d', (0., 5., 10., 15., 20., 25., 30., 35.,
                         40., 45., 50., 75., 100., 250., 500.))

# declare histograms
h_mjj = TH1D('mjj', 'mjj', 150, 0, 150)
h_mjj_uncorrected = TH1D('mjj_uncorrected', 'mjj_uncorrected', 150, 0, 150)
h_truth_mjj = TH1D('truth_mjj', 'truth_mjj', 150, 0, 150)

h_correction_visible = TProfile2D('h_vis', 'h_vis',
                          20, 0, TMath.Pi(), len(arrBins_pT)-1, arrBins_pT,
                          's')
h_correction_equalbins = TProfile2D('h_equalbins', 'h_equalbins',
                          20, 0, TMath.Pi(), 20, 0., 200., 0., 2.,
                          's')

# Histo list for writing to outputs
histos_list = [h_mjj, h_mjj_uncorrected, h_truth_mjj, h_correction_visible, h_correction_equalbins]

for histo in histos_list:
    histo.SetDirectory(0)

#########################

tree = TTree("jets_tree", "jets_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt = array('d', [0])
theta = array('d', [0])
phi = array('d', [0])
Evec = array('d', [0])
pt_nocalib = array('d', [0])
theta_nocalib = array('d', [0])
phi_nocalib = array('d', [0])
Evec_nocalib = array('d', [0])
pt_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
Evec_truth = array('d', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt,  'var/D')
tree.Branch("theta", theta, 'var/D')
tree.Branch("phi", phi, 'var/D')
tree.Branch("E", Evec, 'var/D')
tree.Branch("pT_nocalib",  pt_nocalib,  'var/D')
tree.Branch("theta_nocalib", theta_nocalib, 'var/D')
tree.Branch("phi_nocalib", phi_nocalib, 'var/D')
tree.Branch("E_nocalib", Evec_nocalib, 'var/D')
tree.Branch("pT_truth",  pt_truth,  'var/D')
tree.Branch("theta_truth", theta_truth, 'var/D')
tree.Branch("phi_truth", phi_truth, 'var/D')
tree.Branch("E_truth", Evec_truth, 'var/D')

dijet_tree = TTree("dijet_tree", "dijet_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt1 = array('d', [0])
theta1 = array('d', [0])
phi1 = array('d', [0])
mass1 = array('d', [0])
n_charged1 = array('i', [0])
pt2 = array('d', [0])
theta2 = array('d', [0])
phi2 = array('d', [0])
mass2 = array('d', [0])
n_charged2 = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
dijet_tree.Branch("pT1",  pt1,  'var/D')
dijet_tree.Branch("theta1", theta1, 'var/D')
dijet_tree.Branch("phi1", phi1, 'var/D')
dijet_tree.Branch("mass1", mass1, 'var/D')
dijet_tree.Branch("nCharged1", n_charged1, 'var/I')
dijet_tree.Branch("pT2",  pt2,  'var/D')
dijet_tree.Branch("theta2", theta2, 'var/D')
dijet_tree.Branch("phi2", phi2, 'var/D')
dijet_tree.Branch("mass2", mass2, 'var/D')
dijet_tree.Branch("nCharged2", n_charged2, 'var/I')

calibFile_photons = TFile(options.photonCalibFile, "READ")
calibMap_photons = calibFile_photons.Get("calib_2d")
calibFile_neutrons = TFile(options.neutronCalibFile, "READ")
calibMap_neutrons = calibFile_neutrons.Get("calib_2d")
calibFile_jets = TFile(options.jetCalibFile, "READ")
calibMap_jets = calibFile_jets.Get("h_vis")

def get_calibrated_tlv(particle):
    tlv = TLorentzVector()
    E = particle.getEnergy()
    dp3 = particle.getMomentum()
    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], E)

    theta = tlv.Theta()
    E_cap = E
    correction = 1.0

    if particle.getType() == 22:
        if E_cap > 1000.:
            E_cap = 999.
        correction = calibMap_photons.GetBinContent(calibMap_photons.FindBin(theta, E_cap))
    elif particle.getType() == 2112:
        if E_cap > 250.:
            E_cap = 249.
        correction = calibMap_neutrons.GetBinContent(calibMap_neutrons.FindBin(theta, E_cap))
    else:
        correction = 1.0

    if correction < 0.01:
        correction = 1.0
    if correction > 10.0:
        correction = 1.0

    #compute the transverse momentum from the corrected energy
    tlv.SetPxPyPzE(tlv.Px() * correction, tlv.Py() * correction, tlv.Pz() * correction, tlv.E() * correction)

    return tlv

def get_calibrated_jet(jet_tlv):

    pT = jet_tlv.Perp()
    theta = jet_tlv.Theta()
    correction = 1.0
    if pT > 500.:
        pT = 499.

    correction = calibMap_jets.GetBinContent(calibMap_jets.FindBin(theta, pT))
    if correction < 0.01:
        correction = 1.0

    tlv = TLorentzVector()
    tlv.SetPxPyPzE(jet_tlv.Px() * correction, jet_tlv.Py() * correction, jet_tlv.Pz() * correction, jet_tlv.E() * correction)

    return tlv

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))

    # get the truth jets
    truthJetCollection = event.getCollection('TruthValenciaJetOut')

    if len(truthJetCollection) != 2:
        print("Skipping event " + str(ievt) + " because it does not have exactly two truth jets")
        continue

    tlv_truthJet1 = TLorentzVector()
    dp3 = truthJetCollection[0].getMomentum()
    tlv_truthJet1.SetPxPyPzE(dp3[0], dp3[1], dp3[2], truthJetCollection[0].getEnergy())

    tlv_truthJet2 = TLorentzVector()
    dp3 = truthJetCollection[1].getMomentum()
    tlv_truthJet2.SetPxPyPzE(dp3[0], dp3[1], dp3[2], truthJetCollection[1].getEnergy())

    h_truth_mjj.Fill((tlv_truthJet1 + tlv_truthJet2).M())

    ij1 = -1
    ij2 = -1

    jetCollection = event.getCollection('ValenciaJetOut')
    if len(jetCollection) < 2:
        print("Skipping event " + str(ievt) + " because it does not have at least two jets")
        continue

    deltaR_j = 99
    for ijet, jet in enumerate(jetCollection):
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet1)

        if deltaR < deltaR_j:
            ij1 = ijet
            deltaR_j = deltaR

    deltaR_j = 99
    for ijet, jet in enumerate(jetCollection):
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        deltaR = tlv.DeltaR(tlv_truthJet2)

        if ijet == ij1:
            continue

        if deltaR < deltaR_j:
            ij2 = ijet
            deltaR_j = deltaR
    
    tlv_j1 = TLorentzVector()
    tlv_j2 = TLorentzVector()

    n_charged_constituent1 = 0
    n_charged_constituent2 = 0

    jet1 = jetCollection[ij1]
    for constituent in jet1.getParticles():
        tlv_j1 += get_calibrated_tlv(constituent)
        if fabs(constituent.getCharge()) > 0.:
            n_charged_constituent1 += 1

    jet2 = jetCollection[ij2]
    for constituent in jet2.getParticles():
        tlv_j2 += get_calibrated_tlv(constituent)
        if fabs(constituent.getCharge()) > 0.:
            n_charged_constituent2 += 1

    h_correction_visible.Fill(tlv_j1.Theta(), tlv_j1.Perp(), tlv_truthJet1.Perp() / tlv_j1.Perp())
    h_correction_visible.Fill(tlv_j2.Theta(), tlv_j2.Perp(), tlv_truthJet2.Perp() / tlv_j2.Perp())

    h_correction_equalbins.Fill(tlv_truthJet1.Theta(), tlv_truthJet1.Perp(), tlv_truthJet1.Perp() / tlv_j1.Perp())
    h_correction_equalbins.Fill(tlv_truthJet2.Theta(), tlv_truthJet2.Perp(), tlv_truthJet2.Perp() / tlv_j2.Perp())

    h_mjj_uncorrected.Fill((tlv_j1 + tlv_j2).M())
    h_mjj.Fill((get_calibrated_jet(tlv_j1) + get_calibrated_jet(tlv_j2)).M())

    # fill the tree for first jet
    pt[0] = get_calibrated_jet(tlv_j1).Perp()
    theta[0] = get_calibrated_jet(tlv_j1).Theta()
    phi[0] = get_calibrated_jet(tlv_j1).Phi()
    Evec[0] = get_calibrated_jet(tlv_j1).E()
    pt_nocalib[0] = tlv_j1.Perp()
    theta_nocalib[0] = tlv_j1.Theta()
    phi_nocalib[0] = tlv_j1.Phi()
    Evec_nocalib[0] = tlv_j1.E()
    pt_truth[0] = tlv_truthJet1.Perp()
    theta_truth[0] = tlv_truthJet1.Theta()
    phi_truth[0] = tlv_truthJet1.Phi()
    Evec_truth[0] = tlv_truthJet1.E()
    tree.Fill()
    # fill the tree for second jet
    pt[0] = get_calibrated_jet(tlv_j2).Perp()
    theta[0] = get_calibrated_jet(tlv_j2).Theta()
    phi[0] = get_calibrated_jet(tlv_j2).Phi()
    Evec[0] = get_calibrated_jet(tlv_j2).E()
    pt_nocalib[0] = tlv_j2.Perp()
    theta_nocalib[0] = tlv_j2.Theta()
    phi_nocalib[0] = tlv_j2.Phi()
    Evec_nocalib[0] = tlv_j2.E()
    pt_truth[0] = tlv_truthJet2.Perp()
    theta_truth[0] = tlv_truthJet2.Theta()
    phi_truth[0] = tlv_truthJet2.Phi()
    Evec_truth[0] = tlv_truthJet2.E()
    tree.Fill()

    # fill the dijet tree
    pt1[0] = tlv_j1.Perp()
    theta1[0] = tlv_j1.Theta()
    phi1[0] = tlv_j1.Phi()
    mass1[0] = tlv_j1.M()
    n_charged1[0] = n_charged_constituent1
    pt2[0] = tlv_j2.Perp()
    theta2[0] = tlv_j2.Theta()
    phi2[0] = tlv_j2.Phi()
    mass2[0] = tlv_j2.M()
    n_charged2[0] = n_charged_constituent2
    dijet_tree.Fill()

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
tree.Write()
dijet_tree.Write()
output_file.Close()
