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
parser.add_option('-o', '--outFile', help='--outFile ntup_neutrons.root',
                  type=str, default='ntup_neutrons.root')
(options, args) = parser.parse_args()

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000., 2500., 5000.))

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta, 
               h_matched_theta, h_matched_E]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
neutron_tree = TTree("neutron_tree", "neutron_tree")
E_reco = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
theta_truth = array('d', [0])
neutron_tree.Branch("E",  E_reco,  'var/D')
neutron_tree.Branch("theta", theta, 'var/D')
neutron_tree.Branch("E_truth",  E_truth,  'var/D')
neutron_tree.Branch("theta_truth", theta_truth, 'var/D')

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

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())

        h_truth_E.Fill(mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        E_truth[0] = mcpCollection[0].getEnergy()
        theta_truth[0] = tlv.Theta()

        # Fill the hit-level histos and aggregated energy
        ECAL_sumE = 0.
        ecal_coll = ['EcalBarrelCollectionRec','EcalEndcapCollectionRec']
        he = 0.
        hedr = 99999.

        for icoll, coll in enumerate(ecal_coll):

            try:
                ECALhitCollection = event.getCollection(coll)
                for hit in ECALhitCollection:
                    hitpos = TVector3(hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])
                    deltaR =  tlv.Angle(hitpos)
                    if deltaR<0.4:
                        ECAL_sumE = ECAL_sumE + hit.getEnergy()
            except:
                print("No", coll, "found")

        HCAL_sumE = 0.
        hcal_coll = ['HcalBarrelCollectionRec','HcalEndcapCollectionRec']

        for icoll, coll in enumerate(hcal_coll):

            try:
                HCALhitCollection = event.getCollection(coll)
                for hit in HCALhitCollection:
                    hitpos = TVector3(hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])
                    deltaR =  tlv.Angle(hitpos)
                    if deltaR<0.4:
                        HCAL_sumE = HCAL_sumE + hit.getEnergy()
            except:
                print("No", coll, "found")
        
        E_reco[0] = ECAL_sumE+HCAL_sumE
        theta[0] = tlv.Theta()

        print("MC energy: ", tlv.E())
        print("SumCAL: ", E_reco[0])

        neutron_tree.Fill()

    reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
neutron_tree.Write()
output_file.Close()
