from pyLCIO import IOIMPL, UTIL, EVENT
from ROOT import TLorentzVector
from optparse import OptionParser
from math import *

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
(options, args) = parser.parse_args()

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    mcpCollection = event.getCollection('MCParticle')
    printO = False

    print(" ")

    for mcp in mcpCollection:
        if mcp.getPDG() == 130 :
            if mcp.getGeneratorStatus() == 1:
                if len(mcp.getParents()) == 0:
                    print("Processing event ", ievent)

                    dp3 = mcp.getMomentum()
                    tlv = TLorentzVector()
                    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

                    print("MC part", mcp.getPDG(), mcp.getEnergy(), tlv.Theta(), tlv.Phi())
                    
                    pfoCollection = event.getCollection('PandoraPFOs')
                    print("N pfo ", len(pfoCollection))

                    for pfo in pfoCollection:
                        dp3 = pfo.getMomentum()
                        tlv_pfo = TLorentzVector()
                        tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                        if pfo.getEnergy()>1.:
                            print(pfo.getType(), pfo.getEnergy(), tlv_pfo.Theta(), tlv_pfo.Phi())

                    jetCollection = event.getCollection('JetOut')
            
                    print("N jet ", len(jetCollection))
            
                    for jet in jetCollection:
                        dp3 = jet.getMomentum()
                        tlv_jet = TLorentzVector()
                        tlv_jet.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())
            
                        if jet.getEnergy()>1.:
                            print("Jet ", jet.getEnergy(), tlv_jet.Theta(), tlv_jet.Phi())

                    clusterCollection = event.getCollection('PandoraClusters')
                    for cluster in clusterCollection:

                        px = cluster.getEnergy()*sin(cluster.getITheta())*cos(cluster.getIPhi())
                        py = cluster.getEnergy()*sin(cluster.getITheta())*sin(cluster.getIPhi())
                        pz = cluster.getEnergy()*cos(cluster.getITheta())

                        tlv_clus = TLorentzVector()
                        tlv_clus.SetPxPyPzE(px, py, pz, cluster.getEnergy())

                        dR = tlv_clus.DeltaR(tlv)
                        print("Cluster", cluster.getEnergy(), cluster.getITheta(), cluster.getIPhi())

    HCAL_sumE = 0.
    hcal_coll = ['HcalBarrelsCollectionRec','HcalEndcapsCollectionRec']

    for icoll, coll in enumerate(hcal_coll):

        try:
            HCALhitCollection = event.getCollection(coll)

            for hit in HCALhitCollection:
                HCAL_sumE = HCAL_sumE + hit.getEnergy()
        except:
            print("No", coll, "found")

    print("SumE", HCAL_sumE)

reader.close()
