from ROOT import edm4hep
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree

from podio.root_io import Reader

from optparse import OptionParser
from array import array
import itertools

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.edm4hep.root',
                  type=str, default='Output_REC.edm4hep.root')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()

reader = Reader(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader.get('events')):   

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))
        
    mcpCollection = event.get('MCParticle')
    print("Nparts ", len(mcpCollection))
    good_PID = [1,2,3,4,5,23,25]
    mom_PID = [23,25]

    for mcp in mcpCollection:
        if abs(mcp.getPDG()) in good_PID:
            if abs(mcp.getPDG()) in mom_PID:
                print(" MCP",mcp.getPDG(), mcp.getMass(), mcp.getEnergy())
            else:
                if mcp.getParents()[0].getPDG() in mom_PID:
                    print(" MCP",mcp.getPDG(), mcp.getMass(), mcp.getEnergy())

    # find the jets
    jetCollection = event.get('JetOut')
    print("Njets ", len(jetCollection))

    for ijet, jet in enumerate(jetCollection):
        print(" Jet ", ijet, jet.getEnergy())

    combi = []
    for pair in itertools.combinations(jetCollection,2):
        dp31 = pair[0].getMomentum()
        tlv1 = TLorentzVector()
        tlv1.SetPxPyPzE(dp31[0], dp31[1], dp31[2], pair[0].getEnergy())
        dp32 = pair[0].getMomentum()
        tlv2 = TLorentzVector()
        tlv2.SetPxPyPzE(dp32[0], dp32[1], dp32[2], pair[1].getEnergy())
        combi.append((tlv1+tlv2).M())
    print("Mass combinations", combi)

    # find the PFOs
    pfoCollection = event.get('PandoraPFOs')
    print("Npfo ", len(pfoCollection))

    for ipfo, pfo in enumerate(pfoCollection):
        dp3 = pfo.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

        if tlv.Perp()>20:
            print(" PFO ", ipfo, "(", pfo.getType() ,")", str(tlv.Perp()) + " " + str(len(pfo.getParticles())))

