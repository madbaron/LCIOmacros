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

    print("Processing event ", ievent)

    try:
        relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
        relation = UTIL.LCRelationNavigator(relationCollection)

        # Look at the MC particles
        mcpCollection = event.getCollection('MCParticle')

        for mcp in mcpCollection:

            charge = mcp.getCharge()

            if fabs(charge) > 0:
                if fabs(mcp.getPDG()) == 13:
                    dp3 = mcp.getMomentum()
                    tlv = TLorentzVector()
                    tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

                    if tlv.Perp() > 0.5:
                        print("MCP pT ", tlv.Perp())

                        tracks = relation.getRelatedToObjects(mcp)
                        for track in tracks:
                            print("Track tanLambda ", track.getTanLambda())

    except:
        print("No relation collection!")

reader.close()
