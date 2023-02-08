import os
import sys
import numpy as np
import json

# set work directory to store generated files from this script
workdir = os.path.join(os.environ['ropemodels'], 'scripts', 'JiangModelForKW17S')

# geometrical aspect of the model
design = "bms17S_cmp0d0pct_lf9d5"             # name of the design
nLayers = 3                             # number of layers
nWiresPerLayer = np.array([1, 8, 8])    # number of wires per layer (innermost to outmost)
dWiresExp = np.array([1.204, 0.721, 1.323])  # diameters of wires of each layer from experimental data
rWiresPerLayer = 0.5 * dWiresExp        # radius of wires of each layer
phaseAngle = np.array([0, 0, np.pi/8.0])  # phase angle of starting wire of each layer
wireLayLength = 47.2                    # !unique for all wires
thickness = 0.1                         # thickness of the basic sector
wireLayDirection = [1 for i in range(nLayers)]                    # 1 for right hand lay, -1 for left hand lay

seedSpacing = meshSize = 0.05           # resolution of discretizing elliptical wire profile
print('-'*10)
for (i, d) in enumerate(dWiresExp):
    print('Circumference of wire in layer: ' + str(i) + ' is discretized using ' + str(np.ceil(np.pi*d/seedSpacing)) + '.')
print('-'*10)

# designate orientations of L1 and L2 to get the basic sector
orientOfL1 = 0.0
orientOfL2 = np.pi/8.0

# designate material properties, designate path to folder of wire properties
PathFolderExpData = os.path.join(workdir, 'wireMatProp')
Filei = [
    "D1_1_17S.json",  # wire of 1st layer
    "D2_1_17S.json",  # wire of 2nd layer
    "D3_1_17S.json",  # wire of 3rd layer
         ]
pathsToCalibratedData = [os.path.join(PathFolderExpData, Filei[i]) for i in range(len(Filei))]  # list of str
materialNames = ["MaterialCore", "MaterialLayer1", "MaterialLayer2"]
poissonRatios = [0.3, 0.3, 0.3]
densities = [7.872 * pow(10, -9), 7.872 * pow(10, -9), 7.872 * pow(10, -9)]  # low carbon steel: 7.872 g/cm3

# designate compaction rate and tensile rate
targetRadialStrain = 0.072      # compaction
targetAxialStrain = 0.02        # tension
# -------------------------above are user input------------------

# -------------------------below are automatic operations------------------
# calculate wireLayRadii using geometrical precide model, can be a bit time-consuming!
print('Calculating wireLayRadii for strand...')
from source.utilities.Geometry import getwireLayRadiiInStrand
wireLayRadii = getwireLayRadiiInStrand(nLayers=nLayers,
                                       nWiresPerLayer=nWiresPerLayer,
                                       rWiresPerLayer=rWiresPerLayer,
                                       wireLayLength=wireLayLength,
                                       phaseAngle=phaseAngle,
                                       strandRadius=sum(dWiresExp)*2,     # can be arbitrary large
                                       writeToFile='',
                                       gap=0.)
print('Done with calculating wireLayRadii for strand!')

# ------------------------------------------------------------------------
# calculate elliptical wire profiles
print('Generating elliptical wire cross sections for strand...')
from source.classes.StrandData import StrandData
strand = StrandData(1.0,                    # nLays
                    nWiresPerLayer=nWiresPerLayer,
                    wireLayRadii=wireLayRadii,          # float
                    rWiresPerLayer=rWiresPerLayer,
                    wireLayLength=wireLayLength,        # float! unique for all wires in the strand!
                    phaseAngle=phaseAngle,
                    strandDiameter=wireLayRadii[-1]+rWiresPerLayer[-1],
                    strandLayRadius=0.0,
                    strandLayLength=np.inf,             # !inf for central strand
                    strandLayDirection=1,
                    wireLayDirection=wireLayDirection,                # lay directions of the wires, default right lay
                    onlyConstructData=1,                # not using abaqus modules
                    strandLayer=0,
                    strandIndexInStrandLayer=0,
                    nStrandsInStrandLayer=1,            # default for straight strand
                    maximumInitialPenetrationBeforeAdjustment=0
                    )
strand.projectEllipticalCrossSection(zval=0.0,
                                     seedSpacing=seedSpacing,
                                     noHelicalProj=False)       # time-consuming operation!

# write elliptical wire profiles to json file
data = {}   # empty dict
for iLayer in range(strand.nLayers):
    for iWire in range(strand.nWiresPerLayer[iLayer]):
        wire = strand.wires[iLayer][iWire]
        data[wire.name] = {'cuttingProfile': wire.cuttingProfile,
                           'projectedCenter': wire.projectedCenter}
a_file = open(os.path.join(workdir, 'ellipticalWireProfilesOfUncompactedStrand.json'), "w")    # save to file
json.dump(data, a_file)
a_file.close()

ax = strand.plotCrossSection(plotTitle=None)     # plot if needed
# ------------------------------------------------------------------------
# write all configurations to a json file
settings = {}
settings['workdir'] = workdir
settings['design'] = design
settings['nLayers'] = nLayers
settings['nWiresPerLayer'] = nWiresPerLayer.tolist()
settings['rWiresPerLayer'] = rWiresPerLayer.tolist()
settings['phaseAngle'] = phaseAngle.tolist()
settings['wireLayLength'] = wireLayLength
settings['wireLayRadii'] = wireLayRadii.tolist()
settings['wireLayDirection'] = wireLayDirection
settings['meshSize'] = meshSize
settings['orientOfL1'] = orientOfL1
settings['orientOfL2'] = orientOfL2
settings['thickness'] = thickness
settings['targetRadialStrain'] = targetRadialStrain
settings['targetAxialStrain'] = targetAxialStrain
settings['materialNames'] = materialNames
settings['poissonRatios'] = poissonRatios
settings['densities'] = densities
settings['pathsToCalibratedData'] = pathsToCalibratedData

a_file = open(os.path.join(workdir, 'configuration.json'), "w")    # save to file
json.dump(settings, a_file)
a_file.close()
print('Congfiguration of ' + design + ' saved to ' + os.path.join(workdir, 'configuration.json'))

# ------------------------------------------------------------------------
# below get truncated profiles for wires in the basic sector and designate L1/L2/L1L2 for nodes on the profiles of wires
from source.utilities.utilitiesOfReducedModel import extractBasicSectorFromEllipticalStrandSection, \
    plotProfileAndL1L2Classification
basicSectorDict, L1L2DesignateDict = extractBasicSectorFromEllipticalStrandSection(fullStrandSectionDict=data,
                                                                                   L1Orient=orientOfL1,
                                                                                   L2Orient=orientOfL2,
                                                                                   tol=1e-4)
# plot designation of L1/L2/L1L2
fig = plotProfileAndL1L2Classification(profiles=basicSectorDict, classification=L1L2DesignateDict)
fig.savefig(os.path.join(workdir, 'L1L2DesignationOfProfileNodes.png'), dpi=300, transparent=True)
print('Designation of L1/L2 of profile nodes is plotted and saved to ' +
      os.path.join(workdir, 'L1L2DesignationOfProfileNodes.png'))

# save profile and designation of L1L2 to json files
a_file = open(os.path.join(workdir, 'profileOfBasicSectorBeforeCompaction.json'), 'w')
json.dump(basicSectorDict, a_file)
a_file.close()

b_file = open(os.path.join(workdir, 'setL1L2ForProfileOfBasicSectorBeforeCompaction.json'), 'w')
json.dump(L1L2DesignateDict, b_file)
b_file.close()



