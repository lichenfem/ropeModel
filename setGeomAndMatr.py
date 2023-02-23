import os
import sys
import numpy as np
import json
import warnings
import time
from math import floor, ceil
from source.classes.StrandData import StrandData
from source.utilities.Geometry import getwireLayRadiiInStrand
from source.utilities.calLayRadiiOfHelicalStrandInRope import calLayRadiiOfHelicalStrandInRope

#------------------------------------- below are use input -----------------------------------
# set work directory to store generated files from this script
# workdir = os.path.join(os.environ['ropemodels'], 'scripts', 'ropeModel')
workdir = os.getcwd()
print('workdir set to: ' + str(workdir))

# STEP 2:
# rope level geometrical configs
design = "Rope77"               # name of the design
length = 20.0                   # length of the rope
gap = 0.0                       # artificial gap between wires set by the user

# designate material properties, designate path to folder of wire properties
PathFolderExpData = os.path.join(workdir, 'wireMatProp')
Filei = [
    "d1.json",
    "d23.json",
    "d4.json",
         ]
pathsToCalibratedData = [os.path.join(PathFolderExpData, Filei[i]) for i in range(len(Filei))]  # list of str
materialNames = ["d1", "d23", "d4"]
poissonRatios = [0.3, 0.3, 0.3]
densities = [7.872 * pow(10, -9), 7.872 * pow(10, -9), 7.872 * pow(10, -9)]  # low carbon steel: 7.872 g/cm3

# information of wires in each layer of strand
strandSettings = {}
strandSettings[0] = {'nLays': np.nan,               # number of lays, TBD (to be determined)
                     'nLayers': 2,                  # number of layers
                     'nWiresPerLayer': [1, 6],      # number of wires per layer (innermost to outmost)
                     'wireLayRadii': [0.0, 0.0],    # to be determined
                     'dWiresExp': [1.076, 0.964],   # diameters of wires of each layer from experimental data
                     'phaseAngle': [0, 0],          # phase angle of starting wire of each layer, to be adjusted if required
                     'wireLayLength': 18.10,        # !unique for all wires in the strand
                     'wireLayDirection': [1, 1],    # 1 for right hand lay, -1 for left hand lay
                     'diameter': np.nan,            # strand diameter, TBD
                     'strandLayRadius': 0.0,        # 0 for straight strand, np.nan is the strandLayRadius for the strand is not determined yet.
                     'strandLayLength': np.inf,      # np.inf for straight strand
                     'strandLayDirection': 1,       # 1 for straight or right lay, -1 for left lay
                     'strandLayer': 0,              # counter of the layer of strand in rope, 0 for central strand
                     'nStrandsInStrandLayer': 1,     # number of strand(s) in the strand layer of the rope
                     'angleOffsetOfStartingStrandInStrandLayer': 0.0, # to be determined later
                     'materialNames': [['d1'], ['d23' for i in range(6)]],  # identifiers of material for each wire

                    }   # setting of wires in strand(s) of layer 0
strandSettings[1] = {'nLays': np.nan,               # number of lays, TBD
                     'nLayers': 2,                  # number of layers
                     'nWiresPerLayer': [1, 6],
                     'wireLayRadii': [0.0, 0.0],    # to be determined
                     'dWiresExp': [0.964, 0.889],
                     'phaseAngle': [0, 0],          # init to 0, TBD, [0, 0+0.5236250000000005]
                     'wireLayLength': 22.20,
                     'wireLayDirection': [-1, -1],
                     'diameter': np.nan,                        # strand diameter, TBD
                     'strandLayRadius': np.nan,     # set to np.nan when TBD, 2.6247693695690484
                     'strandLayLength': 40.0,
                     'strandLayDirection': 1,
                     'strandLayer': 1,
                     'nStrandsInStrandLayer': 6,
                     'angleOffsetOfStartingStrandInStrandLayer': 0.0,
                     'materialNames': [['d23'], ['d4' for i in range(6)]],  # identifiers of material for each wire
                    }   # setting of wires in strand(s) of layer 1

# specify loading conditions
targetAxialStrain = 0.02        # tension strain

meshSize = 0.1     # set size of uniform mesh
# prompt info about how many elements will be used to mesh the circumferences of different wires
for iLayerInRope in range(len(strandSettings.keys())):
    for iLayer in range(strandSettings[iLayerInRope]['nLayers']):
        nelem = ceil(np.pi * strandSettings[iLayerInRope]['dWiresExp'][iLayer] / meshSize)
        print('strand layer: ' + str(iLayerInRope) + ', wire layer: ' + str(iLayer) + ': '
              + str(nelem) + ' elements on wire circumference.')

#------------------------------------- above are use input -----------------------------------

#------------------------------------- below are automatic processing -----------------------------------
nLayersOfStrandInRope = len(strandSettings.keys())           # int, number of layers of strands in the rope
strandLayDirections = [strandSettings[i]['strandLayDirection'] for i in range(nLayersOfStrandInRope)]
# list, strandLayDirections[i]=1/-1 means strands in layer i are right/left hand lay
# note that in rope, it is possible for neighbouring layers of strands to have opposite
# lay direction as this is done to make the rope torsion resistent.
nStrandsPerLayerInRope = [strandSettings[i]['nStrandsInStrandLayer'] for i in range(nLayersOfStrandInRope)]    # list of int, number of strands in each layer
strandLayLengthsInRope = [strandSettings[i]['strandLayLength'] for i in range(nLayersOfStrandInRope)]   # list of float, means the lay length of strands in each layer of the rope, inf means straight for central strand


# compute number of lays for strands in each layer of the rope
for iLayerInRope in range(nLayersOfStrandInRope):
    strandLayLength = strandLayLengthsInRope[iLayerInRope]
    if strandLayLength == np.inf:   # in case of central strand
        wireLayLength = strandSettings[iLayerInRope]['wireLayLength']
        nLays = length/wireLayLength
    else:
        nLays = length/strandLayLength
    strandSettings[iLayerInRope]['nLays'] = nLays   # store strand info


# compute wireLayRadii for strand(s) in each layer
for iLayerInRope in range(nLayersOfStrandInRope):
    nLayers = strandSettings[iLayerInRope]['nLayers']
    nWiresPerLayer = np.array(strandSettings[iLayerInRope]['nWiresPerLayer'])
    dWiresExp = np.array(strandSettings[iLayerInRope]['dWiresExp'])
    phaseAngle = np.array(strandSettings[iLayerInRope]['phaseAngle'])
    wireLayLength = strandSettings[iLayerInRope]['wireLayLength']
    wireLayDirection = strandSettings[iLayerInRope]['wireLayDirection']
    rWiresPerLayer = 0.5 * dWiresExp                                    # radius of wires of each layer

    # compute the wireLayRadii for wires in the straight central strand
    wireLayRadii = getwireLayRadiiInStrand(nLayers=nLayers,
                                           nWiresPerLayer=nWiresPerLayer,
                                           rWiresPerLayer=rWiresPerLayer,
                                           wireLayLength=wireLayLength,
                                           phaseAngle=phaseAngle,
                                           strandRadius=sum(rWiresPerLayer) * 2,  # can be arbitrary large
                                           writeToFile='',
                                           gap=gap)
    strandSettings[iLayerInRope]['wireLayRadii'] = wireLayRadii.tolist()                    # store info
    strandSettings[iLayerInRope]['diameter'] = float((wireLayRadii[-1]+rWiresPerLayer[-1])*2.0)    # store strand diameter


# compute strand lay radii as well as the phase angles angleOffsetOfStartingStrandInStrandLayer
if True:
    start = time.time()
    strandLayRadii, strandSettings = calLayRadiiOfHelicalStrandInRope(strandSettings, gap=gap)
    end = time.time()
    print('timelapse for computing strand lay radius = ' + str(end-start) + ' secs.')
else:
    warnings.warn('strandLayRadii provided instead of ccomputed!\n', RuntimeWarning)
    strandLayRadii = [0.0, 2.8666767803099598]
    strandSettings[1]['strandLayRadius'] = 2.8666767803099598
    strandSettings[1]['phaseAngle'] = [0.0, 0.1395625000000001]
phaseAngleOfStrandInRope = [strandSettings[i]['angleOffsetOfStartingStrandInStrandLayer'] for i in range(nLayersOfStrandInRope)]



# creates all strands
strands = [[] for iLayerInRope in range(nLayersOfStrandInRope)]    # init
for iLayerInRope in range(nLayersOfStrandInRope):
    for iStrand in range(nStrandsPerLayerInRope[iLayerInRope]):
        nLays = strandSettings[iLayerInRope]['nLays']
        nWiresPerLayer = np.array(strandSettings[iLayerInRope]['nWiresPerLayer'])
        wireLayRadii = np.array(strandSettings[iLayerInRope]['wireLayRadii'])
        dWiresExp = np.array(strandSettings[iLayerInRope]['dWiresExp'])
        phaseAngle = np.array(strandSettings[iLayerInRope]['phaseAngle'])
        wireLayLength = strandSettings[iLayerInRope]['wireLayLength']
        wireLayDirection = strandSettings[iLayerInRope]['wireLayDirection']
        rWiresPerLayer = 0.5 * dWiresExp                                    # radius of wires of each layer
        strandDiameter = strandSettings[iLayerInRope]['diameter']
        wireLayDirection = strandSettings[iLayerInRope]['wireLayDirection']

        strand = StrandData(nLays=nLays,
                            nWiresPerLayer=nWiresPerLayer,
                            wireLayRadii=wireLayRadii,
                            rWiresPerLayer=rWiresPerLayer,
                            wireLayLength=wireLayLength,
                            phaseAngle=phaseAngle,
                            strandDiameter=strandDiameter,
                            strandLayRadius=strandLayRadii[iLayerInRope],
                            strandLayLength=strandLayLengthsInRope[iLayerInRope],
                            strandLayDirection=strandLayDirections[iLayerInRope],
                            wireLayDirection=wireLayDirection,
                            onlyConstructData=1,                # not using abaqus modules
                            strandLayer=int(iLayerInRope),
                            strandIndexInStrandLayer=int(iStrand),
                            nStrandsInStrandLayer=nStrandsPerLayerInRope[iLayerInRope],
                            angleOffsetOfStartingStrandInStrandLayer=phaseAngleOfStrandInRope[iLayerInRope],
                            maximumInitialPenetrationBeforeAdjustment=0
                            )                                   # construct strand object

        strands[iLayerInRope].append(strand)
        # strand.projectEllipticalCrossSection(zval=0.0, seedSpacing=0.1)     # calculate elliptical profiles of strand


# plot the centroid lines of the wires in the rope if required
if True:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for iLayerInRope in range(len(strands)):
        for strand in strands[iLayerInRope]:
            # ax = strand.plotCrossSection(ax=ax, CircWireSec=True, plotCentroid=True)
            strand.plot(plotLocalBasis=0, ax=ax)
    fig.savefig(os.path.join(workdir, '', 'centroidLinesOfRope.png'), dpi=300, transparent=True)
    plt.close('all')

# plot the rope cross section at z=0
if True:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for iLayerInRope in range(len(strands)):
        for strand in strands[iLayerInRope]:
            # ax = strand.plotCrossSection(ax=ax, CircWireSec=True, plotCentroid=True)
            strand.plotCrossSection(plotTitle=None,
                                    path='',
                                    ax=ax,
                                    color='blue',
                                    CircWireSec=True,       # approximate wire profile using circle
                                    plotCentroid=True       # mark centroid
                                    )
    fig.savefig(os.path.join(workdir, 'crossSectionOfRope.png'), dpi=300, transparent=True)
    plt.close('all')

# write all wire geometry into json files, use to correct Xctr in Wire object
# (this circumvents loading sympy, scipy in abaqus python environment)
folder = os.path.join(workdir, 'jsonOfWireConfigs')
for iLayerInRope in range(len(strands)):
    for iStrand in range(len(strands[iLayerInRope])):
        strand = strands[iLayerInRope][iStrand]
        for iLayer in range(strand.nLayers):
            for iWire in range(strand.nWiresPerLayer[iLayer]):
                wire = strand.wires[iLayer][iWire]
                wire.outputConstructorArgsToJson(jsonFile=os.path.join(folder, wire.name + '.json'))


# write all strand geometry into json files, to be loaded in cae model
folder = os.path.join(workdir, 'jsonOfWireConfigs')
for iLayerInRope in range(len(strands)):
    for iStrand in range(len(strands[iLayerInRope])):
        strand = strands[iLayerInRope][iStrand]
        strand.outputConstructorArgsToJson(jsonFile=os.path.join(folder, strand.name + '.json'))


# write all configurations to a json file
settings = {}
settings['workdir'] = workdir
settings['design'] = design
settings['length'] = length
settings['gap'] = gap
settings['strandLayDirections'] = strandLayDirections
settings['nLayersOfStrandInRope'] = nLayersOfStrandInRope
settings['nStrandsPerLayerInRope'] = nStrandsPerLayerInRope
settings['strandLayLengthsInRope'] = strandLayLengthsInRope
settings['phaseAngleOfStrandInRope'] = phaseAngleOfStrandInRope
settings['strandSettings'] = strandSettings     # dict of dict
settings['meshSize'] = meshSize
settings['targetAxialStrain'] = targetAxialStrain
settings['materialNames'] = materialNames
settings['poissonRatios'] = poissonRatios
settings['densities'] = densities
settings['pathsToCalibratedData'] = pathsToCalibratedData

# save to configuration.json
a_file = open(os.path.join(workdir, 'configuration.json'), "w")    # save to file
json.dump(settings, a_file)
a_file.close()
print('Congfiguration of ' + design + ' saved to ' + os.path.join(workdir, 'configuration.json'))