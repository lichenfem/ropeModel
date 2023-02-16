import os
import sys
import json
import numpy as np
import time
import pickle


def genCAEInAbq(settingJson):
    # loading settings
    a_file = open(settingJson, 'r')
    settings = json.load(a_file)  # settings is a dict
    a_file.close()

    # create necessary folder if required
    try:
        os.mkdir(os.path.join(settings['workdir'], 'tension'))
    except:
        pass
    try:
        os.mkdir(os.path.join(settings['workdir'], 'tension', 'debugDataFiles'))
    except:
        pass


    # create list of Material instances
    from source.classes.Material import createMaterialListsForIsoPlas
    materialList = createMaterialListsForIsoPlas(pathsToCalibratedData=settings['pathsToCalibratedData'],
                                                 materialNames=settings['materialNames'],
                                                 poissonRatios=settings['poissonRatios'],
                                                 densities=settings['densities'])

    # set loading controls
    from source.classes.StaticLoading import StaticLoading
    loading = StaticLoading(typeOfLoadingApplied='tension',
                            minInc=1e-8,
                            maxInc=0.05,
                            initialInc=1e-3,    # preset
                            targetAxialStrain=settings['targetAxialStrain'],
                            )


    # create all strands
    from source.classes.Strand import Strand    # rely on abaqus modules
    strands = [[] for iLayerInRope in range(settings['nLayersOfStrandInRope'])]    # init
    for iLayerInRope in range(settings['nLayersOfStrandInRope']):
        for iStrand in range(settings['nStrandsPerLayerInRope'][iLayerInRope]):
            #
            # read in config of strand
            nLays = settings['strandSettings'][str(iLayerInRope)]['nLays']
            nLayers = settings['strandSettings'][str(iLayerInRope)]['nLayers']
            nWiresPerLayer = np.array(settings['strandSettings'][str(iLayerInRope)]['nWiresPerLayer'])
            wireLayRadii = np.array(settings['strandSettings'][str(iLayerInRope)]['wireLayRadii'])
            dWiresExp = np.array(settings['strandSettings'][str(iLayerInRope)]['dWiresExp'])
            wireLayLength = settings['strandSettings'][str(iLayerInRope)]['wireLayLength']
            phaseAngle = np.array(settings['strandSettings'][str(iLayerInRope)]['phaseAngle'])
            strandDiameter = settings['strandSettings'][str(iLayerInRope)]['diameter']
            strandLayRadius = settings['strandSettings'][str(iLayerInRope)]['strandLayRadius']
            strandLayLength = settings['strandLayLengthsInRope'][iLayerInRope]  # inf is a float
            strandLayDirection = settings['strandLayDirections'][iLayerInRope]
            wireLayDirection = settings['strandSettings'][str(iLayerInRope)]['wireLayDirection']
            nStrandsInStrandLayer = settings['nStrandsPerLayerInRope'][iLayerInRope]
            angleOffsetOfStartingStrandInStrandLayer = settings['phaseAngleOfStrandInRope'][iLayerInRope]
            rWiresPerLayer = 0.5 * dWiresExp  # radius of wires of each layer
            #
            # construct strand object
            strand = Strand(nLays=nLays,
                            nWiresPerLayer=nWiresPerLayer,
                            wireLayRadii=wireLayRadii,
                            rWiresPerLayer=rWiresPerLayer,
                            wireLayLength=wireLayLength,
                            phaseAngle=phaseAngle,
                            strandDiameter=strandDiameter,
                            strandLayRadius=strandLayRadius,
                            strandLayLength=strandLayLength,
                            strandLayDirection=strandLayDirection,
                            wireLayDirection=wireLayDirection,
                            strandLayer=iLayerInRope,
                            strandIndexInStrandLayer=iStrand,                   # id of strand in strand layer
                            nStrandsInStrandLayer=nStrandsInStrandLayer,        # number of strand in strand layer
                            angleOffsetOfStartingStrandInStrandLayer=angleOffsetOfStartingStrandInStrandLayer,
                            )
            #
            strands[iLayerInRope].append(strand)


    # reset xCtr for all wires computed from sympy and scipy
    minWireRadius = np.inf                                  # compute minimum radius among all wires
    for iLayerInRope in range(settings['nLayersOfStrandInRope']):
        for iStrand in range(settings['nStrandsPerLayerInRope'][iLayerInRope]):
            strand = strands[iLayerInRope][iStrand]
            for iLayer in range(strand.nLayers):
                for iWire in range(strand.nWiresPerLayer[iLayer]):
                    wire = strand.wires[iLayer][iWire]
                    # load json file
                    a_file = open(os.path.join(settings['workdir'], 'jsonOfWireConfigs', wire.name + '.json'), "r")
                    wireSettings = json.load(a_file)  # output is a dict
                    a_file.close()
                    wire.xCtr = wireSettings['xCtr']
                    #
                    if wire.radius < minWireRadius:
                        minWireRadius = wire.radius


    # create SimulationSetup object
    launchJob = 0
    useBeams = 0
    elementCode = 'C3D8'
    jobName = loading.typeOfLoadingApplied
    curveName = jobName
    workDir = os.path.join(settings['workdir'], loading.typeOfLoadingApplied)   # must be absolute path
    postProcess = 0
    generateSimImages = 0
    #
    howSurfacesAreRepelled = 'FiniteSlidingContact'
    contactNormalBehaviour = 'Hard'
    normalContactConstraintEnforceMethod = 'Penalty'
    contactTangentialFrictionCoefficient = 0        # frictionless between wires
    numericalStabilization = 0
    interferenceStep = 0
    contactCtrl = 0
    requestRestart = 0
    restartFrequency = 1
    stiffnessScaleFactor = 1
    penaltyType = 'linear'
    constructGeometryMethod = 'sweep'
    dynamic = 'static'
    meshCtrl = {'meshSize': settings['meshSize']}       # set mesh density


    from source.classes.SimulationSetup import SimulationSetup
    setup = SimulationSetup(
        launchJob=launchJob,
        loading=loading,
        materialList=materialList,
        design=settings['design'],
        useBeams=useBeams,
        elementCode=elementCode,
        jobName=jobName,
        curveName=curveName,
        dirSaveImages=workDir,
        workDir=workDir,
        postProcess=postProcess,
        skipJobIfODBPresent=0,
        generateInputFile=1,
        generateImages=generateSimImages,
        nCPU=4,
        toleranceForContactDetection=0.10 * minWireRadius,     # set max gap to define contact pair
        removeContactPairsInvolvingWiresEnds=0,
        howSurfacesAreRepelled=howSurfacesAreRepelled,
        contactNormalBehaviour=contactNormalBehaviour,
        normalContactConstraintEnforceMethod=normalContactConstraintEnforceMethod,
        contactTangentialFrictionCoefficient=contactTangentialFrictionCoefficient,
        numericalStabilization=numericalStabilization,
        interferenceStep=interferenceStep,
        contactCtrl=contactCtrl,
        requestRestart=requestRestart,
        restartFrequency=restartFrequency,
        stiffnessScaleFactor=stiffnessScaleFactor,
        penaltyType=penaltyType,
        constructGeometryMethod=constructGeometryMethod,
        dynamic=dynamic,
        folderToStoreDebugDataFiles=os.path.join(settings['workdir'], loading.typeOfLoadingApplied, 'debugDataFiles'),   # folder to store files generated for debugging purposes
        useBalancedContactBetweenWires=False,                   # asymmetric master-slave contact
        **meshCtrl
    )
    picklefile = open(os.path.join(settings['workdir'], 'pickleSetup'), 'wb')  # pickle setup object for usage in postprocessing
    pickle.dump(setup, picklefile)
    picklefile.close()
    print('setup picked to: ' + os.path.join(settings['workdir'], 'pickleSetup'))


    # map material to strand layers in abaqus model
    from source.classes.MaterialData import MaterialData
    for iLayerInRope in range(settings['nLayersOfStrandInRope']):
        for iStrand in range(settings['nStrandsPerLayerInRope'][iLayerInRope]):
            strand = strands[iLayerInRope][iStrand]
            #
            materialPerWire = [[] for iLayer in range(strand.nLayers)]
            materialNamePerWire = settings['strandSettings'][str(iLayerInRope)]['materialNames']   # list of list of str
            for iLayer in range(strand.nLayers):
                for iWire in range(strand.nWiresPerLayer[iLayer]):
                    materialName = materialNamePerWire[iLayer][iWire]
                    materialPerWire[iLayer].append(materialList[settings['materialNames'].index(materialName)])
            #
            materialData = MaterialData(strand=strand, materialPerWire=materialPerWire)     # the object is not needed



    # construct model and generate inp file
    from source.classes.WireAggregate import WireAggregate
    wireAggregate = WireAggregate(strands)      # instantiate WireAggregate object
    picklefile = open(os.path.join(settings['workdir'], 'pickleWireAggregate'), 'wb')  # pickle setup object for usage in postprocessing
    pickle.dump(wireAggregate, picklefile)
    picklefile.close()
    print('wireAggregate picked to: ' + os.path.join(settings['workdir'], 'pickleWireAggregate'))

    # construct model and generate inp file
    from source.simulationLaunchers.LaunchSimulation import constructModel
    start = time.time()
    constructModel(setup, wireAggregate=wireAggregate)
    end = time.time()
    print('Timelapse of modeling compaction: ' + str(end - start) + 's.')

    return

# ---------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    # settingJson = '/home/lichen/Dropbox/ropemodels/scripts/ropeModel/configuration.json'
    settingJson = sys.argv[-1]
    genCAEInAbq(settingJson)