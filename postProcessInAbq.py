import os
import sys
import pickle

# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior


def postProcessInAbq(pickleSetup, pickleWireAggregate):
    """
    postprocessing work in abaqus
    @param pickleSetup: str, path of file that pickles Setup
    @param pickleWireAggregate: str, path of file that pickles WireAggregate
    @return: none
    """
    # reload setup instance
    picklefile0 = open(pickleSetup, 'rb')
    setup = pickle.load(picklefile0)
    picklefile0.close()

    # reload WireAggregate instance
    picklefile1 = open(pickleWireAggregate, 'rb')
    wireAggregate = pickle.load(picklefile1)
    picklefile1.close()

    # extract result.npz for overall tensile response
    from source.postProcessing.extractInformationFromOdb import extractInformationFromOdb
    extractInformationFromOdb(setup=setup, wireAggregate=wireAggregate)

    # plot contours of stress, strain, displacement magnitude for Z cut, rope ends, rope axial cuts
    from source.simulationLaunchers.LaunchSimulation import generateSimulationImages
    generateSimulationImages(setup, zcut=0.5 * wireAggregate.length, maxS=0)

    # plot contact stress on the surface of ropes
    from source.postProcessing.getPostProcessingImages import (showCPRESSOnSpecifiedInstances,
                                                               showCPRESSOnAllInstancesExceptHiddenOnes,
                                                               getTangentialSlipMagnitude)

    # plot contact stress on central strand
    partInstanceNames = []
    strand = wireAggregate.strandsPerLayer[0][0]
    for iLayer in range(strand.nLayers):
        for iWire in range(strand.nWiresPerLayer[iLayer]):
            partInstanceNames.append(strand.wires[iLayer][iWire].instanceName)
    # show contact pressure contour on the surfaces of desgnated wires
    showCPRESSOnSpecifiedInstances(partInstanceName=tuple(partInstanceNames),
                                   figTitle='CentralStrand',
                                   workDir=str(setup.workDir),
                                   odbFile=str(setup.odbFile))

    # show contact pressure of rope by hiding designated wire
    partInstanceNames = []
    strand = wireAggregate.strandsPerLayer[1][0]
    for iLayer in range(strand.nLayers):
        for iWire in range(strand.nWiresPerLayer[iLayer]):
            partInstanceNames.append(strand.wires[iLayer][iWire].instanceName)
    showCPRESSOnAllInstancesExceptHiddenOnes(hiddenPartInstanceName=tuple(partInstanceNames),
                                             figTitle='removeSideStrand',
                                             workDir=str(setup.workDir),
                                             odbFile=str(setup.odbFile))

    # plot contour of slide between two wires if required
    if False:
        getTangentialSlipMagnitude(
            odbFile=str(setup.odbFile), instanceName0=strand.wires[0][0].instanceName,
            instanceName1=strand.wires[1][0].instanceName, workDir=str(setup.workDir))

    return


# ----------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys

    postProcessInAbq(pickleSetup=sys.argv[-2], pickleWireAggregate=sys.argv[-1])
