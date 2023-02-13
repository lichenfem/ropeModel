import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from source.calibrationExperimentalData.fitMaterialBehaviourFromWireExperiment import (
    nominal2trueStressStrainRelation, fitElastoplasticity, plotNominalAndTrueStressStrainCurvesAndYoungModulus,
    plotFittedAgaistTrueStressStrainCurve)

# path to folder where material properties of wires are stored
pathExp = os.path.join(
    os.environ['ropemodels'], 'scripts', 'ropeModel', 'ropeExperimentalResults', 'Sample Nr.1', '01. Wires')
# pathExp = r'/home/lichen/Dropbox/ropemodels/scripts/ropeModel/ropeExperimentalResults/Sample Nr.1/01. Wires'
pathMatrPropFolder = os.path.join(os.environ['ropemodels'], 'scripts', 'ropeModel', 'wireMatProp')  # destination to store fitted wire material properties

assert os.path.isdir(pathExp)
pathWireType1 = os.path.join(pathExp, 'd1_1.076mm_Strain_Force_Stress.xlsx')
pathWireType23 = os.path.join(pathExp, 'd2, 3_0.964mm_Strain_Force_Stress.xlsx')
pathWireType4 = os.path.join(pathExp, 'd4_0.889mm_Strain_Force_Stress.xlsx')
fileExperimentalWires = [pathWireType1, pathWireType23, pathWireType4]
DiameterWires = [1.076, 0.964, 0.889]       # diameters of wires, used to compute nominal stress if required
labels = ['d1', 'd23', 'd4']
columnIdentifiers = ['A,B,C', 'G,H,I', 'A,B,C']     # which columns of the excel files will be extracted
delimiters = [0.0075, 0.0075, 0.0075]               # manual separators between elastic and elastoplastic loading stages

npzFiles = []
for i in range(len(fileExperimentalWires)):
    fig = plt.figure()
    print("Fitting material wire: " + str(labels[i]))

    df = pd.read_excel(fileExperimentalWires[i], sheet_name='Sheet1', usecols=columnIdentifiers[i])    # load a data frame, nan may present
    headers = df.columns.to_list()
    strnForcStrs = []
    for head in headers:
        column = df[head].dropna().tolist()     # remove trailing nan if required
        strnForcStrs.append(column)
    arr = np.array(strnForcStrs).transpose()
    strn = arr[:, 0]/100.0      # original strain is in the unit of percentage
    forc = arr[:, 1]            # unit of N
    strs = arr[:, 2]            # unit of N/mm^2, MPa

    # double check force = stress * cross sectional area of the wire
    d = [np.sqrt((4*f)/sigma/np.pi) for (f, sigma) in zip(forc, strs)]
    for val in d:
        assert np.allclose(val, DiameterWires[i], rtol=0.05), 'exp: ' + str(DiameterWires[i]) + ', actual: ' + str(val)

    # make strn monotonically increasing
    strn_mono = [strn[0]]
    strs_mono = [strs[0]]
    for pos in range(1, len(strn), 1):
        if strn[pos] <= strn[pos-1]:
            continue
        else:
            strn_mono.append(strn[pos])
            strs_mono.append(strs[pos])
    # plt.plot(strn_mono, strs_mono, '--b', label='monotonic strn')  # plot nominal stress-strain curve

    # neglect decreasing stress
    spacing = 5     # interval of selecting data points
    strn_coarse = strn_mono[::spacing]
    strs_coarse = strs_mono[::spacing]
    #
    negativeSlopeMag = np.tan(5.0/180.0*np.pi)    # slope of -5 degree
    posTrunc = 0
    for pos in range(1, len(strs_mono), 1):
        slope = (strs_coarse[pos] - strs_coarse[pos-1])/(strn_coarse[pos]-strn_coarse[pos-1])
        if slope < 0 and abs(slope) > negativeSlopeMag:
            break
        posTrunc = pos
    strnTrunc = strn_coarse[:posTrunc]
    strsTrunc = strs_coarse[:posTrunc]
    # plt.plot(strnTrunc, strsTrunc, '-b', label='truncated')  # plot nominal stress-strain curve

    # transform nomnal stress-strain to true stress-strain
    strain_true, stress_true = nominal2trueStressStrainRelation(strain_nominal=strnTrunc,
                                                                stress_nominal=strsTrunc)
    ax0 = fig.add_subplot(1, 2, 1)
    ax0 = plotNominalAndTrueStressStrainCurvesAndYoungModulus(strain_nominal=strnTrunc,
                                                             stress_nominal=strsTrunc,
                                                             strain_true=strain_true,
                                                             stress_true=stress_true,
                                                             ax=ax0)    # show difference between true and nominal stress-strain relations

    E, yield_stress, plastic_strain = fitElastoplasticity(materialLabel=labels[i],
                                                          trueStress=stress_true,
                                                          trueStrain=strain_true,
                                                          elasStrainLimit=delimiters[i],       # manual limit to separate elastic loading and elastoplastic loading
                                                          strainFailureThreshold=None,
                                                          pointsBetweenSamples=1,
                                                          save2npz=os.path.join(pathMatrPropFolder, labels[i] + '.npz'))
    npzFiles.append(os.path.join(pathMatrPropFolder, labels[i] + '.npz'))
    ax1 = fig.add_subplot(1, 2, 2)
    ax1 = plotFittedAgaistTrueStressStrainCurve(E, plastic_strain, yield_stress, strain_true, stress_true, ax=ax1)

    ax0.set_title('true vs nominal stress-strain of material: ' + labels[i])
    ax1.set_title('fitted stress-strain of material: ' + labels[i])
    fig.tight_layout()  # Or equivalently,  "plt.tight_layout()"
    fig.savefig(os.path.join(pathMatrPropFolder, labels[i] + '.png'),
                dpi=300, transparent=True)
    plt.close('all')

# ----------------------------------------------------------------------------

def translate(npzFile):
    """
    this function is used to translate the npz files that are generated in the previous step into json files.
    @param npzFile: str, path of npz files
    @return: none
    """
    import json
    data = np.load(npzFile)
    # print(data.files)
    E = data['E'].tolist()
    yield_stress = data['yield_stress'].tolist()
    plastic_strain = data['plastic_strain'].tolist()
    #
    dict = {"E": E, "yield_stress": yield_stress, "plastic_strain": plastic_strain}
    jsonFile = os.path.join(os.path.dirname(npzFile), os.path.basename(npzFile).replace('.npz', '.json'))
    a_file = open(jsonFile, "w")
    json.dump(dict, a_file)
    a_file.close()
    print('Generating ' + jsonFile)
    return

for file in npzFiles:
    translate(npzFile=file)






