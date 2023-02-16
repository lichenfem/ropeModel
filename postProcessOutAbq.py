import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plotFU(pathNumResult, pathExpResult, pngFile):
    """

    @param pathNumResult: str, path of result.npz file
    @param pathExpResult: str, path of excel file containing the experimental results
    @param pngFile: str, path to save the generated png figure
    @return:
    """

    # load numerical prediction
    dataNum = np.load(pathNumResult)
    # dataNum.files   # ['RM', 'RF', 't', 'target']


    # load experimental results
    df = pd.read_excel(pathExpResult,sheet_name='Sheet1', usecols="A,B")
    dataExp = df.iloc[:, 0:2].to_numpy(dtype=float)


    # plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(dataNum['target'], -dataNum['RF'][:, -1], label='numerical')
    ax.plot(dataExp[:, 0]/100, dataExp[:, -1], label='experimental')
    plt.legend()
    ax.set_xlabel(r'$\epsilon_{tensile}$')
    ax.set_ylabel(r'$F_{tensile}$ [N]')
    ax.set_title('Tensile result of rope')

    # save
    fig.savefig(pngFile, dpi=300, transparent=True)
    plt.close('all')
    return

# ------------------------------------------

if __name__ == '__main__':
    import sys
    plotFU(pathNumResult=sys.argv[-3], pathExpResult=sys.argv[-2], pngFile=sys.argv[-1])
