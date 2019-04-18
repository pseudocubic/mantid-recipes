#!/usr/bin/env python

import os
import sys
import numpy as np

mantid_root = "/opt/mantidnightly"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

import mantid.simpleapi as mantidapi
from mantid.simpleapi import *


def outtofile(filename, inten, err, x, y, path=""):
    ""
    fn = filename
    output = []
    f = open(fn, 'w')
    for sig, ne, h, k in zip(inten, err, x, y):
        out = [sig, ne, h, k]
        output.append("\t".join(str(o) for o in out))

    f.write("\n".join(output))
    f.close()

def loadWS(runfiles, inputdir, outputdir, filename='mdWorkspace', save=False):
    if type(runfiles) is list and 'nxspe' in runfiles[0]:
        for file in runfiles:
            ws = mantidapi.LoadNXSPE(file)
        mantidapi.LoadInstrument(ws, InstrumentName='HYS')
    else:
        ws = mantidapi.Load(runfiles)

    mantidapi.SetGoniometer(ws, Axis0="s1,0,1,0,1")
    mantidapi.SetUB(ws, 3.81, 3.81, 6.25, 90, 90, 90, "1,0,0", "0,0,1")

    md = mantidapi.ConvertToMD(ws, QDimensions="Q3D",
                            dEAnalysisMode="Direct",
                            Q3DFrames="HKL",
                            QConversionScales="HKL",
                            PreprocDetectorsWS="",
                            MinValues="-4,-2,-4,-22.5",
                            MaxValues="4,2,4,22.5",
                            MaxRecursionDepth=1)
    mdall = mantidapi.MergeMD(md)

    mantidapi.DeleteWorkspace(ws)
    mantidapi.DeleteWorkspace(md)

    if save: mantidapi.SaveMD(mdall, os.path.join(outputdir, filename + '.nxs'))

    return mdall


def slice_phonon(workspace, outputdir, filename, minE, maxE, stepE, minH, maxH, binsH, minK, maxK, binsK, minL, maxL, binsL):
    start, stop, step = minE, maxE, stepE
    En = np.linspace(start, stop, (stop - start) / step + 1)
    En = np.delete(En, np.where(En == 0.))

    Ens = [str(e - step) + ',' + str(e + step) for e in En]

    for n, deltaE in enumerate(Ens):
        _tempMD = mantidapi.BinMD(workspace, AlignedDim0='[H,0,0],{0},{1},{2}'.format(minH, maxH, binsH),
                                             AlignedDim1='[0,K,0],{0},{1},{2}'.format(minK, maxK, binsK),
                                             AlignedDim2='[0,0,L],{0},{1},{2}'.format(minL, maxL, binsL),
                                             AlignedDim3='DeltaE,' + deltaE + ',1',
                                             Parallel='1')

        signal = _tempMD.getSignalArray()
        numEvents = _tempMD.getNumEventsArray()

        x = np.linspace(minH, maxH, binsH)
        y = np.linspace(minL, maxL, binsL)
        X, Y = np.meshgrid(x, y)

        fn = os.path.join(outputdir, filename + '_' + str(En[n]).replace('.', 'p') + 'mev.iexy')
        outtofile(fn, signal.flatten(), numEvents.flatten(), X.flatten(), Y.flatten())

    mantidapi.DeleteWorkspace(_tempMD)


def slice_HvE(workspace, outputdir='', filename='Material_7K_HK0_HvE'):
    _tempMD = mantidapi.BinMD(workspace,
                              AxisAligned='0',
                              BasisVector0='[H,0,0],in 1.649 A^-1,1,0,0,0',
                              BasisVector1='[0,K,0],in 1.649 A^-1,0,1,0,0',
                              BasisVector2='[0,0,L],in 1.005 A^-1,0,0,1,0',
                              BasisVector3='DeltaE,DeltaE,0,0,0,1',
                              OutputExtents='-4,4,-0.2,0.2,-0.5,0.5,-22.5,22.5',
                              OutputBins='321,1,1,181',
                              Parallel='1',
                              OutputWorkspace='_tempMD')

    signal = _tempMD.getSignalArray()
    numEvents = _tempMD.getNumEventsArray()
    sqErr = _tempMD.getErrorSquaredArray()
    inds = np.where(numEvents > 0)

    inten = signal[inds] / numEvents[inds]
    err = np.sqrt(sqErr[inds]) / numEvents[inds]

    x = np.linspace(-4, 4, 321)
    y = np.linspace(-22.5, 22.5, 181)

    X, Y = np.meshgrid(x, y)
    x, y = X[inds], X[inds]

    fn = os.path.join(outputdir, filename + '.iexy')
    outtofile(fn, signal.flatten(), numEvents.flatten(), X.flatten(), Y.flatten())

    mantidapi.DeleteWorkspace(_tempMD)


if __name__ == '__main__':
    inputdir = ['/SNS/HYS/IPTS-xxxx/shared/h0l_Material_reduced/h0l_I100K_xxxx',
                'SNS/HYS/IPTS-xxxx/shared/h0l_Material_reduced/h0l_I150K_xxxx'
                ]

    outputdir = ['/SNS/HYS/IPTS-xxxx/shared/100K_slices',
                 '/SNS/HYS/IPTS-xxxx/shared/150K_slices'
                 ]

    runfiles1 = str([inputdir[0] + 'h0l_I150K_xxxx_' + str(r) + '.00.nxspe' for r in range(-95, -45)]).strip('[').strip(']').replace("'", "")
    runfiles2 = str([inputdir[1] + 'h0l_I150K_xxxx_' + str(r) + '.00.nxspe' for r in range(-95, -45)]).strip('[').strip(']').replace("'", "")
    runfiles = [runfiles1, runfiles2]

    filenames = ['Material_ei15_100K_H0L_HvL', 
                 'Material_ei15_150K_H0L_HvL'
                 ]

    for indir, outdir, runfile, filename in zip(inputdir, outputdir, runfiles, filenames):
        workspace = loadWS(runfile, indir, outdir)
        slice_phonon(workspace, outdir, filename, -10., 10., 0.5, -4., 4., 321, -0.1, 0.1, 1, -4, 4, 321)
