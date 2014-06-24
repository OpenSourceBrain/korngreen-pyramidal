#
#
#   File to test current configuration of project
#
#   Author: Padraig Gleeson
#
#   This file has been developed as part of the neuroConstruct project
#   This work has been funded by the Wellcome Trust
#
#

import sys
import os

try:
    from java.io import File
except ImportError:
    print "Note: this file should be run using nC.bat -python XXX.py' or 'nC.sh -python XXX.py'"
    print "See http://www.neuroconstruct.org/docs/python.html for more details"
    quit()

sys.path.append(os.environ["NC_HOME"]+"/pythonNeuroML/nCUtils")

import ncutils as nc # Many useful functions such as SimManager.runMultipleSims found here

projFile = File("../KorngreenPyramidal.ncx")


##############  Main settings  ##################

simConfigs = []

simConfigs.append("bac6")

simDt =                 0.025

simulators =            ["NEURON"]

varTimestepNeuron =     False

plotSims =              True
plotVoltageOnly =       True
runInBackground =       True
analyseSims =           True
verbose =               True

#############################################


def testAll(argv=None):
    if argv is None:
        argv = sys.argv

    print "Loading project from "+ projFile.getCanonicalPath()


    simManager = nc.SimulationManager(projFile,
                                      verbose = verbose)

    simManager.runMultipleSims(simConfigs =           simConfigs,
                               simDt =                simDt,
                               simulators =           simulators,
                               runInBackground =      runInBackground,
                               varTimestepNeuron =    varTimestepNeuron)

    simManager.reloadSims(plotVoltageOnly =   plotVoltageOnly,
                          analyseSims =       analyseSims)

             
    spikeTimesToCheck = {'FullCell_nml_0': [403.675, 416.9],
                         'FullCell_nml_0.2674': [405.4, 417.4],
                         'FullCell_nml_0.2720': [405.925]}
    
    spikeTimeAccuracy = 0.01

    report = simManager.checkSims(spikeTimesToCheck = spikeTimesToCheck,
                                  spikeTimeAccuracy = spikeTimeAccuracy,
                                  threshold=-10)

    print report

    return report


if __name__ == "__main__":
    testAll()


