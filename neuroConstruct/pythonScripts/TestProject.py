# -*- coding: utf-8 -*-
#
#
#   File to test project to ensure it is set up correctly for simulations
#
#   To execute this type of file, type '..\..\..\nC.bat -python XXX.py' (Windows)
#   or '../../../nC.sh -python XXX.py' (Linux/Mac). Note: you may have to update the
#   NC_HOME and NC_MAX_MEMORY variables in nC.bat/nC.sh
#
#   Author: Padraig Gleeson
#
#   This file has been developed as part of the neuroConstruct project
#   This work has been funded by the Medical Research Council and the
#   Wellcome Trust
#
#

import sys
import os

try:
    from java.io import File
except ImportError:
    print "Note: this file should be run using ..\\..\\..\\nC.bat -python XXX.py' or '../../../nC.sh -python XXX.py'"
    print "See http://www.neuroconstruct.org/docs/python.html for more details"
    quit()

sys.path.append(os.environ["NC_HOME"]+"/pythonNeuroML/nCUtils")

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.neuron import NeuronSettings
from ucl.physiol.neuroconstruct.utils import Display3DProperties


projFile = File("../KorngreenPyramidal.ncx")

def arrayListsIdentical(arrayList, list):
    assert(arrayList.size() == len(list))
    for ale in arrayList:
        assert(str(ale) in list)


def testAll(argv=None):
    if argv is None:
        argv = sys.argv

    print "Loading project from "+ projFile.getCanonicalPath()
    
    projectManager = ProjectManager()
    project = projectManager.loadProject(projFile)

    assert(len(project.getProjectDescription())>0)

    assert(len(project.cellManager.getAllCells())>=7)

    assert(len(project.cellGroupsInfo.getAllCellGroupNames())>=4)

    assert(project.simulationParameters.getDt()-0.0125<=1e-6)

    assert(project.neuronSettings.isVarTimeStep())

    assert(project.neuronSettings.getDataSaveFormat().equals(NeuronSettings.DataSaveFormat.TEXT_NC))

    assert(project.simulationParameters.getTemperature() == 34)

    '''
    sci = project.simConfigInfo
    cgi = project.cellGroupsInfo

    assert(cgi.getCellType("Channeltestgroup") == "Channeltest_Cell")
    assert(cgi.getCellType("CMLtestGroup") == "CMLtest_Cell")
    assert(cgi.getCellType("pyr_group") == "LarkumPyr")
    assert(cgi.getCellType("pyrCML_group") == "LarkumPyr_NML")

    arrayListsIdentical(sci.getSimConfig("Default Simulation Configuration").getCellGroups(), ["Channeltestgroup", "CMLtestGroup"])
    arrayListsIdentical(sci.getSimConfig("test_IClamp").getCellGroups(), ["pyr_group", "pyrCML_group"])
    arrayListsIdentical(sci.getSimConfig("background activity").getCellGroups(), ["pyr_group"])
    '''

    print "\n**************************************"
    print "    All tests passed!"
    print "**************************************\n"


if __name__ == "__main__":
    testAll()
    exit()