# utility script for setting model parameters for the Korngreen model
# from .params files, via neuroConstruct
import rlcompleter, readline
readline.parse_and_bind('tab: complete')

from java.io import File
from java.util import ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.utils.equation import Variable
from ucl.physiol.neuroconstruct.cell import VariableParameter, VariableMechanism, ChannelMechanism, IonProperties, ParameterisedGroup

from korngreen_utils import AlmogKorngreenPars, Mechanism


class NCProject(object):

    def __init__(self, fname):
        self._project = self.load_project(fname)

    def load_project(self, fname):
        file = File(fname)
        print 'Loading project file: ', file.getAbsolutePath()
        pm = ProjectManager()
        project = pm.loadProject(file)
        print pm.status()
        return project

    def add_cell(self, cell):
        self._project.cellManager.addCellType(k.cell)
    
    def get_cell(self, name):
        return self._project.cellManager.getCell(name)

    def save(self):
        self._project.markProjectAsEdited()
        self._project.saveProject()



class KorngreenCell(object):
    
    def __init__(self, parsname):
        self._pars = AlmogKorngreenPars()
        self._pars.pars_from_file(parsname)


    def create_from_passive(self, pascell, name='test'):
        cell = pascell.clone()
        cell.setInstanceName(name)
        self._cell = cell

    def clear_groups(self):
        for sec in self._cell.allSections:
            groups = sec.getGroups().clone()
            for g in groups:
                sec.removeFromGroup(g)

    def add_groups(self):
        for sec in self._cell.allSections:
            if not 'all' in sec.getGroups():
                sec.addToGroup('all')
            for pat,gname in {'dend':'basal_dend_group', 'apic':'apical_dend_group', 'soma':'soma_group', 'iseg':'iseg_group', 'myelin':'myelin_group', 'hill':'hill_group', 'node':'nodes_group'}.iteritems():
                if pat in sec.getSectionName():
                    sec.addToGroup(gname)
    
    @property
    def groups(self):
       return self._cell.getAllGroupNames()

    @property
    def parametrized_groups(self):
       return self._cell.getParameterisedGroups()

    @property
    def cell(self):
       return self._cell
        
    def create_apical_parametrized_group(self):
        #assuming that there is only one group with "apic" in the name!
        pname = 'p'
        pg = ParameterisedGroup('PathLengthApicalDends', 'apical_dend_group', ParameterisedGroup.Metric.PATH_LENGTH_FROM_ROOT, ParameterisedGroup.ProximalPref.NO_TRANSLATION, ParameterisedGroup.DistalPref.NO_NORMALISATION, pname)
        self._cell.getParameterisedGroups().add(pg);
        return pg


    def add_inhomogeneous_mechanisms(self, mod=False):
        pg = self.create_apical_parametrized_group()
        for mech in self._pars.inhomogeneous_mechs:
            if 'ca' in mech.name:
                par = 'permeability'
            else:
                 par = 'gmax'
            if mod:
                mech.name += '_mod'
            vp = VariableParameter(par, mech.expr, Variable(pg.getVariable()), ArrayList())
            vm = VariableMechanism(mech.name, vp)
            for epname, epval in mech.extra_parameters.iteritems():
                print 'extra par: ', epname, epval
                vm.setExtraParam(epname, epval)
            self._cell.associateParamGroupWithVarMech(pg, vm)
            
    def add_mechanisms_to_group(self, group_name, list_mechs, mod=False):
        for m in list_mechs:
            name = m.name
            val = m.gmax
            if mod:
                name += '_mod'
            cm = ChannelMechanism(name, val)
            for epname, epval in m.extra_parameters.iteritems():
                cm.setExtraParam(epname, epval)
            
            self._cell.associateGroupWithChanMech(group_name, cm)
    
    def add_homogeneous_mechanisms(self, mod=False):
         self.add_mechanisms_to_group('soma_group', self._pars.soma_group.mechanisms, mod)
         self.add_mechanisms_to_group('hill_group', self._pars.hill_group.mechanisms, mod)
         self.add_mechanisms_to_group('iseg_group', self._pars.iseg_group.mechanisms, mod)
         self.add_mechanisms_to_group('nodes_group', self._pars.node_group.mechanisms, mod)
         self.add_mechanisms_to_group('myelin_group', self._pars.myelin_group.mechanisms, mod)
         self.add_mechanisms_to_group('basal_dend_group', self._pars.basal_dend_group.mechanisms, mod)
         self.add_mechanisms_to_group('all', [Mechanism('cad', 0)])
    
    def add_ion_properties(self):
        self._cell.associateGroupWithIonProperties('all', IonProperties('na', self._pars.Ena))
        self._cell.associateGroupWithIonProperties('all', IonProperties('k', self._pars.Ek))
        self._cell.associateGroupWithIonProperties('all', IonProperties('ca', self._pars.cai, self._pars.cao))
        for g in ['apical_dend_group', 'basal_dend_group', 'soma_group']:
            self._cell.associateGroupWithIonProperties(g, IonProperties('h', -33.0))





if __name__ == '__main__':

    proj = NCProject("../neuroConstruct/KorngreenPyramidal.ncx")

    k = KorngreenCell("best.params")
    #k.create_from_passive(proj.get_cell('A140612_pas'), 'A140612_nml')
    k.create_from_passive(proj.get_cell('A140612_pas'), 'test')
    #k.clear_groups()
    k.add_groups()
    k.add_inhomogeneous_mechanisms(mod=False)
    k.add_homogeneous_mechanisms(mod=False)
    k.add_ion_properties()

    proj.add_cell(k.cell)


