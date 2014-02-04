# Utilities for generating homo and heterogeneuos parameters from a ga fit file (best.params originally)
def ineq_to_heaviside(cond, ifTrue, ifFalse):
    import re

    r = re.compile("(.+)([><]+)(.+)")
    s = r.search(cond)
    f1,op,f2 = [g.strip() for g in s.groups()]

    if op == '>':
        harg = f1 + ' - ' + f2
    elif op == '<':
        harg = f2 + ' - ' + f1  

    return "(%s) * H(%s) + (%s) * H(-(%s))" % (ifTrue, harg, ifFalse, harg)
    
def to_condition_str(el1, op, el2):
    assert op in ['<', '>']
    return ' '.join((str(el1), op, str(el2)))

class SectionGroup(object):
    def __init__(self, name, mechs):
        self.name = name
        self.mechanisms = mechs
        
    def __repr__(self):
        r = "section group:" + self.name + "\n"
        for m in self.mechanisms:
            r += "\t" + str(m) + "\n"
        return r

class Mechanism(object):
    def __init__(self, name, gmax=0.0, extrapars={}):
        self.name = name
        self.gmax = gmax
        self.extra_parameters = extrapars
        
    def __repr__(self):
        r = "mechanism:" + self.name + "\n"
        r += "\tgmax:" + str(self.gmax) + "\n"
        r += "\textra pars:" + str(self.extra_parameters) + "\n"
        return r

class VariableMechanism(object):
    def __init__(self, name, expr, extrapars={}):
        self.name = name
        self.expr = expr
        self.extra_parameters = extrapars
        
    def __repr__(self):
        r = "variable mechanism:" + self.name + "\n"
        r += "\texpression:" + self.expr + "\n"
        r += "\textra pars:" + str(self.extra_parameters) + "\n"
        return r
        

class AlmogKorngreenPars(object):
    
    #Distances will be translated by half the length of the soma, as nC measures 
    # from soma(0) and the model used soma(0.5)
    origin = 16.228978
    indepvar = '(p - %g)' % origin

    sigmoid_t = "%(gsoma)g + %(gdend)g/(1+exp(%(slope)g*(%(x)s - %(half)g)))"
    exp_t = "%(gdend)g + %(gsoma)g * exp(%(slope)g * %(x)s)"
    lin_t = "%(gsoma)g + %(x)s * (%(gdend)g - %(gsoma)g) / %(dist)g"

    def __init__(self):
        #This is a rough initial parameter set. It is overwritten by parameters loaded from best.params
        self.cm_myelin = 0.04
        self.g_pas_node = 0.02
        self.v_init = -62.
        self.Ek = -100.
        self.Ena = 60.
        self.Eca = 130.

        #passive
        self.ra = 68.16690
        self.rm = 18017.30000
        self.c_m = 1.43870
        self.epas_sim = -31.40480

        #gih
        self.gih_end = 271.90500*1e-9
        self.gih_x2 = 382.51900
        self.gih_alpha =-0.08090
        self.gih_start = 22.51460*1e-9
        self.ih_q10 = 2.0
        self.Eh = -33.0 # from .mod files
        
        #kslow
        self.gkslow_start = 2.03452*1e-9
        self.gkslow_alpha =-0.00986
        self.gkslow_beta = 127.82800*1e-9
        
        #kfast
        self.gka_start = 20.34520*1e-9
        self.gka_alpha = -0.00297
        self.gka_beta = 320.43100*1e-9

        #na
        self.gna_soma = 128.14300*1e-9
        self.gna_api = 2.69229*1e-9
        self.dist_na = 687.53800
        self.na_shift1 =-5.56301
        self.na_shift2 =-4.52361
        self.na_taum_scale = 1.
        self.na_tauh_scale = 1.

        #Ca
        self.pcah_soma = 1.33e-5
        self.pcah_api = self.pcah_soma
        self.dist_cah = 600.
        self.cah_qm = 4.#changed in accordance with modfiles
        self.cah_qh = 2.
        self.cah_shift = 0.
        self.cah_shifth = 0.
        self.cai = 5e-23
        self.cao = 2e-18 #eca ~140mV

        self.pcar_soma = 0.27e-5
        self.pcar_api = self.pcar_soma
        self.dist_car = 600.
        self.car_qh = 1.
        self.car_qm = 1.
        self.car_shift = 0.
        self.car_shifth = 0.

        #axon parameters
        self.gkslow_node = 1500.*1e-9
        self.gka_node = 1000.*1e-9
        self.gna_node = 30000.*1e-9
        self.shift_na_act_axon = 7.
        self.shift_na_inact_axon = 3.

    
    def pars_from_file(self, fname='best.params'):
        #no numpy in jyton...
        #from numpy import loadtxt
        #p = loadtxt(fname)
        
        import csv
        rdr = csv.reader(open(fname), delimiter=' ')
        p = [float(r[0]) for r in rdr if len(r) > 1]

        self.ra = p[0]
        self.rm = p[1]
        self.c_m = p[2]
        self.epas_sim = p[3]

        self.gih_end  = p[4]*1e-9
        self.gih_x2 = p[5]
        self.gih_alpha = p[6]
        self.gih_start = p[7]*1e-9
        self.ih_q10 = p[8]

        self.gkslow_start = p[9]*1e-9
        self.gkslow_alpha = p[10]
        self.gkslow_beta = p[11]*1e-9

        self.gka_start = p[12]*1e-9
        self.gka_alpha = p[13]
        self.gka_beta = p[14]*1e-9

        self.gna_soma = p[15]*1e-9
        self.gna_api = p[16]*1e-9
        self.dist_na = p[17]
        self.na_shift1 = p[18]
        self.na_shift2 = p[19]

        self.pcah_soma = p[20]*1e-2
        self.pcah_api = p[21]*1e-2
        self.dist_cah = p[22]
        self.cah_shift = p[23]
        self.cah_shifth = p[24]

        self.pcar_soma = p[25]*1e-2
        self.pcar_api = p[26]*1e-2
        self.dist_car = p[27]
        self.car_shift = p[28]
        self.car_shifth = p[29]
        self.car_qm = p[30]	

        self.gsk_soma = p[31]*1e-9
        self.gsk_dend = p[32]*1e-9
        self.dist_sk = p[33]

        self.gbk_soma = p[34]*1e-9
        self.gbk_dend = p[35]*1e-9
        self.dist_bk = p[36]		
        
        
        #global pars
        extra_na = {'vshiftm': self.na_shift1, 'vshifth': self.na_shift2, 'taum_scale': self.na_taum_scale, 'tauh_scale': self.na_tauh_scale}
        extra_iH = {'q10': self.ih_q10}
        extra_cah = {'GHK_permeability': self.pcah_soma, 'shift': self.cah_shift, 'shifth': self.cah_shifth, 'qm':self.cah_qm, 'qh':self.cah_qh}
        extra_car = {'GHK_permeability':self.pcar_soma, 'shift': self.car_shift, 'shifth': self.car_shifth, 'qm':self.car_qm, 'qh':self.car_qh}
        

        #soma pars
        mechsoma = [Mechanism('na', self.gna_soma, extra_na),
                    Mechanism('kslow', self.gkslow_start + self.gkslow_beta),
                    Mechanism('iA', self.gka_start + self.gka_beta),
                    Mechanism('iH', self.gih_start, extra_iH),
                    Mechanism('car', 0., extra_car),
                    Mechanism('cah', 0., extra_cah),
                    Mechanism('sk', self.gsk_soma),
                    Mechanism('bk', self.gbk_soma)]
        self.soma_group = SectionGroup('soma_group', mechsoma)
        

        #basal dends pars
        mechdend = [Mechanism('na', self.gna_soma, extra_na),
                    Mechanism('kslow', self.gkslow_start + self.gkslow_beta),
                    Mechanism('iA', self.gka_start + self.gka_beta),
                    Mechanism('iH', self.gih_start, extra_iH),
                    Mechanism('car', 0., extra_car),
                    Mechanism('cah', 0.,  extra_cah),
                    Mechanism('sk', self.gsk_dend),
                    Mechanism('bk', self.gbk_dend)]
        self.basal_dend_group = SectionGroup('dend_group', mechdend)
        

        #axon pars
        extra_na_axon = extra_na.copy()
        extra_na_axon['vshiftm'] = 7.
        extra_na_axon['vshifth'] = 3.

        extra_cah_axon = extra_cah.copy()
        extra_cah_axon['GHK_permeability'] = 2e-6

        extra_car_axon = extra_car.copy()
        extra_car_axon['GHK_permeability'] = 4e-7

        mechaxon = [Mechanism('na', self.gna_node, extra_na_axon),
                    Mechanism('kslow', self.gkslow_node),
                    Mechanism('iA', self.gka_node),
                    Mechanism('car', 0. , extra_car_axon),
                    Mechanism('cah', 0., extra_cah_axon),
                    Mechanism('cad', 0.),
                    Mechanism('bk', 40*1e-9)]
        self.hill_group = SectionGroup('hill_group', mechaxon)
        self.node_group = SectionGroup('node_group', mechaxon)
        self.iseg_group = SectionGroup('iseg_group', mechaxon)
        

        mechmyelin = [Mechanism('na', self.gna_soma, extra_na),
                      Mechanism('kslow', self.gkslow_beta),
                      Mechanism('iA', self.gka_beta),
                      Mechanism('car', 0., extra_car_axon),
                      Mechanism('cah', 0., extra_cah_axon),
                      Mechanism('cad', 0.),
                      Mechanism('bk', 40*1e-9)]
        self.myelin_group = SectionGroup('myelin_group', mechmyelin)
        
        mechapical = [VariableMechanism('na', self.gna_expr(), extra_na),
                      VariableMechanism('iA', self.giA_expr()),
                      VariableMechanism('kslow', self.gkslow_expr()),
                      VariableMechanism('iH', self.giH_expr(), extra_iH),
                      VariableMechanism('sk', self.gsk_expr()),
                      VariableMechanism('bk', self.gbk_expr()),
                      VariableMechanism('cah', self.pcah_expr(), extra_cah),
                      VariableMechanism('car', self.pcar_expr(), extra_car)]
        self.inhomogeneous_mechs = mechapical


    def giH_expr(self):
        return self.sigmoid_t%{'gsoma':self.gih_start,
                                     'gdend':self.gih_end,
                                     'slope':self.gih_alpha,
                                     'half':self.gih_x2,
                                     'x':self.indepvar}

    def giA_expr(self):
        return self.exp_t%{'gdend':self.gka_start,
                                 'gsoma':self.gka_beta,
                                 'slope':self.gka_alpha,
                                 'x':self.indepvar}

    def gkslow_expr(self):
        return self.exp_t%{'gdend':self.gkslow_start,
                                 'gsoma':self.gkslow_beta,
                                 'slope':self.gkslow_alpha,
                                 'x':self.indepvar}

    def gna_expr(self):
        distr = self.lin_t % {'gsoma':self.gna_soma,
                                  'gdend':self.gna_api,
                                  'dist':self.dist_na,
                                  'x':self.indepvar}
        cond1 = to_condition_str(self.indepvar, '<', self.dist_na)
        condexpr1 = ineq_to_heaviside(cond1, distr, self.gna_api)
        cond2 = to_condition_str(condexpr1, '<', self.gna_api)
        return ineq_to_heaviside(cond2, self.gna_api, condexpr1)

    def gsk_expr(self):
        distr = self.lin_t%{'gsoma':self.gsk_soma,
                                  'gdend':self.gsk_dend,
                                  'dist':self.dist_sk,
                                  'x':self.indepvar}
        cond = to_condition_str(self.indepvar, '<', self.dist_sk)
        condexpr1 =  ineq_to_heaviside(cond, distr, self.gsk_dend)
        cond2 = to_condition_str(condexpr1, '<', self.gsk_dend)
        return ineq_to_heaviside(cond2, self.gsk_dend, condexpr1)

    def gbk_expr(self):
        distr = self.lin_t%{'gsoma':self.gbk_soma,
                                  'gdend':self.gbk_dend,
                                  'dist':self.dist_bk,
                                  'x':self.indepvar}
        cond1 = to_condition_str(self.indepvar, '<', self.dist_bk)
        condexpr1 = ineq_to_heaviside(cond1, distr, self.gbk_dend)
        cond2 = to_condition_str(condexpr1, '<', self.gbk_dend)
        return ineq_to_heaviside(cond2, self.gbk_dend, condexpr1)

    def pcah_expr(self):
        distr = self.lin_t%{'gsoma':self.pcah_soma,
                                  'gdend':self.pcah_api,
                                  'dist':self.dist_cah,
                                  'x':self.indepvar}
        cond = to_condition_str(self.indepvar, '<', self.dist_cah)
        condexpr1 = ineq_to_heaviside(cond, distr, self.pcah_api)
        return condexpr1
        #if self.pcah_soma < self.pcah_api:
        #    cond2 = to_condition_str(condexpr1, '>', self.pcah_api)
        #    ret = ineq_to_heaviside(cond2, self.pcah_api, condexpr1)
        #else:
        #    cond2 = to_condition_str(condexpr1, '<', self.pcah_api)
        #    ret = ineq_to_heaviside(cond2, self.pcah_api, condexpr1)
        #return ret
            
            

    def pcar_expr(self):
        distr = self.lin_t%{'gsoma':self.pcar_soma,
                                  'gdend':self.pcar_api,
                                  'dist':self.dist_car,
                                  'x':self.indepvar}
        cond = to_condition_str(self.indepvar, '<', self.dist_car)
        condexpr1 = ineq_to_heaviside(cond, distr, self.pcar_api)
        return condexpr1
        #if self.pcar_soma < self.pcar_api:
        #    cond2 = to_condition_str(condexpr1, '>', self.pcar_api)
        #    ret = ineq_to_heaviside(cond2, self.pcar_api, condexpr1)
        #else:
        #    cond2 = to_condition_str(condexpr1, '<', self.pcar_api)
        #    ret = ineq_to_heaviside(cond2, self.pcar_api, condexpr1)
        #return ret
    

if __name__ == '__main__':
    p = AlmogKorngreenPars()
    p.pars_from_file('best.params')

    
    print p.soma_group
    print p.basal_dend_group
    print p.iseg_group
    print p.hill_group
    print p.node_group
    print p.myelin_group


    print 'giH:', p.giH_expr()
    print 'giA:', p.giA_expr()
    print 'gkslow:', p.gkslow_expr()
    print 'gna:', p.gna_expr()
    print 'gsk:', p.gsk_expr()
    print 'gbk:', p.gbk_expr()
    print 'pcah:', p.pcah_expr()
    print 'pcar:', p.pcar_expr()


