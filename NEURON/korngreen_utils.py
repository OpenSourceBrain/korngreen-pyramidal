
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

class AlmogKorngreenPars(object):
    
    #Distances will be translated by half the length of the soma, as nC measures 
    # from soma(0) and the model used soma(0.5)
    origin = 16.228978
    indepvar = '(p - %g)' % origin

    sigmoid_t = "{gsoma} + {gdend}/(1+exp({slope}*({x} - {half})))"
    exp_t = "{gdend} + {gsoma} * exp({slope} * {x})"
    lin_t = "{gsoma} + {x} * ({gdend} - {gsoma}) / {dist}"

    def __init__(self):
        #This is a rough initial parameter set. It is overwritten by parameters loaded from best.params
        self.cm_myelin = 0.04
        self.g_pas_node = 0.02
        self.v_init = -62
        self.Ek = -100
        self.Ena = 60
        self.Eca = 130

        #passive
        self.ra = 68.16690
        self.rm = 18017.30000
        self.c_m = 1.43870
        self.epas_sim = -31.40480

        #gih
        self.gih_end = 271.90500*1e-5
        self.gih_x2 = 382.51900
        self.gih_alpha =-0.08090
        self.gih_start = 22.51460*1e-5
        self.ih_q10 = 2
        
        #kslow
        self.gkslow_start = 2.03452*1e-5
        self.gkslow_alpha =-0.00986
        self.gkslow_beta = 127.82800*1e-5
        
        #kfast
        self.gka_start = 20.34520*1e-5
        self.gka_alpha = -0.00297
        self.gka_beta = 320.43100*1e-5

        #na
        self.gna_soma = 128.14300*1e-5
        self.gna_api = 2.69229*1e-5
        self.dist_na = 687.53800
        self.na_shift1 =-5.56301
        self.na_shift2 =-4.52361
        self.na_taum_scale = 1
        self.na_tauh_scale = 1

        #Ca
        self.pcah_soma = 1.33e-3
        self.pcah_api = self.pcah_soma
        self.dist_cah = 600
        self.cah_qm = 2
        self.ca_qh = 2
        self.cah_shift = 0
        self.cah_shifth = 0

        self.pcar_soma = 0.27e-3
        self.pcar_api = self.pcar_soma
        self.dist_car = 600
        self.car_qh = 1
        self.car_qm = 1
        self.car_shift = 0
        self.car_shifth = 0

        #axon parameters
        self.gkslow_node = 1500
        self.gka_node = 1000
        self.gna_node = 30000
        self.shift_na_act_axon = 7
        self.shift_na_inact_axon = 3

    
    def pars_from_file(self, fname):
        from numpy import loadtxt
        p = loadtxt('best.params')

        self.ra = p[0]
        self.rm = p[1]
        self.c_m = p[2]
        self.epas_sim = p[3]

        self.gih_end  = p[4]*1e-5
        self.gih_x2 = p[5]
        self.gih_alpha = p[6]
        self.gih_start = p[7]*1e-5
        self.ih_q10 = p[8]

        self.gkslow_start = p[9]*1e-5
        self.gkslow_alpha = p[10]
        self.gkslow_beta = p[11]*1e-5

        self.gka_start = p[12]*1e-5
        self.gka_alpha = p[13]
        self.gka_beta = p[14]*1e-5

        self.gna_soma = p[15]*1e-5
        self.gna_api = p[16]*1e-5
        self.dist_na = p[17]
        self.na_shift1 = p[18]
        self.na_shift2 = p[19]

        self.pcah_soma = p[20]*1e-5
        self.pcah_api = p[21]*1e-5
        self.dist_cah = p[22]
        self.cah_shift = p[23]
        self.cah_shifth = p[24]

        self.pcar_soma = p[25]*1e-5
        self.pcar_api = p[26]*1e-5
        self.dist_car = p[27]
        self.car_shift = p[28]
        self.car_shifth = p[29]
        self.car_qm = p[30]	

        self.gsk_soma = p[31]*1e-5
        self.gsk_dend = p[32]*1e-5
        self.dist_sk = p[33]

        self.gbk_soma = p[34]*1e-5
        self.gbk_dend = p[35]*1e-5
        self.dist_bk = p[36]		

    def forall(self):
        return dict(Ra = self.ra, cm = self.c_m, g_pas = 1./self.rm, e_pas = self.epas_sim, vshiftm_na= self.na_shift1, vshifth_na= self.na_shift2, taum_scale_na= self.na_taum_scale, tauh_scale_na= self.na_tauh_scale, q10_iH= self.ih_q10, shift_cah= self.cah_shift, shifth_cah= self.cah_shifth, shift_car= self.car_shift, shifth_car= self.car_shifth, qm_car= self.car_qm)

        

    def giH_expr(self):
        return self.sigmoid_t.format(gsoma=self.gih_start,
                                     gdend=self.gih_end,
                                     slope=self.gih_alpha,
                                     half=self.gih_x2,
                                     x=self.indepvar)

    def giA_expr(self):
        return self.exp_t.format(gdend=self.gka_start,
                                 gsoma=self.gka_beta,
                                 slope=self.gka_alpha,
                                 x=self.indepvar)

    def gkslow_expr(self):
        return self.exp_t.format(gdend=self.gkslow_start,
                                 gsoma=self.gkslow_beta,
                                 slope=self.gkslow_alpha,
                                 x=self.indepvar)

    def gna_expr(self):
        distr = self.lin_t.format(gsoma=self.gna_soma,
                                  gdend=self.gna_api,
                                  dist=self.dist_na,
                                  x=self.indepvar)
        cond = ' '.join((self.indepvar, '<', str(self.dist_na)))
        return ineq_to_heaviside(cond, distr, self.gna_api)

    def gsk_expr(self):
        distr = self.lin_t.format(gsoma=self.gsk_soma,
                                  gdend=self.gsk_dend,
                                  dist=self.dist_sk,
                                  x=self.indepvar)
        cond = ' '.join((self.indepvar, '<', str(self.dist_sk)))
        return ineq_to_heaviside(cond, distr, self.gsk_dend)

    def gbk_expr(self):
        distr = self.lin_t.format(gsoma=self.gbk_soma,
                                  gdend=self.gbk_dend,
                                  dist=self.dist_bk,
                                  x=self.indepvar)
        cond = ' '.join((self.indepvar, '<', str(self.dist_bk)))
        return ineq_to_heaviside(cond, distr, self.gbk_dend)

    def pcah_expr(self):
        distr = self.lin_t.format(gsoma=self.pcah_soma,
                                  gdend=self.pcah_api,
                                  dist=self.dist_cah,
                                  x=self.indepvar)
        cond = ' '.join((self.indepvar, '<', str(self.dist_cah)))
        return ineq_to_heaviside(cond, distr, self.pcah_api)

    def pcar_expr(self):
        distr = self.lin_t.format(gsoma=self.pcar_soma,
                                  gdend=self.pcar_api,
                                  dist=self.dist_car,
                                  x=self.indepvar)
        cond = ' '.join((self.indepvar, '<', str(self.dist_car)))
        return ineq_to_heaviside(cond, distr, self.pcar_api)

if __name__ == '__main__':
    p = AlmogKorngreenPars()
    p.pars_from_file('best.params')
    print 'giH:', p.giH_expr()
    print 'giA:', p.giA_expr()
    print 'gkslow:', p.gkslow_expr()
    print 'gna:', p.gna_expr()
    print 'gsk:', p.gsk_expr()
    print 'gbk:', p.gbk_expr()
    print "attention, the ones below require extra checks\n\t (see last 'forsec ApicalDendSectionName' in model.hoc)"
    print 'pcah:', p.pcah_expr()
    print 'pcar:', p.pcar_expr()

