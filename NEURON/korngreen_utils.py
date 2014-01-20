class AlmogKorngreenPars(object):

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

    def limit_value(self, expr, cond, lim):
        harg = "-".join((expr,lim))
        if cond == '>':
            f1 = expr
            f2 = lim
        elif cond == '<':
                f1 = expr
                f2 = lim
        return "%s*H(%s)+%s*H(-(%s))"%(f1,harg,f2,harg)

    sigmoid_t = "%f+%f/(1+exp(%f*(p-%f)))"
    exp_t = "%f+%f*exp(%f*p)"
    lin_t = "%f+p*(%f-%f)/%f"

    def giH_expr(self):
        p = (self.gih_start,self.gih_end,self.gih_alpha,self.gih_x2)
        return self.sigmoid_t % p

    def giA_expr(self):
        p = (self.gka_start,self.gka_beta,self.gka_alpha)
        return self.exp_t % p

    def gkslow_expr(self):
        p = (self.gkslow_start,self.gkslow_beta,self.gkslow_alpha)
        return self.exp_t % p

    def gna_expr(self):
        p = (self.gna_soma,self.gna_api,self.gna_soma, self.dist_na)
        return self.lin_t % p

    def gsk_expr(self):
        p = (self.gsk_soma,self.gsk_dend,self.gsk_soma, self.dist_sk)
        return self.lin_t % p

    def gbk_expr(self):
        p = (self.gbk_soma,self.gbk_dend,self.gbk_soma, self.dist_bk)
        return self.lin_t % p

if __name__ == '__main__':
    p = AlmogKorngreenPars()
    p.pars_from_file('best.params')
    print p.giH_expr()
    print p.giA_expr()
    print p.gkslow_expr()
    print p.gna_expr()
    print p.gsk_expr()
    print p.gbk_expr()
