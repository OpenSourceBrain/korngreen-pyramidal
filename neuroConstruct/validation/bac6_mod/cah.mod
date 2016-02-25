TITLE Mod file for component: Component(id=cah type=MOD_cah)

COMMENT

    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.4.4
         org.neuroml.model   v1.4.4
         jLEMS               v0.9.8.4

ENDCOMMENT

NEURON {
    SUFFIX cah
    USEION ca READ cai, cao WRITE ica VALENCE 2 ? Assuming valence = 2 (Ca ion); TODO check this!!
    
    RANGE gion                           
    RANGE permeability                      : Will be changed when ion channel mechanism placed on cell!
    RANGE cai
    RANGE cao
    RANGE conductance                       : parameter
    
    RANGE g                                 : exposure
    
    RANGE fopen                             : exposure
    RANGE m_instances                       : parameter
    
    RANGE m_tau                             : exposure
    
    RANGE m_inf                             : exposure
    
    RANGE m_rateScale                       : exposure
    
    RANGE m_fcond                           : exposure
    RANGE m_timeCourse_TIME_SCALE           : parameter
    RANGE m_timeCourse_VOLT_SCALE           : parameter
    RANGE m_timeCourse_shift                : parameter
    RANGE m_timeCourse_shifth               : parameter
    RANGE m_timeCourse_qm                   : parameter
    RANGE m_timeCourse_qh                   : parameter
    
    RANGE m_timeCourse_t                    : exposure
    RANGE m_steadyState_TIME_SCALE          : parameter
    RANGE m_steadyState_VOLT_SCALE          : parameter
    RANGE m_steadyState_shift               : parameter
    RANGE m_steadyState_shifth              : parameter
    RANGE m_steadyState_qm                  : parameter
    RANGE m_steadyState_qh                  : parameter
    
    RANGE m_steadyState_x                   : exposure
    RANGE m_q10Settings_q10Factor           : parameter
    RANGE m_q10Settings_experimentalTemp    : parameter
    RANGE m_q10Settings_TENDEGREES          : parameter
    
    RANGE m_q10Settings_q10                 : exposure
    RANGE h_instances                       : parameter
    
    RANGE h_tau                             : exposure
    
    RANGE h_inf                             : exposure
    
    RANGE h_rateScale                       : exposure
    
    RANGE h_fcond                           : exposure
    RANGE h_timeCourse_TIME_SCALE           : parameter
    RANGE h_timeCourse_VOLT_SCALE           : parameter
    RANGE h_timeCourse_shift                : parameter
    RANGE h_timeCourse_shifth               : parameter
    RANGE h_timeCourse_qm                   : parameter
    RANGE h_timeCourse_qh                   : parameter
    
    RANGE h_timeCourse_t                    : exposure
    RANGE h_steadyState_TIME_SCALE          : parameter
    RANGE h_steadyState_VOLT_SCALE          : parameter
    RANGE h_steadyState_shift               : parameter
    RANGE h_steadyState_shifth              : parameter
    RANGE h_steadyState_qm                  : parameter
    RANGE h_steadyState_qh                  : parameter
    
    RANGE h_steadyState_x                   : exposure
    RANGE h_q10Settings_q10Factor           : parameter
    RANGE h_q10Settings_experimentalTemp    : parameter
    RANGE h_q10Settings_TENDEGREES          : parameter
    
    RANGE h_q10Settings_q10                 : exposure
    RANGE m_timeCourse_V                    : derived variable
    RANGE m_steadyState_V                   : derived variable
    RANGE m_tauUnscaled                     : derived variable
    RANGE h_timeCourse_V                    : derived variable
    RANGE h_steadyState_V                   : derived variable
    RANGE h_tauUnscaled                     : derived variable
    RANGE conductanceScale                  : derived variable
    RANGE fopenHHrates                      : derived variable
    RANGE fopenHHtauInf                     : derived variable
    RANGE fopenHHratesTau                   : derived variable
    RANGE fopenHHratesInf                   : derived variable
    RANGE fopenHHratesTauInf                : derived variable
    RANGE fopenHHInstantaneous              : derived variable
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    : bypass nrn default faraday const
    FARADAY = 96485.3 (coulomb)
    R = (k-mole) (joule/degC)
    
}

PARAMETER {
    
    permeability = 0  (cm/s)                       : Will be changed when ion channel mechanism placed on cell!
    
    conductance = 1.0E-5 (uS)
    m_instances = 2 
    m_timeCourse_TIME_SCALE = 1 (ms)
    m_timeCourse_VOLT_SCALE = 1 (mV)
    m_timeCourse_shift = -4.49601 
    m_timeCourse_shifth = -7.11157 
    m_timeCourse_qm = 4 
    m_timeCourse_qh = 2 
    m_steadyState_TIME_SCALE = 1 (ms)
    m_steadyState_VOLT_SCALE = 1 (mV)
    m_steadyState_shift = -4.49601 
    m_steadyState_shifth = -7.11157 
    m_steadyState_qm = 4 
    m_steadyState_qh = 2 
    m_q10Settings_q10Factor = 4 
    m_q10Settings_experimentalTemp = 297.15 (K)
    m_q10Settings_TENDEGREES = 10 (K)
    h_instances = 1 
    h_timeCourse_TIME_SCALE = 1 (ms)
    h_timeCourse_VOLT_SCALE = 1 (mV)
    h_timeCourse_shift = -4.49601 
    h_timeCourse_shifth = -7.11157 
    h_timeCourse_qm = 4 
    h_timeCourse_qh = 2 
    h_steadyState_TIME_SCALE = 1 (ms)
    h_steadyState_VOLT_SCALE = 1 (mV)
    h_steadyState_shift = -4.49601 
    h_steadyState_shifth = -7.11157 
    h_steadyState_qm = 4 
    h_steadyState_qh = 2 
    h_q10Settings_q10Factor = 2 
    h_q10Settings_experimentalTemp = 297.15 (K)
    h_q10Settings_TENDEGREES = 10 (K)
}

ASSIGNED {
    cai (mM)
    cao (mM)
    
    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel
    v (mV)
    celsius (degC)
    temperature (K)
    eca (mV)
    ica (mA/cm2)
    
    
    m_timeCourse_V                         : derived variable
    
    m_timeCourse_t (ms)                    : derived variable
    
    m_steadyState_V                        : derived variable
    
    m_steadyState_x                        : derived variable
    
    m_q10Settings_q10                      : derived variable
    
    m_rateScale                            : derived variable
    
    m_fcond                                : derived variable
    
    m_inf                                  : derived variable
    
    m_tauUnscaled (ms)                     : derived variable
    
    m_tau (ms)                             : derived variable
    
    h_timeCourse_V                         : derived variable
    
    h_timeCourse_t (ms)                    : derived variable
    
    h_steadyState_V                        : derived variable
    
    h_steadyState_x                        : derived variable
    
    h_q10Settings_q10                      : derived variable
    
    h_rateScale                            : derived variable
    
    h_fcond                                : derived variable
    
    h_inf                                  : derived variable
    
    h_tauUnscaled (ms)                     : derived variable
    
    h_tau (ms)                             : derived variable
    
    conductanceScale                       : derived variable
    
    fopenHHrates                           : derived variable
    
    fopenHHtauInf                          : derived variable
    
    fopenHHratesTau                        : derived variable
    
    fopenHHratesInf                        : derived variable
    
    fopenHHratesTauInf                     : derived variable
    
    fopenHHInstantaneous                   : derived variable
    
    fopen                                  : derived variable
    
    g (uS)                                 : derived variable
    rate_m_q (/ms)
    rate_h_q (/ms)
    
}

STATE {
    m_q  
    h_q  
    
}

INITIAL {
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
    m_q = m_inf
    
    h_q = h_inf
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=cah type=MOD_cah), from conductanceScaling; null
    ? Path not present in component, using factor: 1
    
    conductanceScale = 1 
    
    ? DerivedVariable is based on path: gatesHHrates[*]/fcond, on: Component(id=cah type=MOD_cah), from gatesHHrates; null
    ? Path not present in component, using factor: 1
    
    fopenHHrates = 1 
    
    ? DerivedVariable is based on path: gatesHHtauInf[*]/fcond, on: Component(id=cah type=MOD_cah), from gatesHHtauInf; Component(id=m type=gateHHtauInf)
    ? multiply applied to all instances of fcond in: <gatesHHtauInf> ([Component(id=m type=gateHHtauInf), Component(id=h type=gateHHtauInf)]))
    fopenHHtauInf = m_fcond * h_fcond ? path based
    
    ? DerivedVariable is based on path: gatesHHratesTau[*]/fcond, on: Component(id=cah type=MOD_cah), from gatesHHratesTau; null
    ? Path not present in component, using factor: 1
    
    fopenHHratesTau = 1 
    
    ? DerivedVariable is based on path: gatesHHratesInf[*]/fcond, on: Component(id=cah type=MOD_cah), from gatesHHratesInf; null
    ? Path not present in component, using factor: 1
    
    fopenHHratesInf = 1 
    
    ? DerivedVariable is based on path: gatesHHratesTauInf[*]/fcond, on: Component(id=cah type=MOD_cah), from gatesHHratesTauInf; null
    ? Path not present in component, using factor: 1
    
    fopenHHratesTauInf = 1 
    
    ? DerivedVariable is based on path: gateHHInstantaneous[*]/fcond, on: Component(id=cah type=MOD_cah), from gateHHInstantaneous; null
    ? Path not present in component, using factor: 1
    
    fopenHHInstantaneous = 1 
    
    fopen = conductanceScale  *  fopenHHrates  *  fopenHHtauInf  *  fopenHHratesTau  *  fopenHHratesInf  *  fopenHHratesTauInf  *  fopenHHInstantaneous ? evaluable
    g = conductance  *  fopen ? evaluable
    gion = permeability * fopen 
    ica = gion * ghk(v, cai, cao)
    
}

DERIVATIVE states {
    rates()
    m_q' = rate_m_q 
    h_q' = rate_h_q 
    
}

PROCEDURE rates() {
    
    m_timeCourse_V = v /  m_timeCourse_VOLT_SCALE ? evaluable
    m_timeCourse_t = (0.97/(cosh(0.032*( m_timeCourse_V + m_timeCourse_shift +26.31)))) *  m_timeCourse_TIME_SCALE ? evaluable
    m_steadyState_V = v /  m_steadyState_VOLT_SCALE ? evaluable
    m_steadyState_x = 1.092/(1+exp(-( m_steadyState_V + m_steadyState_shift +14.17)/9.76)) ? evaluable
    m_q10Settings_q10 = m_q10Settings_q10Factor ^((temperature -  m_q10Settings_experimentalTemp )/ m_q10Settings_TENDEGREES ) ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=m type=gateHHtauInf), from q10Settings; Component(id=null type=q10ExpTemp)
    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10ExpTemp)]))
    m_rateScale = m_q10Settings_q10 ? path based
    
    m_fcond = m_q ^ m_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=m type=gateHHtauInf), from steadyState; Component(id=null type=cah_m_inf_inf)
    m_inf = m_steadyState_x ? path based
    
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=m type=gateHHtauInf), from timeCourse; Component(id=null type=cah_m_tau_tau)
    m_tauUnscaled = m_timeCourse_t ? path based
    
    m_tau = m_tauUnscaled  /  m_rateScale ? evaluable
    h_timeCourse_V = v /  h_timeCourse_VOLT_SCALE ? evaluable
    h_timeCourse_t = (70/(cosh(0.047*( h_timeCourse_V + h_timeCourse_shifth -19.73)))) *  h_timeCourse_TIME_SCALE ? evaluable
    h_steadyState_V = v /  h_steadyState_VOLT_SCALE ? evaluable
    h_steadyState_x = 0.75/(1+exp(( h_steadyState_V + h_steadyState_shifth +22.63)/6.6)) ? evaluable
    h_q10Settings_q10 = h_q10Settings_q10Factor ^((temperature -  h_q10Settings_experimentalTemp )/ h_q10Settings_TENDEGREES ) ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=h type=gateHHtauInf), from q10Settings; Component(id=null type=q10ExpTemp)
    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10ExpTemp)]))
    h_rateScale = h_q10Settings_q10 ? path based
    
    h_fcond = h_q ^ h_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=h type=gateHHtauInf), from steadyState; Component(id=null type=cah_h_inf_inf)
    h_inf = h_steadyState_x ? path based
    
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=h type=gateHHtauInf), from timeCourse; Component(id=null type=cah_h_tau_tau)
    h_tauUnscaled = h_timeCourse_t ? path based
    
    h_tau = h_tauUnscaled  /  h_rateScale ? evaluable
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    rate_m_q = ( m_inf  -  m_q ) /  m_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    rate_h_q = ( h_inf  -  h_q ) /  h_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
        LOCAL z, eci, eco
        z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
        eco = co*efun(z)
        eci = ci*efun(-z)
        :high cao charge moves inward
        :negative potential charge moves inward
        ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
        if (fabs(z) < 1e-4) {
                efun = 1 - z/2
        }else{
                efun = z/(exp(z) - 1)
        }
}
