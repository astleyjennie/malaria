# Times ####
times <- seq(0, 10000, 1) # approx 27 years

# Initial conditions ####
# Based on 2010 numbers from World Malaria Report 2021

inits <- read_excel("populations_and_prevalence.xlsx", sheet="Initial_values_Thailand", range="B7:C27", col_names=TRUE, col_types=c("text","numeric"))

Testf_0 <- as.numeric(inits[1,2])
Outf_0 <- as.numeric(inits[2,2])
Sf_0 <- as.numeric(inits[3,2])
Ef_0 <- as.numeric(inits[4,2])
Af_0 <- as.numeric(inits[5,2])
Cf_0 <- as.numeric(inits[6,2])
Tf_0 <- as.numeric(inits[7,2])
Rf_0 <- as.numeric(inits[8,2])
E2f_0 <- as.numeric(inits[9,2])

Testv_0 <- as.numeric(inits[10,2])
Outv_0 <- as.numeric(inits[11,2])
Sv_0 <- as.numeric(inits[12,2])
Ev_0 <- as.numeric(inits[13,2])
Av_0 <- as.numeric(inits[14,2])
Cv_0 <- as.numeric(inits[15,2])
TPv_0 <- as.numeric(inits[16,2]) 
TAv_0 <- as.numeric(inits[17,2])
Lv_0 <- as.numeric(inits[18,2])
Rv_0 <- as.numeric(inits[19,2])
E2v_0 <- as.numeric(inits[20,2])

CIncf_0 <- 0
CIncv_0 <- 0
CTrtf_0 <- 0
CPrim_0 <- 0
CACT_0 <- 0
CCv_0 <- 0
CTestf_0 <-0
CTestv_0 <- 0

istate <- c(Testf=Testf_0, Outf=Outf_0, Sf = Sf_0, Ef = Ef_0, Af = Af_0, Cf = Cf_0, Tf=Tf_0, Rf = Rf_0, E2f=E2f_0,
            Testv=Testv_0, Outv=Outv_0, Sv = Sv_0, Ev = Ev_0, Av = Av_0, Cv = Cv_0, TPv=TPv_0, TAv=TAv_0, Lv=Lv_0, Rv = Rv_0, E2v=E2v_0,
            CIncf=CIncf_0, CIncv=CIncv_0, CTrtf=CTrtf_0, CPrim=CPrim_0, CACT=CACT_0, CCv=CCv_0, CTestf=CTestf_0, CTestv=CTestv_0 )

# Parameters ####

params <- as_tibble(as.data.frame(read_excel("parameters.xlsx", range="A1:B80", col_names=TRUE, col_types=c("text","numeric")))) # Read excel file with all parameters

parameters <- rep(0,length(params$PARAMETER)) # create empty parameter vector

for (i in 1:length(params$PARAMETER)){ # populate parameter vector with parameter names and values
  parameters[i] <- params$VALUE[i]
  names(parameters) <- params$PARAMETER
}

#Read in reduction data
net_data <- read_excel("populations_and_prevalence.xlsx", sheet="Reduction_parameter", range="A1:D30", col_names=TRUE)


# Define model function ####

thailand_model<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # Populations
         
         PThaif <- (Sf + Ef + Af + Cf + Tf + Rf + E2f) # Total population of Thailand falciparum transmission system
         PThaiv <- (Sv + Ev + Av + Cv + TPv + TAv + Lv + Rv + E2v) # Total population of Thailand vivax transmission system
         
         Tv <- (TPv + TAv) # Total treated for vivax: ACT and primaquine
         
         mu1 <- g * mu2 # Birth rate is higher than death rate
         
         nets<-approx(net_data$Day, net_data$nets, t)$y # Reading artificial decreasing function from excel
         #nets2<-approx(net_data$Day, net_data$nets2, t)$y # Reading artificial decreasing function from excel
         
         # Infection variables
         
         Infectiousf <- Cf + zeta_af * Af + zeta_tf * Tf # Relative infectiousness for asymptomatic and treated considered
         Infectiousv <- Cv + zeta_av * Av + zeta_tv * Tv
         
         seas <- 1 + amp * cos(2 * pi * (t / 365 - phi)) ^ peak # Seasonality
         
         # Force of infection equations
         
         lambdaf <- nets * seas*(a^2*b*c*m*Infectiousf/PThaif)/(a*c*Infectiousf/PThaif+mu_m)*(gam_mf/(gam_mf+mu_m))
         lambdav <- relv * nets * seas*(a^2*b*c*m*Infectiousv/PThaiv)/(a*c*Infectiousv/PThaiv+mu_m)*(gam_mv/(gam_mv+mu_m))
         
         
         # Falciparum
         
         dTestf <- ( t1 * nu3 * Outf # Individuals being tested at the border
                     - t6 * nu1 * Testf # Individuals entering Thailand transmission cycle with border testing on
                     - (1 - t6) * nu1 * Testf ## Individuals entering Thailand transmission cycle with border testing off
         )
         
         dOutf <- (  t1 * nu2 * PThaif # Leaving Thailand and reentering 'outside
                     + (1 - t8) * t6 * ppf * nu1 * Testf # Imported cases turned away at border
                     - t1 * nu3 * Outf ) # Entering test compartment
         
         dSf <- (  mu1 * PThaif # Birth
                   + t6 * nu1 * (1 - ppf) * Testf # Negative border screening result entering susceptible
                   + (1 - t6) * (1 - pfal) * nu1 * Testf # No screening, disease-free individuals enter susceptible 
                   - lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Sf # Individuals becoming exposed with cross immunity
                   + rhof * Rf # Loss of immunity
                   - t1 * nu2 * Sf # Leaving Thailand
                   - mu2 * Sf ) # Natural death
         
         dEf <- (  lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Sf # Individuals becoming infected with cross immunity
                   - gam_hf * Ef # Infection developing into asymptomatic or clinical
                   - t1 * nu2 * Ef # Leaving Thailand
                   - mu2 * Ef ) # Natural death
         
         dAf <- (  pa * gam_hf * Ef # Asymptomatic disease
                   + pa2 * gam_hf * E2f # Asymptomatic disease from secondary infection
                   + omegaf * Cf # Natural improvement from clinical to asymptomatic disease
                   - phif * Af # Natural recovery from asymptomatic disease
                   - t1 * nu2 * Af # Leaving Thailand
                   - mu2 * Af ) # Natural death
         
         dCf <- (  (1 - pa) * gam_hf * Ef # Exposed individuals developing clinical disease
                   + (1 - pa2) * gam_hf * E2f # Exposed individuals developing clinical disease from secondary infection
                   + (1 - t6) * pfal * nu1 * Testf # With no border screening, infected individuals from border are clinical
                   - t3 * rv * (Tv / PThaiv) * Cf # Treatment for vivax treats falciparum
                   - omegaf * Cf # Improvement from clinical to asymptomatic disease
                   - tauf * Cf # Rate of receiving ACT treatment
                   - t1 * nu2 * Cf # Leaving Thailand
                   - mu2 * Cf ) # Natural death
         
         dTf <- (  tauf * Cf # Receiving ACT treatment
                   + t8 * t6 * ppf * nu1 * Testf # Individuals screened at the border with a positive result receive treatment
                   - rf * Tf # Recovering after ACT treatment
                   - t1 * nu2 * Tf # Leaving Thailand
                   - mu2 * Tf ) # Natural death
         
         dRf <- (  rf * Tf # Recovering after ACT treatment
                   + phif * Af # Natural recovery from asymptomatic disease
                   + t3 * rv * (Tv / PThaiv) * Cf # Treatment for vivax treats falciparum
                   - lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Rf # Reinfection with cross immunity
                   - rhof * Rf # Loss of immunity
                   - t1 * nu2 * Rf # Leaving Thailand
                   - mu2 * Rf ) # Natural death
         
         dE2f <- (  lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Rf # Reinfection with cross immunity
                    - gam_hf * E2f # Developing asymptomatic or clinical disease from secondary infection
                    - t1 * nu2 * E2f # Leaving Thailand
                    - mu2 * E2f ) # Natural death
         
         # Vivax
         
         dTestv <-( t1 * nu3 * Outv # Individuals from outside Thailand arriving at border screening
                    - (1 - t6) * nu1 * Testv # Border screening off, individuals entering Sv or Ev depending on disease status
                    - t6 * nu1 * (1 - ppv) * Testv # Border screening on, disease free individuals entering susceptible
                    - t7 * psi * t6 * nu1 * ppv * Testv # Border screening and G6PD screening on, deficient vivax-positive individuals receiving ACT treatment
                    - t7 * (1 - psi) * t6 * nu1 * ppv * Testv # Border screening and G6P screening on, non-deficient vivax-positive individuals receive primaquine
                    - (1 - t7) * t6 * nu1 * ppv * Testv ) # Border screening, no G6PD. All vivax-positive individuals receive primaquine
         
         dOutv <- (t1 * nu2 * PThaiv # Individuals leaving Thailand and reentering 'outside'
                   + (1 - t8) * t6 * nu1 * ppv * Testv # Individuals testing positive at the border are turned away
                   - t1 * nu3 * Outv) # Individuals entering border testing
         
         dSv <- (  mu1 * PThaiv # Birth
                   - lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Sv # Infection with cross immunity
                   + (1 - t6) * nu1 * (1 - pviv) * Testv # No border screening, vivax-negative individuals enter susceptible population
                   + t6 * nu1 * (1 - ppv) * Testv # Border screening on, individuals with negative-vivax result enter susceptible population
                   + rhov * Rv # Loss of immunity
                   + dhyp * Lv # Death of hypnozoites
                   - t1 * nu2 * Sv # Leaving Thailand
                   - mu2 * Sv ) # Natural death
         
         dEv <- ( lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Sv # Infection with cross immunity
                  - gam_hv * Ev # Developing asymptomatic or clinical disease
                  - t1 * nu2 * Ev  # Leaving Thailand
                  - mu2 * Ev ) # Natural death
         
         dAv <- (  pa * gam_hv * Ev # Developing asymptomatic disease
                   + pa2 * gam_hv * E2v # Developing asymptomatic disease from secondary infection
                   + omegav * Cv # Natural improvement from clinical to asymptomatic
                   - phiv * Av # Natural recovery from asymptomatic disease
                   - t1 * nu2 * Av # Leaving Thailand
                   - mu2 * Av ) # Natural death
         
         dCv <- (  (1 - pa) * gam_hv * Ev # Developing clinical disease
                   + (1 - pa2) * gam_hv * E2v # Developing clinical disease from secondary infection
                   + (1 - t6) * nu1 * pviv * Testv # No border screening, individuals with vivax enter clinical compartment
                   - t3 * rf * (Tf / PThaif) * Cv # Treatment for falciparum treats vivax
                   - omegav * Cv # Natural improvement from clinical to asymptomatic
                   - (1 - t7) * tauv * Cv # All vivax cases treated with primaquine if G6PD testing is switched off
                   - t7 * tauv * (1 - psi) * Cv # Only non-deficient treated with primaquine if G6PD testing is switched on
                   - t7 * tauf * psi * Cv # Deficient individuals treated with ACT when G6PD testing is switched on
                   - t1 * nu2 * Cv # Leaving Thailand
                   - mu2 * Cv ) # Natural death
         
         dTPv <- ( (1 - t7) * tauv * Cv # All vivax cases treated with primaquine if G6PD testing is switched off
                   + t7 * tauv * (1 - psi) * Cv # Only non-deficient treated with primaquine if G6PD testing is switched on
                   - rp * TPv # Recovering with or without hypnozoites with primaquine
                   + t8 * (1 - t7) * t6 * nu1 * ppv * Testv # No G6PD screening. Border-vivax-positive individuals all receive primaquine
                   + t8 * t7 * (1 - psi) * t6 * nu1 * ppv * Testv # Non-deficient border-vivax-positive individuals receive primaquine
                   - t1 * nu2 * TPv # Leaving Thailand
                   - mu2 * TPv ) # Natural death
         
         dTAv <- (  t7 * tauf * psi * Cv # Deficient individuals treated with ACT when G6PD testing is switched on)
                    - t7 * rv * TAv # Recovering with or without hypnozoites after ACT treatment (no primaquine))
                    + t8 * t7 * psi * t6 * nu1 * ppv * Testv # Border and G6PD screening on, vivax-positive  deficient individuals receive ACT
                    - t1 * nu2 * TAv # Leaving Thailand
                    - mu2 * TAv ) # Natural death
         
         dLv <-  (  prelp * rp * TPv # Recovery with hypnozoites after primaquine treatment
                    + t7 * prel * rv * TAv # Recovery with hypnozoites after ACT treatment  
                    + t3 * prel * rf * (Tf / PThaif) * Cv # Dual treatment: treatment for falciparum treats vivax
                    - rel * Lv # Relapse to secondary infection
                    - t2 * increl * (Tf / PThaif) * Lv # Falciparum infection triggering vivax relapse
                    - dhyp * Lv # Death of hypnozoites
                    - t1 * nu2 * Lv # Leaving Thailand
                    - mu2 * Lv ) # Natural death
         
         dRv <- (   (1 - prelp) * rp * TPv # Recovery without hypnozoites after primaquine treatment
                    + t7 * (1 - prel) * rv * TAv # Recovery without hypnozoites after ACT treatment
                    + phiv * Av # Recovery without hypnozoites from asymptomatic infection
                    + t3 * (1 - prel) * rf * (Tf / PThaif) * Cv # Dual treatment: treatment for falciparum treats vivax
                    - rhov * Rv # Loss of immunity
                    - lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Rv # Secondary infection with cross-immunity
                    - t1 * nu2 * Rv # Leaving Thailand
                    - mu2 * Rv ) # Natural death
         
         dE2v <- (  lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Rv # Reinfection with cross immunity
                    + rel * Lv # Relapse due to hypnozoites
                    + t2 * increl * (Tf / PThaif) * Lv # Falciparum infection triggering vivax relapse 
                    - gam_hv * E2v # Developing asymptomatic or clinical disease from secondary infection
                    - t1 * nu2 * E2v # Leaving Thailand
                    - mu2 * E2v ) # Natural death
         
         
         # Counters
         dCIncf <- (  lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Sf
                      + lambdaf * (1 - t4 * alpha21 * Infectiousv / PThaiv) * Rf )
         
         dCIncv <- (  lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Sv
                      + lambdav * (1 - t4 * alpha12 * Infectiousf / PThaif) * Rv
                      + rel * Lv
                      + t2 * increl * Tf / PThaif * Lv )
         
         dCTrtf <- (  tauf * Cf
                      + t3 * rv * ((TAv+TPv) / PThaiv) * Cf
                      + t8 * t6 * ppf * nu1 * Testf)
         
         # Cumulative cost of primaquine treatments
         
         dCPrim <- ( c_prim * ((1 - t7) * tauv * Cv # All vivax cases treated with primaquine if G6PD testing is switched off
                               + t7 * tauv * (1 - psi) * Cv # Only non-deficient individuals treated with primaquine if G6PD testing is switched on
                               + (1 - t7) * t6 * nu1 * ppv * t8 * Testv
                               + t7 * (1 - psi) * t6 * nu1 * ppv * t8 * Testv) ) 
         
         # Cumulative cost of ACT treatments
         
         dCACT <- ( c_ACT * ( t7 * tauf * psi * Cv # G6PD deficient individuals treated with ACT for vivax (if G6PD screening is on)
                              + tauf * Cf # All individuals treated with ACT for falciparum
                              + t3 * rf * (Tf / PThaif) * Cv # Dual treatment
                              + t7 * psi * t6 * nu1 * ppv * t8 * Testv)  ) 
         
         # Cumulative cost of G6PD testing
         dCCv <- ( t7 * c_G6PD * (  (1 - pa) * gam_hv * Ev # Cumulative clinical cases
                                    + (1 - pa2) * gam_hv * E2v
                                    + (1 - t6) * nu1 * pviv * Testv) )
         
         # Cumulative cost of RDTs. Code function check: the two below should be the same
         
         dCTestf <- t6 * c_RDT * (t1 * nu3 * Outf) # Cumulative tests
         
         dCTestv <- t6 * c_RDT * (t1 * nu3 * Outv) # Cumulative tests
         
         # return the rate of change
         list(c(dTestf, dOutf, dSf, dEf, dAf, dCf, dTf, dRf, dE2f,
                dTestv, dOutv, dSv, dEv, dAv, dCv, dTPv, dTAv, dLv, dRv, dE2v,
                dCIncf, dCIncv, dCTrtf, dCPrim, dCACT, dCCv, dCTestf, dCTestv ))
       }
  ) 
  
}