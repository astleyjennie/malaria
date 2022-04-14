# Load packages ####
library(pacman)
p_load(deSolve, tidyverse, doParallel, manipulate, readxl)

# Plot real case data ####

data <- as_tibble(as.data.frame(read_excel("populations_and_prevalence.xlsx", sheet="Case_data_Thailand", range="B20:D31", col_names=TRUE)))
cases <- reshape2::melt(data, id.var='Year')

(ggplot(cases, aes(x=Year, y=value, col=variable)) %>% 
    + geom_line() %>% 
    + xlim(2009,2021)
  + ylim(0,1500000)
  + labs(title="Assumed cases by species", x ="Year", y = "Cases"))

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
# Run the model ####
output <- ode(times = times, y = istate, func = thailand_model, parms = parameters)

# Manipulate data ####

df1<-as_tibble(as.data.frame(output)) %>% 
  mutate(PThaif = (Sf+Ef+Af+Cf+Tf+Rf+E2f),
         PThaiv = (Sv+Ev+Av+Cv+TPv+TAv+Lv+Rv+E2v),
         Pf=(Outf+Testf+PThaif),
         Pv=(Outv+Testv+PThaiv),
         If = (Ef+Af+Cf+Tf),
         Iv = (Ev+Av+Cv+TPv+TAv),
         P = PThaif+PThaiv,
         CTrtv=CACT+CPrim,
         Incf = c(0, diff(CIncf)),
         Incv = c(0, diff(CIncv)),
         Trtf = c(0, diff(CTrtf)),
         Trtv = c(0, diff(CTrtv)),
         CostPrim = c(0, diff(CPrim)),
         CostACT = c(0,diff(CACT)),
         CostG6PD = c(0, diff(CCv)),
         CostRDTf = c(0, diff(CTestf)),
         CostRDTv = c(0, diff(CTestv))) %>% 
  pivot_longer(names_to = "variable", cols = !1) %>%
  mutate(SP = ifelse(str_ends(variable, "f"), "Pf", "Pv")
  )

# Yearly incidence from model to compare with data

predicted_cases <- as.data.frame(seq(from=2010,to=2037,by=1))
names(predicted_cases)[names(predicted_cases) == colnames(predicted_cases)[1]] <- "Year"

df_falc <- as.data.frame(df1 %>%  filter(variable %in% c("Incf")))$value
falc_by_year <- unname(tapply(df_falc, (seq_along(df_falc)-1) %/% 365, sum))
names(falc_by_year)[names(falc_by_year) == colnames(falc_by_year)[1]] <- "Predicted_Falciparum"

df_viv <- as.data.frame(df1 %>%  filter(variable %in% c("Incv")))$value
viv_by_year <- unname(tapply(df_viv, (seq_along(df_viv)-1) %/% 365, sum))
names(viv_by_year)[names(viv_by_year) == colnames(viv_by_year)[1]] <- "Predicted_Vivax"

predicted_cases$Predicted_Falciparum <- as.numeric(falc_by_year)
predicted_cases$Predicted_Vivax <- as.numeric(viv_by_year)

data_compare <- as.data.frame(read_excel("populations_and_prevalence.xlsx", sheet="Case_data_Thailand", range="B20:D48", col_names=TRUE))

cases_compare <- merge(data_compare,predicted_cases,by="Year")

all_cases1<-cases_compare %>% 
  pivot_longer(names_to = "variable", cols = !1) %>%
  mutate(species = ifelse(str_ends(variable, "um"), "Falciparum", "Vivax")
  )


# Prediction vs data ####

ggplot(all_cases1, aes(x=Year, y=value, group=variable)) +
  geom_line(aes(color=variable)) +
  facet_wrap(~species)

# Population checks ####
# Population of Thailand only
df1 %>% 
  filter(variable %in% c("PThaif", "PThaiv")) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value))+
  theme_minimal() +
  labs(title = "Populations", y =("population")) +
  facet_wrap(~SP)
# Population of total system
df1 %>% 
  filter(variable %in% c("Pf", "Pv")) %>% 
  ggplot()+
  geom_line(aes(x = time, y=value))+
  theme_minimal() +
  labs(title = "Populations", y =("population")) +
  facet_wrap(~SP)

tail(output)


# Primaquine costing ####
yearly_cost_prim <- as.data.frame(seq(from=2010,to=2037,by=1))
names(yearly_cost_prim)[names(yearly_cost_prim) == colnames(yearly_cost_prim)[1]] <- "Year"

df_cprim <- as.data.frame(df1 %>%  filter(variable %in% c("CostPrim")))$value
yearly_prim_cost <- unname(tapply(df_cprim, (seq_along(df_cprim)-1) %/% 365, sum))
names(yearly_prim_cost)[names(yearly_prim_cost) == colnames(yearly_prim_cost)[1]] <- "Primaquine_Cost"

yearly_cost_prim$Primaquine_Cost <- as.numeric(yearly_prim_cost)

ggplot(data=yearly_cost_prim, aes(x=Year, y=Primaquine_Cost, group=1)) +
  geom_line(colour="darkorchid") +
  labs(title = "Yearly primaquine cost", y =("USD"))

# ACT costing ####

yearly_cost_ACT <- as.data.frame(seq(from=2010,to=2037,by=1))
names(yearly_cost_ACT)[names(yearly_cost_ACT) == colnames(yearly_cost_ACT)[1]] <- "Year"

df_cACT <- as.data.frame(df1 %>%  filter(variable %in% c("CostACT")))$value
yearly_ACT_cost <- unname(tapply(df_cACT, (seq_along(df_cACT)-1) %/% 365, sum))
names(yearly_ACT_cost)[names(yearly_ACT_cost) == colnames(yearly_ACT_cost)[1]] <- "ACT_Cost"

yearly_cost_ACT$ACT_Cost <- as.numeric(yearly_ACT_cost)

ggplot(data=yearly_cost_ACT, aes(x=Year, y=ACT_Cost, group=1)) +
  geom_line(colour="deeppink") +
  labs(title = "Yearly ACT cost", y =("USD"))

# G6PD costing ####

yearly_cost_G6PD <- as.data.frame(seq(from=2010,to=2037,by=1))
names(yearly_cost_G6PD)[names(yearly_cost_G6PD) == colnames(yearly_cost_G6PD)[1]] <- "Year"

df_cG6PD <- as.data.frame(df1 %>%  filter(variable %in% c("CostG6PD")))$value
yearly_G6PD_cost <- unname(tapply(df_cG6PD, (seq_along(df_cG6PD)-1) %/% 365, sum))
names(yearly_G6PD_cost)[names(yearly_G6PD_cost) == colnames(yearly_G6PD_cost)[1]] <- "G6PD_Cost"

yearly_cost_G6PD$G6PD_Cost <- as.numeric(yearly_G6PD_cost)

ggplot(data=yearly_cost_G6PD, aes(x=Year, y=G6PD_Cost, group=1)) +
  geom_line(colour="orange") +
  labs(title = "Yearly G6PD cost", y =("USD"))

# RDT costing ####

yearly_cost_RDT <- as.data.frame(seq(from=2010,to=2037,by=1))
names(yearly_cost_RDT)[names(yearly_cost_RDT) == colnames(yearly_cost_RDT)[1]] <- "Year"

df_cRDT <- as.data.frame(df1 %>%  filter(variable %in% c("CostRDTf")))$value
yearly_RDT_cost <- unname(tapply(df_cRDT, (seq_along(df_cRDT)-1) %/% 365, sum))
names(yearly_RDT_cost)[names(yearly_RDT_cost) == colnames(yearly_RDT_cost)[1]] <- "RDT_Cost"

yearly_cost_RDT$RDT_Cost <- as.numeric(yearly_RDT_cost)

ggplot(data=yearly_cost_RDT, aes(x=Year, y=RDT_Cost, group=1)) +
  geom_line(colour="red") +
  labs(title = "Yearly RDT cost", y =("USD"))








# Human Compartments ####
`%nin%` = Negate(`%in%`)

df1 %>% 
  filter(variable %nin% c("PThaif", "PThaiv","P", "CIncf", "CIncv", "Incf", "Incv", "Trtv", "Trtf", "CTrtf", "CTrtv")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Human Compartments", y =("population"), colour="Compartment")+
  facet_wrap(~SP)
# Daily incidence ####

df1 %>% 
  filter(variable %in% c("Incf", "Incv")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Incidence", y =("population"), colour="Compartment")+
  facet_wrap(~SP)

# Treatment ####
df1 %>% 
  filter(variable %in% c("Trtf", "Trtv")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Treated cases", y =("population"), colour="Compartment")+
  facet_wrap(~SP)

# Cases ####
# Plot real case data

data <- as_tibble(as.data.frame(read_excel("populations_and_prevalence.xlsx", sheet="Case_data_Thailand", range="C20:E31", col_names=TRUE)))
comparison <- merge(data, df1, by="time")
comparison <- reshape2::melt(data, id.var='Days')

(ggplot(cases, aes(x=Year, y=value, col=variable)) %>% 
    + geom_line() %>% 
    + xlim(2009,2021)
  + ylim(0,1500000)
  + labs(title="Assumed cases by species", x ="Year", y = "Cases"))



df1 %>% 
  filter(variable %in% c("If", "Iv")) %>% 
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = time, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Cases", y =("population"), colour="Compartment")+
  facet_wrap(~SP)

tail(output)