anderson<-function(X)
{
  
  ###############################
  # Parameters dispatching
  ###############################
  
  ex_D1in = X[1]    # POC export, mg C m-2 d-1
  ex_DOCin = X[2]    # DOC export: direct, mg C m-2 d-1
  y_B = X[3] # Psi
  w_att = X[4] # PGE attached
  w_fl = X[5] # PGE Free Living
  alpha = X[6] 
  cf_att = X[7] # Leu conversion factor for attached bacteria
  cf_fl = X[8] # Leu conversion factor for Free Living bacteria
  
  vatt_DOC = 0.05
  vatt_beta = 0.72
  vatt_npe = 0.44
  vfl_DOC = 0.05
  vfl_beta = 0.72
  vfl_npe = 0.44
  z_DOC = 0.05
  z_beta = 0.66
  z_D2 = 0.15
  z_npe = 0.39
  h_DOC = 0.05
  h_beta = 0.6
  h_D2 = 0.3
  h_npe = 0.39
  zi = 0.48 
  zi2 = 0.48 
  
  # All the DOC is accounted by ex_DOCin here
  ex_act = 0
  
  ###############################
  # Initialise state variables
  ###############################
  
  D1 <- 1.0    # sinking POC
  D2 <- 1.0    # suspended POC
  DOC <- 1.0   # DOC
  Batt <- 1.0  # attached prokaryotes: D1
  Bf <- 1.0    # free-living prokaryotes
  Vatt <- 1.0  # attached prokaryote consumers: D1
  Vf  <- 1.0   # free-living prokaryote consumers
  H <- 1.0     # detritivores
  Z <- rep(1.0, times=6)   #6 trophic levels of Z
  Batt2 <- 0.5 # attached prokaryotes: D2
  Vatt2 <- 0.5 # attached prokaryote consumers: D2
  
  # 2D array to hold fluxes 14 cols for the different state vars (incl. 6Z)
  flux <- matrix(rep(0.0, times=160), nrow=10, ncol=16)
  
  tstep <- 0.1    # time step (day)
  tfin <- 100     # run duration (days)
  tloop <- tfin/10
  rate <- 0.5     # coeff for flows: arbitrary
  
  # Store intermediary outputs
  #RZ.store = rep(NA, tfin * 10)
  #RB.store = rep(NA, tfin * 10)
  #Production_NonSinking.store = rep(NA, tfin * 10)
  #Production_Sinking.store = rep(NA, tfin * 10)
  
  for (t in seq(1,10)) {
    for (t2 in seq(tstep,tloop, by=tstep)) { #1 time loop
      
      # calculate rates of change
      
      HgrazingD <- rate*D1*(1.0-y_B)                         # grazing by detritivores on detritus
      HgrazingBatt <- rate*(1.0-y_B)*zi*Batt                 # grazing by detritivores on attached prokaryotes
      HgrazingVatt <- rate*(1.0-y_B)*zi2*Vatt                # grazing by detritivores on attached prokaryote consumers
      Hgrazing <- HgrazingD+HgrazingBatt+HgrazingVatt        # grazing by detritivores: total
      ZgrazingVf <- rate*Vf                                  # grazing by carnivores on free-living prokaryote consumers
      ZgrazingH <- rate*H                                    # grazing by carnivores on detritivores
      ZgrazingZ <- rate*(Z[1]+Z[2]+Z[3]+Z[4]+Z[5])           # grazing by carnivores on carnivores
      Zgrazing <- ZgrazingVf+ZgrazingH+ZgrazingZ             # grazing by carnivores
      Vattgrazing <- rate*(1.0-(1.0-y_B)*zi)*Batt            # grazing by attached prokaryote consumers
      Vflgrazing <- rate*Bf                                  # grazing by free-living prokaryote consumers
      grazingZ2 <- rate*Z[1]                                 # grazing by Z2 on Z1
      grazingZ3 <- rate*Z[2]                                 # grazing by Z3 on Z2
      grazingZ4 <- rate*Z[3]                                 # grazing by Z4 on Z3
      grazingZ5 <- rate*Z[4]                                 # grazing by Z5 on Z4
      grazingZ6 <- rate*Z[5]                                 # grazing by Z6 on Z5
      Vatt2grazing <- rate*Batt2                              # grazing by attached prokaryote consumers: D2
      
      D1toDOCsolub <- rate*D1*y_B*alpha                      # solubilization D1 to DOC
      D1toBatt <- rate*D1*y_B*(1.0-alpha)                    # uptake D1 by attached prokaryotes
      D2toDOCsolub <- rate*D2*alpha                          # solubilization D2 to DOC
      D2toBatt2 <- rate*D2*(1.0-alpha)                        # uptake D2 by attached prokaryotes
      DOCuptakeBfl <- rate*DOC                               # uptake DOC by free-living prokaryotes
      
      # Detritus D1 (sinking)
      flux[1,1] <- ex_D1in                    # detritus export from euphotic zone
      flux[2,1] <- -D1toDOCsolub              # solubilization to DOC
      flux[3,1] <- -D1toBatt                  # uptake by attached B
      flux[4,1] <- -HgrazingD                 # grazing by detritivores
      flux[5,1] <- Hgrazing*(1.0-h_D2-h_DOC)*(1.0-h_beta)  # faecal pellets: detritivores 
      flux[6,1] <- Zgrazing*(1.0-z_D2-z_DOC)*(1.0-z_beta)  # faecal pellets: carnivores
      
      # Detritus D2 (colloidal)
      flux[1,2] <- -D2toDOCsolub              # solubilisation to DOC (no consumption by detritivores)
      flux[2,2] <- -D2toBatt2                  # uptake by attached B
      flux[3,2] <- Hgrazing*h_D2              # from detritivores 
      flux[4,2] <- Zgrazing*z_D2              # from Z
      flux[5,2] <- (Vattgrazing+Vatt2grazing)*(1.0-vatt_beta)*(1.0-vatt_DOC)    # faecal pellets: attached prokaryote consumers 
      flux[6,2] <- Vflgrazing*(1.0-vfl_beta)*(1.0-vfl_DOC)       # faecal pellets: free-living prokaryote consumers
      
      # DOC
      flux[1,3] <- ex_DOCin                  # DOC from euphotic zone (direct flux)
      flux[2,3] <- ex_act                    # DOC from eupthotic zone (active flux)
      flux[3,3] <- D1toDOCsolub              # solubilization D1 to DOC
      flux[4,3] <- D2toDOCsolub              # solubilization D2 to DOC
      flux[5,3] <- Hgrazing*h_DOC            # from detritivores 
      flux[6,3] <- Zgrazing*z_DOC            # from Z
      flux[7,3] <- (Vattgrazing+Vatt2grazing)*vatt_DOC      # from attached prokaryote consumers 
      flux[8,3] <- Vflgrazing*vfl_DOC        # from free living prokaryote consumers
      flux[9,3] <- -DOCuptakeBfl             # uptake by free living prokaryotes
      
      # Attached prokaryotes
      flux[1,4] <- D1toBatt                  # utilisation of D1 and D2 (= carbon demand)
      flux[2,4] <- -flux[1,4]*(1.0-w_att)    # respiration
      flux[3,4] <- -HgrazingBatt             # loss to detritivores
      flux[4,4] <- -Vattgrazing              # loss to attached prokaryote consumers
      
      # Free living prokaryotes
      flux[1,5] <- DOCuptakeBfl              # consumption (= carbon demand)
      flux[2,5] <- -flux[1,5]*(1.0-w_fl)     # respiration
      flux[3,5] <- -Vflgrazing               # loss to free living bacterivores
      
      # Attached prokaryote consumers
      flux[1,6] <- Vattgrazing               # consumption of attached prokaryotes (= carbon demand)
      flux[2,6] <- -Vattgrazing*(1.0-vatt_DOC)*vatt_beta*(1.0-vatt_npe) # respiration
      flux[3,6] <- -Vattgrazing*vatt_DOC     # excretion of DOC
      flux[4,6] <- -Vattgrazing*(1.0-vatt_beta)*(1.0-vatt_DOC)          # faecal pellets to D2
      flux[5,6] <- -HgrazingVatt             # grazing by detritivores
      flux[6,6] <- -rate*(1.0-(1.0-y_B)*zi2)*Vatt      # closure respiration
      
      # Free living prokaryote consumers
      flux[1,7] <- Vflgrazing                # consumption of free-living prokaryotes (= carbon demand)
      flux[2,7] <- -flux[1,7]*(1.0-vfl_DOC)*vfl_beta*(1.0-vfl_npe)      # respiration
      flux[3,7] <- -Vflgrazing*vfl_DOC       # excretion of DOC
      flux[4,7] <- -Vflgrazing*(1.0-vfl_beta)*(1.0-vfl_DOC)             # faecal pellets to D2
      flux[5,7] <- -ZgrazingVf               # grazing by Z1
      
      # Detritivores
      flux[1,8] <- Hgrazing                  # grazing (= carbon demand)
      flux[2,8] <- -Hgrazing*(1.0-h_DOC-h_D2)*h_beta*(1.0-h_npe)        # respiration
      flux[3,8] <- -Hgrazing*h_DOC           # excretion of DOC
      flux[4,8] <- -Hgrazing*(1.0-h_D2-h_DOC)*(1.0-h_beta)              # faecal pellets to D1
      flux[5,8] <- -Hgrazing*h_D2            # losses to D2
      flux[6,8] <- -ZgrazingH                # grazing by Z1
      
      # Zooplankton 1   (1st in chain of 6)
      flux[1,9] <- ZgrazingVf+ZgrazingH      # consumption of free-living prokaryote consumers and detritivores (= carbon demand)
      flux[2,9] <- -(ZgrazingVf+ZgrazingH)*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[3,9] <- -(ZgrazingVf+ZgrazingH)*z_DOC   # excretion of DOC
      flux[4,9] <- -(ZgrazingVf+ZgrazingH)*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,9] <- -(ZgrazingVf+ZgrazingH)*z_D2    # faecal pellets to D2
      flux[6,9] <- -grazingZ2                      # grazing by Z2
      
      # Zooplankton 2
      flux[1,10] <- grazingZ2    # consumption (= carbon demand)
      flux[2,10] <- -grazingZ2*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[3,10] <- -grazingZ2*z_DOC   # excretion DOC
      flux[4,10] <- -grazingZ2*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,10] <- -grazingZ2*z_D2    # losses to D2
      flux[6,10] <- -grazingZ3   # grazing by Z3
      
      # Zooplankton 3
      flux[1,11] <- grazingZ3   # consumption (= carbon demand)
      flux[2,11] <- -grazingZ3*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[3,11] <- -grazingZ3*z_DOC   # excretion of DOC
      flux[4,11] <- -grazingZ3*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,11] <- -grazingZ3*z_D2    # losses to D2
      flux[6,11] <- -grazingZ4   # grazing by Z4
      
      # Zooplankton 4
      flux[1,12] <- grazingZ4   # consumption (= carbon demand)
      flux[2,12] <- -grazingZ4*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[3,12] <- -grazingZ4*z_DOC   # excretion of DOC
      flux[4,12] <- -grazingZ4*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,12] <- -grazingZ4*z_D2    # losses to D2
      flux[6,12] <- -grazingZ5   # grazing by Z5
      
      # Zooplankton 5
      flux[1,13] <- grazingZ5   # consumption (= carbon demand)
      flux[2,13] <- -grazingZ5*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[3,13] <- -grazingZ5*z_DOC   # excretion of DOC
      flux[4,13] <- -grazingZ5*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,13] <- -grazingZ5*z_D2    # flosses to D2
      flux[6,13] <- -grazingZ6   # grazing by Z5
      
      # Zooplankton 6   (last in chain)
      flux[1,14] <- grazingZ6   # consumption (= carbon demand)
      flux[2,14] <- -grazingZ6*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)-0.5*Z[6]   # respiration
      flux[3,14] <- -grazingZ6*z_DOC   # excretion of DOC
      flux[4,14] <- -grazingZ6*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[5,14] <- -grazingZ6*z_D2    # flosses to D2
      
      # Attached prokaryotes: D2
      flux[1,15] <- D2toBatt2                  # utilisation of D1 and D2 (= carbon demand)
      flux[2,15] <- -flux[1,15]*(1.0-w_fl)    # respiration. ! Changed compared to the original version which was -flux[1,15]*(1.0-w_att)
      flux[3,15] <- 0.0                        # loss to detritivores
      flux[4,15] <- -Vatt2grazing              # loss to attached prokaryote consumers
      
      # Attached prokaryote consumers: D2
      flux[1,16] <- Vatt2grazing               # consumption of attached prokaryotes (= carbon demand)
      flux[2,16] <- -Vatt2grazing*(1.0-vatt_DOC)*vatt_beta*(1.0-vatt_npe) # respiration
      flux[3,16] <- -Vatt2grazing*vatt_DOC     # excretion of DOC
      flux[4,16] <- -Vatt2grazing*(1.0-vatt_beta)*(1.0-vatt_DOC)          # faecal pellets to D2
      flux[5,16] <- 0.0                       # grazing by detritivores
      flux[6,16] <- -rate*Vatt2               # closure respiration
      
      # update state variables
      D1 <- D1 + sum(flux[1:10])*tstep          # detritus D1
      D2 <- D2 + sum(flux[11:20])*tstep         # detritus D2
      DOC <- DOC + sum(flux[21:30])*tstep       # DOC
      Batt <- Batt + sum(flux[31:40])*tstep     # attached prokaryotes
      Bf <- Bf + sum(flux[41:50])*tstep         # free-living prokaryotes
      Vatt <- Vatt + sum(flux[51:60])*tstep     # attached prokaryote consumers
      Vf <- Vf + sum(flux[61:70])*tstep         # free-living prokaryote consumers
      H <- H + sum(flux[71:80])*tstep           # detritivores
      Z[1] <- Z[1] + sum(flux[81:90])*tstep     # Z1
      Z[2] <- Z[2] + sum(flux[91:100])*tstep    # Z2
      Z[3] <- Z[3] + sum(flux[101:110])*tstep   # Z3
      Z[4] <- Z[4] + sum(flux[111:120])*tstep   # Z4
      Z[5] <- Z[5] + sum(flux[121:130])*tstep   # Z5
      Z[6] <- Z[6] + sum(flux[131:140])*tstep   # Z6
      Batt2 <- Batt2 + sum(flux[141:150])*tstep # attached prokaryotes: D2
      Vatt2 <- Vatt2 + sum(flux[151:160])*tstep # attached prokaryote consumers: D2
      
      # Store intermediate output
      #Prod_Bfl <- DOCuptakeBfl*w_fl            # Free living prokaryotes         
      #Prod_Batt <- (D1toBatt+D2toBatt2)*w_att  # attached prokaryotes
      
      #Production_NonSinking.store[10 * ((t-1) * 10 + t2)] = (Prod_Bfl + (D2toBatt2 * w_fl)) * cf_fl * 1E6  # production of the Free living prokaryotes + attached to suspended particles prokaryotes 
      #Production_Sinking.store[10 * ((t-1) * 10 + t2)] = D1toBatt * w_att * cf_att * 1E6 # Production of prokaryotes attached to sinking particles 
      
      # respiration
      #RZ <- flux[2,8]+flux[2,9]+flux[2,10]+flux[2,11]+flux[2,12]+flux[2,13]+flux[2,14]  # detritivores and carnivores
      #RB <- flux[2,4]+flux[2,5]+flux[2,15]      # attached and free-living prokaryotes
      
      # Store intermediate output
      #RZ.store[10 * ((t-1) * 10 + t2)] = RZ
      #RB.store[10 * ((t-1) * 10 + t2)] = RB
      
      
    }  # time loop
    
    # respiration
    #RZ <- flux[2,8]+flux[2,9]+flux[2,10]+flux[2,11]+flux[2,12]+flux[2,13]+flux[2,14]  # detritivores and carnivores
    #RB <- flux[2,4]+flux[2,5]+flux[2,15]      # attached and free-living prokaryotes
    
    # Store intermediate output
    #RZ.store[t] = RZ
    #RB.store[t] = RB
  }  # time loop 1-10
  
  
  # Inputs of C to mesopelagic zone
  n <- flux[1,1]+flux[1,3]+flux[2,3]                # total C input (POC plus DOC)
  
  # Respiration
  ZR <- -(flux[2,9]+flux[2,10]+flux[2,11]+flux[2,12]+flux[2,13]+flux[2,14])   # carnivores
  BattR <- -flux[2,4]-flux[2,15]
  VattR <- -flux[2,6]-flux[6,6]-flux[2,16]-flux[6,16]  #4 attached prokaryote consumers (includes closure respiration)
  
  # Carbon demand
  ZCD <-flux[1,9]+flux[1,10]+flux[1,11]+flux[1,12]+flux[1,13]+flux[1,14]      # carnivores
  BattCD <- flux[1,4]+flux[1,15]
  VattCD <- flux[1,6]+flux[1,16]
  
  # Production: 
  Prod_Bfl <- DOCuptakeBfl*w_fl                     # free-living prokaryotes
  Prod_Batt <- (D1toBatt+D2toBatt2)*w_att            # attached prokaryotes
  Prod_Vfl <- Vflgrazing*(1.0-vfl_DOC)*vfl_beta*vfl_npe      # free-living prokaryote consumers
  Prod_Vatt <- (Vattgrazing+Vatt2grazing)*(1.0-vatt_DOC)*vatt_beta*vatt_npe # attached prokaryote consumers
  Prod_H <- Hgrazing*(1.0-h_DOC-h_D2)*h_beta*h_npe  # detritivores
  Prod_Z <- Zgrazing*(1.0-z_DOC-z_D2)*z_beta*z_npe  # carnivores
  
  # D1 sources
  ex_actD1 <- 0
  
  # D2 sources
  ex_actD2 <- 0
  
  # DOC sources: ex,act,sol,Vfl,Vatt,H,Z
  DOC_sol <- flux[3,3]+flux[4,3]
  
  #closure respi
  closure = -(flux[6,6]+flux[6,16])
  D1_Batt <- -flux[2,1]-flux[3,1] 
  D2_Batt <- -flux[1,2]-flux[2,2] 
  
  #################################################
  # Plot the quantities of interest
  #################################################
  #plot(Production_NonSinking.store, xlab = 'Iteration Loop',
       #ylab = 'Value', main = 'Production Non Sinking through model run')
  #plot(Production_Sinking.store, xlab = 'Iteration Loop',
       #ylab = 'Value', main = 'Production Sinking through model run')
  #plot(-RB.store, xlab = 'Iteration Loop',
       #ylab = 'Value', main = 'Respiration Prokaryotes through model run')
  #plot(-RZ.store, xlab = 'Iteration Loop',
       #ylab = 'Value', main = 'Respiration Zoo. through model run')
  
  
  ##################################################
  # Returns the 4 in situ measurements counterparts 
  ##################################################
  
  Respi_attached = -flux[2,4]-flux[2,15]
  Respi_zooplancton = (-flux[2,7]) + (-flux[2,6]-flux[6,6]-flux[2,16]-flux[6,16]) + (-flux[2,8]) + (-(flux[2,9]+flux[2,10]+flux[2,11]+flux[2,12]+flux[2,13]+flux[2,14]))
  
  Production_Free_living = Prod_Bfl 
  Production_attached = Prod_Batt
  
  # The production flows are converted using dedicated conversion factors (CF)
  # 1 CF for the attached prokaryotes and 1 CF for the Free Living prokaryotes
  Production_NonSinking = ((Prod_Bfl + (D2toBatt2 * w_fl)) / cf_fl) * 1E6 # production of the Free living prokaryotes + attached to suspended particles prokaryotes 
  Production_Sinking = ((D1toBatt * w_att) / cf_att) * 1E6 # Production of prokaryotes attached to sinking particles 
  
  Respiration_Sinking  = -flux[2,4]
  Respiration_NonSinking = -flux[2,5]-flux[2,15]
  
  Pordzoo = Prod_Z+Prod_Vfl+Prod_Vatt+Prod_H

  # Check that the model is at equilibrium
  (ex_D1in + ex_DOCin) - (Respiration_Sinking + Respiration_NonSinking + Respi_zooplancton) < 1E-6
  
  # Return the four outflows of interest
  c(Production_NonSinking, Production_Sinking, Respiration_Sinking, Respi_zooplancton)
  
}
