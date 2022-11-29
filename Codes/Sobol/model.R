anderson<-function(X, flux_out_nb = NULL)
{
  # Parameters dispatching
  y_B = X[,1]
  w_att = X[,2]
  w_fl = X[,3]
  alpha = X[,4]
  vatt_DOC = X[,5]
  vatt_beta = X[,6]
  vatt_npe = X[,7]
  vfl_DOC = X[,8]
  vfl_beta = X[,9]
  vfl_npe = X[,10]
  z_DOC = X[,11]
  z_beta = X[,12]
  z_D2 = X[,13]
  z_npe = X[,14]
  h_DOC = X[,15]
  h_beta = X[,16]
  h_D2 = X[,17]
  h_npe = X[,18]
  zi = X[,19]
  zi2 = X[,20]

  
  n_bootstrap = dim(X)[1]
   
  ex_D1in     <- 74     # POC export, mg C m-2 d-1
  ex_DOCin    <- 15     # DOC export: direct, mg C m-2 d-1
  ex_act      <- 3     # DOC export: via active flux, mg C m-2 d-1
  
  ###############################le debut de mon modele anderson et tang 2010  
  # Initialise state variables
  D1 <- rep(1.0, n_bootstrap)    # sinking POC
  D2 <- rep(1.0, n_bootstrap)     # suspended POC
  DOC <- rep(1.0, n_bootstrap)     # DOC
  Batt <- rep(1.0, n_bootstrap)    # attached prokaryotes: D1
  Bf <- rep(1.0, n_bootstrap)     # free-living prokaryotes
  Vatt <- rep(1.0, n_bootstrap)    # attached prokaryote consumers: D1
  Vf  <- rep(1.0, n_bootstrap)     # free-living prokaryote consumers
  H <- rep(1.0, n_bootstrap)       # detritivores
  Z <- matrix(1.0, n_bootstrap, 6)   #6 trophic levels of Z
  Batt2 <- rep(0.5, n_bootstrap)   # attached prokaryotes: D2
  Vatt2 <- rep(0.5, n_bootstrap) # attached prokaryote consumers: D2
  
  # 3D array to hold fluxes 14 cols for the different state vars (incl. 6Z)
  flux = array(0.0, dim=c(n_bootstrap, 10, 16))
  
  tstep <- 0.1    # time step (day)
  tfin <- 100     # run duration (days)
  tloop <- tfin/10
  
  rate <- 0.5     # coeff for flows: arbitrary
  
  
  for (t in seq(1,10)) {
    for (t2 in seq(tstep,tloop, by=tstep)) { #1 time loop
      
      # calculate rates of change
      
      HgrazingD <- rate*D1*(1.0-y_B)                         # grazing by detritivores on detritus
      HgrazingBatt <- rate*(1.0-y_B)*zi*Batt                 # grazing by detritivores on attached prokaryotes
      HgrazingVatt <- rate*(1.0-y_B)*zi2*Vatt                # grazing by detritivores on attached prokaryote consumers
      Hgrazing <- HgrazingD+HgrazingBatt+HgrazingVatt        # grazing by detritivores: total
      ZgrazingVf <- rate*Vf                                  # grazing by carnivores on free-living prokaryote consumers
      ZgrazingH <- rate*H                                    # grazing by carnivores on detritivores
      ZgrazingZ <- rate*(Z[,1]+Z[,2]+Z[,3]+Z[,4]+Z[,5])           # grazing by carnivores on carnivores
      Zgrazing <- ZgrazingVf+ZgrazingH+ZgrazingZ             # grazing by carnivores
      Vattgrazing <- rate*(1.0-(1.0-y_B)*zi)*Batt            # grazing by attached prokaryote consumers
      Vflgrazing <- rate*Bf                                  # grazing by free-living prokaryote consumers
      grazingZ2 <- rate*Z[,1]                                 # grazing by Z2 on Z1
      grazingZ3 <- rate*Z[,2]                                 # grazing by Z3 on Z2
      grazingZ4 <- rate*Z[,3]                                 # grazing by Z4 on Z3
      grazingZ5 <- rate*Z[,4]                                 # grazing by Z5 on Z4
      grazingZ6 <- rate*Z[,5]                                 # grazing by Z6 on Z5
      Vatt2grazing <- rate*Batt2                              # grazing by attached prokaryote consumers: D2
      
      D1toDOCsolub <- rate*D1*y_B*alpha                      # solubilization D1 to DOC
      D1toBatt <- rate*D1*y_B*(1.0-alpha)                    # uptake D1 by attached prokaryotes
      D2toDOCsolub <- rate*D2*alpha                          # solubilization D2 to DOC
      D2toBatt2 <- rate*D2*(1.0-alpha)                        # uptake D2 by attached prokaryotes
      DOCuptakeBfl <- rate*DOC   
      
      # Detritus D1 (sinking)
      flux[,1,1] <- ex_D1in                    # detritus export from euphotic zone
      flux[,2,1] <- -D1toDOCsolub              # solubilization to DOC
      flux[,3,1] <- -D1toBatt                  # uptake by attached B
      flux[,4,1] <- -HgrazingD                 # grazing by detritivores
      flux[,5,1] <- Hgrazing*(1.0-h_D2-h_DOC)*(1.0-h_beta)  # faecal pellets: detritivores 
      flux[,6,1] <- Zgrazing*(1.0-z_D2-z_DOC)*(1.0-z_beta)  # faecal pellets: carnivores
      
      # Detritus D2 (colloidal)
      flux[,1,2] <- -D2toDOCsolub              # solubilisation to DOC (no consumption by detritivores)
      flux[,2,2] <- -D2toBatt2                  # uptake by attached B
      flux[,3,2] <- Hgrazing*h_D2              # from detritivores 
      flux[,4,2] <- Zgrazing*z_D2              # from Z
      flux[,5,2] <- (Vattgrazing+Vatt2grazing)*(1.0-vatt_beta)*(1.0-vatt_DOC)    # faecal pellets: attached prokaryote consumers 
      flux[,6,2] <- Vflgrazing*(1.0-vfl_beta)*(1.0-vfl_DOC)       # faecal pellets: free-living prokaryote consumers
      
      # DOC
      flux[,1,3] <- ex_DOCin                  # DOC from euphotic zone (direct flux)
      flux[,2,3] <- ex_act                    # DOC from eupthotic zone (active flux)
      flux[,3,3] <- D1toDOCsolub              # solubilization D1 to DOC
      flux[,4,3] <- D2toDOCsolub              # solubilization D2 to DOC
      flux[,5,3] <- Hgrazing*h_DOC            # from detritivores 
      flux[,6,3] <- Zgrazing*z_DOC            # from Z
      flux[,7,3] <- (Vattgrazing+Vatt2grazing)*vatt_DOC      # from attached prokaryote consumers 
      flux[,8,3] <- Vflgrazing*vfl_DOC        # from free living prokaryote consumers
      flux[,9,3] <- -DOCuptakeBfl             # uptake by free living prokaryotes
      
      # Attached prokaryotes
      flux[,1,4] <- D1toBatt                  # utilisation of D1 and D2 (= carbon demand)
      flux[,2,4] <- -flux[,1,4]*(1.0-w_att)    # respiration
      flux[,3,4] <- -HgrazingBatt             # loss to detritivores
      flux[,4,4] <- -Vattgrazing              # loss to attached prokaryote consumers
      
      # Free living prokaryotes
      flux[,1,5] <- DOCuptakeBfl              # consumption (= carbon demand)
      flux[,2,5] <- -flux[,1,5]*(1.0-w_fl)     # respiration
      flux[,3,5] <- -Vflgrazing               # loss to free living bacterivores
      
      # Attached prokaryote consumers
      flux[,1,6] <- Vattgrazing               # consumption of attached prokaryotes (= carbon demand)
      flux[,2,6] <- -Vattgrazing*(1.0-vatt_DOC)*vatt_beta*(1.0-vatt_npe) # respiration
      flux[,3,6] <- -Vattgrazing*vatt_DOC     # excretion of DOC
      flux[,4,6] <- -Vattgrazing*(1.0-vatt_beta)*(1.0-vatt_DOC)          # faecal pellets to D2
      flux[,5,6] <- -HgrazingVatt             # grazing by detritivores
      flux[,6,6] <- -rate*(1.0-(1.0-y_B)*zi2)*Vatt      # closure respiration
      
      # Free living prokaryote consumers
      flux[,1,7] <- Vflgrazing                # consumption of free-living prokaryotes (= carbon demand)
      flux[,2,7] <- -flux[,1,7]*(1.0-vfl_DOC)*vfl_beta*(1.0-vfl_npe)      # respiration
      flux[,3,7] <- -Vflgrazing*vfl_DOC       # excretion of DOC
      flux[,4,7] <- -Vflgrazing*(1.0-vfl_beta)*(1.0-vfl_DOC)             # faecal pellets to D2
      flux[,5,7] <- -ZgrazingVf               # grazing by Z1
      
      # Detritivores
      flux[,1,8] <- Hgrazing                  # grazing (= carbon demand)
      flux[,2,8] <- -Hgrazing*(1.0-h_DOC-h_D2)*h_beta*(1.0-h_npe)        # respiration
      flux[,3,8] <- -Hgrazing*h_DOC           # excretion of DOC
      flux[,4,8] <- -Hgrazing*(1.0-h_D2-h_DOC)*(1.0-h_beta)              # faecal pellets to D1
      flux[,5,8] <- -Hgrazing*h_D2            # losses to D2
      flux[,6,8] <- -ZgrazingH                # grazing by Z1
      
      # Zooplankton 1   (1st in chain of 6)
      flux[,1,9] <- ZgrazingVf+ZgrazingH      # consumption of free-living prokaryote consumers and detritivores (= carbon demand)
      flux[,2,9] <- -(ZgrazingVf+ZgrazingH)*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[,3,9] <- -(ZgrazingVf+ZgrazingH)*z_DOC   # excretion of DOC
      flux[,4,9] <- -(ZgrazingVf+ZgrazingH)*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,9] <- -(ZgrazingVf+ZgrazingH)*z_D2    # faecal pellets to D2
      flux[,6,9] <- -grazingZ2                      # grazing by Z2
      
      # Zooplankton 2
      flux[,1,10] <- grazingZ2    # consumption (= carbon demand)
      flux[,2,10] <- -grazingZ2*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[,3,10] <- -grazingZ2*z_DOC   # excretion DOC
      flux[,4,10] <- -grazingZ2*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,10] <- -grazingZ2*z_D2    # losses to D2
      flux[,6,10] <- -grazingZ3   # grazing by Z3
      
      # Zooplankton 3
      flux[,1,11] <- grazingZ3   # consumption (= carbon demand)
      flux[,2,11] <- -grazingZ3*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[,3,11] <- -grazingZ3*z_DOC   # excretion of DOC
      flux[,4,11] <- -grazingZ3*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,11] <- -grazingZ3*z_D2    # losses to D2
      flux[,6,11] <- -grazingZ4   # grazing by Z4
      
      # Zooplankton 4
      flux[,1,12] <- grazingZ4   # consumption (= carbon demand)
      flux[,2,12] <- -grazingZ4*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[,3,12] <- -grazingZ4*z_DOC   # excretion of DOC
      flux[,4,12] <- -grazingZ4*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,12] <- -grazingZ4*z_D2    # losses to D2
      flux[,6,12] <- -grazingZ5   # grazing by Z5
      
      # Zooplankton 5
      flux[,1,13] <- grazingZ5   # consumption (= carbon demand)
      flux[,2,13] <- -grazingZ5*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)   # respiration
      flux[,3,13] <- -grazingZ5*z_DOC   # excretion of DOC
      flux[,4,13] <- -grazingZ5*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,13] <- -grazingZ5*z_D2    # flosses to D2
      flux[,6,13] <- -grazingZ6   # grazing by Z5
      
      # Zooplankton 6   (last in chain)
      flux[,1,14] <- grazingZ6   # consumption (= carbon demand)
      flux[,2,14] <- -grazingZ6*(1.0-z_DOC-z_D2)*z_beta*(1.0-z_npe)-0.5*Z[,6]   # respiration
      flux[,3,14] <- -grazingZ6*z_DOC   # excretion of DOC
      flux[,4,14] <- -grazingZ6*(1.0-z_D2-z_DOC)*(1.0-z_beta)         # faecal pellets to D1
      flux[,5,14] <- -grazingZ6*z_D2    # flosses to D2
      
      # Attached prokaryotes: D2
      flux[,1,15] <- D2toBatt2                  # utilisation of D1 and D2 (= carbon demand)
      flux[,2,15] <- -flux[,1,15]*(1.0-w_fl) # respiration
      #flux[,2,15] <- -flux[,1,15]*(1.0-w_att)    
      flux[,3,15] <- 0.0                        # loss to detritivores
      flux[,4,15] <- -Vatt2grazing              # loss to attached prokaryote consumers
      
      # Attached prokaryote consumers: D2
      flux[,1,16] <- Vatt2grazing               # consumption of attached prokaryotes (= carbon demand)
      flux[,2,16] <- -Vatt2grazing*(1.0-vatt_DOC)*vatt_beta*(1.0-vatt_npe) # respiration
      flux[,3,16] <- -Vatt2grazing*vatt_DOC     # excretion of DOC
      flux[,4,16] <- -Vatt2grazing*(1.0-vatt_beta)*(1.0-vatt_DOC)          # faecal pellets to D2
      flux[,5,16] <- 0.0                       # grazing by detritivores
      flux[,6,16] <- -rate*Vatt2               # closure respiration
      
      flux.df = as.data.frame(flux)
      
      # update state variables
      D1 <- D1 + rowSums(flux.df[,1:10])*tstep          # detritus D1
      D2 <- D2 + rowSums(flux.df[,11:20])*tstep         # detritus D2
      DOC <- DOC + rowSums(flux.df[,21:30])*tstep       # DOC
      Batt <- Batt + rowSums(flux.df[,31:40])*tstep     # attached prokaryotes
      Bf <- Bf + rowSums(flux.df[,41:50])*tstep         # free-living prokaryotes
      Vatt <- Vatt + rowSums(flux.df[,51:60])*tstep     # attached prokaryote consumers
      Vf <- Vf + rowSums(flux.df[,61:70])*tstep         # free-living prokaryote consumers
      H <- H + rowSums(flux.df[,71:80])*tstep           # detritivores
      Z[,1] <- Z[,1] + rowSums(flux.df[,81:90])*tstep     # Z1
      Z[,2] <- Z[,2] + rowSums(flux.df[,91:100])*tstep    # Z2
      Z[,3] <- Z[,3] + rowSums(flux.df[,101:110])*tstep   # Z3
      Z[,4] <- Z[,4] + rowSums(flux.df[,111:120])*tstep   # Z4
      Z[,5] <- Z[,5] + rowSums(flux.df[,121:130])*tstep   # Z5
      Z[,6] <- Z[,6] + rowSums(flux.df[,131:140])*tstep   # Z6
      Batt2 <- Batt2 + rowSums(flux.df[,141:150])*tstep # attached prokaryotes: D2
      Vatt2 <- Vatt2 + rowSums(flux.df[,151:160])*tstep # attached prokaryote consumers: D2
      
    }  # time loop
    # respiration
    RZ <- flux[,2,8]+flux[,2,9]+flux[,2,10]+flux[,2,11]+flux[,2,12]+flux[,2,13]+flux[,2,14]  # detritivores and carnivores
    RB <- flux[,2,4]+flux[,2,5]+flux[,2,15]      # attached and free-living prokaryotes
    
  }  # time loop 1-10
  
  
  # Inputs of C to mesopelagic zone
  n <- flux[,1,1]+flux[,1,3]+flux[,2,3]                # total C input (POC plus DOC)
  
  # Respiration
  ZR <- -(flux[,2,9]+flux[,2,10]+flux[,2,11]+flux[,2,12]+flux[,2,13]+flux[,2,14])   # carnivores
  BattR <- -flux[,2,4]-flux[,2,15]
  VattR <- -flux[,2,6]-flux[,6,6]-flux[,2,16]-flux[,6,16]  #4 attached prokaryote consumers (includes closure respiration)
  
  
  # Carbon demand
  ZCD <-flux[,1,9]+flux[,1,10]+flux[,1,11]+flux[,1,12]+flux[,1,13]+flux[,1,14]      # carnivores
  BattCD <- flux[,1,4]+flux[,1,15]
  VattCD <- flux[,1,6]+flux[,1,16]
  
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
  DOC_sol <- flux[,3,3]+flux[,4,3]
  
  #closure respi
  closure = -(flux[,6,6]+flux[,6,16])
  D1_Batt <- -flux[,2,1]-flux[,3,1] 
  D2_Batt <- -flux[,1,2]-flux[,2,2] 
  
  #########################SORTIE POUR ANALYSE ABC :
  
  Respi_attached = -flux[,2,4]-flux[,2,15]
  Respi_zooplancton = (-flux[,2,7]) + (-flux[,2,6]-flux[,6,6]-flux[,2,16]-flux[,6,16]) + (-flux[,2,8]) + (-(flux[,2,9]+flux[,2,10]+flux[,2,11]+flux[,2,12]+flux[,2,13]+flux[,2,14]))
  Production_Free_living = Prod_Bfl
  Production_attached = Prod_Batt
  
  #Production_NonSinking =  Prod_Bfl + (D2toBatt2 * w_fl ) #production free living + production bacteries attachÃ©es aux particules suspendues
  #Production_Sinking = D1toBatt*w_att   #production bacterie attachÃ©es aux particules qui chutent (peu importe la vitesse)
  Production_NonSinking = (Prod_Bfl + (D2toBatt2 * w_fl)) #production free living + production bacteries attachées aux particules suspendues
  Production_Sinking = D1toBatt * w_att  #production bacterie attachées aux particules qui chutent (peu importe la vitesse)
  
  
  Respiration_Sinking  = -flux[,2,4]
  
  Respiration_zoo = Respi_zooplancton
  
  
  res <- c( Production_NonSinking, Production_Sinking, Respiration_Sinking, Respiration_zoo)
  res <- matrix(data = res, n_bootstrap, 4)
 
  if (!is.null(flux_out_nb)){
    return(res[,flux_out_nb])
  }
  else{
    return(res)
  } 
  
}