#' funzione spettro di potenza
P <- function(K, lambda=0.15) {
  cost <- 1/( (6*pi)/(4*lambda^2) )
  return(sqrt(cost * (exp(-K*lambda)*(1-exp(-K*lambda)))))
}


#' noise sulla CT
create.noise.ct <- function(ct, lambda=0.15, sigma=0.025, sigmaHU=NULL, WGN=FALSE, return.noise=FALSE) {
  
  cat('adding noise to CT...\n')
  
  new.ct <- ct
  
  Nx <- ct[['Nx']]
  Ny <- ct[['Ny']]
  Nz <- ct[['Nz']]
  x <- ct[['x']]
  y <- ct[['y']]
  z <- ct[['z']]
  
  # se lambda è 0, usa direttamente un white gaussian noise
  if(lambda==0) {WGN=TRUE}

  # Calcola coordinate  
  dx <- mean(diff(ct[['x']]))
  dy <- mean(diff(ct[['y']]))
  k.x <- 1/dx
  k.y <- 1/dy
  
  # cordinate numeri d'onda
  u <- seq(0, (Nx-1)*k.x, by=k.x)
  u <- u - mean(u)
  v <- seq(0, (Ny-1)*k.y, by=k.y)
  v <- v - mean(v)
  
  # mesh
  U <- matrix(data=u, nrow=Nx, ncol=Ny)
  V <- matrix(data=v, nrow=Nx, ncol=Ny, byrow=TRUE)
  
  # calcola spettro di potenza
  K <- sqrt(U^2 + V^2)
  Pm <- P(K, lambda=lambda)
  
  if(is.null(sigmaHU)) {
    # calcola sigma white gaussian noise in HU complessiva (da controllare se ha senso)
    S <- abs(mean(ct$values) + min(ct$values)) * sigma
  } else {
    S <- sigmaHU
  }
  cat('Sigma:', S, '\n')
  
  # ciclo sulle slices
  for(i in 1:Nz) {
  
    # estrai slice
    my.slice <- ct[['values']][,,i]
    
    # calcola white gaussian noise
    #S <- abs(mean(my.slice) + min(my.slice)) * sigma
    Ng <- array(rnorm(Nx*Ny, mean=0, sd=S), dim=c(Nx,Ny))
    
    # trasformata FFT del noise gaussiano
    Ng.f <- fft(Ng)
    
    # shift FFT
    Nx2 <- floor(Nx/2)
    Ny2 <- floor(Ny/2)
    Ng.f.s1 <- Ng.f.s <- Ng.f * 0
    # swap x
    Ng.f.s1[ 1:(Nx-Nx2) , ] <- Ng.f[ (Nx2+1):Nx , ]
    Ng.f.s1[ (Nx2+1):Nx , ] <- Ng.f[ 1:(Nx-Nx2) , ]
    # swap y
    Ng.f.s[ , 1:(Ny-Ny2) ] <- Ng.f.s1[ , (Ny2+1):Ny ]
    Ng.f.s[ , (Ny2+1):Ny ] <- Ng.f.s1[ , 1:(Ny-Ny2) ]
    
    # calcola noise sullo spettro di potenza
    Np.f.s <- Ng.f.s * Pm
    
    # shift inverso
    Np.f.s1 <- Np.f <- Np.f.s * 0
    # swap x
    Np.f.s1[ (Nx2+1):Nx , ] <- Np.f.s[ 1:(Nx-Nx2) , ]
    Np.f.s1[ 1:(Nx-Nx2) , ] <- Np.f.s[ (Nx2+1):Nx , ]
    # swap y
    Np.f[ , (Ny2+1):Ny ] <- Np.f.s1[ , 1:(Ny-Ny2) ]
    Np.f[ , 1:(Ny-Ny2) ] <- Np.f.s1[ , (Ny2+1):Ny ]
    
    # trasformata inversa
    Np <- Re(fft(Np.f, inverse=TRUE))/(Nx*Ny)
    
    # normalizza noise
    #sd.p <- apply(Np, 2, sd)
    #sd.g <- apply(Ng, 2, sd)
    #sd.p <- sd(Np)
    #sd.g <- sd(Ng)
    #Np <- Np * (sd.g/sd.p)
      
    # CT con rumore
    if(WGN){
      my.slice.n <- my.slice + Ng
    } else {
      my.slice.n <- my.slice + Np
    }
    new.ct[['values']][,,i] <- my.slice.n
    
  }
  
  # normalizza fluttuazioni noise complessivo
  noise.ct <- new.ct
  noise.ct$values <- new.ct$values - ct$values
  noise.ct$values <- noise.ct$values/sd(noise.ct$values)*S
  new.ct$values <- ct$values + noise.ct$values
  
  
  if(return.noise) {
    return(list(ct=new.ct, noise=noise.ct))
  } else {
    return(new.ct)
  }
}
