Volumes <- function(totalVol, rCount, pCount, kRratio, PtoR=2.5, kR=74.3,      
                    measVol=0.1) {
  # rCount:   n rotifers ml^-1 in stock
  # pCount:   n protists ml^-1 in stock
  # kRratio:  targeted fraction of rotifer carrying capacity (k)
  # PtoR:     ratio of protists to rotifers (default 2.5)
  # kR:       carrying capacity of average rotifer clone (constant)
  
  Rtarget <- kR * kRratio
  Ptarget <- Rtarget * PtoR
  
  rotifer_volume <- ((Rtarget * measVol) / rCount) * (totalVol / measVol)
  protist_volume <- ((Ptarget * measVol) / pCount) * (totalVol / measVol)
  
  output <- data.frame(rotifer_volume, protist_volume)
  return(output)
}
