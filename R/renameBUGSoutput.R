
# rename BUGSoutput object for further use
renameBUGSoutput <- function(BUGSoutput,
                             thetaUnique,
                             model=c("traitMPT", "betaMPT")){

  rownames(BUGSoutput$median$theta) <-
    rownames(BUGSoutput$mean$theta) <-
    rownames(BUGSoutput$sd$theta) <-
    names(BUGSoutput$mean$mean) <-
    names(BUGSoutput$sd$mean) <-
    names(BUGSoutput$median$mean) <-
    thetaUnique

  if(model == "traitMPT"){
    rownames(BUGSoutput$median$rho) <-
      rownames(BUGSoutput$mean$rho) <-
      rownames(BUGSoutput$sd$rho) <-
      colnames(BUGSoutput$median$rho) <-
      colnames(BUGSoutput$mean$rho) <-
      colnames(BUGSoutput$sd$rho) <-
      names(BUGSoutput$mean$mu) <-
      names(BUGSoutput$sd$mu) <-
      names(BUGSoutput$median$mu) <-
      names(BUGSoutput$mean$sigma) <-
      names(BUGSoutput$sd$sigma) <-
      names(BUGSoutput$median$sigma) <-
      thetaUnique
  }else if (model == "betaMPT"){
    names(BUGSoutput$mean$alph) <-
      names(BUGSoutput$mean$bet) <-
      names(BUGSoutput$sd$alph) <-
      names(BUGSoutput$sd$bet) <-
      names(BUGSoutput$median$alph) <-
      names(BUGSoutput$median$bet) <-
      names(BUGSoutput$mean$sd) <-
      names(BUGSoutput$sd$sd) <-
      names(BUGSoutput$median$sd) <-
      thetaUnique
  }
  BUGSoutput
}
