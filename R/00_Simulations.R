SPRobj <- function(h=0.75,Dobj=0.4){
  Dobj*((1-h)+Dobj*(5*h-1))/(4*h*Dobj)
}

stepness <- function(Dobj=0.55,SPRobj=0.6){
  (1-Dobj)/(1-5*Dobj+4*SPRobj)
}

stepness(Dobj=0.4,SPRobj=0.45)
SPRobj(h=0.7,Dobj=0.4)
