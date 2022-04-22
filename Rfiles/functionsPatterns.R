#helper functions later in use
###################################
#calculates between a given start point and and end point the points in between for a given total number of points
calc.X.values<-function(x.B, x.E, time.points){
  
  i<-0:(time.points-1)
  x<-x.B+(x.E-x.B)/(time.points-1)*i
  return(x)
}


#calculates for given start Fidelity and end Fidelity the within Fidelity 
#for a given number of timepoints
#using a logarithmus function, 
#where with the par.slope Stretching and compression of the function can be influenced
#par.slope has to be greater than 0, near zero then the values are almost constant, the greater then the values are going against the linear line
find.Fidelity.log<-function(time.points, Fid.End, Fid.T1, par.slope=1){
  
  #using Logarithmusfunktion, Basis ist e
  
  i<-0:(time.points-1)
  ##Reference Logarithmus-Funktion
  #Start-X-Wert
  ref.x1<-0.005#ändert sich in allen log. Kurven
  #End-X-Wert
  ref.x2<-1#bleibt bei allen Kurven (ist deren Schnittpunkt)
  #merke Referenz-Y-Startwert für andere Funktionen
  ref.y1<-log(ref.x1)
  
  #Streckungs/Stauchungsparameter der log Funktions
  m<-par.slope
  #Berechnung neuen Start-X-Wert, hat den selben Funktionswert wie in der Referenzfunktion
  x.B.neu<-(exp(1)^(ref.y1))^(1/m)
  #Berechnung für alle Zeitpunkte x-Werte
  x<-calc.X.values(x.B=x.B.neu, x.E=ref.x2, time.points)
  #Berechne Funktionswerte der Funktion zu x-Werten
  f<-m*log(x)
  #NormierungsFaktor für Funktionswerte zum Mappen auf [Fid_B,Fid_E]
  norm<-(Fid.End-Fid.T1)/(m*log(x[time.points])-m*log(x[1]))
  #Normierung der Funktionswerte zum Mappen in neuen Funktionsbereich
  f.norm<-Fid.End+norm*f
  res<-cbind(time=1:time.points, Fidelity.Prozent=round(f.norm*100))
  return(res)
}

#calculates for given start Fidelity and end Fidelity the within Fidelity 
#for a given number of timepoints
#using a linear function
find.Fidelity.linear<-function(time.points, Fid.End, Fid.T1){
  
  m<-(Fid.T1-Fid.End)/(1-time.points)
  n<-(Fid.End-time.points*Fid.T1)/(1-time.points)
  x<-1:time.points
  
  f<-m*x+n
  res<-cbind(time=1:time.points, Fidelity.Prozent=round((m*x+n)*100)
  )
  return(res)
}

#calculates for given start Fidelity and end Fidelity the within Fidelity 
#for a given number of timepoints
#using an exponential function, 
#where with the par.slope Stretching and compression of the function can be influenced
#wherby here the exponential function is derived as the reflection of the logarithm function on the straight line
#par.slope has to be greater than 0, near zero then the values are almost constant, the greater then the values are going against the linear line
find.Fidelity.exp<-function(time.points, Fid.End, Fid.T1, par.slope=1){
  
  #Berechnung Logarithmusfunktion
  res.log<-find.Fidelity.log(time.points, Fid.End, Fid.T1, par.slope=par.slope)
  #Berechnung Spiegelgeraden
  res.linear<-find.Fidelity.linear(time.points, Fid.End, Fid.T1)
  #Differenzberechnung
  diff.exp<-rev(res.linear[1,"Fidelity.Prozent"]+
                  (res.linear[time.points,"Fidelity.Prozent"]-res.log[,"Fidelity.Prozent"]))
  res<-cbind(time=1:time.points, 
             Fidelity.Prozent= diff.exp)
  return(res)
}
