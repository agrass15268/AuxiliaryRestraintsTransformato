#Specify Colors as required by the University of Vienna

uv_blue="#0063A6"
uv_grey="#666666"
uv_wine="#A71C49"
uv_orangered="#DD4814"
uv_lightgreen="#94C154"
uv_mintgreen="#11897A"


calculate_kjoule_kcal<-function(dataframe){
  dataframe["kjvalue"]<-dataframe["ktvalue"]*303.15/1000
  dataframe["kcalvalue"]<-dataframe["kjvalue"]*0.239006
  returnValue(dataframe)
}

LJ_potential<-function(xvalue,lambda){
  #xvalue=distance between atoms r
  y=(1-lambda)*((1/xvalue^12)-(1/xvalue^6))
}

softcore_potential<-function(xvalue,lambda){
  delta=1
  y=(1-lambda)*((1/((xvalue^2+lambda*delta)^6))-(1/((xvalue^6+lambda*delta)^3)))
}

xvals=seq(-0.5,5,0.02)
yvals=LJ_potential(xvals,0.5)
yvalssmalllambda=LJ_potential(xvals,0.9)
yvalssoftcore=softcore_potential(xvals,0.7)
yvalssoftcore2=softcore_potential(xvals,0.75)
ljvals=data.frame(xvals,yvals,yvalssoftcore)

ljplot=ggplot(ljvals)+geom_line(aes(xvals,yvals,color="Lennard-Jones-potential with λ=0.5"))+
  geom_line(aes(xvals,yvalssoftcore,color="Soft-core-potential with λ=0.70"))+
  geom_line(aes(xvals,yvalssoftcore2,color="Soft-core-potential with λ=0.75"))+
  coord_cartesian(xlim=c(0,2),ylim=c(-0.5,0.5))+xlab("Distance r")+ylab("Potential Energy")+
  geom_line(aes(xvals,yvalssmalllambda,color="Lennard-Jones-potential with λ=0.9"))+
  theme_light() + labs(color="Type of Potential")+
  scale_color_manual(name="Potential:",values=c("Soft-core-potential with λ=0.75"=uv_orangered,
                                                "Soft-core-potential with λ=0.70"=uv_blue,
                                                "Lennard-Jones-potential with λ=0.5"="black",
                                                "Lennard-Jones-potential with λ=0.9"=uv_wine))
