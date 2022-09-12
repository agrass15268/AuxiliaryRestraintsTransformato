require("ggplot2")
require("dplyr")
require("tidyr")
require("reshape2")

### trafo-functions.r
# Script to transform data and generate the plots used in the thesis
# to work, the working directory must be one folder below the data source (as it is in the repository)
# Plots will be generated in the current working directory


## constant setup

# set plot generation settings

PWIDTH=1500
PHEIGHT=600
PDPI=800

# specify Colors as required by the University of Vienna

uv_blue="#0063A6"
uv_grey="#666666"
uv_wine="#A71C49"
uv_orangered="#DD4814"
uv_lightgreen="#94C154"
uv_mintgreen="#11897A"

## data functions
calculate_kjoule_kcal<-function(dataframe){
  dataframe["kjvalue"]<-dataframe["ktvalue"]*303.15/1000
  dataframe["kcalvalue"]<-dataframe["kjvalue"]*0.239006
  returnValue(dataframe)
}

form_boxplot_data<-function(plotdata){
  
  # Replace kvalue - Strings with sortable ints
  plotdata["kvalue"][plotdata["kvalue"]=="k10000"]<-10000
  plotdata["kvalue"][plotdata["kvalue"]=="k400"]<-400
  plotdata["kvalue"][plotdata["kvalue"]=="k100"]<-100
  plotdata["kvalue"][plotdata["kvalue"]=="k10"]<-10
  plotdata["kvalue"][plotdata["kvalue"]=="k3"]<-3
  plotdata["kvalue"][plotdata["kvalue"]==""]<-0
  
  plotdata$kvalue=as.numeric(plotdata$kvalue)
  plotdata$kvalue=as.factor(plotdata$kvalue)
  # Reassign norestraint data to enable visual comparision
  
  nrd=filter(plotdata,restraint == "norestraint")
  plotdata["restraint"][plotdata["restraint"]=="norestraint"]<-"simple"
  nrd["restraint"]<-"ex3"
  plotdata=bind_rows(plotdata,nrd)
  
  plotdata=arrange(plotdata,kvalue)
  
  
  returnValue(plotdata)
}

form_RelativeDistance_data<-function(plotdata){
  
  
  
  
  returnValue(plotdata)
}
quicksave<-function(name){
  ggsave(name,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE)
}

## plotting functions

quickplot_rd <- function(plotdata,type="simple"){
  
  
  dfbond1=melt(select(plotdata,contains("bond_1")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
  if (type!="simple"){
    dfbond2=melt(select(plotdata,contains("bond_2")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
    dfbond3=melt(select(plotdata,contains("bond_3")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
  }
  plot<-ggplot(plotdata,aes(x=time))+
    xlab("Simulation time [ps]")+
    ylab("Relative Distance to starting configuration [nm]")+
    #geom_line(size=1,aes(y=intst1.bond_1))+
    theme_light()+
    geom_smooth(data=dfbond1,aes(x=time,y=reldist),span=0.3,color=uv_blue)
    #geom_point(data=dfbond1,aes(x=time,y=reldist))
  
  if (type!="simple"){
    plot<-plot+
      geom_smooth(data=dfbond2,aes(x=time,y=reldist),span=0.3,color=uv_wine)+
      geom_smooth(data=dfbond3,aes(x=time,y=reldist),span=0.3,,color=uv_lightgreen)
  }
    returnValue(plot)
}
## Lennard - Jones visualisation plot
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

## data read-in: RBFE

ma_taapdb=read.csv("../ma_taapdb.csv")
ma_tablit24to25=read.csv("../ma_tablit24to25.csv")
ma_tablit24to26=read.csv("../ma_tablit24to26.csv")
ma_zn222=read.csv("../ma_zn222.csv")
ma_zn223a=read.csv("../ma_zn223a.csv")
ma_zn223b=read.csv("../ma_zn223b.csv")

## data read-in: relative distances

rd_tablit2425_k3ex3=read.csv("../relative_distances/tablit_24to25/5ns-k3ex3-1.csv")
rd_tablit2425_k100ex3=read.csv("../relative_distances/tablit_24to25/5ns-k100ex3-1.csv")
rd_tablit2425_flatbottom=read.csv("../relative_distances/tablit_24to25/5ns-k10000ex3flatbottom-1.csv")

rd_taapdb_k400ex3=read.csv("../relative_distances/taapdb24to25/1.25ns-k400ex3-noscaling-2.csv")
rd_taapdb_k10ex3=read.csv("../relative_distances/taapdb24to25/1.25ns-k10ex3-noscaling-2.csv")
rd_taapdb_k0ex3=read.csv("../relative_distances/taapdb24to25/1.25ns-norestraints-1.csv")

rd_zn222_k3simple=read.csv("../relative_distances/zn222tozn148/5ns-k3simple-1.csv")
rd_zn222_k100simple=read.csv("../relative_distances/zn222tozn148/5ns-k100simple-1.csv")

rd_zn222_k0ex3=read.csv("../relative_distances/zn222tozn148/5ns-norestraints-1.csv")
rd_zn222_k3ex3=read.csv("../relative_distances/zn222tozn148/5ns-k3ex3-1.csv")
rd_zn222_k100ex3=read.csv("../relative_distances/zn222tozn148/5ns-k100ex3-1.csv")

## create plots: relative distances zn 222

pdata=form_RelativeDistance_data(rd_zn222_k3simple)
plot_rd=quickplot_rd(pdata)

quicksave("reldistplot_zn222_k3simple.png")

pdata=form_RelativeDistance_data(rd_zn222_k100simple)
plot_rd=quickplot_rd(pdata)

quicksave("reldistplot_zn222_k100simple.png")

pdata=form_RelativeDistance_data(rd_zn222_k0ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_zn222_k0ex3.png")

pdata=form_RelativeDistance_data(rd_zn222_k3ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_zn222_k3ex3.png")

pdata=form_RelativeDistance_data(rd_zn222_k100ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_zn222_k100ex3.png")

## create plots: relative distances taapdb 24to25

pdata=form_RelativeDistance_data(rd_taapdb_k10ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_taapdb2425_k10ex3.png")

pdata=form_RelativeDistance_data(rd_taapdb_k400ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_taapdb2425_k400ex3.png")

pdata=form_RelativeDistance_data(rd_taapdb_k0ex3)
plot_rd=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_taapdb2425_k0ex3.png")

## create plot: relative distances tablit2425

pdata=form_RelativeDistance_data(rd_tablit2425_k3ex3)
plot_rd_k3ex3=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_tablit2425_k3ex3.png")

pdata=form_RelativeDistance_data(rd_tablit2425_k100ex3)
plot_rd_k100ex3=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_tablit2425_k100ex3.png")

pdata=form_RelativeDistance_data(rd_tablit2425_flatbottom)
plot_rd_flatbottom=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_tablit2425_k100ex3flatbottom.png")
## create plot: boxplot taapdb

pdata=form_boxplot_data(ma_taapdb)



plot_taapdb=ggplot(filter(pdata,timescale=="1.25ns"),aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_taapdb_1.25ns.png")

## create plot: boxplot tablit24to25

pdata=form_boxplot_data(ma_tablit24to25)



plot_taapdb=ggplot(filter(pdata,timescale=="1.25ns" | kvalue==0),aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_tablit24to25_1.25ns.png")

plot_taapdb=ggplot(filter(pdata,timescale=="5ns" | kvalue==0),aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_tablit24to25_5ns.png")
## create plot: boxplot tablit24to26

pdata=form_boxplot_data(ma_tablit24to26)



plot_taapdb=ggplot(pdata,aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_tablit24to26_1.25ns.png")

## create plot: boxplot zn222

pdata=form_boxplot_data(ma_zn222)



plot_taapdb=ggplot(pdata,aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_zn222_5ns.png")

## create plot: boxplot zn223a

pdata=form_boxplot_data(ma_zn223a)



plot_taapdb=ggplot(pdata,aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_zn223a_5ns.png")

## create plot: boxplot zn223a

pdata=form_boxplot_data(ma_zn223b)



plot_taapdb=ggplot(pdata,aes(x=kvalue,y=kcalvalue))+
  geom_boxplot(aes(fill=kvalue))+
  stat_boxplot(geom="errorbar",width=0.15)+
  facet_grid(rows=vars(restraint))+
  xlab("Value of force constant k in energy expression")+
  ylab("RBFE [kcal]")+
  theme_light()+
  theme(legend.position="none")

quicksave("boxplot_rbfe_zn223b_5ns.png")