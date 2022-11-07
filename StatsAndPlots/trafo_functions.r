require("ggplot2")
require("dplyr")
require("tidyr")
require("reshape2")
require("Cairo")
require("patchwork")
require("egg")
library("ggpubr")


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
quicksave<-function(name,plottosave){
  CairoPNG()
  ggsave(name,plottosave,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE,type="cairo-png")
  dev.off()
  
}

## plotting functions

quickplot_rd <- function(plotdata_or,type="simple"){
  
  plotdata<-plotdata_or
  dfbond1=melt(select(plotdata,contains("bond_1")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
  if (type!="simple"){
    dfbond2=melt(select(plotdata,contains("bond_2")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
    dfbond3=melt(select(plotdata,contains("bond_3")|contains("time")), id.vars=c("time"), variable.name="name",value.name = "reldist")
  }
  plot<-ggplot(plotdata,aes(x=time))+
    xlab("Simulation time [ps]")+
    ylab("Relative Distance [臸")+
    #geom_smooth(size=1,aes(y=intst1.bond_1))+
    theme_classic()+
    theme(text=element_text(size=30))+
    geom_smooth(data=dfbond1,aes(x=time,y=reldist),span=0.3,color=uv_blue,size=1.5)
    #geom_point(data=dfbond1,aes(x=time,y=reldist))
  
  if (type!="simple"){
    plot<-plot+
      geom_smooth(data=dfbond2,aes(x=time,y=reldist),span=0.3,color=uv_wine,size=1.5)+
      geom_smooth(data=dfbond3,aes(x=time,y=reldist),span=0.3,color=uv_lightgreen,size=1.5)
  }
    returnValue(plot)
}

quickplot_bp<-function(plotdata){
  
  plot=ggplot(plotdata,aes(x=kvalue,y=kcalvalue))+
    geom_boxplot(aes(fill=kvalue))+
    stat_boxplot(geom="errorbar",width=0.15)+
    facet_grid(rows=vars(restraint))+
    xlab("Value of force constant k in energy expression")+
    ylab("RBFE [kcal]")+
    theme_classic()+
    theme(legend.position="none",text=element_text(size=30))
  
  
  
  
  
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

ljplot=ggplot(ljvals)+geom_line(aes(xvals,yvals,color="Lennard-Jones-potential with 位=0.5"))+
  geom_line(aes(xvals,yvalssoftcore,color="Soft-core-potential with 位=0.70"))+
  geom_line(aes(xvals,yvalssoftcore2,color="Soft-core-potential with 位=0.75"))+
  coord_cartesian(xlim=c(0,2),ylim=c(-0.5,0.5))+xlab("Distance r")+ylab("Potential Energy")+
  geom_line(aes(xvals,yvalssmalllambda,color="Lennard-Jones-potential with 位=0.9"))+
  theme_light() + labs(color="Type of Potential")+
  scale_color_manual(name="Potential:",values=c("Soft-core-potential with 位=0.75"=uv_orangered,
                                                "Soft-core-potential with 位=0.70"=uv_blue,
                                                "Lennard-Jones-potential with 位=0.5"="black",
                                                "Lennard-Jones-potential with 位=0.9"=uv_wine))

## data read-in: RBFE

ma_taapdb=read.csv("../ma_taapdb.csv")
ma_tablit24to25=read.csv("../ma_tablit24to25.csv")
ma_tablit24to26=read.csv("../ma_tablit24to26.csv")
ma_zn222=read.csv("../ma_zn222.csv")
ma_zn223a=read.csv("../ma_zn223a.csv")
ma_zn223b=read.csv("../ma_zn223b.csv")


## data read-in: binding site relative distances

bs_taapdb=read.csv("../bindingsitedynamics/taBlit_run2.csv")
bs_taapdbn=read.csv("../bindingsitedynamics/tablit new class test_run3.csv")

## data read-in: restraints relative distances 

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

## create plots: bindings site relative distances

bsplot1=ggplot(bs_taapdb)+theme_classic()+theme(text=element_text(size=30))+geom_smooth(aes(frames,d_bond1),color=uv_blue)+geom_smooth(aes(frames,d_bond2),color=uv_wine)+geom_smooth(aes(frames,d_bond3),color=uv_lightgreen)+ylab("Absolute distance [臸")+xlab("Simulation Time [ps]")

bsplot2=ggplot(bs_taapdbn)+theme_classic()+theme(text=element_text(size=30),axis.ticks.y=element_blank(),axis.text.y=element_blank())+geom_smooth(aes(frames,d_1),color=uv_blue)+geom_smooth(aes(frames,d_2),color=uv_wine)+geom_smooth(aes(frames,d_3),color=uv_lightgreen)+rremove("ylab")

bsplotcombi=bsplot1+bsplot2+ylim(1,7)+xlab("Simulation Time [ps]")

quicksave("bsplotcombi.png",bsplotcombi)
## create plots: relative distances zn 222

pdata=form_RelativeDistance_data(rd_zn222_k3simple)
plot_rd1=quickplot_rd(pdata)+theme(plot.margin=margin(0,1,0,0, 'cm'))

quicksave("reldistplot_zn222_k3simple.png",plot_rd1)

pdata=form_RelativeDistance_data(rd_zn222_k100simple)
plot_rd2=quickplot_rd(pdata)

quicksave("reldistplot_zn222_k100simple.png",plot_rd2)

plot_rd2<-plot_rd2+rremove("ylab")+
  
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )

combiplot1=plot_rd1+ ylim(0,3)+plot_rd2+theme(plot.margin=unit(c(0,1,0,0), 'cm'))
ggsave("rdcombi_zn222_k3_k100_simple.png",combiplot1,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE,type="cairo-png")

pdatard3=form_RelativeDistance_data(rd_zn222_k0ex3)
plot_rd3<-quickplot_rd(pdatard3,"ex3")+theme(plot.margin=margin(0,1,0,0, 'cm'))



quicksave("reldistplot_zn222_k0ex3.png",plot_rd3)




pdatard5=form_RelativeDistance_data(rd_zn222_k100ex3)
plot_rd5<-quickplot_rd(pdatard5,"ex3")



quicksave("reldistplot_zn222_k100ex3.png",plot_rd5)




#pdatard5=form_RelativeDistance_data(rd_zn222_k100ex3)
#plot_rd5<-quickplot_rd(pdatard5,"ex3")

#plot_rd3<-plot_rd3+ylim(0,7)
plot_rd5<-plot_rd5+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
combiplot2=plot_rd3+plot_rd5+ylim(0,9)+theme(plot.margin=unit(c(0,1,0,0), 'cm'))
                       


ggsave("rdcombi_zn222_k0_k100_ex3.png",combiplot2,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE,type="cairo-png")


pdata=form_RelativeDistance_data(rd_zn222_k3ex3)
plot_rd4=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_zn222_k3ex3.png",plot_rd4)
## create plots: relative distances taapdb 24to25

pdata=form_RelativeDistance_data(rd_taapdb_k10ex3)
plot_rd6=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_taapdb2425_k10ex3.png",plot_rd6)

pdata=form_RelativeDistance_data(rd_taapdb_k400ex3)
plot_rd7=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_taapdb2425_k400ex3.png",plot_rd7)

pdata=form_RelativeDistance_data(rd_taapdb_k0ex3)
plot_rd8=quickplot_rd(pdata,"ex3")+theme(plot.margin=margin(0,1,0,0, 'cm'))

quicksave("reldistplot_taapdb2425_k0ex3.png",plot_rd8)

plot_rd7<-plot_rd7+rremove("ylab")+
  ylim(0,1)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
combiplot4<-plot_rd8+plot_rd7+theme(plot.margin=unit(c(0,1,0,0), 'cm'))
ggsave("rdcombi_taapdb2425_k0_k400_ex3.png",combiplot4,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE,type="cairo-png")

## create plot: relative distances tablit2425

pdata=form_RelativeDistance_data(rd_tablit2425_k3ex3)
plot_rd9=quickplot_rd(pdata,"ex3")+theme(plot.margin=margin(0,1,0,0, 'cm'))

quicksave("reldistplot_tablit2425_k3ex3.png",plot_rd9)

pdata=form_RelativeDistance_data(rd_tablit2425_k100ex3)
plot_rd10=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_tablit2425_k100ex3.png",plot_rd10)

pdata=form_RelativeDistance_data(rd_tablit2425_flatbottom)
plot_rd11=quickplot_rd(pdata,"ex3")

quicksave("reldistplot_tablit2425_k100ex3flatbottom.png",plot_rd11)

plot_rd10<-plot_rd10+rremove("ylab")+
  
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
combiplot3<-plot_rd9+plot_rd10+ylim(0,1.5)+theme(plot.margin=unit(c(0,1,0,0), 'cm'))
                     
ggsave("rdcombi_tablit2425_k3_k100_ex3.png",combiplot3,device=png,width=PWIDTH,height=PHEIGHT,dpi=PDPI,limitsize=FALSE,type="cairo-png")


## create plot: boxplot taapdb

pdata=form_boxplot_data(ma_taapdb)



plot_taapdb=quickplot_bp(filter(pdata,timescale=="1.25ns"))

quicksave("boxplot_rbfe_taapdb_1.25ns.png",plot_taapdb)

## create plot: boxplot tablit24to25

pdata=form_boxplot_data(ma_tablit24to25)



plot_taapdb=quickplot_bp(filter(pdata,timescale=="1.25ns" | kvalue==0))

quicksave("boxplot_rbfe_tablit24to25_1.25ns.png",plot_taapdb)

plot_taapdb=quickplot_bp(filter(pdata,timescale=="5ns" | kvalue==0))

quicksave("boxplot_rbfe_tablit24to25_5ns.png",plot_taapdb)
## create plot: boxplot tablit24to26

pdata=form_boxplot_data(ma_tablit24to26)



plot_taapdb=quickplot_bp(pdata)

quicksave("boxplot_rbfe_tablit24to26_1.25ns.png",plot_taapdb)

## create plot: boxplot zn222

pdata=form_boxplot_data(ma_zn222)



plot_taapdb=quickplot_bp(pdata)
quicksave("boxplot_rbfe_zn222_5ns.png",plot_taapdb)

## create plot: boxplot zn223a

pdata=form_boxplot_data(ma_zn223a)



plot_taapdb=quickplot_bp(pdata)

quicksave("boxplot_rbfe_zn223a_5ns.png",plot_taapdb)

## create plot: boxplot zn223a

pdata=form_boxplot_data(ma_zn223b)



plot_taapdb=quickplot_bp(pdata)

quicksave("boxplot_rbfe_zn223b_5ns.png",plot_taapdb)