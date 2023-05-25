library(dplR)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(tidyr)
library(treeclim) # climatic signal
library(ncdf4) # nc loading for climadata
library(raster)
library(rgdal)
library(multcompView)

# set working directory
setwd("C:/Users/jirka/Desktop/Juniper/")

serie <- read.rwl("Data/RWL/KGJ.rwl")
serie$GJ1<- NULL ### WEIRDLY THICK!

### series ploting
spag.plot(serie)

GG_serie<-serie 
GG_serie$YEAR<- rownames(GG_serie)
GG_serie<-gather(GG_serie, "LOC", "VAL", 1:41)
GG_serie<- na.omit(GG_serie)

ggplot(GG_serie, aes(x=as.numeric(YEAR), y=VAL, color=LOC)) + geom_line()

### detrend a chronology building

Det<- detrend(serie, method = "Spline", nyrs = 50)

Chrono<- chron(Det)
crn.plot(Chrono)

### Statistics
Stat_RWL<- rwl.stats(serie)
Stat_RWI<-rwi.stats.running(Det, running.window = TRUE, window.length = 30, window.overlap = 15)

## age
mean(Stat_RWL$year)
sd(Stat_RWL$year)

## TRW
mean(Stat_RWL$mean)
sd(Stat_RWL$mean)

## EPS and Rbar
rwi.stats(Det)

## diameter
diameter<- as.data.frame(t(serie))
diameter$diameter<- rowSums(diameter, na.rm = T)

mean(diameter$diameter)
sd(diameter$diameter)

## cor with Chron

for (i in c(1:ncol(serie))) {
  
  diameter[i, "COR"]<- cor(Chrono$std, serie[,i], use = "complete.obs")
  
}

mean(diameter$COR)
sd(diameter$COR)

###############################
#         Fig. 2
###############################

ylim.prim<- c(0.4, 1.8)
ylim.sec<- c(0, 42)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

F2A<- ggplot(Chrono)+
  geom_area(aes(x=as.numeric(rownames(Chrono)),y=a+samp.depth*b), fill="#acb6c2")+
  geom_line(aes(x=as.numeric(rownames(Chrono)), y=std))+
  scale_y_continuous("Ring-width index", sec.axis = sec_axis(~ (. - a)/b, name = "Sample depth"))+
  labs(x= "Time (years)", y="")+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom")+ ## umístění legendy dolu
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ylim.prim1<- c(-0.2, 0.35)
ylim.sec1<- c(-0.5, 0.9)

b1 <- diff(ylim.prim1)/diff(ylim.sec1)
a1 <- b1*(ylim.prim1[1] - ylim.sec1[1])

F2B<-  ggplot(Stat_RWI)+
  geom_line(aes(x=mid.year,y=a1+eps*b1), color="#156cd4")+
  geom_line(aes(x=mid.year, y=rbar.tot), color="#e00d14")+
  scale_y_continuous("Rbar", sec.axis = sec_axis(~ (. - a1)/b1, name = "EPS"))+
  labs(x= "Middle year of running window", y="")+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom")+ ## umístění legendy dolu
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

plot_grid(F2A, F2B, nrow=2, ncol=1, labels = c("A", "B", "C"), label_size = 12, label_x = 0.08)

ggsave("Graphs/Juniper_Troll/Fig. 2 (Chrono_info).tiff", height = 120, width = 120, units = "mm", dpi = 300)

###############################################################################
###                             climatic signal
###############################################################################

temp <- nc_open("Data/Clima_data/CRU_temp.nc")
prec <- nc_open("Data/Clima_data/CRU_prec.nc")

STORAGE_T <- data.frame(TIMESTAMP = seq(from = as.Date("19010101", tryFormats = c("%Y%m%d")), to = as.Date("20211201", tryFormats = c("%Y%m%d")), by = "month"))
STORAGE_P <- data.frame(TIMESTAMP = seq(from = as.Date("19010101", tryFormats = c("%Y%m%d")), to = as.Date("20211201", tryFormats = c("%Y%m%d")), by = "month"))

## zavedení souřadnic dané lokality
site1x <- 30.848 ; site1y <- 69.713

### Urceni pixelu podle souradnice zajmoveho bodu
Order.X <- data.frame(ORDER = c(1:temp$dim$lon$len), GRID = ncvar_get(nc = temp, varid = "lon"))
Order.Y <- data.frame(ORDER = c(1:temp$dim$lat$len), GRID = ncvar_get(nc = temp, varid = "lat"))
Order.X$DIFFERENCE <- abs(Order.X$GRID - site1x)
Order.Y$DIFFERENCE <- abs(Order.Y$GRID - site1y)
site1x.order <- Order.X[Order.X$DIFFERENCE == min(Order.X$DIFFERENCE),"ORDER"] # Poradi pixelu, pro ktery je rozdil pozice stredu pixelu od zadaneho bodu nejmensi
site1y.order <- Order.Y[Order.Y$DIFFERENCE == min(Order.Y$DIFFERENCE),"ORDER"]

### Extrakce dat z pixelu
STORAGE_T[,"Temp"] <- ncvar_get(nc = temp, varid = "tmp", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 
STORAGE_P[,"Prec"] <- ncvar_get(nc = prec, varid = "pre", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 

STORAGE_T$YEAR<- lapply(strsplit(as.character(STORAGE_T$TIMESTAMP), "\\-"), "[",1)
STORAGE_T$MONTH<- lapply(strsplit(as.character(STORAGE_T$TIMESTAMP), "\\-"), "[",2)
STORAGE_T$TIMESTAMP<- NULL
Temp<- pivot_wider(STORAGE_T, names_from = MONTH, values_from = Temp)
Temp <- as.data.frame(lapply(Temp, unlist))

STORAGE_P$YEAR<- lapply(strsplit(as.character(STORAGE_P$TIMESTAMP), "\\-"), "[",1)
STORAGE_P$MONTH<- lapply(strsplit(as.character(STORAGE_P$TIMESTAMP), "\\-"), "[",2)
STORAGE_P$TIMESTAMP<- NULL
Prec<- pivot_wider(STORAGE_P, names_from = MONTH, values_from = Prec)
Prec <- as.data.frame(lapply(Prec, unlist))

clima <- list(temperature=Temp, precipitation=Prec) # Making a list from separate dataframes

MCor<- dcc(Chrono, clima, selection = -6:9, method = "correlation", dynamic = "moving", win_size = 35, win_offset = 1)

SIG<- MCor$coef$significant
SIG$ROWS<- rownames(SIG)
SIG<- gather(SIG, "WINDOWS", "SIG", 1:ncol(SIG)-1)

MCor<- MCor$coef$coef
MCor$ROWS<- rownames(MCor)
MCor<- gather(MCor, "WINDOWS", "COR", 1:ncol(MCor)-1)
MCor$VAR<- lapply(strsplit(as.character(MCor$ROWS), "\\."), "[",1)
MCor$YEAR<- lapply(strsplit(as.character(MCor$ROWS), "\\."), "[",2)
MCor$MONTH<- lapply(strsplit(as.character(MCor$ROWS), "\\."), "[",3)
MCor$CODE<- paste(MCor$YEAR, MCor$MONTH, sep = " ")
MCor <- as.data.frame(lapply(MCor, unlist))
MCor$MID_W<- rep(c(1919:1996),each=32)
MCor$SIG<- SIG$SIG

order<- c("curr sep", "curr aug", "curr jul", "curr jun", "curr may", "curr apr", "curr mar", "curr feb", "curr jan", "prev dec", "prev nov", "prev oct", "prev sep", "prev aug", "prev jul", "prev jun")

ggplot(MCor, aes(x=MID_W, y=factor(CODE,level=order), fill=COR, alpha=SIG=="TRUE"))+geom_tile()+
  facet_wrap2(~VAR, axes = "all")+
  guides(alpha = "none")+
  scale_fill_gradient2(low = "#fc0303", mid= "#f2e355", high = "#3bd923",  midpoint = 0, limit = c(-1,1), name="Correlation")+
  labs(x= "Midle year of running window", y="")+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "right", legend.margin=margin(t=-10))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave("Graphs/Juniper_Troll/Fig. 3 (Clima_cor).tiff", height = 100, width = 250, units = "mm", dpi = 300)


##############################################################################################################
##                        Juniper Troll
##############################################################################################################

Juniper<- as.data.frame(KGJ$GJ10F); rownames(Juniper)<- rownames(KGJ); colnames(Juniper)<- colnames(KGJ)[2]

MCor_Junip<- dcc(Juniper, clima, selection = -6:9, method = "correlation", dynamic = "moving", win_size = 35, win_offset = 1)

SIG<- MCor_Junip$coef$significant
SIG$ROWS<- rownames(SIG)
SIG<- gather(SIG, "WINDOWS", "SIG", 1:ncol(SIG)-1)

MCor_Junip<- MCor_Junip$coef$coef
MCor_Junip$ROWS<- rownames(MCor_Junip)
MCor_Junip<- gather(MCor_Junip, "WINDOWS", "COR", 1:ncol(MCor_Junip)-1)
MCor_Junip$VAR<- lapply(strsplit(as.character(MCor_Junip$ROWS), "\\."), "[",1)
MCor_Junip$YEAR<- lapply(strsplit(as.character(MCor_Junip$ROWS), "\\."), "[",2)
MCor_Junip$MONTH<- lapply(strsplit(as.character(MCor_Junip$ROWS), "\\."), "[",3)
MCor_Junip$CODE<- paste(MCor_Junip$YEAR, MCor_Junip$MONTH, sep = " ")
MCor_Junip <- as.data.frame(lapply(MCor_Junip, unlist))

MCor_Junip$SIG<- SIG$SIG

order<- c("curr sep", "curr aug", "curr jul", "curr jun", "curr may", "curr apr", "curr mar", "curr feb", "curr jan", "prev dec", "prev nov", "prev oct", "prev sep", "prev aug", "prev jul", "prev jun")

ggplot(MCor_Junip, aes(x=WINDOWS, y=factor(CODE,level=order), fill=COR, alpha=SIG=="TRUE"))+geom_tile()+
  facet_wrap2(~VAR, axes = "all")+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  labs(x= "", y= "")+
  guides(alpha = "none")+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Graphs/OLD_juniper_cor.tiff", height = 150, width = 400, units = "mm", dpi = 300)

ggplot(Juniper, aes(x=as.numeric(rownames(Juniper)), y=GJ10F))+geom_line()+
  labs(x= "", y="TRW (mm)")+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom")+ ## umístění legendy dolu
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))


ggsave("Graphs/OLD_juniper.tiff", height = 100, width = 200, units = "mm", dpi = 300)

