getlegends <- function(col.tab){
  lgd.posx = c(1.1)
  lgd.posy = c(1.7)
  colors <- col.tab
  names(colors) <- names(col.tab)
  ###legend is plotted in second column ## see the layout function call in the script.
  par(mar = c(2, 0, 2, 2) + 0.1)
  #plot(-1:1,-1:1,type='n', axes = FALSE, ann = FALSE)
  par(lend=0.2); legend(lgd.posx,lgd.posy, inset=.02, fill=colors, border=gray(.3), legend=names(colors), cex=0.9)
}

readfiles <- function(filinName){
  alagudatatable <- fread(filinName, header = TRUE, sep ="\t", na.string = "NA",colClasses=c("numeric", "character","numeric","numeric","character"))
  return (alagudatatable)
}

makepfamabundaceTable <- function(alagudatatables, sample122, sample123, sample124, sample125){
  # Summing the RPKM by each pfam
  ST122 <- alagudatatable[sampleid %in% sample122 ,list(alagurpkmST122=sum(RPKM)), keyby='pfam']
  ST123 <- alagudatatable[sampleid %in% sample123 ,list(alagurpkmST123=sum(RPKM)), keyby='pfam']
  ST124 <- alagudatatable[sampleid %in% sample124 ,list(alagurpkmST124=sum(RPKM)), keyby='pfam']
  ST125 <- alagudatatable[sampleid %in% sample125 ,list(alagurpkmST125=sum(RPKM)), keyby='pfam']
  Tablelist <- list(ST123, ST124, ST125)
  All4Table <- ST122
  for ( table in Tablelist ) {
    All4Table <-merge(All4Table,table,by.x="pfam", all=T)
  }
  All4Table[is.na(All4Table)] <- 0
  return (All4Table)
}


makeorganismtable <- function(alagudatatable, sampleids, All4Table){
  alagudatatableonsampleid <- alagudatatable[sampleid %in% sampleids]
  organismtables = c()
  AllorganismTable <- alagudatatableonsampleid[pfam==All4Table$pfam[1] ,list(A=sum(RPKM)), keyby='organism']
  toreplacedefaultname=All4Table$pfam[1]
  names(AllorganismTable)[2] = paste(toreplacedefaultname)
  for (i in 2:length(unique(alagudatatable$pfam))){
    toreplacedefaultname=All4Table$pfam[i]
    organismtables <- alagudatatableonsampleid[pfam==All4Table$pfam[i] ,list(replaceme=sum(RPKM)), keyby='organism']
    names(organismtables)[2] = paste(toreplacedefaultname)
    AllorganismTable <-merge(AllorganismTable,organismtables,by.x=c("organism"), all=T)
  }
  AllorganismTable[is.na(AllorganismTable)] <- 0
  return (AllorganismTable)
}

getgeometry <- function(){
  par(lend=0.5);text(-1.04,0, label=c("ST122"), font=2, cex =1.2, srt=90 ) 
  par(lend=11);text(0,1.04, label=c("ST123"),font=2, cex =1.2, srt=0) 
  par(lend=12);text(1.04,0, label=c("ST124"),font=2, cex =1.2, srt=-90) 
  par(lend=13);text(0,-1.04, label=c("ST125"),font=2, cex =1.2, srt=0)
}

addgridlines <- function(){
  ## white line to overwrite the triangle base
  par(lend=10);segments(x0=c(-0.98), y0=c(0), x1=c(0.98), y1=c(0), col = "white",lwd =4 )
  ## Add grid lines to the plots
  par(lend=8);segments(x0=c(-0.166,-0.333,-0.500,-0.666,-0.833), y0=c(-0.833,-0.666,-0.500,-0.333,-0.166 ), y1=c(0.166,0.333,0.500,0.666,0.833), x1=c(0.833, 0.666, 0.500, 0.333, 0.166), col = "grey50", lty=2) 
  par(lend=9);segments(x0=c(-0.833,-0.666,-0.500,-0.333,-0.166), y0=c(0.166,0.333,0.500,0.666,0.833), y1=c(-0.833,-0.666,-0.500,-0.333,-0.166 ), x1=c(0.166,0.333,0.500,0.666,0.833),col = "grey50", lty=2) 
  ### line to extend all four sides
  par(lend=11);segments(x0=c(1,-1,0,0), y0=c(0,0,1,-1), y1=c(0,0,2,-2), x1=c(2,-2,0,0), col = "grey50", lty=2)
  #par(lend=12);segments(x0=c(-1), y0=c(0), y1=c(0), x1=c(-2), col = "grey50", lty=2 lwd=2)
}

proportionCalc <- function(f122,f123,f124,f125){
  top <- 1
  #      two
  #        -
  #      -----
  #  one---------three
  #       ---
  #        -
  #      four
  
  ##### upper triangle 
  xpoint1=0
  xpoint2=0
  ypoint1=0
  ypoint2=0
  xpoint3=0
  xpoint4=0
  ypoint3=0
  ypoint4=0
  ###upper
  st1 = f122 + f123 + f124
  xpoint1 = -1*(1-2*((f124/st1) + ((f123/st1)/2)))
  ypoint1 = (f123/st1) * top
  ##down
  st2 = f122+ f124+f125
  xpoint2 = -1*(1-2*((f124/st2) + ((f125/st2)/2)))
  ypoint = (f125/st2) * top
  ypoint2 = -1 * ypoint
  
  ##left
  st3 = f122+ f123+f125
  xpoint = -1*(1-2*((f123/st3) + ((f122/st3)/2)))
  ypoint3 =xpoint
  ypoint = (f122/st3) * top
  xpoint3 = -1 * ypoint
  ## right
  st4 = f123+ f124+f125
  xpoint = 1*(1-2*((f125/st4) + ((f124/st4)/2)))
  ypoint = (f124/st4) * top
  xpoint4 = ypoint
  ypoint4 = xpoint
  ##### point of intersection of the 4 triangles ######## 
  p0_x = xpoint1; p0_y = ypoint1; p1_x =xpoint2; p1_y = ypoint2
  p2_x = xpoint3; p2_y = ypoint3; p3_x =xpoint4; p3_y = ypoint4
  s1_x = p1_x - p0_x; s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;s2_y = p3_y - p2_y;

  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
  if (s >= 0 && s <= 1 && t >= 0 && t <= 1){
    midx = p0_x + (t * s1_x)
    midy = p0_y + (t * s1_y)
  }
  
  

  #######To check whether the points are at really intersecting ###
#   if (xpoint1 !=0 & xpoint2 != 0 & xpoint3 !=0 & xpoint4 != 0 & ypoint1 !=0 & ypoint2 != 0 & ypoint3 !=0 & ypoint4 != 0){
#   segments(x0=xpoint1, y0=ypoint1, x1=xpoint2, y1=ypoint2, col = "red",lwd =4 )
#   segments(x0=xpoint3, y0=ypoint3, x1=xpoint4, y1=ypoint4, col = "blue",lwd =4 )
#   }
#   
  print (midx)
  print (midy)
  ######OUTER### for top and down
  #if ((midx > - 0.1 & midx < 0.1) & (0.9*f123 > f122) & (0.9*f123 > f124) & (0.9*f125 > f122) & (0.9*f125 > f124) ){
  if ((0.9*f123 > f122) & (0.9*f123 > f124) & (0.9*f125 > f122) & (0.9*f125 > f124) ){
  
    if (midy <= 0){
      midy = (midy/2)-1.1
    }
    if (midy > 0){
      midy = (midy/2)+1.1
    }
  }
  ######OUTER#### for left and right
  if ((0.9*f122 > f123) & (0.9*f122 > f125) & (0.9*f124 > f123) & (0.9*f124 > f125)){
    if (midx <= 0) {
      midx = (midx/2)-1.1}
    if (midx > 0){
      midx = (midx/2)+1.1
    }
  }
  ##################################

  #### special case as it is NaN, NaN ??
  if (f122 > 0 & f123 > 0 & f124 == 0 & f125 == 0){midx = xpoint1; midy = ypoint1}
  if (f122==0 & f123 > 0 & f124 > 0 & f125 == 0){midx = xpoint1; midy = ypoint1}
  if (f122 == 0 & f123 == 0 & f124 > 0 & f125 > 0){midx = xpoint2; midy = ypoint2}
  if (f122 > 0 & f123 ==0 & f124 == 0 & f125 > 0){midx = xpoint2; midy = ypoint2}

  if (f122==0 & f123 ==0 & f124 ==0){midx = xpoint2; midy = ypoint2}
  if (f123==0 & f124 ==0 & f125 ==0){midx = xpoint3; midy = ypoint3}
  if (f122==0 & f124 ==0 & f125 ==0){midx = xpoint1; midy = ypoint1}
  if (f122==0 & f123 ==0 & f125 ==0){midx = xpoint1; midy = ypoint1}
  
  midxmidy<-c(midx,midy)
  return(midxmidy)
}

getnewpieradii <- function(get_pie_radii){
  icount = 1
  get_newpie_list = c()
  for (oldradii in get_pie_radii){
    if(oldradii >= 1){get_newpie_list[[icount]] = c(5)}
    else if(oldradii >= 0.5 & oldradii < 1){get_newpie_list[icount] = c(4)}
    else if(oldradii >= 0.1 & oldradii < 0.5){get_newpie_list[icount] = c(3.5)}
    else if(oldradii > 0.01 & oldradii < 0.1){get_newpie_list[icount] = c(3)}
    else if(oldradii > 0 & oldradii < 0.01){get_newpie_list[icount] = c(2)}
    icount = icount +1
  }
  return (get_newpie_list)
}

##############pie parameters ###############
ploteachpie <- function(piex, piey, pie_totalsum, AllorganismTable,pie_names){ 
  pie_x = piex
  pie_y = piey
  ### get the abundance for pie size
  get_pie_radii =c()
  get_pie_list =c()
  AllorganismTableorgdel <- data.table(AllorganismTable)
  AllorganismTableorgdel[,organism:=NULL]
  AllorganismTableorgdelmatrix = as.matrix(AllorganismTableorgdel)
  maxall = max(sapply(colnames(AllorganismTableorgdelmatrix), function(x) sum(AllorganismTableorgdelmatrix[,x])))
  for (i in 1:ncol(AllorganismTableorgdelmatrix)){
    get_pie_radii[[i]] <- (sum(AllorganismTableorgdelmatrix[,i])*5)/maxall
    vectorx<- nv(AllorganismTableorgdelmatrix[,i], AllorganismTable$organism)
    assign(paste("vector",i,sep=""),vectorx) 
    get_pie_list[[i]] <- vectorx
  }
 
  ### function to adjust the size of pie to make it visible ###
  get_newpie_list <- getnewpieradii(get_pie_radii)
  #### Forcing default color selected to a selected ones ####
  col.tab =c("blue", "green","red","violet","darkgreen","darkorange1","darkkhaki","deepskyblue","darkturquoise","gold","lightpink","lightseagreen","midnightblue","steelblue4","maroon","grey85")
  colors <- nv(leghead(sstable(data.frame(AllorganismTable), idx.clmns=c('organism'), ct.clmns=colnames(AllorganismTable)[-1]), n=15), 'color')
  names(col.tab) <-  names(colors)
  par(new=TRUE); pies(get_pie_list, x0=pie_x, y0= pie_y, radii= get_newpie_list, border=c('black'), color.table = col.tab)
  text(pie_x, pie_y, labels=pie_names, cex=0.75, pos=3, offset = 0.8)
  return (col.tab)
}

###############Program starts here ########
library(caroline)
library("grid")
library(data.table)
### layout for two columns one for diamond and the other for legend
#layout(matrix(c(1,1), nrow = 1), widths = c(0.7, 0.3))
layout(matrix(c(1,1), nrow = 1), widths = c(1))
par(mar = c(2, 2, 2, 0) + 0.1)
pointdesc <- seq(-1.6, 1.6, by = 0.1)
plot(pointdesc,pointdesc, type='n', axes = TRUE, ann = FALSE)
#plot(-1:1,-1:1,type='n', axes = FALSE, ann = FALSE)
###polygon and grid lines #######
top <- 1
size=unit(1, "lines")
polygon(c(-1, 0, 1), c(0, top, 0), lwd = 2)
polygon(c(-1, 0, 1), c(0, -top, 0), lwd = 2)
getgridlines <- addgridlines()
#### Declaare variables ####
piex = c()
piey = c()
pie_totalsum = c()
###organism table ### ## you need to use ff library to read large files see that later
filinName = "/Users/velucham/ENS/TARA_STEFI/STEFI_IronPlot/All4stations/AllgenesHits_Taxadded_243553_numberedFormatnew.txt"
alagudatatable <- readfiles(filinName)
###All samples
#sample122<-c(782,784,786,788,790,2054,2264,2265,2266,2267);sample123<-c(806,808,810,812,822,2059,2061,2062,2063,2274);sample124<-c(826,828,830,832,834,2065,2066,2067,2287,2288);sample125<-c(842,844,846,848,850,2073,2074,2075,2076,2295);
###0.8
#sample122<-c(782, 2264);sample123<-c(806,2059);sample124<-c(826,2065);sample125<-c(842,2073);
#sample S5
#sample122<-c(784,2265);sample123<-c(808,2061);sample124<-c(828,2287);sample125<-c(844,2074);
#sample S20
#sample122<-c(786,2266);sample123<-c(810,2062);sample124<-c(830,2288);sample125<-c(846,2295);
#sample 180
sample122<-c(788,2267);sample123<-c(812,2274);sample124<-c(832,2066);sample125<-c(848,2075);

sampleids <- c(sample122,sample123,sample124,sample125)
######Two main data table creating function ### both makes pfam abundance table
All4Table <- makepfamabundaceTable(alagudatatable, sample122,sample123,sample124,sample125)
###Table of organism with pfam
AllorganismTable <- makeorganismtable(alagudatatable, sampleids, All4Table)
## Vertex naming ## if you want to automate use a vector here and pass to it.

### Calculate the cartesian co-ordinate ######
for (i in 1:nrow(All4Table)){
  totalst = All4Table$alagurpkmST122[i]+All4Table$alagurpkmST123[i]+All4Table$alagurpkmST124[i]+All4Table$alagurpkmST125[i]
  pie_totalsum[[i]] <- totalst
  print (totalst)
  xandy <- proportionCalc(All4Table$alagurpkmST122[i],All4Table$alagurpkmST123[i],All4Table$alagurpkmST124[i],All4Table$alagurpkmST125[i])
  piex[[i]] <-xandy[1]
  piey[[i]] <-xandy[2]
}
#### ranking according to abundance and then color selection
col.tab <- ploteachpie(piex,piey, pie_totalsum, AllorganismTable,All4Table$pfam)
plotlegend <- getlegends(col.tab)
getthegeometry <- getgeometry()
##################FINISHED CODING##################################
