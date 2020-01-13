# Utilites need as part of R plotting
# From Nancy Kiang

#utils.R
#source("entdiag995fn.R")
#source("mvn.R")
#source("plot_fit.R")
library(sp)
library(spam) 
library(maps) 
library(fields) 
library(maptools) 
library(rworldmap) 
library(SDMTools) 
library(plotrix)
library("RNetCDF")
#library("animation") #For using im.convert. Not working due OS 10.14 bug with gs.

Kelvin=273.15

rho_H2O.kg.m3 = 997 #kg-H2O m-3

mod = function(n, m) {
    return(n - m * floor(n/m))
}

mod.mid = function(n,m) {
    #Rounds to nearest multiple of m (instead of floor)
    diffn = mod(n,m)
    if (diffn >= m/2) {
        return(m*floor(n/m) + m)        
    } else {
        return(m*floor(n/m))
    }
}

div0.array2 =  function(num, div, undefin = -1e30, undefout=0) {
    #For 2D array
    dims = dim(num)
    ndim = length(dims)
    
    divresult = array(0, dim=dims)
    for (i in dims[1]) {
        for (j in dims[2]) {
            if (num[i,j]==undefin | div[i,j]==undefin) {
                divresult[i,j] = undefout
            } else if (div[i,j]==0) {
                divresult[i,j] = 0
            } else {
                divresult[i,j] = num[i,j]/div[i,j]
            }
        }
    }
    return(divresult)
}


resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


plot.blank = function(x=0, y=0) {
    plot(x,y,type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
}

ma = function(y,window) {
    #Moving average.  Window should be odd length centered.
    movingaverage = NULL
    mid = trunc(window/2)
    for (i in (mid+1):(length(y)-mid)) {
        movingaverage = c(movingaverage, mean(y[(i-mid):(i+mid)]))
        }
}

my.ma = function(y,period=3) {
    y.ma = NULL
    for (i in 1:(length(y)-(period-1))) {
        y.ma = c(y.ma,na.mean(y[i:(i+period-1)]))
        }
    endNA = rep(NA,(period-1)/2)
    return(c(endNA, y.ma,endNA))
    }
    
my.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


na.mode = function(v) {
    index = !is.infinite(v) & !is.na(v)
    if (sum(index)==0) {
        return(NA) 
    } else {
        return(mode(na.omit(v[index])))
    }
}

na.mean = function(v) {
    index = !is.infinite(v) & !is.na(v)
    if (sum(index)==0) {
        return(NA) 
    } else {    
        return(mean(na.omit(v[index])))
    }
    }
    
na.max = function(v) {
    index = !is.infinite(v) & !is.na(v)
    if (sum(index)==0) {
        return(NA)
    } else {
        return(max(na.omit(v[index])))
    }
    }
    
na.min = function(v) {
    index = !is.infinite(v) & !is.na(v)
    if (sum(index)==0) {
        return(NA)
    } else {
        return(min(na.omit(v[index])))
    }
}

na.sum = function(v) {
    index = !is.infinite(v)
    return(sum(na.omit(v[index])))
    }
    
na.var = function(v) {
    index = !is.infinite(v)
    return(var(na.omit(v[index])))
    }

plot.ma = function(t,y,period=3, type="l",xlab="day",ylab=NULL) {
    y.ma = my.ma(y,period=period)
    t.ma = my.ma(t,period=period)
    plot(t,y,type=type,,xlab=xlab,ylab=ylab)
    lines(t.ma,y.ma,col=4)
    title(paste(ylab))
    }
    

prf = function(dat, titleouter="", type="l") {
    
    for (i in 1:ncol(dat)) {
        plot(dat[,i],type=type, ask=TRUE)
        title(names(dat)[i])
        mtext(outer=TRUE,paste(titleouter),line=-1.5)
    }
}


prfd = function(dat, xlab="time", titleouter="", type="p", line=-1.5, if.NA=FALSE, legx=0) {
    namesnull = NULL
    for (i in 2:ncol(dat)) {
        if (sum(is.na(dat[,i]))==nrow(dat)) {
            namesnull = c(namesnull, names(dat)[i])
            print(paste("All NA:",names(dat)[i]))
        } else {
            plot(dat[,1], dat[,i],pch=".",xlab=xlab, ylab=names(dat)[i], type=type)
            if (if.NA) {
                index = is.na(dat[,i])
                points(dat[index,1],rep(min(na.omit(dat[,i])),sum(index)),pch=".",col=2)
            }           
            title(names(dat)[i])
        }
        mtext(outer=TRUE,paste(titleouter),line=line)
    }
    return(namesnull)
}

MONTH = c("January","February","March","April","May","June","July","August","September","October", "November","December")
MONcap = c('JAN','FEB','MAR','APR','MAY','JUN'
        ,'JUL','AUG','SEP','OCT','NOV','DEC' )
MON = c('Jan','Feb','Mar','Apr','May','Jun'
        ,'Jul','Aug','Sep','Oct','Nov','Dec' )
        
SEASON4cap = c('JFM', 'AMJ', 'JAS', 'OND')

monthdays.reg = c(31,28,31,30,31,30,31,31,30,31,30,31)
monthdays.leap = c(31,29,31,30,31,30,31,31,30,31,30,31)

get.month = function(jday, year) {
         if (((year-1988)/4)!=trunc((year-1988)/4)) {
                monthdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
        } else { #is leap year
                monthdays = c(31,29,31,30,31,30,31,31,30,31,30,31)
        }
        jcount = 0
        m = 0
        while (m<12 & jcount<jday) {
            m = m+1
            jcount = jcount+monthdays[m]
        }
        return(list(m, MONcap[m]))
}

jday.month = function(jday, is.leap) {
         if (!is.leap) {
                monthdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
        } else { #is leap year
                monthdays = c(31,29,31,30,31,30,31,31,30,31,30,31)
        }
        jcount = 0
        m = 0
        while (m<12 & jcount<jday) {
            m = m+1
            jcount = jcount+monthdays[m]
        }
        return(list(m, MONcap[m]))
}

month.char2num = function(MONchar) {
    return((1:12)[MONcap==MONchar])
}

get.jday = function(monthdayyear) {
        month = monthdayyear[1]
        day = monthdayyear[2]
        year = monthdayyear[3]
#        if (mod(year-1988,4)!=0) {
        if (((year-1988)/4)!=trunc((year-1988)/4)) {
                monthdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
        } else { #is leap year
                monthdays = c(31,29,31,30,31,30,31,31,30,31,30,31)
        }
        return(sum(monthdays[0:(month-1)]) + day)
}

is.leap = function(year) {
        if (((year-1988)/4)!=trunc((year-1988)/4)) {
                return(FALSE)
        } else { #is leap year
                return(TRUE)
        }
    }

year.days = function(year) {
    #year must be a vector
    year = as.vector(year)
    isleap = sapply(year,is.leap)
    ydays = rep(365,length(year))
    if (sum(isleap)>0) {
        ydays[isleap]=366
        }
    #ydays = NULL
    #for (i in 1:length(year)) {
    #   if (is.leap(year[i])) {
    #       ydays = c(ydays,366)
    #   } else {
    #       ydays = c(ydays,365)    
    #   } 
    #}
    return(ydays)
}

get.yeartime.ymdhm = function(ymdhm=as.data.frame(cbind(month=1,day=1,year=1995))) {
    #Return yeartime as year.fraction for given data frame of
    # ymdhm = [year month day hour minute]
    year=1
    month=2
    day=3
    hour=4
    minute=5
    #jday = apply(ymdhm[,c(month,day,year)], 1, get.jday)-1
    #return(ymdhm[,year]
    #   +(jday + (ymdhm[,hour] + ymdhm[,minute]/60)/24)
    #   /year.days(ymdhm[,year]))

    ymdh = ymdhm[,1:4]
    ymdh[,4] = ymdh[,hour]+ymdhm[,minute]/60
    return(apply(ymdh,1,get.yeartime))
    }

get.yeartime = function(ymdh) {
        #Returns times in years and fractions down to hour
        #Use as apply(YMDHmatrix, 1, get.yeartime)
        #ymd = c(year, month, dayofmonth, hour), or nx4 matrix of these values     
        return(ymdh[1] + 
                (get.jday(ymdh[c(2,3,1)])-1 + ymdh[4]/24)/year.days(ymdh[1]))
}

get.yearhalfhr = function(mdy1=c(1,1,2001),mdy2=NULL) {
    #Make yeartime vector that is half-hour time points starting with
    #monthdayyear1 and ending with last half hour of monthdayyear2
    y1 = mdy1[3]
    jd1 = get.jday(mdy1)
    yd = year.days(y1)
    
    #Do first year piece
    yrhalfhr = y1 + (((jd1-1)*48+1):(yd*48)-1)/(yd*48)

    if (!is.null(mdy2)) {
        y2 = mdy2[3]
        jd2 = get.jday(mdy2)
        yd = year.days(y2)
        #Do middle full years
        if (y2>(y1+1)) { 
        for (y in (y1+1):(y2-1)) {
            yd = year.days(y)
            yrhalfhr = c(yrhalfhr, y+((1:(yd*48))-1)/(yd*48))
            }
        }
        #Do end year piece
        y = y2
        yrhalfhr = c(yrhalfhr, y + (0:(jd2*48-1))/(year.days(y)*48))
    }
    return(yrhalfhr)
}

#---------------------  
    
plot.cat = function(x,y,categ, xlab="x",ylab="y",if.legend=TRUE, legx=NULL,legy=NULL,xlim=NULL,ylim=NULL) {
        ycat = unique(categ)
        if (is.null(xlim)) xlim=c(na.min(x),na.max(x))
        if (is.null(ylim)) ylim=c(na.min(y),na.max(y))
        index = categ==ycat[1]
        plot(x[index],y[index],xlab=xlab,ylab=ylab,xlim=xlim,
            ylim=ylim,pch=16, cex=.5)
        for (i in 2:length(ycat)) {
            index = categ==ycat[i]
            points(x[index],y[index],pch=16,cex=.5,col=i)
        }
        if (if.legend) {
            legend(legx,legy,legend=ycat,col=1:length(ycat),pch=16)
        }
}

#--------------------------------------------------------------
plot.multi =
function(x="", ydat, xlab="", ylab="",ymin=0.0, legx=0, type="l",if.sepleg=FALSE, if.leg=TRUE) {
    if (x=="") {
        ymax = max(ydat,na.rm=TRUE)
        plot(ydat[,1],type=type, ylab=ylab, xlab=xlab, ylim=c(ymin,ymax))
        for (j in 2:ncol(ydat)) {
            if (type=="l") {
                lines(ydat[,j],col=j,lty=j)
            } else {  
                points(ydat[,j],col=j)
            }
        }
        if (if.leg) {
            if (if.sepleg) {
                plot(0,0,type="n", axes=F, xlab="",ylab="")
                legx = -0.01
                legy = 0.01
            } else {
                legx = legx
                legy = ymax*0.9
            }
            legend(legx,legy,legend=paste(names(ydat)), col=c(1:ncol(ydat)),lty=1:ncol(ydat))
        }
    } else {
        ymax = max(ydat,na.rm=TRUE)
        plot(x, ydat[,1],type=type, ylab=ylab, xlab=xlab, ylim=c(ymin,ymax))
        for (j in 2:ncol(ydat)) {
            if (type=="l") {
                lines(x, ydat[,j],col=j,lty=j)
            } else {
                points(x, ydat[,j],col=j)
            }
        }
        if (if.leg) {
            if (if.sepleg) {
                plot(c(0,1),c(0,1),type="n", axes=F, xlab="",ylab="")
                legx = 0.5
                legy = 1
            } else {
                legx = legx
                legy = ymax*0.9
            }
            legend(legx,legy,legend=paste(names(ydat)), col=c(1:ncol(ydat)),lty=1:ncol(ydat))
        }
    }
}
#------------------------------------------------------------
dsum = function(dhalfhr,y) {
    jday = floor(dhalfhr)
    yd = 60*60*24 * tapply(y,jday, FUN="na.mean")
    return( cbind(jday,yd))
    
    }
#------------------------------------------------------------
prfa = function(dtime, dat, xlab="time", ylabin=NULL, titleouter="", type="p", line=-1.5, if.NA=FALSE, legx=0) {
    
    ddat = NULL
    jday = floor(dtime)

    for (i in 1:ncol(dat)) {
        #Plot given data
        if (is.null(ylabin)) {
            ylab= names(dat)[i]
        }else{
            ylab=ylabin
        }
        plot(dtime, dat[,i],pch=".",xlab=xlab, ylab=names(dat)[i], type=type)
        #Plot daily means (assume dtime is in jday.fraction)
        ddat=tapply(dat[,i], jday, FUN="na.mean")
        lines(unique(jday),ddat,col=4)
        
        #Plot points where NA are
        if (if.NA) {
            index = is.na(dat[,i])
            points(dat[index,1],rep(min(na.omit(dat[,i])),sum(index)),pch=".",col=2)
        }
        title(names(dat)[i])
        mtext(outer=TRUE,paste(titleouter),line=line)
    }
}

#------------------------------------------------------------
prfp = function(dat, titleouter="") {
    
    for (i in 1:ncol(dat)) {
        plot(dat[,i],pch=".", cex=3)
        title(names(dat)[i])
        mtext(outer=TRUE,paste(titleouter),line=-1.5)
    }
}

#------------------------------------------------------------
gridxy = function(lat=0,lon=0, res="4x5") {
    #lat = latitude in degrees
    #lon = longitude in degrees
    #res = grid solution in degrees
    #   options: 1x1, 2x2.5, 4x5
    # Returns x=lon,y=lat of grid cell
    if (res=="1x1") {
        x = max(1,ceiling(lon+180))
        y = max(1,ceiling((lat+90)))
    } else if (res=="2x2.5") {
        x = max(1,ceiling((lon+180)/2.5))
        y = max(1,ceiling((lat+90)/2))
    } else if (res=="4x5") {
        x = max(1,ceiling((lon+180)/5))
        y = ceiling((lat+90+2)/4)
    }
    return(c(x,y))
}

#-------------------
grid.lon.lat = function(res) {
        if (res=="5Mx5M") {
        #1/12 degree
        i = (-12*180):(12*180 -1) + 0.5
        j = (-12*90):(12*90 -1) + 0.5
    } else if (res=="QxQ" | res=="qxq"| res=="1440x720") {
        #0.25 degree
        i = ((-4*180):(4*180 -1) + 0.5)/(4)
        j = ((-4*90):(4*90 -1) + 0.5    )/(4)
    } else if (res=="HxH" | res=="HXH") {
        #0.5 degree
        i = -360:359 + 0.5
        i = i*0.5
        j = -180:179 + 0.5
        j = j*0.5
    } else if   (res=="1x1") {
        i = (-180:179) + 0.5
        j = (-90:89) + 0.5
    } else if (res=="2x2.5") {
        i = (-72:71)*2.5 + 1.25
        i = i*2.5
        j = (-45:44)*2 + 1
        j = j*2
        i = i/2.5
        j = j/2
    } else if (res=="4x5") {
        di = 5
        i = ((-180/di):((180-di)/di))*di + di/2
        dj = 4
        j = ((-90/dj):((90)/dj))*dj
        j[1] = -89 #cannot have point at pole
        j[length(j)] = 89
    }
    return(list(i,j))
}

#-------------------
plot.grid.categorical = function(mapz, res="1x1", zlim=NULL, colors=terrain.colors(40), xlab="longitude", ylab="latitude", ADD=FALSE) {
    #Plot map of categorical values.  
    #NOTE:  image function stretches z values over entire colors range.  Therefore, pass in colors array that
    #       matches range of your mapz values.
    #Grid centers (x,y)
    gridij = grid.lon.lat(res)
    i = gridij[[1]]
    j = gridij[[2]]
    
    image(x=i,y=j,mapz, xlab=xlab, ylab=ylab, col=colors, add=ADD)
}

#-----------------
plot.grid.continuous = function(mapz, res="1x1", colors=terrain.colors(40), legend.lab=NULL, xlab="longitude", ylab="latitude", titletext="", zlim=NULL, ADD=FALSE, if.fill=TRUE, if.coasts=FALSE, ask=TRUE) {
    #Plot map of continuous values, gridded, and fill in z extremes with colors limits.
    #Grid centers (x,y)
    #mapz is continuous values
    if (res=="5Mx5M") {
        #1/12 degree
        i = (-12*180):(12*180 -1) + 0.5
        j = (-12*90):(12*90 -1) + 0.5
    } else if (res=="QxQ" | res=="qxq" | res=="1440x720") {
        #0.25 degree
        i = ((-4*180):(4*180 -1) + 0.5)/(4)
        j = ((-4*90):(4*90 -1) + 0.5    )/(4)
    } else if (res=="HxH" | res=="HXH") {
        #0.5 degree
        i = -360:359 + 0.5
        i = i*0.5
        j = -180:179 + 0.5
        j = j*0.5
    } else if   (res=="1x1") {
        i = (-180:179) + 0.5
        j = (-90:89) + 0.5
    } else if (res=="2x2.5" | res=="2HX2" ) {
        i = (-72:71)*2.5 + 1.25
        i = i*2.5
        j = (-45:44)*2 + 1
        j = j*2
        i = i/2.5
        j = j/2
    } else if (res=="4x5") {
        di = 5
        i = ((-180/di):((180-di)/di))*di + di/2
        dj = 4
        j = ((-90/dj):((90)/dj))*dj
    }
    #print(res)
    #print(i)
    #print(j)
    if (sum(!is.na(mapz))==0) {
        mapzlim = matrix(0, dim(mapz)[1], dim(mapz)[2]) 
    } else {
        mapzlim = mapz
    }
    if (if.fill & !is.null(zlim)) {  #Color extreme values with the zlim colors
        mapzlim[mapzlim<zlim[1]] = zlim[1]
        mapzlim[mapzlim>zlim[2]] = zlim[2]
    }
    image.plot(x=i,y=j,mapzlim, xlab=xlab, ylab=ylab, zlim=zlim, col=colors, xaxt="n", yaxt="n", add=ADD, ask=ask)
    if (if.coasts) {
            plot(coastsCoarse, add=TRUE)
    }
    title(titletext)
}


#-----------------
plot.grid3 = function(mapz, res="1x1", colors=terrain.colors(40), legend.lab=NULL, xlab="longitude", ylab="latitude", zlim=NULL, ADD=FALSE) {
    #Plot map of continuous values contoured
    #Grid centers (x,y)
    if (res=="30arcsec") {
        #30 arc sec ~ 1 km
        i = ((-120*180):(120*180 - 1) + 0.5)/120
        j = ((-120*90):(120*90 - 1) + 0.5)/120
    } else  if (res=="1x1") {
        i = (-180:179) + 0.5
        j = (-90:89) + 0.5
    } else if (res=="2x2.5") {
        i = (-72:71)*2.5 + 1.25
        j = (-45:44)*2 + 1
    }
    filled.contour(x=i, y=j, z=t(mapz), xlab=xlab, ylab=ylab, col=colors, zlim=zlim)
}

#----------------

giss.palette = colorRampPalette(c("dark blue", "blue","light blue", "cyan","white", "yellow", "orange",  "red","dark red"),    space = "rgb")
#giss.palette.nowhite = colorRampPalette(c("dark blue", "blue", "light blue", "cyan", "yellow", "orange", "red","dark red"),    space = "rgb")
giss.palette.nowhite = colorRampPalette(c("blue", "light blue", "cyan", "yellow", "orange", "red"),    space = "rgb")

drywet = colorRampPalette(c("tan", "yellow", "green", "dark green", "dark blue"), spac="rgb")

giss.band = c(300,770,860,1250,1500,2200,4000)
modis.band.visnir = c(300, 700, 5000)
#----------------

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11,  x=0, y0=0, title='', if.add=TRUE) {
    scale = (length(lut)-1)/(max-min)
	ticks=seq(min, max, len=nticks)
    if (if.add) {
	    dev.new(width=1.75, height=5)
	    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	}
    axis(2, ticks, las=1)
    if (is.null(x)) {
    	x = 0
    }
    if (is.null(y0)) {
    	y0=0
    }
    for (i in 1:(length(lut)-1)) {
    	 y = (i-1)/scale + min
    	 rect(x,y,x+10,y+1/scale, col=lut[i], border=NA)
  	}
}

#----------------
albedo.avg.fn = function(bandlims, wavin, albin) {
    #bandlims:  bound nm boundaries, therefore length number of bands + 1
    nbands = length(bandlims)-1
    alband = array(0, nbands)
    for (b in 1:nbands) {
        alband[b] = mean(albin[wavin>=bandlims[b] & wavin<bandlims[b+1]])
    }
    return(alband)
}
lines.bands = function(bandlims, abands, col="black", lty=1, lwd=1) {
    #bandlims:  bound nm boundaries, therefore length number of bands + 1
    #abands:  should be length nbands, albedo of each band  
    nbands = length(bandlims)-1
    lines( array(t(cbind(bandlims, bandlims)),  2*length(bandlims))[2:(2*length(bandlims)-1)],
       array(t(cbind(abands, abands)),  2*nbands)[1:(2*nbands)], col=col, lty=lty, lwd=lwd)
}

giss.albedo.avg.fn = function(wavin, albin) {
    alband = array(0,6)
    for (b in 1:6) {
        alband[b] = mean(albin[wavin>=giss.band[b] & wavin<giss.band[b+1]])
    }
    return(alband)
}

lines.gissbands = function(agiss, col="black", lty=1, lwd=1) {
    lines( array(t(cbind(giss.band, giss.band)),  14)[2:13],
       array(t(cbind(agiss, agiss)),  12)[1:12], col=col, lty=lty, lwd=lwd)
}

#---------------- Ent specific functions ----------------------------------------------------
entgvsd_pagetitles = function(paths, fname, dy=0, cex=1) {
    mtext(outer=TRUE, paths, line=1.5+dy, cex=cex)
    mtext(outer=TRUE, fname, line=0+dy, cex=cex)
    #mtext(outer=TRUE, date(), line=-1,cex=0.5, adj=1)
    mtext(outer=TRUE, date(), line=3+dy,cex=0.5, adj=1)
}

Ent_Type13 = c("evergreen broadleaf trees", "evergreen needleleaf trees", "cold deciduous broadleaf trees", "drought deciduous broadleaf", "deciduous needleleaf", "cold-adapted shrub", "arid-adapted shrub", "C3 grass perennial", "C4 grass", "C3 grass annual", "arctic C3 grass", "C4 herb crops", "tree crops")

EntGVSD_PFTs = c("ever_br_early",
         "ever_br_late ",
         "ever_nd_early",
         "ever_nd_late ",
         "cold_br_early",
         "cold_br_late ",
         "drought_br   ",
         "decid_nd     ",
         "cold_shrub   ",
         "arid_shrub   ",
         "c3_grass_per ",
         "c4_grass     ",
         "c3_grass_ann ",
         "c3_grass_arct",
         "crops_herb   ",
         "crops_woody  ",
         "bare_bright  ",
         "bare_dark    ")
  
EntGVSD_PFT13 = EntGVSD_PFTs[c(2,4,6,7:15)]
EntGVSD_COV13 = EntGVSD_PFTs[c(2,4,6,7:18)]

EntGVSD_COVER20 = c("ever_br_early",
         "ever_br_late ",
         "ever_nd_early",
         "ever_nd_late ",
         "cold_br_early",
         "cold_br_late ",
         "drought_br   ",
         "decid_nd     ",
         "cold_shrub   ",
         "arid_shrub   ",
         "c3_grass_per ",
         "c4_grass     ",
         "c3_grass_ann ",
         "c3_grass_arct",
         "crops_c3_herb",
         "crops_c4_herb",
         "crops_woody  ",
         "permanent_ice", #"snow_ice     ",
         "bare_sparse  ",
         "water        ")
EntGVSD_COV17 = EntGVSD_COVER20[c(2,4,6,7:16,18:20)]         

Ent_diags = c("vf", "Anet", "Atot", "Rd", "GCANOPY", "TRANS_SW", "LAI", "Resp_fol",
"Resp_sw", "Resp_lab", "Resp_root", "Resp_maint", "Resp_growth_1", "Resp_growth",
"GPP", "R_auto", "C_total", "C_lab", "C_fol", "C_sw", "C_hw", "C_froot", "C_croot", "C_soil",
"Resp_soil", "phenofactor", "betad")

Ent_diags_LUT = data.frame(entdiagname=Ent_diags, units=c("cover fraction", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "m s-1", "fraction", "m^2 m-2", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 s-1", "kgC m-2 gnd", "kgC m-2 gnd", "kgC m-2 gnd", "kgC m-2 gnd", "kgC m-2 gnd", "kgC m-2 gnd", "kgC m-2 gnd","kgC m-2 gnd", "kgC m-2 s-1", "index", "fraction") )


EntPFTchar3 = c(paste("00",1:9,sep=""), paste("0",10:16,sep=""))
EntGCMdiags27 = c(paste("ra", "00",1:9, sep=""), paste("ra" ,"0",10:27, sep=""))
temp = t(matrix(rep(EntGCMdiags27, 16), 27,16))
EntGCMdiags27pft = paste(temp, EntPFTchar3, sep="")

Entcolors17 = as.data.frame(t(array(dim=c(5, 1+1+length(EntGVSD_COVER20)),
c(
0.7, 0.7, 0.7, "grey",       "water        ",
0.0, 0.3, 0.0, "dark green", "ever_br_early",
0.0, 0.3, 0.0, "dark green", "ever_br_late ",
0.0, 0.4, 0.4, "dark blue-green", "ever_nd_early",
0.0, 0.4, 0.4, "dark blue-green", "ever_nd_late ",
0.0, 0.6, 0.0, "medium green", "cold_br_early",
0.0, 0.6, 0.0, "medium green", "cold_br_late ",
0.4, 0.5, 0.3, "taupe", "drought_br   ",
0.5, 0.1, 0.1, "dark red", "decid_nd     ",
0.0, 0.6, 0.6, "blue-green", "cold_shrub   ",
0.95, 0.85, 0.65, "tan", "arid_shrub   ",
0.0, 0.8, 0.1, "bright green", "c3_grass_per ",
0.7, 0.8, 0.0, "leaf green", "c4_grass     ",
1.0, 1.0, 0.0, "yellow", "c3_grass_ann ",
0.1, 0.8, 0.8, "blue-grey", "c3_grass_arct",
0.8, 0.7, 0.1, "gold", "crops_c3_herb",
0.8, 0.4, 0.0, "orange", "crops_c4_herb",
0.8, 0.0, 0.0, "red", "crops_woody  ",
0.9, 1.0, 1.0, "blue-white", "snow_ice     ",
1.0, 0.95, 0.8, "pale tan", "bare_sparse  ",
0.5, 0.5, 0.5, "grey",      "water        ",
1.0, 1.0, 1.0, "white",     "undef        "
#0.0, 0.3, 0.0, "dark green", "dummy R plotting ever_br_early",
#0.0, 0.4, 0.4, "dark blue-green", "ever_nd_early",
#0.0, 0.6, 0.0, "medium green", "cold_br_early"
))))

names(Entcolors17) = c("r","g","b", "color", "lc_type")
Entcolors17 = as.data.frame(cbind(num=0:(nrow(Entcolors17)-1), Entcolors17))
Entcolors17[,c("r")] = as.numeric(as.character(Entcolors17[,c("r")]))
Entcolors17[,c("g")] = as.numeric(as.character(Entcolors17[,c("g")]))
Entcolors17[,c("b")] = as.numeric(as.character(Entcolors17[,c("b")]))


Entrgbhex = function(Entcolors=Entcolors17) {
	#Colors in hex for Ent 17 PFTs + 3 non-veg cover
	rgbhex = NULL
	for (i in 1:nrow(Entcolors)) {
		rgbhex = c(rgbhex, rgb(Entcolors[i,"r"],Entcolors[i,"g"], Entcolors[i,"b"]))
	}
	return(rgbhex)
}

Ent17rgbhex = Entrgbhex(Entcolors=Entcolors17)

Ent17legend = function(colors=Ent17rgbhex(Entcolors=Entcolors17), newwindow=FALSE) {
    #Quick check of map colors
    if (newwindow) {
        quartz(width=6,height=4)
    }
    cex = 0.3
    plot(0,0, xlim=c(-180,180),ylim=c(-90,90), type="n", bty="n", xaxt="n", yaxt="n",xlab="",ylab="")
    par(xpd=TRUE)
    colors=Ent17rgbhex  
    legend(-180, 150, legend=Entcolors17[1:20,"num"], col=colors[1:20], pch=15, cex=cex)
    legend(-150, 150, legend=Entcolors17[1:20,"KGcode"], col=colors[1:20], pch=15, cex=cex)
    legend(-110,150, legend=Entcolors17[1:20,"color"], col=colors[1:20], pch=15, cex=cex)

    legend(0,150, legend=Entcolors17[21:39, "num"], col=colors[21:39], pch=15, cex=cex)
    legend(30,150, legend=Entcolors17[21:39, "KGcode"], col=colors[21:39], pch=15, cex=cex)
    legend(70,150, legend=Entcolors17[21:39, "color"], col=colors[21:39], pch=15, cex=cex)
    
}

#Version with water last.
Entcolors17b = as.data.frame(t(array(dim=c(5, 1+1+length(EntGVSD_COVER20)),
c(
0.0, 0.3, 0.0, "dark green", "ever_br_early",
0.0, 0.3, 0.0, "dark green", "ever_br_late ",
0.0, 0.4, 0.4, "dark blue-green", "ever_nd_early",
0.0, 0.4, 0.4, "dark blue-green", "ever_nd_late ",
0.0, 0.6, 0.0, "medium green", "cold_br_early",
0.0, 0.6, 0.0, "medium green", "cold_br_late ",
0.4, 0.5, 0.3, "taupe", "drought_br   ",
0.5, 0.1, 0.1, "dark red", "decid_nd     ",
0.0, 0.6, 0.6, "blue-green", "cold_shrub   ",
0.95, 0.85, 0.65, "tan", "arid_shrub   ",
0.0, 0.8, 0.1, "bright green", "c3_grass_per ",
0.7, 0.8, 0.0, "leaf green", "c4_grass     ",
1.0, 1.0, 0.0, "yellow", "c3_grass_ann ",
0.1, 0.8, 0.8, "blue-grey", "c3_grass_arct",
0.8, 0.7, 0.1, "gold", "crops_c3_herb",
0.8, 0.4, 0.0, "orange", "crops_c4_herb",
0.8, 0.0, 0.0, "red", "crops_woody  ",
0.9, 1.0, 1.0, "blue-white", "snow_ice     ",
1.0, 0.95, 0.8, "pale tan", "bare_sparse  ",
0.5, 0.5, 0.5, "grey",      "water        ",
1.0, 1.0, 1.0, "white",     "undef        ",
0.7, 0.7, 0.7, "grey",       "water        "
#0.0, 0.3, 0.0, "dark green", "dummy R plotting ever_br_early",
#0.0, 0.4, 0.4, "dark blue-green", "ever_nd_early",
#0.0, 0.6, 0.0, "medium green", "cold_br_early"
))))

names(Entcolors17b) = c("r","g","b", "color", "lc_type")
Entcolors17b = as.data.frame(cbind(num=0:(nrow(Entcolors17b)-1), Entcolors17b))
Entcolors17b[,c("r")] = as.numeric(as.character(Entcolors17b[,c("r")]))
Entcolors17b[,c("g")] = as.numeric(as.character(Entcolors17b[,c("g")]))
Entcolors17b[,c("b")] = as.numeric(as.character(Entcolors17b[,c("b")]))


EntcolorsbareGISS = as.data.frame(t(array(dim=c(6, 2),
c(17, 1.0, 0.95, 0.8, "pale tan",   "bare_bright  ",
  18, 0.5, 0.45, 0.4, "dark brown", "bare_dark    " ))))
names(EntcolorsbareGISS) = names(Entcolors17)
Entcolorsundef = 
as.data.frame(t(array(dim=c(6, 1),
c(0, .9, .9, .9, "white",        "undef        "))))
names(Entcolorsundef) = names(Entcolors17)

  
Entcolors16 = rbind( Entcolors17[match(c(21, 1:15,17), Entcolors17[,"num"]),], EntcolorsbareGISS)
Entcolors16[pmatch("undef", Entcolors16[,"lc_type"]), "num"] = 0

#-------------------------------------------------------------------------------------------------
map.entgvsd.steps = function(entlclaidir, res, enttyp=enttyp, varname, trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = TRUE, pathplot="", do.checksum=TRUE) { 
    # Plot maps of lc, laimax, lai, or hgt
    # entlclaidir:  Directory containing subdirectories of trim options Ent GVSD files
    # res:          Grid resolution of data file.  "ent17" and "pure" V1km are plotted at qzq.  Trimmed files are at HXH.
    # enttyp:         ent17: 1:20;  pure, trimmmed...nocrops:  1:18
    # varname:    lc, laimax, hgt, lai
    # trimopt:    ent17, pure, trimmed, trimmed_scaled, trimmed_scaled_nocrops, trimmed_scaled_crops_ext1
    # filepre:    File prefix (Vres_EntGVSD<PFTs>_<LAILdata>, e.g. V1km_EntGVSD17G_BNUM, VHXH_EntGVSD16G_BNUM 
    # datatime:   Time point of data, e.g.:  2004
    # version:    EntGVSD version, e.g.:  v0.1, v1.1
    # filesuf:    Extra suffix for miscellaneous, e.g. "_forplot" for scaled up V1km for plotting at coarser resolution.
    # do.pdf = TRUE  Logical if to output pdf instead of to the screen.
    # pathplot=""     If do.pdf, then give output path for plots.

   #quartz(width=11,height=6) #Open plotting window
   print(varname)
    
   if (varname == "lc") {
      zlim = c(0,1)
      restime = "_ann"
      colors = giss.palette.nowhite(40)
   } else if (varname == "laimax" ) {
      zlim = c(0,7)
      restime = "_ann"
      colors=drywet(40)
   } else if (varname=="lai" | varname=="lclai_checksum") {
      zlim = c(0,7)
      restime = ""
      colors=drywet(40)   
   } else if (varname=="laimax_err" ) {
      zlim = c(-1,1)
      restime = "_ann"
      colors=giss.palette(40)
   } else if (varname == "hgt") {
      zlim = c(0,40)
      restime = "_ann"
      colors=drywet(40)
   } else {
    print(paste(varname, "Whoops, only does lc or laimax"))
    return()
   }
        
  for (opt in trimopt) {
    #if (if.trim) { #TRUE only for Ent 16 PFTs
    if (opt=="ent17") {  #Ent 17 PFTs
#        fname = paste(filepre, "_",varname,"_",datatime,restime, "_", opt, "_",version, filesuf,".nc", sep="")
        fname = paste(filepre, "_",version,"_", icov, "_", idat, "_", varname,"_",datatime,restime, "_", opt, filesuf,".nc", sep="")
    } else { #Ent 16 PFTs pure, trimmed, trimmed_scaled, trimmed_scaled_nocrops
        #fname = paste(filepre, "_lc_max_", opt, "_", version, filesuf, ".nc", sep="")
	#        fname = paste(filepre, "_",varname,"_", datatime, restime,"_", opt, "_",  version, filesuf, ".nc", sep="")
        fname = paste(filepre, "_",version,"_", icov, "_", idat, "_", varname,"_",datatime,restime, "_", opt, filesuf,".nc", sep="")	
    }
    
    #fnamelc = paste(filepre, "_","lc","_",datatime,restime, "_", opt, "_",version, filesuf,".nc", sep="")
    fnamelc = paste(filepre, "_",version,"_", icov, "_", idat, "_", "lc","_",datatime,restime, "_", opt, filesuf,".nc", sep="")
    
    if (do.pdf) { 
    	filepdf = paste(pathplot, fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7) 
    } else if (add.new) {
            quartz(width=11,height=7) #Open plotting window
    } else {
        print("Printing to screen.")
    }
    #filelc = paste(entlclaidir, opt, "/", fnamelc, sep="")
    filevar = paste(entlclaidir, opt, "/", fname, sep="")
    #print(filelc)
    print(filevar)
    #par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
    par(mfrow=c(4,5), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
    if (do.pdf) {
        par(mfrow=c(4,5), omi=c(0,0.0,1.0,0.5), mar=c(1,1,2,2)+0.1)
    }
    titletop = paste(entlclaidir, sep="")
    map.EntGVSD.v1.1(file=filevar, res=res, zlim=zlim, varname=varname, lcnum=enttyp, colors=colors)
    
    entgvsd_pagetitles(titletop, fname)

    # Checksum
    if (do.checksum) {
            print(varname)
            #for (checksum in c("lai_allmonths_checksum", "lai_checksum",  "laimax_checksum", "lc_checksum", "lc_modis_checksum")) {
#       for (checksum in c("lc", "lclai", "lclaimax", "lchgt")) { # c( "lclai_allmonth")) {
        if (varname == "lc") {
            checksum = "lc"
        } else {
            checksum = paste("lc", varname, sep="")
        }
        varnamecheck = paste(checksum, "_checksum", sep="")
        print(paste("netcdf varnamecheck: ", varnamecheck))         
        if (opt=="ent17") {  #Ent 17 PFTs
            fname = paste(filepre, "_",version, "_",icov, "_",idat, "_",varnamecheck, "_",datatime,restime, "_",opt, filesuf,".nc", sep="")
        
        } else { #Ent 16 PFTs pure, trimmed, trimmed_scaled, trimmed_scaled_nocrops
            #fname = paste(filepre, "_lc_max_", opt, "_", version, filesuf, ".nc", sep="")
            fname = paste(filepre, "_",version, "_",icov, "_",idat, "_",varnamecheck, "_",datatime,restime,"_", opt, filesuf, ".nc", sep="")
        }
        
        filevar = paste(entlclaidir, opt, "/", fname, sep="")   
        print(filevar)
        map.GCM(file=filevar, varname=varnamecheck, res=res,colors=colors,  zlim=zlim, if.zeroNA=TRUE, titletype=1) 

   } #if do.checksum
 
   if (do.pdf) {
        dev.off()
        #im.convert(filepdf, output=paste(filepdf, ".png", sep=""), extra.opts="-density 150") #Bug in Mac gs, otherwise should work
    }   
  } #for opt
}


map.entgvsd.check.misc = function(entlclaidir, res, enttyp=enttyp, varnamecheck, trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = TRUE, pathplot="") { 
    # Plot maps of npftgrid, lc_dompft, lc_dompftlc, and add any other misc here

	if (trimopt=="ent17") {
		npft=17
	} else {
		npft=16
	}

	if (varnamecheck ==  "lc_dompft" | varnamecheck == "lc_domlc") {
		zlim = c(min(enttyp), max(enttyp))
		if (trimopt=="ent17") {
			color = Entrgbhex(Entcolors17b[1:20,])
			leg=Entcolors17b[1:20, "lc_type"]
		} else {
			color = Entrgbhex(Entcolors16)
			leg=Entcolors16[,"lc_type"]
		}
		restime = "_ann"
		if.cat=TRUE
		
        } else if (varnamecheck == "lc_npftgrid" ) {
                if (FALSE) { #categorical colors do not work
	   	zlim = c(0,17)
	   	color = giss.palette.nowhite(18)
	   	color[1] = rgb(0.5, 0.5, 0.5)
	   	leg=zlim[1]:zlim[2]
	   	restime="_ann"
		if.cat=TRUE
                }
                zlim = c(0,17)
                color = giss.palette.nowhite(18)
                color[1] = rgb(0.5, 0.5, 0.5)
                leg=zlim[1]:zlim[2]
                restime="_ann"
                if.cat=FALSE                
	} else if (varnamecheck == "lc_dompftlc" | varnamecheck == "lc_checksum") {
		zlim = c(0,1)
		color = giss.palette.nowhite(11)
		color[1] = rgb(0.5, 0.5, 0.5)
		leg=(0:10)/10
		restime = "_ann"
		if.cat=FALSE
	} else if (varnamecheck == "laimax_err" | varnamecheck == "lclaimax_err"  ) {
		zlim = c(-1,1)
		color = giss.palette(40)
		leg = NULL
		restime = "_ann"
		if.cat=FALSE
	} else if (varnamecheck == "lclai_err" |varnamecheck == "lclai_checksum_diff") {
		zlim = c(-1,1)
		color = giss.palette(40)
		leg = NULL
		restime = ""
		if.cat=FALSE
	} else if (varnamecheck == "lclai_checksum" ) {
		if (length(grep('diff', filesuf))>0) {
			zlim = c(-1,1)
		} else {
			zlim = c(0,7)
		}
		color = drywet(40)
		leg = NULL
		restime = ""
		if.cat=FALSE
	} else if (varnamecheck == "lclaimax_checksum") {
		if (length(grep('diff', filesuf))>0) {
			zlim = c(-1,1)
		} else {
			zlim = c(0,7)
		}
		color = drywet(40)
		leg = NULL
		restime = "_ann"
		if.cat=FALSE		
    } else if (varnamecheck == "lchgt_checksum" ) {
    	if (length(grep('diff', filesuf))>0) {
    		zlim=c(-15,15)
    	} else {
			zlim = c(0,40)
		}
		color = drywet(40)
		leg = NULL
		restime = "_ann"
		if.cat = FALSE
    } else if (varnamecheck == "hgt_err" | varnamecheck == "lchgt_err" | varnamecheck=="lchgt_checksum_diff") {
		zlim = c(-15,15)
		color = giss.palette(40)
		leg = NULL
		restime = "_ann"
		if.cat = FALSE
	} else if (varnamecheck == "bs_brightratio") {
		zlim=NULL
		color = giss.palette(40)
		leg = NULL
		restime = ""
		if.cat = FALSE
    } else { #if (varnamecheck ==  ) {
#		zlim = c(
#		color = 
#		restime =	} else {
		print(paste(varnamecheck, "Whoops, not set up for that variable, yet."))
		return()	
	}
	
    fname = paste(filepre, "_",version, "_",icov, "_",idat, "_",varnamecheck, "_",datatime,restime,"_", trimopt, filesuf, ".nc", sep="")
	if (length(grep('diff', filesuf))>0) {
		if (varnamecheck=="lclaimax_checksum" ) {  #HACK because varname got cut off in file
	    		 varnamecheck = paste(varnamecheck, '_diff', sep="")
		} else {	 
		   varnamecheck = paste(varnamecheck, '_diff', sep="")
		}
	}		

	
    file = paste(entlclaidir, trimopt, "/", fname, sep="")
	print(file)
	
	#ncid <- open.nc(con=file, write=FALSE)
	#x = var.get.nc(ncid, varnamecheck)
	#plot.grid.categorical(x, res, colors=color)
	
	if (do.pdf) {
		pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=8, height=5) 
	} else if(!add.new) {
		quartz(width=8, height=5)
	}

   map.GCM (file, varname=varnamecheck, res=res,colors=color,  zlim=zlim, if.cat=if.cat, if.zeroNA=FALSE, titletype=1) 
   
   if (if.cat) { #Do legend for categorical variable
   par(xpd=NA)  #Let legend be outside the plot area
	#legend(200, 90, legend=zlim[1]:zlim[2], col=color, pch=15, cex=0.9)	
	legend(180, 90, legend=leg, col=color, pch=15, cex=.7, bty="n")
	
	#Alternative color bar or legend
	#par(xpd=NA)
	#legend.gradient(cbind(x = c(200,210,210,200), y = c(80,80,-80,-80)), 
    #             cols = rgbhex, title = "", limits = zlim)
	}
	
    #Ent_domlc_plot(lctype, numpft=numpft, res="qxq", legend.cex=0.6, Entcolors=Entcolors, if.new=FALSE) 
  	
  	cex=0.8
  	mtext(outer=TRUE, entlclaidir, line=-1, cex=cex)
  	mtext(outer=TRUE, fname, line=-2, cex=cex)
  	mtext(outer=TRUE, date(), line=-1, cex=0.5, adj=1)
#    entgvsd_pagetitles(paths=entlclaidir, fname, dy=-2.5, cex=0.9)

  	if (do.pdf) {
	  	dev.off()
        #im.convert(file, output=paste(file, ".png", sep=""), extra.opts="-density 150") #Bug in Mac OS10.14 gs, otherwise should work
	}
}

Ent_calc_npftgrid = function(pathin, fname, pathout, npft=17) {
	file = paste(pathin, fname, sep="")
	print(file)
	
	ncid <- open.nc(con=file, write=FALSE)
    lcin = var.get.nc(ncid, "lc")

	npftgrid = array(0, dim=dim(lcin)[1:2])
	
	for (p in 1:npft) {
		for (i in 1:dim(lcin)[1]) {
			for (j in 1:dim(lcin)[2]) {
				if (!is.na(lcin[i,j,p])) {
					if (lcin[i,j,p]>0.0 ) {
						npftgrid[i,j] = npftgrid[i,j] + 1
					}
				}
			}
		}
	}
	
	zlim = c(0,npft)
   	color = giss.palette.nowhite(max(npftgrid)-min(npftgrid)+1)
   	color[1] = rgb(0.5, 0.5, 0.5)
   	leg=zlim[1]:zlim[2]

   	par(mar=c(5,4,4,5)+0.1)
   	plot.grid.categorical(mapz=npftgrid, res=res, colors=color	, xlab="", ylab="")
   	#plot.grid.continuous(mapz=npftgrid, res=res, colors=color	, xlab="", ylab="")
    plot(coastsCoarse, add=TRUE)
	par(xpd=NA)
    legend(180, 90, legend=0:max(npftgrid), col=color, pch=15, cex=.7, bty="n")
    mtext(paste("number of PFTs per grid cell"))
 
 	fileout = paste(pathout, fname,"_npftgrid.nc", sep="")
	create.map.template.nc(res, varname="npftgrid", longname="number_of_PFTs", units="fraction", undef=-1e30, description="number of PFTs in grid cell", fileout, contact="Nancy.Y.Kiang@nasa.gov")
	close.nc(ncid)
	ncid = open.nc(fileout, write=TRUE)
	var.put.nc(ncid, 'npftgrid', npftgrid)
   
    return(npftgrid.list=list(npftgrid=npftgrid, file.nc=fileout) )
}
	
Ent_calc_lc_checksum = function(pathin, fname, pathout, enttyp=1:20) {
	
	file = paste(pathin, fname, sep="")
	print(file)
	ncid <- open.nc(con=file, write=FALSE)
    lcin = var.get.nc(ncid, "lc")
	
	lcchecksum = apply(lcin, c(1,2), na.sum)

	
	fileout = paste(pathout, fname,"_checksum.nc", sep="")
	create.map.template.nc(res, varname="lc_checksum", longname="cover", units="fraction", undef=-1e30, description="lc checksum", fileout, contact="Nancy.Y.Kiang@nasa.gov")
	close.nc(ncid)
	ncid = open.nc(fileout, write=TRUE)
	var.put.nc(ncid, 'lc_checksum', lcchecksum)
	return(lcchecksum.list = list(lcchecksum=lcchecksum, file.nc=fileout))
}


Ent_calc_domlc = function(file, enttyp=1:20) {
	
	print(file)
	ncid <- open.nc(con=file, write=FALSE)
    lcin = var.get.nc(ncid, "lc")
    lc = lcin
	lc[is.na(lc)] = 0
	
	domlc = array(NA, dim(lc)[1:2])
	for (p in enttyp) {
		maxlc = apply(lc, c(1,2), max)
		for (i in 1:dim(lc)[1]) {
			for (j in 1:dim(lc)[2]) {
				if (lc[i,j,p]>0 & lc[i,j,p]==maxlc[i,j]) {
					domlc[i,j] = p
				}
			}
		}
	}
	#domlc[is.na(domlc)] = 0
	#domlc[lc==0] = NA
	return(domlc)
}

Ent_domlc_plot = function(lctype, numpft=17, res="HXH", legend.cex=0.6, Entcolors=Entcolors17[1:20], if.new=FALSE) {
	#Plot maps of Ent GVSD dominant cover types with nice color scheme.
	#Works for any number of Ent PFTs.  MUST carefully specify the Entcolors table according to order of cover types in lctype.!!!
	
	rgbhex = Entrgbhex(Entcolors)
	ncov = max(Entcolors[,"num"]) #dim(Entcolors)[1]
	
	#Screwy R skipping over drought-broad if water is last.  Annoying
	#lctype[lctype==20] = 0 #Put water first at zero to plot in R.  
	
	if (if.new) {
		quartz(width=9.6, height=6)
	}
	#quartz(width=10.6, height=6)
	#par(omi=c(0,0,0,1)) #(bottom, left, top, right)
	par(omi=c(0,0,0,0), oma=c(0,0,0,4)) #(bottom, left, top, right) 	#Use for single

	plot.grid.categorical(lctype, res=res,color=rgbhex)
	par(xpd=NA)
	#legend.gradient(cbind(x = c(200,210,210,200), y = c(80,80,-80,-80)), 
    #             cols = rgbhex, title = "", limits = c(0,ncov))
    legend(180, 90, legend=Entcolors[,"lc_type"], col=rgbhex, pch=15, cex=.7, bty="n")

	#par(xpd=TRUE)
	#pt.cex=1.5
	#legend(-180, 124, legend=Entcolors[1:7,"lc_type"],col=rgbhex[1:7], pt.cex=pt.cex, pch=15, cex=legend.cex, horiz=TRUE, bty="n")
	#legend(-180, 113, legend=Entcolors[8:14,"lc_type"],col=rgbhex[8:14], pt.cex=pt.cex, pch=15, cex=legend.cex, horiz=TRUE, bty="n")
	#legend(-180, 102, legend=Entcolors[15:ncov,"lc_type"],col=rgbhex[15:ncov], pt.cex=pt.cex, pch=15, cex=legend.cex, horiz=TRUE, bty="n")
		
	plot(coastsCoarse, add=TRUE)
}




make.GCM.Ent.diag.name = function(varname="vf", p="ever_br_late") {
    #Make GCM diagnostic name for Ent, raVVVPPP are names of diagnostics Ent_diags VVV for EntGVSD_PFTs PPP
    # in 3-character digits.  E.g ra001001 is vf ever_br_early
    #Default example should produce ra001002
    vindex = which(Ent_diags %in% varname)
    pindex = which(trim(EntGVSD_PFTs) %in% p)
    
    if (vindex<10) {
        v3 = paste("00",vindex, sep="")
    } else {
        v3 = paste("0", vindex, sep="")
    }
    if(pindex <10) {
        p3 = paste("00", pindex, sep="")
    } else {
        p3 = paste("0", pindex, sep="")
    }
    
    return(paste("ra", v3, p3, sep=""))
}

get.GCM.Ent.diag.units = function(varname="vf") {
    #Return units for GCM Ent diagnostic ra00#.
    #Default returns for "vf" cover fraction.
    vindex = which(Ent_diags %in% varname)
    
    return(Ent_diags_LUT[vindex,"units"])
}

#------------
map.EntGVSD <- function(filelc=NULL, file, res="2x2.5", varpre="", varlist=EntGVSD_PFT13, colors=giss.palette.nowhite(40), type="any", zlim=c(0,1), if.zeroNA=TRUE, titletype=1) {
    #file = netcdf file path and name
    #type = 
    # any-plot all layers continuous
    # lc-land cover all layers continuous
    # lcdom - dominant PFT categorical
    # lai-lai all layers continuous;
    # laiall - average lai per grid cell
    # height - height all layers continuous
    # heightall - cover-weighted average height per grid cell
    # heighttree - height of tree PFTs only
    
    ncid <- open.nc(con=file, write=FALSE)
    if (type=="any") {
        for (p in varlist) {
            p = paste(varpre, trim(p), sep="")
            print(paste(p))
            x = var.get.nc(ncid, p)
            if (if.zeroNA) {x[x==0]=NA}
            plot.grid.continuous(mapz=x, res=res,colors=colors,         
                    xlab="", ylab="", 
                    zlim=zlim)
            plot(coastsCoarse, add=TRUE, lwd=0.25)
            if (titletype==1) {
                long_name = att.get.nc(ncid, p, attribute="long_name")
                mtext(paste(long_name), cex=0.6)
            } else {
                units = att.get.nc(ncid, p, attribute="units")
                title(paste(p, " (",units,")", sep=""))
            }
        }
    } 
}

#------------
map.EntGVSD.v1.1 <- function(filelc=NULL, file, res="2x2.5", varpre="", varname="lc", lcnum=c(2,4,6,7:15), colors=giss.palette.nowhite(40), type="any", zlim=c(0,1), if.zeroNA=TRUE, titletype=1) {
    #This version reads files formatted as 3D arrays of var[IM,JM,layers]
    #file = netcdf file path and name
    #type = 
    # any-plot all layers continuous
    # lc-land cover all layers continuous
    # lcdom - dominant PFT categorical
    # lai-lai all layers continuous;
    # laiall - average lai per grid cell
    # hgt - height all layers continuous
    # heightall - cover-weighted average height per grid cell
    # heighttree - height of tree PFTs only
    
    ncid <- open.nc(con=file, write=FALSE)
    if (type=="any") {
        x = var.get.nc(ncid, varname)
        if (grepl("err", varname)) {
        	print("err variable")
        	x[is.na(x)] = 0
        }
        lcnames = var.get.nc(ncid, "lctype")
        for (p in lcnum) {
            print(paste(p))
            if (if.zeroNA) {x[x==0]=NA}
            plot.grid.continuous(mapz=x[,,p], res=res,colors=colors,        
                    xlab="", ylab="", 
                    zlim=zlim)
            plot(coastsCoarse, add=TRUE, lwd=0.25)
            if (titletype==1) {
                long_name = att.get.nc(ncid, varname, attribute="long_name")
                mtext(paste(long_name), cex=0.6, line=1)
                mtext(paste(trim(lcnames[p]), "(", round(na.min(x[,,p]),2), round(na.mean(x[,,p]),2), round(na.max(x[,,p]),2), ")"), cex=0.6)
            } else {
                units = att.get.nc(ncid, varname, attribute="units")
                title(paste(lcname[p], " (",units,")", sep=""))
            }
        }
    } 
    close.nc(con=ncid)
}

#------------

map.EntGVSD.3Darray <- function(file, res="2x2.5", varnc="lc", varlist=1:20, colors=giss.palette.nowhite(40), zlim=c(0,1), if.zeroNA=TRUE, titletype=1, titletop="", fname="") {
    #file = netcdf file path and name
    #type = 
    # any-plot all layers continuous
    # lc-land cover all layers continuous
    # lcdom - dominant PFT categorical
    # lai-lai all layers continuous;
    # laiall - average lai per grid cell
    # height - height all layers continuous
    # heightall - cover-weighted average height per grid cell
    # heighttree - height of tree PFTs only
    
    ncid <- open.nc(con=file, write=FALSE)
    varin <- var.get.nc(ncid, varnc)
    layers <- var.get.nc(ncid, "layers")
    print(paste(dim(varin)))
    for (p in varlist) {
        print(paste(varnc, p))
        x = varin[,,p]
        if (if.zeroNA) {x[x==0]=NA}
        plot.grid.continuous(mapz=x, res=res,colors=colors,         
                xlab="", ylab="", 
                zlim=zlim)
        plot(coastsCoarse, add=TRUE, lwd=0.25)
        if (titletype==1) {
            long_name = att.get.nc(ncid, varnc, attribute="long_name")
            mtext(paste(varnc, layers[p]), cex=0.6)
        } else {
            units = att.get.nc(ncid, varnc, attribute="units")
            title(paste(varnc,  " (",units,")   ", layers[p], sep=""))
        }
        entgvsd_pagetitles(titletop, fname) #file
    } 
}


#----------------
map.GCM <- function(file, varname="tsurf", res="2x2.5",colors=giss.palette(40),  zlim=NULL, if.cat=FALSE, if.zeroNA=TRUE, titletype=1) {
    #Map a single GCM diagnostic with varname
    #file = netcdf file path and name
    #varname = netcdf variable name
    
    print(varname)
    ncid <- open.nc(con=file, write=FALSE)
    x = var.get.nc(ncid, varname)
    if (if.zeroNA) {x[x==0]=NA}
    if (!if.cat) { 
	    plot.grid.continuous(mapz=x, res=res,colors=colors,         
    	    xlab="", ylab="", 
        	zlim=zlim)
    } else {
    	par(mar=c(5,4,4,5)+0.1)
    	plot.grid.categorical(mapz=x, res=res, colors=colors	, xlab="", ylab="")
    }
    plot(coastsCoarse, add=TRUE)
    if (titletype==1) {
        long_name = att.get.nc(ncid, varname, attribute="long_name")
        units = att.get.nc(ncid, varname, attribute="units")
        mtext(paste(long_name," (",units,")", sep=""), cex=0.6, line=1)
        mtext(paste("min=",round(na.min(x),1), ", mean=", round(na.mean(x),1),
            "  max=", round(na.max(x))), cex=0.5)
    } else if (titletype==2){
        units = att.get.nc(ncid, varname, attribute="units")
        title(paste(varname, " (",units,")", sep=""))
    } else if (titletype==3) {
        title(paste(varname))
    }
} 

#----------------
map.GCM.Ent <- function(filelc=NULL, res="2x2.5", file, varname="vf", pftlist=EntGVSD_PFT13, colors=giss.palette.nowhite(40), type="any", zlim=NULL, unitstype=1, if.zeroNA=TRUE, titletype=1) {
    #Map a GCM Ent diagnostic with varname, for all PFTs in pftlist (netcdf names)
    #file = netcdf file path and name
    #type = 
    # any-plot all layers continuous
    # lc-land cover all layers continuous
    # lcdom - dominant PFT categorical
    # lai-lai all layers continuous;
    # laiall - average lai per grid cell
    # height - height all layers continuous
    # heightall - cover-weighted average height per grid cell
    # heighttree - height of tree PFTs only
    
    ncid <- open.nc(con=file, write=FALSE)
    if (type=="any") {
        for (p in pftlist) {
            p = trim(p)
            diagname = make.GCM.Ent.diag.name(varname, p)
            diagunits = get.GCM.Ent.diag.units(varname)
            print(diagname)
            x = var.get.nc(ncid, diagname)
            if (if.zeroNA) {x[x==0]=NA}
            plot.grid.continuous(mapz=x, res=res,colors=colors,         
                    xlab="", ylab="", 
                    zlim=zlim)
            plot(coastsCoarse, add=TRUE)
            if (titletype==1) {
                long_name = att.get.nc(ncid, p, attribute="long_name")
                mtext(paste(long_name), cex=0.6, line=1)
            } else if (titletype==2){
                units = att.get.nc(ncid, p, attribute="units")
                mtext(paste(p, " (",units,")", sep=""), line=1)
            } else if (titletype==3) {
                mtext(paste(varname, p), line=1)
            }
            mtext(paste(diagunits, "(", round(na.min(x[,]),2), round(na.mean(x[,]),2), round(na.max(x[,]),2), ")"), cex=0.6)

        }
    } 
}

#------------
create.map.template.nc = function(res, varname, longname, units, undef=-1e30, description, fileout, contact="Nancy.Y.Kiang@nasa.gov") {
    lon.lat = grid.lon.lat(res)
    IM = length(lon.lat[[1]])
    JM = length(lon.lat[[2]])

    ncid <- create.nc(filename=fileout)
    dim.def.nc(ncid, "lon", IM)
    dim.def.nc(ncid, "lat", JM)

    var.def.nc(ncid, 'lon', 'NC_FLOAT', 'lon')
    var.def.nc(ncid, 'lat', 'NC_FLOAT', 'lat')
    var.def.nc(ncid, varname, 'NC_FLOAT', dimensions=c('lon','lat'))
    
    att.put.nc(ncid, 'lon', 'long_name', 'NC_CHAR', 'longitude degrees east')
    att.put.nc(ncid, 'lat', 'long_name', 'NC_CHAR', 'latitude degrees north')
    att.put.nc(ncid, varname, 'long_name', 'NC_CHAR', longname)

    att.put.nc(ncid, varname, 'units', 'NC_FLOAT', units)
    att.put.nc(ncid, varname, '_FillValue', 'NC_FLOAT', undef)

    att.put.nc(ncid, 'NC_GLOBAL', 'Description', 'NC_CHAR', description)
    att.put.nc(ncid, 'NC_GLOBAL', 'Contact','NC_CHAR', contact)
    att.put.nc(ncid, 'NC_GLOBAL', 'Date created','NC_CHAR', paste(date()))

    close.nc(ncid)
    open.nc(fileout, write=TRUE)
    var.put.nc(ncid, 'lon', lon.lat[[1]])
    var.put.nc(ncid, 'lat', lon.lat[[2]])
    var.put.nc(ncid, varname, matrix(undef, IM,JM), start=c(1,1), count=c(IM,JM))
    close.nc(ncid)
} 


#------------
plot.obs.diff2 <- function(x0, x1, x2
        , txt0="OBS", txt1="MODEL1", txt2="MODEL2", txt0b="", txt1b="", txt2b=""
    , var0, var1,var2
        , zlim0=NULL, zlim0n=NULL, zlim12=NULL
    ,res="2x2.5", colors=giss.palette(40))
{#Compare two (2) runs x1 and x2 to obs x0.

   #x0
    plot.grid.continuous(mapz=x0,res=res,colors=colors,
                           xlab="", ylab="", zlim=zlim0)
    plot(coastsCoarse, add=TRUE)
    mtext(paste(txt0), line=1.5, cex=0.7)
    mtext(paste(var0,txt0b), cex=0.7)

   #x1-x0
        plot.grid.continuous(mapz=x1-x0,res=res,colors=colors,
                           xlab="", ylab="", zlim=zlim0n)
        plot(coastsCoarse, add=TRUE)
        mtext(paste(txt1,"- OBS"), line=1.5, cex=0.7)
        mtext(paste(var1, txt1b, "  ",round(na.mean(x1-x0),2)), cex=0.7)

   #x2-x0
        plot.grid.continuous(mapz=x2-x0,res=res,colors=colors,
                           xlab="", ylab="", zlim=zlim0n)
        plot(coastsCoarse, add=TRUE)
        mtext(paste(txt2,"- OBS"), line=1.5, cex=0.7)
        mtext(paste(var2, txt2b, "  ",round(na.mean(x2-x0),2)), cex=0.7)

   #x2-x1
        plot.grid.continuous(mapz=x2-x1,res=res,colors=colors,
                           xlab="", ylab="", zlim=zlim12)
        plot(coastsCoarse, add=TRUE)
        mtext(paste(txt2,"-",txt1), line=1.5, cex=0.7)
        mtext(paste(var2, txt2b, "  ",round(na.mean(x2-x1),2)), cex=0.7)

}

make.newx = function(rmod, nout=20) {
newx = NULL
rdat = rmod$model
factorset = NULL
varset = NULL
for (j in 2:ncol(rdat)) {
    print(names(rdat)[j])
    if (names(rdat)[j]!="(weights)") {
        if (!is.factor(rdat[,j])) {
            newx = cbind(newx, seq(min(rdat[,j]), max(rdat[,j]), length.out=nout))
            varset = c(varset, names(rdat)[j])
        } else {
        factorset = c(factorset, names(rdat)[j])
        }
    }
}
newx= data.frame(newx)
names(newx) = varset
print(paste(varset))

print(factorset)

if (length(which(factorset=="TL"))>0) {
    TL = as.factor(newx[,"rotation.XEarth.day"] == newx[,"orb.period.Earthdays"])
    newx = cbind(newx, TL)
}
if (length(which(factorset=="E0"))>0) {
    E0 = as.factor(newx[,"eccentricity"] == 0 )
    newx = cbind(newx, E0)
}
if (length(which(factorset=="OBL20"))>0) {
    #OBL0 = ensemble[,"obliquity.degrees"] == 0
    OBL20 = as.factor(newx[,"obliquity.degrees"] > 20 )
    newx = cbind(newx, OBL20)
}
if (length(which(factorset=="AQUA"))>0) {
    AQUA = as.factor(newx[,"Land.fraction"] == 0 )
    newx = cbind(newx, AQUA)
}
if (length(which(factorset=="mars"))>0) {
    mars = as.factor(newx["Experiment"] == '"Ancient Mars"' )
    newx = cbind(newx, mars)
}
if (length(which(factorset=="Stellar.type"))>0) {   
    Stellar.type = as.factor(newx[,"Tstar.K"]<3500 )
    index = Stellar.type=="M"
    newx = cbind(newx, Stellar.type)
}
if (length(which(factorset=="S0Xlt1"))>0) { 
    if (length(which(varset=="S0X"))>0) {
        S0Xlt1 = as.factor(newx[,"S0X"] < 1 )
    } else {
        S0Xlt1 = as.factor(newx[,"logS0X"] < 0)
    }
    newx = cbind(newx, S0Xlt1)
}
if (length(which(factorset=="Venus"))>0) {  
    Venus = as.factor(newx[,"Experiment"]=="Venus 1.5" | newx[,"Experiment"]=="Venus 1.9" | newx[,"Experiment"]=="Venus 2.4" )
    newx = cbind(newx, Venus)
}
if (length(which(factorset=="CO2mbge0.2"))>0) { 
    CO2mbge0.2 = as.factor(ensemble[,"CO2.mb"] >= 0.2)
}

#names(newx) = names(rmod$model)[2:ncol(rmod$model)]
#names(newx) = c(varset, factorset)
return(newx)
}

simultaneous_CBs <- function(linear_model, newdata, level = 0.95, if.lines=TRUE){
    #Source: https://stats.stackexchange.com/questions/231632/how-to-plot-simultaneous-and-pointwise-confidence-bands-for-linear-regression-wi
    # Working-Hotelling 1 – α confidence bands for the model linear_model
    # at points newdata with α = 1 - level
    #** THIS JUST GIVES THE SAME THING AS tcrit(alpha=0.05, df) * se.mean, or same as
    #**      dt(0.05/2, df=n-p) * predict(lmmod, se.fit)$se

    # estimate of residual standard error
    lm_summary <- summary(linear_model)
    # degrees of freedom 
    p <- lm_summary$df[1] - 1
    # residual degrees of freedom
    nmp <-lm_summary$df[2]
    # F-distribution
    Fvalue <- qf(level,p,nmp)
    # multiplier
    W <- sqrt(p*Fvalue)
    # confidence intervals for the mean response at the new points
    CI <- predict(linear_model, newdata, se.fit = TRUE, interval = "confidence",  level = level)
    # mean value at new points
    if (class(linear_model)[1] == "glm") {
            Y_h <- CI$fit
    } else {
        Y_h <- CI$fit[,"fit"]       
    }
    # Working-Hotelling 1 – α confidence bands
    LB <- Y_h - W*CI$se.fit
    UB <- Y_h + W*CI$se.fit
    
    sim_CB <- data.frame(LowerBound = LB, Mean = Y_h, UpperBound = UB)
}

#Was lines.confint
#my.confint.simultaneous 
lines.confint.simultaneous = function(lmod, varx="", newx=lmod$model, level=0.95, if.response = FALSE, if.mean=TRUE) {
    #Plot CI around the predicted value.
    #lmod is a linear model object from lm
    cf = simultaneous_CBs(lmod, newdata=newx, level = level)
    #ord = order(lmod$model[,varx])
    #dlo = cf$Mean - cf$Lower
    #dup = cf$Upper - cf$Mean

    if (if.response) {
        ord = order(cf$Mean)
        xvec = cf$Mean[ord] #lmod$model[ord,varx]
        #lines(c(min(xvec),max(xvec)), c(min(xvec), max(xvec))) #1:1, predicted range
    } else {
        ord = order(newx[,varx])
        xvec = newx[ord,varx]
    }
    if (if.mean) {
        #lines(lmod$model[ord,varx], cf$Mean[ord], lty=2)
        lines(xvec, cf$Mean[ord])
    }       
    lines(xvec, cf$Lower[ord], lty=3)
    lines(xvec, cf$Upper[ord], lty=3)
    #lines(xvec, xvec - dlo[ord], lty=3)
    #lines(xvec, xvec + dup[ord], lty=3)
    return(cf)
}

lines.confint.simultaneous2 = function(lmod, varx="", level=0.95, if.response = FALSE, if.mean=TRUE) {
    #Plot CI around observed response.
    #lmod is object from an lm
    
    newdat = lmod$model
    cf = simultaneous_CBs(lmod, newdata=newdat, level = level)
    #ord = order(lmod$model[,varx])
    #dlo = cf$Mean - cf$Lower
    #dup = cf$Upper - cf$Mean
    obs = lmod$model[,1]

    if (if.response) {
        ord = order(obs)
        #xvec = obs[ord] #lmod$model[ord,varx]
        xvec = cf$Mean[ord] #lmod$model[ord,varx]
        lines(c(min(obs),max(obs)), c(min(obs), max(obs))) #1:1, observed range

    } else {
        ord = order(newdat[,varx])
        xvec = newdat[ord,varx]
    }
    
    if (if.mean) {
        #lines(lmod$model[ord,varx], cf$Mean[ord], lty=2)
        lines(xvec, cf$Mean[ord])
    }       
    lines(xvec, cf$Lower[ord], lty=3)
    lines(xvec, cf$Upper[ord], lty=3)
    #lines(xvec, xvec - dlo[ord], lty=3)
    #lines(xvec, xvec + dup[ord], lty=3)
    return(cf)
}

#Was my.confint.se 
lines.confint.se <- function(pred, x, numse=1.96) {
    fit = pred$fit
    se.fit = pred$se.fit
    ord = order(x)
    fit = fit[ord]
    se.fit = se.fit[ord]
    xo = x[ord]
    #lines(x, exp(fit+numse*se.fit)/(1+exp(fit+numse*se.fit)), lty=2)
    #lines(x, exp(fit-numse*se.fit)/(1+exp(fit-numse*se.fit)), lty=2)
    #lines(xo, fit)
    lines(c(min(c(xo,fit)), max(c(xo,fit))), c(min(c(xo,fit)), max(c(xo,fit))))
    lines(xo, fit+numse*se.fit, lty=2)
    lines(xo, fit-numse*se.fit, lty=2)
}

my.confint.band <- function(rmod, newx, namex, if.response=FALSE, if.wide=FALSE, sel=NULL) {
    #For rmod of class glm
    print(paste(class(rmod)))
    if (class(rmod)[1] == "glm") {
        y = rmod$y
        pred = predict(rmod, newdata=newx, type="link", se.fit=TRUE)
        linkinv <- family(rmod)$linkinv
    } else { #Assume lm
        y = rmod$model[,1]
        pred = predict(rmod, newdata=newx, se.fit=TRUE)
        linkinv <- family(rmod)$linkinv
    }
    
    
    pred$pred <- linkinv(pred$fit)
    pred$LC <- linkinv(pred$fit - 1.96*pred$se.fit)
    pred$UC <- linkinv(pred$fit + 1.96*pred$se.fit)

    pred$pred <- pred$fit
        
    #points(newx[,namex], pred$pred, col=2)
    #ord = order(pred$pred)
    fit = pred$pred
    LC = pred$LC
    UC = pred$UC
    
    if (if.response) {
        xo = y
        
        plot(y, rmod$fit, xlab="", ylab="fitted", pch=16
        #, xlim=c(min(c(rmod$y, rmod$fit, fit)), max(c(rmod$y, rmod$fit, fit)))
        #,ylim=c(min(c(rmod$y, rmod$fit, fit)), max(c(rmod$y, rmod$fit, fit)))
        )
        lines(c(min(c(xo,fit)), max(c(xo,fit))), 
            c(min(c(xo,fit)), max(c(xo,fit))))
        if (!is.null(sel)) {
            points(y[sel], rmod$fit[sel], pch=16, col=2, cex=0.8)
        }
    } else {
        xo = newx[,namex]
        plot(rmod$model[,namex], rmod$fit, xlab="", ylab="fitted", pch=16
        , ylim=c(min(c(rmod$y, rmod$fit, fit)), max(c(rmod$y, rmod$fit, fit)))) #obs
        lines(xo, fit)
        points(rmod$model[,namex], rmod$y, pch=16, cex=0.8) #fit
        #points(xo, fit, pch=16, col=3, cex=0.6) #fit new
        if (!is.null(sel)) {
            points(rmod$model[sel,namex], rmod$y[sel], pch=16, cex=0.8)
        }
    }

    mtext(side=1, namex, line=3.5, cex=0.7)
    #lines(x, exp(fit+numse*se.fit)/(1+exp(fit+numse*se.fit)), lty=2)
    #lines(x, exp(fit-numse*se.fit)/(1+exp(fit-numse*se.fit)), lty=2)
    #lines(xo, fit)

    
    #plot(xo, fit, pch=16, col=3)
    #mtext(side=1, namex, line=2.5, cex=0.8)
    ord = order(xo)
    lines(xo[ord], LC[ord], lty=2)
    lines(xo[ord], UC[ord], lty=2)
    
    return(data.frame(fit, LC, UC))
}

my.plot.mod.conf = function(rmod, response.name="surface.albedo.planetary") {
par(mfrow=c(4,3))
newx = make.newx(rmod)
print(response.name)
my.confint.band(rmod=rmod, newx=newx, namex=response.name, if.response=TRUE)
for (j in 1:ncol(newx)) {
    print(names(newx)[j])
    if (names(newx)[j]!="(weights)" & !is.factor(newx[,j])) 
    my.confint.band(rmod=rmod, newx=newx, namex=names(newx)[j]) 
}
}

my.confint.lm = function(lmmod, x) {
newx = NULL
for (j in 2:ncol(lmmod$model)) {
    lmdat = lmmod$model
    newx = cbind(newx, seq(min(lmdat[,j]), max(lmdat[,j]), length.out=20))
}
newx= data.frame(newx)
names(newx) = names(lmmod$model)[2:ncol(lmmod$model)]

pred = predict(lmmod, newdata=newx, interval="confidence")
par(mfrow=c(4,3))
plot(lmmod$y, pred[,"fit"])
polygon(c(rev(pred[,1]), pred[,1]), c(rev(pred[ ,3]), pred[ ,2]), col = 'grey80', border = NA)
# intervals
dat = cbind(y=pred[,"fit"], newx)
par(mfrow=c(4,3), mar=c(3,4,3,3))
for (j in 1:ncol(dat)) {
    #plot(dat[,j], pred[,"fit"], xlab=names(dat)[j], ylab="y")
    plot(lmmod$model[,j], lmmod$model[,1], pch=16, cex=0.7, xlab=names(lmmod$model)[j], ylab="y", ylim=c(-0.3,0.4))
    mtext(names(lmmod$model)[j], side=1, line=2.5, cex=0.8)
    if (j==1) {
        polygon(c(rev(dat[,j]), dat[,j]), c(rev(pred[ ,3]), pred[ ,2]), col = 'grey80', border = NA)
        abline(lmmod)
        ord = order(dat[,j])
        lines(dat[ord,j], pred[ ord,3], lty = 'dashed')#, col = 'red')
        lines(dat[ord,j], pred[ ord,2], lty = 'dashed')#, col = 'red')
        points(dat[,j], lmmod$y, pch=16)
        points(dat[,j], pred[,"fit"], pch=16)
    }
    points(lmmod$model[,j], lmmod$model[,1], pch=16, cex=0.7) #obs
    points(lmmod$model[,j], lmmod$fit, pch=16, col=2, cex=0.7) #fit
    points(dat[,j], pred[,"fit"], pch=16, cex=0.5, col=3) #fit new
}

}


mean.pred.intervals <- function(lm.model, pred.x, y, pred.y) {
  n <- length(y) # Find sample size
  #lm.model <- lm(y ~ x) # Fit linear model
  y.fitted <- lm.model$fitted.values # Extract the fitted values of y
  
  # Coefficients of the linear model, beta0 and beta1
  b0 <- lm.model$coefficients[1]
  b1 <- lm.model$coefficients[2:length(lm.model$coefficients)]
  
  pred.y <- b1 * pred.x + b0 # Predict y at the given value of x (argument pred.x)
  
  # Find SSE and MSE
  sse <- sum((y - y.fitted)^2)
  mse <- sse / (n - 2)
  
  t.val <- qt(0.975, n - 2) # Critical value of t
  
  mean.se.fit <- (1 / n + ((pred.x - apply(pred.x,2,mean))^2) / (apply((pred.x - apply(pred.x,2,mean))^2, 2, sum))) # Standard error of the mean estimate
  pred.se.fit <- (1 + (1 / n) + (pred.x - apply(pred.x, 2, mean)^2) / (apply((pred.x - apply(pred.x, 2, mean)^2), 2, sum))) # Standard error of the prediction
  
  # Mean Estimate Upper and Lower Confidence limits at 95% Confidence
  mean.conf.upper <- pred.y + t.val * sqrt(mse * mean.se.fit)
  mean.conf.lower <- pred.y - t.val * sqrt(mse * mean.se.fit)
  
  # Prediction Upper and Lower Confidence limits at 95% Confidence
  pred.conf.upper <- pred.y + t.val * sqrt(mse * pred.se.fit)
  pred.conf.lower <- pred.y - t.val * sqrt(mse * pred.se.fit)
  
  # Beta 1 Upper and Lower Confidence limits at 95% Confidence
 # b1.conf.upper <- b1 + t.val * sqrt(mse) / sqrt(sum((x - mean(x))^2))
 # b1.conf.lower <- b1 - t.val * sqrt(mse) / sqrt(sum((x - mean(x))^2))
  b1.conf.upper <- b1 + t.val * sqrt(mse) / sqrt(sum((pred.x - apply(pred.x, 2, mean))^2))
  b1.conf.lower <- b1 - t.val * sqrt(mse) / sqrt(sum((pred.x - apply(x,2,mean))^2))
  
  # Build data.frame of upper and lower limits calculated above, as well as the predicted y and beta 1 values
  upper <- data.frame(rbind(round(mean.conf.upper, 2), round(pred.conf.upper, 2), round(b1.conf.upper, 2)))
  lower <- data.frame(rbind(round(mean.conf.lower, 2), round(pred.conf.lower, 2), round(b1.conf.lower, 2)))
  fit <- data.frame(rbind(round(pred.y, 2), round(pred.y, 2), round(b1, 2)))
  
  # Collect all into data.frame and rename columns
  results <- data.frame(cbind(lower, upper, fit), row.names = c('Mean', 'Prediction', 'Coefficient'))
  colnames(results) <- c('Lower', 'Upper', 'Fit')
  
  return(results)
}

#----------
xnorm = function(data, dim=2) {
    #Returns normalized data set, default applied to columns
    xmax = t(array(rep(apply(data, dim, FUN="max"), nrow(data)), dim=c(ncol(data), nrow(data))))
    xmin = t(array(rep(apply(data, dim, FUN="min"), nrow(data)), dim=c(ncol(data), nrow(data))))
    xmean = t(array(rep(apply(data, dim, FUN="mean"), nrow(data)), dim=c(ncol(data), nrow(data))))
    xnorm = (data - xmean)/(xmax - xmin)
    return(xnorm)
}

xstd = function(data, dim=2) {
    #Returns standardized data set, subtract mean and divide by standard deviation, default applied to columns
    xmean = t(array(rep(apply(data, dim, FUN="mean"), nrow(data)), dim=c(ncol(data), nrow(data))))
    xsdev = t(array(rep(apply(data, dim, FUN="var"), nrow(data)), dim=c(ncol(data), nrow(data))))   
    xsdev = sqrt(xsdev)
    xstd = (data - xmean)/xsdev
    return(xstd)
}


my.ace = function(x, y, yname='y', wt=NULL, lin=NULL, if.lm=TRUE, if.plot=FALSE, titletext="", col=1, line=-1) {
    #Converts any columns in dat that are factor to as.numeric and provides ACE a vector of columns to treat as categorical. Return ACE model and lm of ACE transforms.
    #x is data frame of the x's
    
    if (!is.null(wt)) {
        print("my.ace not set up to take weights, yet.  Assuming all 1's")
        #return(NA)
    }
    
    acex = NULL
    acexnames = NULL
    categ = NULL
    categnames = NULL
    ic = 0
    for (i in 1:ncol(x)) {
        print(paste(i, names(x)[i]))
        if (is.numeric(x[,i])) {
            ic = ic + 1
            acex = cbind(acex, x[,i])
            acexnames = c(acexnames, names(x)[i])
        } else if (is.factor(x[,i]) | is.logical(x[,i])) {
            if (length(unique(x[,i]))>1) {
                acex = cbind(acex, as.numeric(x[,i]))
                ic = ic+1
                categ = c(categ, ic)
                categnames = c(categnames, names(x)[i])
                acexnames = c(acexnames, names(x)[i])
            } else {  #Leave out factors with only one value, must have 2 or more levels
                print(paste("factor", names(x)[i], "has less than 2 levels; leaving it out."))
                # Do not increment ic.
            }
        } else {
            print ("It's something other than numeric, factor, or logical")
            return(NA)
        }
    }
    acex = data.frame(acex)
    print("dimensions of x:")
    print(dim(acex))
    names(acex) = acexnames
    print("categ")
    print(categ)
    print("names of ace x")
    print(names(acex))
    mod.ace = ace(x=acex, y=y, cat=categ, lin=lin)
    
    lm.dat = data.frame(mod.ace$tx)
    #Make lm.dat have factors as original factors, not as.numeric
    #for (i in 1:ncol(x)) {
        #if (is.factor(x[,i])) {
        #   lm.dat[,i] = x[,i]
        #   lm.dat[,i] = as.factor(lm.dat[,i])
        #} else if (is.logical(x[,i])) {
        #   lm.dat[,i] = x[,i]
        #}
    for (xname in names(lm.dat)) {
        if (!is.na(match(xname, categnames))) {
            lm.dat[,xname] = x[,xname]
        }
    }
    lm.dat = data.frame(cbind(ty=mod.ace$ty,lm.dat))
    tyname = paste('t', yname, sep="")
    names(lm.dat)[1] = tyname
    print(names(lm.dat))
    
    lm.form = make.formula(xdat=lm.dat[,-1], yname=tyname)
    print("lm.form")
    print(lm.form)
    
    acet.lm = lm(as.formula(lm.form), data=lm.dat)
    
    if (if.plot) {
        plot.ace(main=yname, aceobj=mod.ace, categ=categ, titletext=titletext, col=col, line=line)
    }
    
    return(list(yname=yname, mod.ace=mod.ace, acet.lm=acet.lm, lm.formula=lm.form, categ=categ, categnames=categnames))
}


make.formula = function(xdat, yname) {
    #Make a formula from a data set and yname.
    #Return as just text, not formula.
    
     lm.form = NULL 
     xnames = names(xdat)
    for (i in 1:ncol(xdat)) {
        #print(xnames[i])
        if (i == 1) {
            lm.form = paste(lm.form, xnames[i])
        } else {
            lm.form = paste(lm.form, '+', xnames[i])
        }

    }
    lm.form = paste(yname, ' ~', lm.form, sep="")
    print(lm.form)
    return(lm.form)
}

plot.lm = function(lmobj) {
    varnames = names(lmobj$model)
    y = lmobj$model[,1]
    x = lmobj$model[,2:ncol(lmobj$model)]
    
    plot(y, lmobj$fit, xlab=varnames[1], ylab="predicted", pch=16)
    lines( c(min(c(y, lmobj$fit)), max(c(y,lmobj$fit))), c(min(c(y, lmobj$fit)), max(c(y,lmobj$fit))), lty=2)
    
    for (i in 1:ncol(x)) {
        plot(x[,i], lmobj$fit, xlab=varnames[i+1], ylab="predicted", pch=16)
        points(x[,i], lmobj$fit, pch=16)
        points(x[,i], y, pch=1, col=2)
    }
    
    #Plot residuals
    for (i in 1:ncol(x)) {
        plot(x[,i], lmobj$resid,  xlab=varnames[i+1], ylab="residuals", pch=16)
    }
    hist( lmobj$resid, main="residuals")
}

plot.ace = function(main='y', aceobj, categ=NA, titletext="", line=-1, mfrow=c(3,3), if.plotdata=TRUE, col=1) {

    
    ylab = paste("acet(", main, ")", sep="")

    ord=order(aceobj$y) 

    xa = as.data.frame(t(aceobj$x))

    #Data
    if (if.plotdata) {
        for (i in 1:ncol(xa)) {
            if (is.na(match(i, categ))) {
                plot(xa[,i], aceobj$y, pch=16, col=col, xlab=names(xa)[i], ylab=main )
            } else {
                plot(as.factor(xa[,i]), aceobj$y, pch=16, xlab=names(xa)[i], ylab=main, pch=16)     
                #points(xa[,i], aceobj$y, pch=16)
            }
            mtext(outer=TRUE, titletext, line=line)
        }
    }
    
    #plot.new() #new plot
    par(mfrow=mfrow)
    #y vs. ty
    plot(aceobj$y[ord], aceobj$ty[ord], pch=16, type="c", xlab=main, ylab=ylab,     
        main=main, col=1,
        xlim=c(min(aceobj$y), max(aceobj$y)),
        ylim=c(min(aceobj$ty), max(aceobj$ty)))
    points(aceobj$y[ord], aceobj$ty[ord], pch=16,  col=col[ord])
    #lines(aceobj$y[ord], aceobj$ty[ord])
    mtext(outer=TRUE, titletext, line=line)
    
    #x's vs. tx
    for (i in 1:ncol(xa)) {
        xlab=names(xa)[i]
        ord = order(xa[i])
        if (ncol(xa)==1) { #Argh! Stupid R can't deal with singles.
            vari=1
        } else {
            vari = names(xa)[i]
        }
        if (is.na(match(i, categ))) {
            plot(xa[ord,i], aceobj$tx[ord,vari], pch=16,type="c", 
                xlab=xlab, ylab=paste("acet(",xlab,")",sep=""), main="", col=1, 
                xlim=c(min(aceobj$x[vari,]), max(aceobj$x[vari,])),
                ylim=c(min(aceobj$tx[,vari]), max(aceobj$tx[,vari])))
            points(xa[,i], aceobj$tx[,vari],  pch=16, col=col)
        } else {
            plot(as.factor(xa[,i]), aceobj$tx[,i],xlab=xlab, 
                ylab=paste("acet(",xlab,")"),  pch=16)
            #points(xa[,i], aceobj$tx[,vari], pch=16)
        }
        mtext(outer=TRUE, titletext, line=line)
    }
    
    #Predicted ty   
    sumtx = apply(aceobj$tx, 1, sum)
    plot(sumtx, aceobj$ty, xlab="sum of acet(x's)", ylab=ylab, pch=16, col=col)
    lines( c(min(c(sumtx, aceobj$ty)), max(c(sumtx, aceobj$ty))),  
            c(min(c(sumtx, aceobj$ty)), max(c(sumtx, aceobj$ty))), lty=3)
    title(paste("Rsq =", round(aceobj$rsq,3)))
    mtext(outer=TRUE, titletext)
    
    #histogram
    hist(sumtx - aceobj$ty, main="", xlab=paste("Error in acet(",main,")", sep=""))
}

my.glm.summary = function(glmm) {
    print(summary(glmm))
    print(anova(glmm, test="Chisq"))
    lmmod = lm(glmm$formula, data=glmm$model)
    print(summary(lmmod))
    return(lmmod)
}

PCA.contrib = function(p) {
    #p is a PCA object from prcomp
    #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
    rot = p$rotation  #variable rows x PC columns
    npc = ncol(rot)
    nvar = nrow(rot)
    sdev.mat = t(matrix(rep(p$sdev, nvar), npc, nvar))     #sdev is length number of PC's
    
    var.coord = rot * sdev.mat
    var.cos2 = var.coord^2
    comp.cos2 = apply(var.cos2, 2, sum)
    contrib <- function(var.cos2, comp.cos2){var.cos2/comp.cos2}
    var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2)) #Contrib of vars to each PC, PC column sums to 1.
    
    varfrac = summary(p)$importance["Proportion of Variance", ] #Rel importance of each PC
    contrib.varfrac = var.contrib * t(matrix(rep(varfrac, nvar), npc, nvar)) #Sums to 1 over all PC's and vars
    
    
    print(round(var.contrib,3)*100)
    print(round(contrib.varfrac,3)*100)
    print(summary(p)$import)
    
    return(list(var.contrib=var.contrib, contrib.varfrac=contrib.varfrac, 
        varfrac.sum=apply(contrib.varfrac, 2, )))
}


my.R2 = function(y, pred) {
 return(1 - sum((pred - y)^2) / sum((y - mean(y))^2)  )
}

my.adjusted.R2 = function(y, pred, p) {
    #p = 0 gives sames as regular R^2
     n = length(y)
     return(1 - (1/(n-p-1)*sum((pred - y)^2)) / (1/(n-1)*sum((y - mean(y))^2)) )
}

lm_predict <- function (lmObject, newdata, diag = TRUE) {
  ## From https://stackoverflow.com/questions/39337862/linear-model-with-lm-how-to-get-prediction-variance-of-sum-of-predicted-value/39338512#39338512
  ## input checking
  if (!inherits(lmObject, "lm")) stop("'lmObject' is not a valid 'lm' object!")
  ## extract "terms" object from the fitted model, but delete response variable
  tm <- delete.response(terms(lmObject))      
  ## linear predictor matrix
  Xp <- model.matrix(tm, newdata)
  ## predicted values by direct matrix-vector multiplication
  pred <- c(Xp %*% coef(lmObject))
  ## efficiently form the complete variance-covariance matrix
  QR <- lmObject$qr   ## qr object of fitted model
  piv <- QR$pivot     ## pivoting index
  r <- QR$rank        ## model rank / numeric rank
  if (is.unsorted(piv)) {
    ## pivoting has been done
    B <- forwardsolve(t(QR$qr), t(Xp[, piv]), r)
    } else {
    ## no pivoting is done
    B <- forwardsolve(t(QR$qr), t(Xp), r)
    }
  ## residual variance
  sig2 <- c(crossprod(residuals(lmObject))) / df.residual(lmObject)
  if (diag) {
    ## return point-wise prediction variance
    VCOV <- colSums(B ^ 2) * sig2
    } else {
    ## return full variance-covariance matrix of predicted values
    VCOV <- crossprod(B) * sig2
    }
  list(fit = pred, var.fit = VCOV, df = lmObject$df.residual, residual.var = sig2)
  }
  
  
