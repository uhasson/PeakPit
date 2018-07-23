## June 2917. Uri Hasson, uri@hasson.org.
## This code derives visibility graphs for temperature in a given year,derives a distribution and compare top vs. bottom distributions
## The core visibility graph algorithm created by Thomas Jagger (function 'get visibility').
## Looping over one year take a few min minutes. Input file "air.2m.gauss.2015.nc" included here for completion

library(ncdf4)
library(moments)
overallcors <- c(); # stores correlations between DIST and marginal quantities

source("getVisibilityFunction.r")

allyears <- c(2015)
for (i in allyears) {
year <-  i
print(year);

inputfile=paste("air.2m.gauss.",year,".nc", sep="")
outfile=paste("results.",year,".topdown.txt", sep="")

ncin <- nc_open(inputfile)
print(ncin)

# some diagnostics on longitude and lattitude for sanity
lonvals <- ncvar_get(ncin, "lon")
nlon <- dim(lonvals)
head(lonvals)

latvals <- ncvar_get(ncin, "lat", verbose = F)
nlat <- dim(latvals)
head(latvals)

# some more diagnostics if need (uncomment next two lines)
# attributes(ncin$var)$names
# ncatt_get(ncin, attributes(ncin$var)$names[1])

# all relevant data here
ncatt_get(ncin, attributes(ncin$var)$names[2])

# get the core air temperature data
airtmpmat <- ncvar_get(ncin, attributes(ncin$var)$names[2])

dim(airtmpmat); # diagnostic: should be lon 192, lat94 time365

results <- c(); # stores the result for given year.

# the following 2 loops go over each long/lat combination, get the airtemperature time series and analyze it
for (lon in 1:192) {
     for (lat in 1:94) {
     # print(c(lon, lat))
    
# lon=104; lat=10; short circuit example to instantiate specific loation
testTS <- airtmpmat[lon,lat,]; # great example of badly named variable
revTestTS <- -1 * testTS; # this is the inverse time series

net1 <- get.visibility(testTS) ; # get the visibility of original
net2 <- get.visibility(revTestTS); # get visibility for inverse
# get visibility function by J elsner and THomas Jagger  http://myweb.fsu.edu/jelsner/Book/Chap10/Chapter10.html

updist <- net1$pk
downdist <- net2$pk

# need to generate two degree distributions with the same max bin.
# this is essential for comparign the two distributions
if (max(updist$k) > max(downdist$k)) {overmax <- max(updist$k)} else {overmax <- max(downdist$k)}
# prepare two empty historgrams matched in max bin
uptemp <- cbind(2:overmax, rep(0, overmax-1))
downtemp <- cbind(2:overmax, rep(0, overmax-1))

#  fill in empty histograms with as many degree slots as max of up or down
for (i in 1:length(updist$k)) {
  currentK <- updist$k[i]
  currentDeg <- updist$degree[i] 
  uptemp[(which(uptemp[,1]==currentK)),2] <- currentDeg
  }
  
for (i in 1:length(downdist$k)) {
  currentK <- downdist$k[i]
  currentDeg <- downdist$degree[i] 
  downtemp[(which(downtemp[,1]==currentK)),2] <- currentDeg
  }

updistfinal <- data.frame(uptemp)
downdistfinal <- data.frame(downtemp)
colnames(updistfinal) <- c("k", "degree")
colnames(downdistfinal) <- c("k", "degree")

# now that we have updistfinal and downdistfinal we can get a 'distance' value
currentDist <-dist(rbind(updistfinal$degree, downdistfinal$degree), method = "manhattan")
currentRes <- cbind(lon, lat, currentDist)
results <- rbind(results, currentRes)
 }
} ; # end loop over all long/lat. at this point #results conains all results.

write.table(results, file=outfile)

## Still within year loop; calculate correlations
avgTmp1 <- rowMeans(airtmpmat, dim=2)
avgTmp <- apply(airtmpmat, c(1,2), max)
sdTmp <- apply(airtmpmat, c(1,2), sd)
skewnessTmp <- apply(airtmpmat, c(1,2), skewness)
kurtosisTmp <- apply(airtmpmat, c(1,2), kurtosis)

avgtmpres <- c()
for (i in 1:192) {
     for (j in 1:94) {
     curval <- avgTmp[i, j]
     curres <- c(lon[i], lat[j], curval)
     avgtmpres <- rbind(avgtmpres, curres)
    }
}

sdtmpres <- c()
for (i in 1:192) {
     for (j in 1:94) {
     curval <- sdTmp[i, j]
     curres <- c(lon[i], lat[j], curval)
     sdtmpres <- rbind(sdtmpres, curres)
    }
}

skewtmpres <- c()
for (i in 1:192) {
     for (j in 1:94) {
     curval <- skewnessTmp[i, j]
     curres <- c(lon[i], lat[j], curval)
     skewtmpres <- rbind(skewtmpres, curres)
    }
}

kurttmpres <- c()
for (i in 1:192) {
     for (j in 1:94) {
     curval <- kurtosisTmp[i, j]
     curres <- c(lon[i], lat[j], curval)
     kurttmpres <- rbind(kurttmpres, curres)
    }
}


# bind the current cors to all cors
tmpcor <- cor.test(results[,3], avgtmpres[,3])
sdcor <- cor.test(results[,3], sdtmpres[,3])
skewcor <- cor.test(results[,3], skewtmpres[,3])
kurtcor <- cor.test(results[,3], kurttmpres[,3])
allCorrelRes <- c(year, tmpcor$estimate, sdcor$estimate,skewcor$estimate, kurtcor$estimate)
overallcors <- rbind(overallcors, allCorrelRes)
print(overallcors); #not saved to file in this version

} ; # end year loop
date()




   
