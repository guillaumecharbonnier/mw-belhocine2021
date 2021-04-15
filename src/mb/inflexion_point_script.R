m= read.delim("Peaks_Broad.txt" ,sep="," , header=TRUE)
m = m[order(m$Size_Peak),]
m$ord=1:nrow(m)
x=m$ord
y=m$Size_Peak
lo <- loess(y~x)
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
out = predict(lo,xl)
infl <- c(FALSE, diff( diff(diff(out))>0)!=0)
xl[infl]
out[infl]
