pro covtripgal

; This makes the covariance matrix of the results file. 
 
readcol, 'mcorbitker2gal.txt', x, y, z, vx, vy, vz, masample, mbsample, mproximasample, format='d,d,d,d,d,d,d,d,d'

n = n_elements(x) 
aarr=dblarr(n)
earr=dblarr(n)
iarr=dblarr(n)
loarr=dblarr(n)
boarr=dblarr(n)
farr=dblarr(n)

for i=0l,n-1 do begin
   
   m = (2*!dpi)^2*(masample(i)+mbsample(i)+mproximasample(i))
   elom = statetokep(x(i), y(i), z(i), vx(i), vy(i), vz(i), m=m)
   a = elom(0)
   e = elom(1)
   inc = elom(2)
   lo = elom(3)
   bo = elom(4)
   f = elom(5)
   
   aarr(i)=a
   earr(i)=e
   iarr(i)=inc
   loarr(i)=lo
   boarr(i)=bo
   farr(i)=f

endfor 

print, mean(x), sigma(x) 
print, mean(y), sigma(y) 
print, mean(z), sigma(z) 
print, mean(vx), sigma(vx) 
print, mean(vy), sigma(vy) 
print, mean(vz), sigma(vz) 

print, '--------'
print, VARIANCE(x), CORRELATE(x, y, /COVARIANCE), CORRELATE(x, z, /COVARIANCE), CORRELATE(x, vx, /COVARIANCE), CORRELATE(x, vy, /COVARIANCE), CORRELATE(x, vz, /COVARIANCE)
print, CORRELATE(x, y, /COVARIANCE), VARIANCE(y), CORRELATE(y, z, /COVARIANCE), CORRELATE(y, vx, /COVARIANCE), CORRELATE(y, vy, /COVARIANCE), CORRELATE(y, vz, /COVARIANCE)
print, CORRELATE(x, z, /COVARIANCE), CORRELATE(y, z, /COVARIANCE), VARIANCE(z), CORRELATE(z, vx, /COVARIANCE), CORRELATE(z, vy, /COVARIANCE), CORRELATE(z, vz, /COVARIANCE)
print, CORRELATE(x, vx, /COVARIANCE), CORRELATE(y, vx, /COVARIANCE), CORRELATE(z, vx, /COVARIANCE), VARIANCE(vx), CORRELATE(vx, vy, /COVARIANCE), CORRELATE(vx, vz, /COVARIANCE)
print, CORRELATE(x, vy, /COVARIANCE), CORRELATE(y, vy, /COVARIANCE), CORRELATE(z, vy, /COVARIANCE), CORRELATE(vx, vy, /COVARIANCE), VARIANCE(vy), CORRELATE(vy, vz, /COVARIANCE)
print, CORRELATE(x, vz, /COVARIANCE), CORRELATE(y, vz, /COVARIANCE), CORRELATE(z, vz, /COVARIANCE), CORRELATE(vx, vz, /COVARIANCE), CORRELATE(vy, vz, /COVARIANCE), VARIANCE(vz)

print, '--------'
print, VARIANCE(x), CORRELATE(x, y), CORRELATE(x, z), CORRELATE(x, vx), CORRELATE(x, vy), CORRELATE(x, vz)
print, CORRELATE(x, y), VARIANCE(y), CORRELATE(y, z), CORRELATE(y, vx), CORRELATE(y, vy), CORRELATE(y, vz)
print, CORRELATE(x, z), CORRELATE(y, z), VARIANCE(z), CORRELATE(z, vx), CORRELATE(z, vy), CORRELATE(z, vz)
print, CORRELATE(x, vx), CORRELATE(y, vx), CORRELATE(z, vx), VARIANCE(vx), CORRELATE(vx, vy), CORRELATE(vx, vz)
print, CORRELATE(x, vy), CORRELATE(y, vy), CORRELATE(z, vy), CORRELATE(vx, vy), VARIANCE(vy), CORRELATE(vy, vz)
print, CORRELATE(x, vz), CORRELATE(y, vz), CORRELATE(z, vz), CORRELATE(vx, vz), CORRELATE(vy, vz), VARIANCE(vz)
print, '--------'
print, '$1$ & $', CORRELATE(x, y), '$&$', CORRELATE(x, z), '$&$', CORRELATE(x, vx), '$&$', CORRELATE(x, vy), '$&$', CORRELATE(x, vz), '$ \\'
print, '$', CORRELATE(x, y), '$, $1$ & $', CORRELATE(y, z), '$&$', CORRELATE(y, vx), '$&$', CORRELATE(y, vy), '$&$', CORRELATE(y, vz), '$\\'
print, '$', CORRELATE(x, z), '$&$', CORRELATE(y, z), '$, $1$ & $',  CORRELATE(z, vx), '$&$', CORRELATE(z, vy), '$&$', CORRELATE(z, vz), '$\\'
print, '$', CORRELATE(x, vx), '$&$', CORRELATE(y, vx), '$&$', CORRELATE(z, vx), '$, $1$ & $', CORRELATE(vx, vy), '$&$', CORRELATE(vx, vz), '$\\'
print, '$', CORRELATE(x, vy), '$&$', CORRELATE(y, vy), '$&$', CORRELATE(z, vy), '$&$', CORRELATE(vx, vy),'$, $1$ & $',  CORRELATE(vy, vz), '$ \\'
print, '$', CORRELATE(x, vz), '$&$', CORRELATE(y, vz), '$&$', CORRELATE(z, vz), '$&$', CORRELATE(vx, vz), '$&$', CORRELATE(vy, vz), '$, $1$ \\'

print, '--------'
print, mean(aarr), sigma(aarr) 
print, mean(earr), sigma(earr) 
print, mean(iarr), sigma(iarr) 
print, mean(loarr), sigma(loarr) 
print, mean(boarr), sigma(boarr) 
print, mean(farr), sigma(farr) 

cs=2.5
!p.multi=[0,3,2]

histogram_ez, aarr/1000, binsize=0.1, xtitle='a [kAU]', ytitle='N', xrange=[6,12], /xstyle, charsize=cs
histogram_ez, earr, binsize=0.02, xtitle='e', ytitle='N', charsize=cs, xrange=[0.1,0.9], /xstyle
histogram_ez, pushpi(iarr*180/!dpi,/deg), binsize=2, xtitle='i [deg]', ytitle='N', charsize=cs, xrange=[50,110], /xstyle
histogram_ez, (pushpi(loarr*180/!dpi,/deg)+360) mod 360, binsize=2, xtitle='omega [deg]', ytitle='N', charsize=cs, xrange=[160,230], /xstyle
histogram_ez, pushpi(boarr*180/!dpi,/deg), binsize=0.5, xtitle='Omega [deg]', ytitle='N', charsize=cs
histogram_ez, (pushpi(farr*180/!dpi,/deg)+360) mod 360, binsize=2, xtitle='f [deg]', ytitle='N', charsize=cs, xrange=[150,220], /xstyle

print, VARIANCE(aarr), CORRELATE(aarr, earr, /COVARIANCE), CORRELATE(aarr, iarr, /COVARIANCE), CORRELATE(aarr, loarr, /COVARIANCE), CORRELATE(aarr, boarr, /COVARIANCE), CORRELATE(aarr, farr, /COVARIANCE)
print, CORRELATE(aarr, earr, /COVARIANCE), VARIANCE(earr), CORRELATE(earr, iarr, /COVARIANCE), CORRELATE(earr, loarr, /COVARIANCE), CORRELATE(earr, boarr, /COVARIANCE), CORRELATE(earr, farr, /COVARIANCE)
print, CORRELATE(aarr, iarr, /COVARIANCE), CORRELATE(earr, iarr, /COVARIANCE), VARIANCE(iarr), CORRELATE(iarr, loarr, /COVARIANCE), CORRELATE(iarr, boarr, /COVARIANCE), CORRELATE(iarr, farr, /COVARIANCE)
print, CORRELATE(aarr, loarr, /COVARIANCE), CORRELATE(earr, loarr, /COVARIANCE), CORRELATE(iarr, loarr, /COVARIANCE), VARIANCE(loarr), CORRELATE(loarr, boarr, /COVARIANCE), CORRELATE(loarr, farr, /COVARIANCE)
print, CORRELATE(aarr, boarr, /COVARIANCE), CORRELATE(earr, boarr, /COVARIANCE), CORRELATE(iarr, boarr, /COVARIANCE), CORRELATE(loarr, boarr, /COVARIANCE), VARIANCE(boarr), CORRELATE(boarr, farr, /COVARIANCE)
print, CORRELATE(aarr, farr, /COVARIANCE), CORRELATE(earr, farr, /COVARIANCE), CORRELATE(iarr, farr, /COVARIANCE), CORRELATE(loarr, farr, /COVARIANCE), CORRELATE(boarr, farr, /COVARIANCE), VARIANCE(farr)
print, '--------'
print, VARIANCE(aarr), CORRELATE(aarr, earr), CORRELATE(aarr, iarr), CORRELATE(aarr, loarr), CORRELATE(aarr, boarr), CORRELATE(aarr, farr)
print, CORRELATE(aarr, earr), VARIANCE(earr), CORRELATE(earr, iarr), CORRELATE(earr, loarr), CORRELATE(earr, boarr), CORRELATE(earr, farr)
print, CORRELATE(aarr, iarr), CORRELATE(earr, iarr), VARIANCE(iarr), CORRELATE(iarr, loarr), CORRELATE(iarr, boarr), CORRELATE(iarr, farr)
print, CORRELATE(aarr, loarr), CORRELATE(earr, loarr), CORRELATE(iarr, loarr), VARIANCE(loarr), CORRELATE(loarr, boarr), CORRELATE(loarr, farr)
print, CORRELATE(aarr, boarr), CORRELATE(earr, boarr), CORRELATE(iarr, boarr), CORRELATE(loarr, boarr), VARIANCE(boarr), CORRELATE(boarr, farr)
print, CORRELATE(aarr, farr), CORRELATE(earr, farr), CORRELATE(iarr, farr), CORRELATE(loarr, farr), CORRELATE(boarr, farr), VARIANCE(farr)
print, '--------'
print, '$1$ & $', CORRELATE(aarr, earr), '$&$', CORRELATE(aarr, iarr), '$&$', CORRELATE(aarr, loarr), '$&$', CORRELATE(aarr, boarr), '$&$', CORRELATE(aarr, farr), '$ \\'
print, '$', CORRELATE(aarr, earr), '$, $1$ & $', CORRELATE(earr, iarr), '$&$', CORRELATE(earr, loarr), '$&$', CORRELATE(earr, boarr), '$&$', CORRELATE(earr, farr), '$\\'
print, '$', CORRELATE(aarr, iarr), '$&$', CORRELATE(earr, iarr), '$, $1$ & $',  CORRELATE(iarr, loarr), '$&$', CORRELATE(iarr, boarr), '$&$', CORRELATE(iarr, farr), '$\\'
print, '$', CORRELATE(aarr, loarr), '$&$', CORRELATE(earr, loarr), '$&$', CORRELATE(iarr, loarr), '$, $1$ & $', CORRELATE(loarr, boarr), '$&$', CORRELATE(loarr, farr), '$\\'
print, '$', CORRELATE(aarr, boarr), '$&$', CORRELATE(earr, boarr), '$&$', CORRELATE(iarr, boarr), '$&$', CORRELATE(loarr, boarr),'$, $1$ & $',  CORRELATE(boarr, farr), '$ \\'
print, '$', CORRELATE(aarr, farr), '$&$', CORRELATE(earr, farr), '$&$', CORRELATE(iarr, farr), '$&$', CORRELATE(loarr, farr), '$&$', CORRELATE(boarr, farr), '$, $1$ \\'


end
