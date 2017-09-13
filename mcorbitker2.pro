pro mcorbitker2, seed=seed, ps=ps, gal=gal

; using Kervella's published orbit. 
; IDL> mcorbitker2, /gal, /ps

if not keyword_set(seed) then seed=215
  
n = 100000l ; number of samples we want. 

;PROXIMA CENTARURI
alpha = 217.44894750d0  ; checked, for published K17.
delta = -62.68135207d0
parallax = 768.77d0
mualpha = -3773.8d0
mudelta = 770.5d0
vr = -22.204d0   ; checked, for published K17.
mproxima = 0.1221d0

dalpha = 3.6d-7
ddelta = 4.2d-7
dparallax = 0.37d0
dmualpha = 0.4d0
dmudelta = 2.0d0
dvr = 0.032d0
dmproxima = 0.0022d0

alphasample = dblarr(n)
deltasample = dblarr(n)
parallaxsample = dblarr(n)
mualphasample = dblarr(n)
mudeltasample = dblarr(n)
vrsample = dblarr(n)
mproximasample = dblarr(n)

n_unbound = 0l
n_bound = 0l
n_small_value = 0l
n_dynamics = 0l

xarr=dblarr(n)
yarr=dblarr(n)
zarr=dblarr(n)
vxarr=dblarr(n)
vyarr=dblarr(n)
vzarr=dblarr(n)
x1arr=dblarr(n)
y1arr=dblarr(n)
z1arr=dblarr(n)
vx1arr=dblarr(n)
vy1arr=dblarr(n)
vz1arr=dblarr(n)

aarrk=dblarr(n)
earrk=dblarr(n)
iarrk=dblarr(n)
loarrk=dblarr(n)
boarrk=dblarr(n)
farrk=dblarr(n)

aarr=dblarr(n)
earr=dblarr(n)
iarr=dblarr(n)
loarr=dblarr(n)
boarr=dblarr(n)
farr=dblarr(n)

;ALPHA CENTAURI

alpha2 = 219.917528d0  ; checked, for published K17.
delta2 = -60.837128d0
parallax2 = 747.17d0
mualpha2 = -3619.9d0
mudelta2 = 693.8d0
vr2 = -22.332d0  ; checked, for published K17.
ma = 1.1055d0
mb = 0.9373d0
pab =79.91

dalpha2 = 1.4d-5 
ddelta2 = 1.0d-5
dparallax2 = 0.61d0
dmualpha2 = 3.9d0
dmudelta2 = 3.9d0
dvr2 = 0.005d0
dma = 0.0039d0
dmb = 0.0033d0

alpha2sample = dblarr(n)
delta2sample = dblarr(n)
parallax2sample = dblarr(n)
mualpha2sample = dblarr(n)
mudelta2sample = dblarr(n)
vr2sample = dblarr(n)
masample = dblarr(n)
mbsample = dblarr(n)

openw, out, 'mcorbitker2.txt', /get_lun
openw, out2, 'mcorbitunker2.txt', /get_lun
for i=0l,n-1 do begin 

   alphasample(i) = alpha + dalpha*randomn(seed)
   deltasample(i) = delta + ddelta*randomn(seed)
   parallaxsample(i) = parallax + dparallax*randomn(seed)
   mualphasample(i) = mualpha + dmualpha*randomn(seed)
   mudeltasample(i) = mudelta + dmudelta*randomn(seed)
   vrsample(i) = vr + dvr*randomn(seed)
   mproximasample(i) = mproxima + dmproxima*randomn(seed)

   alpha2sample(i) = alpha2 + dalpha2*randomn(seed)
   delta2sample(i) = delta2 + ddelta2*randomn(seed)
   parallax2sample(i) = parallax2 + dparallax2*randomn(seed)
   mualpha2sample(i) = mualpha2 + dmualpha2*randomn(seed)
   mudelta2sample(i) = mudelta2 + dmudelta2*randomn(seed)
   vr2sample(i) = vr2 + dvr2*randomn(seed)
   masample(i) = ma + dma*randomn(seed)
   mbsample(i) = mb + dmb*randomn(seed)

 ; ...
  ; then compute X,Y,Z,VX,VY,VZ to start dynamical simulation.

valpha = (mualphasample(i)/1000)*(1000/parallaxsample(i))
;*(cos(deltasample(i)*!dpi/180))  -- not needed. 
;print, valpha

vdelta = (mudeltasample(i)/1000)*(1000/parallaxsample(i))
;print, vdelta

vrsample(i) *= 6.6846d-9*3.15576d7  ; converting km/s -> to AU/d
;print, vr

valpha2 = (mualpha2sample(i)/1000)*(1000/parallax2sample(i))
;*(cos(delta2sample(i)*!dpi/180)) -- not needed.
;print, valpha2

vdelta2 = (mudelta2sample(i)/1000)*(1000/parallax2sample(i))
;print, vdelta2

vr2sample(i) *= 6.6846d-9*3.15576d7  ; converting km/s -> to AU/d
;print, vr2


;Difference Between X , Y , & Z for Proxima Minus Alpha Centauri


rproxima = 1/(parallaxsample(i)/1000)*206264.806d0

xproxima = (rproxima)*(cos(deltasample(i)*!dpi/180))*(cos(alphasample(i)*!dpi/180))
yproxima = (rproxima)*(cos(deltasample(i)*!dpi/180))*(sin(alphasample(i)*!dpi/180))
zproxima = (rproxima)*(sin(deltasample(i)*!dpi/180))

ralphacentauri = 1/(parallax2sample(i)/1000)*206264.806d0

xalphacentauri = (ralphacentauri)*(cos(delta2sample(i)*!dpi/180))*(cos(alpha2sample(i)*!dpi/180))
yalphacentauri = (ralphacentauri)*(cos(delta2sample(i)*!dpi/180))*(sin(alpha2sample(i)*!dpi/180))
zalphacentauri = (ralphacentauri)*(sin(delta2sample(i)*!dpi/180))

xdiff = (xproxima) - (xalphacentauri)
print, xdiff

ydiff = (yproxima) - (yalphacentauri)
print, ydiff

zdiff = (zproxima) - (zalphacentauri)
print, zdiff

print, sqrt(xdiff^2 + ydiff^2 + zdiff^2)


;Differences in Velocity


valphadiff = (valpha) - (valpha2)
print, "valphadiff", valphadiff, valpha, valpha2

vdeltadiff = (vdelta) - (vdelta2)
print, "vdeltadiff", vdeltadiff

vrdiff = (vrsample(i)) - (vr2sample(i))
print, "vrdiff", vrdiff


; VX, VY, & VZ for Proxima


vx = (-(sin(deltasample(i)*!dpi/180))*(cos(alphasample(i)*!dpi/180)))*(vdelta) + (-(sin(alphasample(i)*!dpi/180)))*(valpha) + ((cos(deltasample(i)*!dpi/180))*(cos(alphasample(i)*!dpi/180)))*(vrsample(i))

vy = (-((sin(deltasample(i)*!dpi/180)))*(sin(alphasample(i)*!dpi/180)))*(vdelta) + (cos(alphasample(i)*!dpi/180))*(valpha) + ((cos(deltasample(i)*!dpi/180))*(sin(alphasample(i)*!dpi/180)))*(vrsample(i))

vz = (cos(deltasample(i)*!dpi/180))*(vdelta) + 0*valpha + (sin(deltasample(i)*!dpi/180))*(vrsample(i))


; VX, VY, & VZ for Alpha Centauri


vx2 = (-(sin(delta2sample(i)*!dpi/180))*(cos(alpha2sample(i)*!dpi/180)))*(vdelta2) + (-(sin(alpha2sample(i)*!dpi/180)))*(valpha2) + ((cos(delta2sample(i)*!dpi/180))*(cos(alpha2sample(i)*!dpi/180)))*(vr2sample(i))

vy2 = (-((sin(delta2sample(i)*!dpi/180)))*(sin(alpha2sample(i)*!dpi/180)))*(vdelta2) + (cos(alpha2sample(i)*!dpi/180))*(valpha2) + ((cos(delta2sample(i)*!dpi/180))*(sin(alpha2sample(i)*!dpi/180)))*(vr2sample(i))

vz2 = (cos(delta2sample(i)*!dpi/180))*(vdelta2) + 0*valpha2 + (sin(delta2sample(i)*!dpi/180))*(vr2sample(i))



;Checking Velocities

;print, "vprox" , valpha, vdelta, vrsample(i)
;print, "valph" , valpha2, vdelta2, vr2sample(i)


;Difference Between VX, VY, & VZ for Proxima Minus Alpha Centauri

vxdiff = (vx) - (vx2) 
;print, vxdiff;, vx, vx2

vydiff = (vy) - (vy2)
;print, vydiff;, vy, vy2

vzdiff = (vz) - (vz2)
;print,vzdiff;, vz, vz2

vtot =  sqrt(vxdiff^2 + vydiff^2 + vzdiff^2)
;print, vtot

if keyword_set(gal) then begin  ; casper.berkeley.edu/astrobaki/index.php/Coordinates
  print, "my coords proxima: ", xproxima/206264.806d0, yproxima/206264.806d0, zproxima/206264.806d0
  xgal = -0.054876d0 * xproxima - 0.873437d0 * yproxima - 0.483835d0 * zproxima 
  ygal =  0.494109d0 * xproxima - 0.444830d0 * yproxima + 0.746982d0 * zproxima 
  zgal = -0.867666d0 * xproxima - 0.198076d0 * yproxima + 0.455984d0 * zproxima 
;  print, "gal coords proxima: ", xgal/206264.806d0, ygal/206264.806d0, zgal/206264.806d0

  print, "my coords alphacentauri: ", xalphacentauri/206264.806d0, yalphacentauri/206264.806d0, zalphacentauri/206264.806d0
  xgala = -0.054876d0 * xalphacentauri - 0.873437d0 * yalphacentauri - 0.483835d0 * zalphacentauri 
  ygala =  0.494109d0 * xalphacentauri - 0.444830d0 * yalphacentauri + 0.746982d0 * zalphacentauri 
  zgala = -0.867666d0 * xalphacentauri - 0.198076d0 * yalphacentauri + 0.455984d0 * zalphacentauri 
;  print, "gal coords alphacentauri: ", xgala/206264.806d0, ygala/206264.806d0, zgala/206264.806d0

  vxgal = -0.054876d0 * vx - 0.873437d0 * vy - 0.483835d0 * vz 
  vygal =  0.494109d0 * vx - 0.444830d0 * vy + 0.746982d0 * vz 
  vzgal = -0.867666d0 * vx - 0.198076d0 * vy + 0.455984d0 * vz   
  print, "gal vel prox: ", [vxgal, vygal, vzgal]*4.7406

  vx2gal = -0.054876d0 * vx2 - 0.873437d0 * vy2 - 0.483835d0 * vz2 
  vy2gal =  0.494109d0 * vx2 - 0.444830d0 * vy2 + 0.746982d0 * vz2 
  vz2gal = -0.867666d0 * vx2 - 0.198076d0 * vy2 + 0.455984d0 * vz2   
  print, "gal vel alphacentari: ", [vx2gal, vy2gal, vz2gal]*4.7406

  xgal = -0.054876d0 * xdiff - 0.873437d0 * ydiff - 0.483835d0 * zdiff 
  ygal =  0.494109d0 * xdiff - 0.444830d0 * ydiff + 0.746982d0 * zdiff 
  zgal = -0.867666d0 * xdiff - 0.198076d0 * ydiff + 0.455984d0 * zdiff 
  vxgal = -0.054876d0 * vxdiff - 0.873437d0 * vydiff - 0.483835d0 * vzdiff 
  vygal =  0.494109d0 * vxdiff - 0.444830d0 * vydiff + 0.746982d0 * vzdiff 
  vzgal = -0.867666d0 * vxdiff - 0.198076d0 * vydiff + 0.455984d0 * vzdiff 
;  xdiff = xgal
;  ydiff = ygal
;  zdiff = zgal
;  vxdiff = vxgal
;  vydiff = vygal
;  vzdiff = vzgal
  print, "gal coords prox-alpha: ", [xgal, ygal, zgal]/206264.806d0
  print, "gal vels prox-alpha: ", [vxgal, vygal, vzgal]*4.7406

   print, xdiff, ydiff, zdiff, vxdiff, vydiff, vzdiff
   xdiff1 = (-0.05622d0 + 0.00089d0*randomn(seed))*206264.806d0
   ydiff1 = (-0.00198d0 + 0.00089d0*randomn(seed))*206264.806d0
   zdiff1 = (-0.02785d0 + 0.00002d0*randomn(seed))*206264.806d0
   vxdiff1 = (-0.099d0 + 0.038d0*randomn(seed))/4.7406
   vydiff1 = (0.173d0 + 0.027d0*randomn(seed))/4.7406
   vzdiff1 = (0.187d0 + 0.016d0*randomn(seed))/4.7406
   print, xdiff1, ydiff1, zdiff1, vxdiff1, vydiff1, vzdiff1
   xarr(i)=xdiff
   yarr(i)=ydiff
   zarr(i)=zdiff
   vxarr(i)=vxdiff
   vyarr(i)=vydiff
   vzarr(i)=vzdiff
   x1arr(i)=xdiff1
   y1arr(i)=ydiff1
   z1arr(i)=zdiff1
   vx1arr(i)=vxdiff1
   vy1arr(i)=vydiff1
   vz1arr(i)=vzdiff1

   ; let's see what happens when we inject the Kervella answers: 
   if 1 eq 0 then begin
      xdiff = xdiff1
      ydiff = ydiff1
      zdiff = zdiff1
      vxdiff = vxdiff1
      vydiff = vydiff1
      vzdiff = vzdiff1
   endif
endif 

;STATETOKEP

m = (2*!dpi)^2*(masample(i)+mbsample(i)+mproximasample(i))
elom = statetokep(xdiff, ydiff, zdiff, vxdiff, vydiff, vzdiff, m=m)
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

; elom = statetokep(ydiff, -xdiff, zdiff, vydiff, -vxdiff, vzdiff, m=m)
; elom = statetokep(ydiff, zdiff, xdiff, vydiff, vzdiff, vxdiff, m=m)
elom = statetokep(ydiff, zdiff, xdiff, -vydiff, vzdiff, vxdiff, m=m)

a = elom(0)
e = elom(1)
inc = elom(2)
lo = elom(3)
bo = elom(4)
f = elom(5)

aarrk(i)=a
earrk(i)=e
iarrk(i)=inc
loarrk(i)=lo
boarrk(i)=bo
farrk(i)=f

;print, "a", a
;print, "e", e
;print, "i", inc
;print, "lo", lo
;print, "bo", bo
;print, "f", f

;AOUT & EOUT

value = a*sqrt(1-e^2)
;print, value

tkozai = (2/(3*!pi*mproximasample(i)*pab)*value^3)/1e9

;print, "i=", i, "   gives ", (2/(3*!pi*mproxima*79)*value^3)/1e9, "gyr*kozaitimescale"

print, ".............................................................."

if(value le 2000 and value gt 0) then n_small_value += 1

if(value gt 0) then begin
printf, out, xdiff, ydiff, zdiff, vxdiff, vydiff, vzdiff, masample(i), mbsample(i), mproximasample(i), format='(9e16.8)'

n_bound += 1
endif

if(a lt 0) then begin
 n_unbound += 1
printf, out2, xdiff, ydiff, zdiff, vxdiff, vydiff, vzdiff, masample(i), mbsample(i), mproximasample(i), format='(9e16.8)'
endif

if(a gt 0 and tkozai lt 6) then n_dynamics += 1

endfor

free_lun, out 
free_lun, out2


if keyword_set(ps) then begin
   set_plot, 'ps'
   device, filename='outfile.eps', /encapsulated, /inches, xsize=5, ysize=4
endif else begin
   cs=2.5
   window, xsize=900, ysize=600
endelse

!p.multi=[0,3,2]

histogram_ez, aarr/1000, binsize=0.1, xtitle='a [kAU]', ytitle='N', xrange=[6,12], /xstyle, charsize=cs
histogram_ez, earr, binsize=0.02, xtitle='e', ytitle='N', charsize=cs, xrange=[0.1,0.9], /xstyle
histogram_ez, pushpi(iarr*180/!dpi,/deg), binsize=2, xtitle='i [deg]', ytitle='N', charsize=cs, xrange=[50,110], /xstyle
histogram_ez, (pushpi(loarr*180/!dpi,/deg)+360) mod 360, binsize=2, xtitle='omega [deg]', ytitle='N', charsize=cs, xrange=[160,230], /xstyle
histogram_ez, pushpi(boarr*180/!dpi,/deg), binsize=0.5, xtitle='Omega [deg]', ytitle='N', charsize=cs
histogram_ez, (pushpi(farr*180/!dpi,/deg)+360) mod 360, binsize=2, xtitle='f [deg]', ytitle='N', charsize=cs, xrange=[150,220], /xstyle

if keyword_set(ps) then begin
    device, /close_file
    set_plot, 'x'
endif


if keyword_set(gal) then begin  

if keyword_set(ps) then begin
   set_plot, 'ps'
   device, filename='outfilegal.eps', /encapsulated, /inches, xsize=5, ysize=4
endif else begin
   cs=2.5
   window, xsize=900, ysize=600, 3
endelse

!p.multi=[0,3,2]

histogram_ez, aarrk/1000, binsize=0.1, xtitle='a [kAU]', ytitle='N', xrange=[6,12], /xstyle, charsize=cs
histogram_ez, earrk, binsize=0.02, xtitle='e', ytitle='N', charsize=cs, xrange=[0.1,0.9], /xstyle
histogram_ez, pushpi(iarrk*180/!dpi,/deg), binsize=1, xtitle='i [deg]', ytitle='N', charsize=cs, /xstyle
histogram_ez, (pushpi(loarrk*180/!dpi,/deg)+540) mod 360 + 180, binsize=1, xtitle='omega [deg]', ytitle='N', charsize=cs, /xstyle
histogram_ez, (pushpi(boarrk*180/!dpi,/deg) + 360) mod 360, binsize=1, xtitle='Omega [deg]', ytitle='N', charsize=cs
histogram_ez, (pushpi(farrk*180/!dpi,/deg)+360) mod 360, binsize=2, xtitle='f [deg]', ytitle='N', charsize=cs, /xstyle

if keyword_set(ps) then begin
    device, /close_file
    set_plot, 'x'
endif

endif

nxyz=1000
xarr=dblarr(nxyz)
yarr=dblarr(nxyz)
zarr=dblarr(nxyz)

window, 1, xsize=1000, ysize=600
!p.multi=[0,3,2]

; plot the orbit and the current position. 
if 1 eq 1 then begin
; first, compare xyzuvw between observations and derived from orbital
; elements. 
	a = 8.7d3	
	e = 0.5d0	
	inc=107.6*!dpi/180
	lo=72.3*!dpi/180
	bo=126*!dpi/180
	t0=284d3
	p=547d3
	bigm=2*!dpi*(1-t0/p)
	f=getf(bigm+lo+bo, e, lo+bo) 

	m = (2*!dpi)^2*(1.1055+0.9373+0.1221)
	xyzs = keptostate(a, e, inc, lo, bo, f, m=m)	

	print, "paper's xyz from aei: ", xyzs(0:2)/206265, xyzs(3:5)*4.7406
	cafe_pause

	if keyword_set(gal) then begin  

	for j=0, nxyz-1 do begin
  	  xyzs = keptostate(aarr(i-1), earr(i-1), iarr(i-1), loarr(i-1), boarr(i-1), 2*!pi*j/nxyz, m=m)
  	  xarr(j) = xyzs(0)
  	  yarr(j) = xyzs(1)
 	  zarr(j) = xyzs(2)
	endfor
endif

plot, xarr/1000, yarr/1000, /isotropic, xtitle='X', ytitle='Y', charsize=2
oplot, [xgal]/1000, [ygal]/1000, psym=7
oplot, xgal/1000+[0,vxgal*100], ygal/1000+[0,vygal*100]

plot, yarr/1000, zarr/1000, /isotropic, xtitle='Y', ytitle='Z', charsize=2
oplot, [ygal]/1000, [zgal]/1000, psym=7
oplot, ygal/1000+[0,vygal*100], zgal/1000+[0,vzgal*100]

plot, xarr/1000, zarr/1000, /isotropic, xtitle= 'X', ytitle='Z', charsize=2
oplot, [xgal]/1000, [zgal]/1000, psym=7
print, xgal, ygal, zgal
oplot, xgal/1000+[0,vxgal*100], zgal/1000+[0,vzgal*100]

print, "orbit parameters = "
print, mean(aarr), sigma(aarr)
print, mean(earr), sigma(earr)

niaark = pushpi(iarrk*180/!dpi,/deg)
print, mean(niaark), sigma(niaark)
nloarrk = (pushpi(loarrk*180/!dpi,/deg)+540) mod 360 + 180
print, mean(nloarrk), sigma(nloarrk)
nboarrk = (pushpi(boarrk*180/!dpi,/deg) + 360) mod 360
print, mean(nboarrk), sigma(nboarrk)
nfarrk = (pushpi(farrk*180/!dpi,/deg)+360) mod 360
print, mean(nfarrk), sigma(nfarrk)

print, n_small_value, n_bound, n_unbound, n_dynamics

endif

end
