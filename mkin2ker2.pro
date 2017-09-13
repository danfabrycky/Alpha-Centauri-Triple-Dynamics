pro mkin2ker2

readcol, '../mcorbitker2.txt', xp, yp, zp, vxp, vyp, vzp, maarr, mbarr, mcarr, format='d,d,d,d,d,d,d,d,d'

;ma=1.1055d0
;mb=0.9373d0   ; 
parallax=742.12d0
n=100
incs =dblarr(n)

for i=0,n-1 do begin

ma = maarr(i)
mb = mbarr(i)
mc = mcarr(i)
a = keptostate(17.66*1000/parallax, 0.524, 79.32*!dpi/180, 204.75*!dpi/180, 232.3*!dpi/180, getf((2016-1955.66)*(2*!dpi)/79.91, 0.524, 0), m = (4*!dpi^2)*(ma+mb))
print, a
xs=a(0)
ys=a(1)
zs=a(2)
vxs=a(3)
vys=a(4)
vzs=a(5)

; rotate state into X,Y,Z plane
dec = 219.91753275d0
ra = -60.83712790d0
xint = -sin(dec)*xs - cos(dec)*zs
yint = ys
zint = -cos(dec)*xs + sin(dec)*zs
x = - (cos(ra)*xint - sin(ra)*yint)
y = sin(ra)*xint + cos(ra)*yint
z = zint
vxint = -sin(dec)*vxs - cos(dec)*vzs
vyint = vys
vzint = -cos(dec)*vxs + sin(dec)*vzs
vx = - (cos(ra)*vxint - sin(ra)*vyint)
vy = sin(ra)*vxint + cos(ra)*vyint
vz = vzint

; mc=0.1221
mtot = ma+mb+mc

; check the mutual inclination
lxb = ys*vzs - zs*vys
lyb = zs*vxs - xs*vzs
lzb = xs*vys - ys*vxs
lbtot = sqrt(lxb^2+lyb^2+lzb^2)

lxc = yp(i)*vzp(i) - zp(i)*vyp(i)
lyc = zp(i)*vxp(i) - xp(i)*vzp(i)
lzc = xp(i)*vyp(i) - yp(i)*vxp(i)
lctot = sqrt(lxc^2+lyc^2+lzc^2)

cosi = (lxb*lxc + lyb*lyc + lzb*lzc)/lbtot/lctot
incs(i) = acos(cosi)*180/!pi

; make it primary-centric: 
xc = xp(i) + mb/(ma+mb) * x
yc = yp(i) + mb/(ma+mb) * y
zc = zp(i) + mb/(ma+mb) * z
vxc = vxp(i) + mb/(ma+mb) * vx
vyc = vyp(i) + mb/(ma+mb) * vy
vzc = vzp(i) + mb/(ma+mb) * vz

   print, "writing ", 'insker2/integration2.'+string(i,format='(i02)')+'.in'
   openw, outf, 'insker2/integration2.'+string(i,format='(i02)')+'.in', /get_lun
   printf, outf, 3
   printf, outf, 6d9
   printf, outf, format='(7(e22.15, " "))', ma, 0,0,0,0,0,0
   printf, outf, format='(7(e22.15, " "))', mb, x, y, z, vx, vy, vz
   printf, outf, format='(7(e22.15, " "))', mc, xc, yc, zc, vxc, vyc, vzc
   free_lun, outf

endfor

histogram_ez, incs, binsize=1, xtitle='mutual inclination (deg)', ytitle='N', charsize=1.5
print, "inc = ", mean(incs), "+/-", sigma(incs)

end

