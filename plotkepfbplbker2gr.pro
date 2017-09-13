pro plotkepfbplbker2gr, ps=ps

;   This version plots both forwards and backwards in time. 
  red = [0,255,0,20,255]
  green = [0,0,220,20,255]
  blue = [0,0,0,255,255]

  red =   [0,255,  0,  0, 21, 42, 63, 84,105,126,147,168,189,210,232,253,255]
  green = [0,  0,205,  0, 21, 42, 63, 84,105,126,147,168,189,210,232,253,255]
  blue =  [0,  0,  0,253,232,210,189,168,147,126,105, 84, 63, 42, 21,  0, 255]
  
  if keyword_set(ps) then begin
     set_plot, 'ps'
     !p.font=0
     device, filename='kozstatsfbplb2gr.eps', /color, bits_per_pixel=8, /encapsulated, /inches, xsize=7, ysize=5, /isolatin, /helvetica
     cs=1.2
     TVLCT, red, green, blue
     mycolors=indgen(17)
  endif else begin
     cs=2
     mycolors=red+256l*green+256l*256l*blue 
  endelse

  plot, [-6,3.5], [0.3,60], /ylog, /ystyle, xtitle='t (Gyr from present)', ytitle='distances as labeled (AU)', charsize=cs, /nodata, xthick=2, ythick=2, charthick=2, pos=[0.11,0.11,0.98,0.99], /xstyle

  nint=100
  narr=100
  tarr = -5.99e9+findgen(narr)*0.999e8
;  narr=80
;  tarr = -4.0e9+findgen(narr)*.999e8
  aBcarr = dblarr(narr,nint)
  aAcarr = dblarr(narr,nint)

  stableones = []
  
  for j=0,99 do begin
     readcol, 'insker2/integration2.'+string(j,format='(i02)')+'.in', masses, x, y, z
     ma=masses(0)
     mb=masses(1)
     mc=masses(2)
  
     file='outsker2/integration2.'+string(j,format='(i02)')+'.out'
     
     readcol, file, t, x, y, z, vx, vy, vz, xc, yc, zc, vxc, vyc, vzc, /silent
;     readcol, file, t, x, y, z, vx, vy, vz, xc, yc, zc, vxc, vyc, vzc

     a=[]
     e=[]
     
     for i=0,n_elements(x)-1 do begin
        kep = statetokep(x(i), y(i), z(i), vx(i), vy(i), vz(i), m = (4*!dpi^2)*(ma+mb))
        a = [a,kep(0)]
        e = [e,kep(1)]
     endfor
     
     qle67 = where(a*(1-e) le 4, nq)
     if nq ge 1 then print, "j = ", j, " has a low-dipper."

     rc = sqrt(xc(-1)^2+yc(-1)^2+zc(-1)^2)
     if rc ge 1e5 then begin
        print, "j = ", j, " ejected Proxima."
     endif else begin
        stableones = [stableones, j]
     endelse
     
     ttot = t
     atot = a
     etot = e
     
     st = sort(ttot)
;     t = ttot(st)
;     a = atot(st)
;     e = etot(st)

     if rc le 1e5 then begin 
       oplot, (t(st)-6e9)/1e9, a(st)
       oplot, (t(st)-6e9)/1e9, a(st)*(1-e(st))
       oplot, (t(st)-6e9)/1e9, a(st)*(1+e(st))
     endif
     
     aAcrit = 6.87 - 8.58*e(st) + 1.39*e(st)^2
     aBcrit = 6.17 - 7.51*e(st) + 1.03*e(st)^2

     aAcritcum = aAcrit
     aBcritcum = aBcrit
     for k=1,n_elements(aAcrit)-1 do begin
          if aAcritcum(k) ge aAcritcum(k-1) then aAcritcum(k)=aAcritcum(k-1)
          if aBcritcum(k) ge aBcritcum(k-1) then aBcritcum(k)=aBcritcum(k-1)
     endfor

     ; plot the fully time-resolved critical radii
;     oplot, (ttot-6e9)/1e9, aAcrit, color=mycolors(2)
;     oplot, (ttot-6e9)/1e9, aBcrit, color=mycolors(2), linestyle=2
;     oplot, (t-6e9)/1e9, aAcritcum, color=mycolors(2)
;     oplot, (t-6e9)/1e9, aBcritcum, color=mycolors(2), linestyle=2

;     aAcarr(*,j) = interpol(aAcritcum, t-6e9, tarr)
     aBcarr(*,j) = interpol(aBcritcum, t(st)-6e9, tarr)

     ; check that the interpolation works. 
;     oplot, tarr/1e9, aAcarr(*,j), psym=7
;     oplot, tarr/1e9, aBcarr(*,j), psym=7
     
     print, j
;     cafe_pause

                                ; ones that abruptly stopped cycling:
                                ; 4, 7, 13, 15, 22, 27, 30, 47, 57, 83
               ; one that went real low without halting: 94
  endfor
  
  if 1 eq 1 then begin

  ; paint one of them blue-to-yellow. 
     j=51
      ; 42 goes a bit wild
                                ; 51 and 84 remain regular yet have a
                                ; pretty low periastron in the past.
  file='outsker2/integration2.'+string(j,format='(i02)')+'.out'
  readcol, file, t, x, y, z, vx, vy, vz, xc, yc, zc, vxc, vyc, vzc
  a=[]
  e=[]
  for i=0,n_elements(x)-1 do begin
     kep = statetokep(x(i), y(i), z(i), vx(i), vy(i), vz(i), m = (4*!dpi^2)*(ma+mb))
     a = [a,kep(0)]
     e = [e,kep(1)]
  endfor

  st = sort(t)

  for i=0,9 do begin
    dothese = st(i*100:(i*100+101))
    oplot, (t(dothese)-6e9)/1e9, a(dothese)
    oplot, (t(dothese)-6e9)/1e9, a(dothese)*(1-e(dothese)), thick=3, color=mycolors(3+i)
    oplot, (t(dothese)-6e9)/1e9, a(dothese)*(1+e(dothese)), thick=3, color=mycolors(3+i)
 endfor
  
  end
     
  ; now work out the levels you want to represent.
  aBcp2s = dblarr(narr)  ; acrit around A, plus 2 sigma
  aBcp1s = dblarr(narr)
  aBcmed = dblarr(narr)  ; acrit around B, median
  aBcm1s = dblarr(narr)  ; acrit around B, minus 2 sigma
  aBcm2s = dblarr(narr)
  
  for k=0,narr-1 do begin
      ; for this point in time, rank order the runs, assign levels.
     saBc = sort( aBcarr(k,stableones) )
;     aBcp2s(k) = aBcarr(k,saBc(97))
;     aBcp1s(k) = aBcarr(k,saBc(84))
;     aBcmed(k) = aBcarr(k,saBc(50))
;     aBcm1s(k) = aBcarr(k,saBc(16))
;     aBcm2s(k) = aBcarr(k,saBc(2))
     aBcp2s(k) = aBcarr(k,stableones(saBc(87)))
     aBcp1s(k) = aBcarr(k,stableones(saBc(75)))
     aBcmed(k) = aBcarr(k,stableones(saBc(45)))
     aBcm1s(k) = aBcarr(k,stableones(saBc(14)))
     aBcm2s(k) = aBcarr(k,stableones(saBc(2)))
     if k eq 12 then print, "ones making it go unstable: ", stableones(saBc(0:2))
  endfor

  print, "aBcp2s = ", aBcp2s
  oplot, tarr/1e9, aBcp2s, color=mycolors(1)
  oplot, tarr/1e9, aBcp1s, color=mycolors(1)  
  oplot, tarr/1e9, aBcmed, color=mycolors(1), thick=2
  oplot, tarr/1e9, aBcm1s, color=mycolors(1)
  oplot, tarr/1e9, aBcm2s, color=mycolors(1)

  oplot, [0.00000, 1.85722e+09, 4.22675e+09, 7.77691e+09, 9.52701e+09]/1e9-6, [0.76037669, 0.76715070, 0.82168440, 0.93059693, 0.99998024], color=mycolors(2)
  oplot, [0.00000, 1.85722e+09, 4.22675e+09, 7.77691e+09, 9.52701e+09]/1e9-6, [1.3536179, 1.3663466, 1.4598883, 1.6489019, 1.7713600], color=mycolors(2)

  xyouts, -0.2, 40, 'Q=a(1+e)', charthick=2 
  xyouts, -0.2, 25, 'a', charthick=2 
  xyouts, -0.2, 9.5, 'q=a(1-e)', charthick=2
  xyouts, -2.2, 0.97, 'Habitable Zone for star B', color=mycolors(2), charthick=2 
  xyouts, -1.0, 2.8, 'probability levels of destabilization', color=mycolors(1), charthick=2 
  
  tnow=0
  print, tarr(tnow)/1e9, "Gyr: " ,  aBcp2s(tnow), aBcp1s(tnow), aBcmed(tnow), aBcm1s(tnow), aBcm2s(tnow)
  tnow=60
  print, tarr(tnow)/1e9, "Gyr: " ,  aBcp2s(tnow), aBcp1s(tnow), aBcmed(tnow), aBcm1s(tnow), aBcm2s(tnow)
  
if keyword_set(ps) then begin
    device, /close_file
    set_plot, 'x'
endif

end
