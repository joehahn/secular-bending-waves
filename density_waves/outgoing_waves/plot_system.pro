;plot_system.pro

;set flags
pflag=0

;read model output
restore,'h_over_da=1/results.dat'
;restore,'results.dat'
rings=where(idj eq 'r')
planets=where(idj ne 'r')
planet0=planets(0)
rings_in=where((idj eq 'r') and (aj lt aj(planet0)))
rings_out=where((idj eq 'r') and (aj gt aj(planet0)))
sz=size(e)
Ntimes=sz(1)
Ntotal=sz(2)
Nrings=n_elements(rings)
Nplanets=n_elements(planets)
mr=mj(rings)
ar=aj(rings)
hr0=h_over_da*(ar(1)-ar(0))

;periapse longitude relative to planet0
twopi=2.0*!pi
w_rel=fltarr(Ntimes, Ntotal)
for p=0,Ntotal-1 do begin $
  wrel=w(*,p)-w(*,planet0) & $
  j=where(wrel gt !pi) & $
  if (j(0) ne -1) then wrel(j)=wrel(j)-2.0*!pi & $
  j=where(wrel lt -!pi) & $
  if (j(0) ne -1) then wrel(j)=wrel(j)+2.0*!pi & $
  w_rel(0,p)=wrel & $
endfor

;node longitude relative to planet0
O_rel=fltarr(Ntimes, Ntotal)
for p=0,Ntotal-1 do begin $
  Orel=O(*,p)-O(*,planet0) & $
  j=where(Orel gt !pi) & $
  if (j(0) ne -1) then Orel(j)=Orel(j)-2.0*!pi & $
  j=where(Orel lt -!pi) & $
  if (j(0) ne -1) then Orel(j)=Orel(j)+2.0*!pi & $
  O_rel(0,p)=Orel & $
endfor

;reconstructed longitudes w_recon,
;for calculating KD=wavenumber*delta
smth=1
w_recon=fltarr(Ntimes, Nrings)
KDw=fltarr(Ntimes, Nrings)
for tm=0,Ntimes-1 do begin $
  w_recon(tm,0)=w(tm,rings) & $
  for p=1,Nrings-1 do begin & $
    diff=w_recon(tm,p)-w_recon(tm,p-1) & $
    if (diff gt !pi) then w_recon(tm,p)=w_recon(tm,p:*)-twopi & $
    if (diff lt -!pi) then w_recon(tm,p)=w_recon(tm,p:*)+twopi & $
  endfor & $
  w_recon(tm,0)=smooth(w_recon(tm,*),smth) & $
  KDw(tm,0)=-deriv(ar, w_recon(tm,*))*delta & $
endfor

;reconstructed longitudes O_recon,
;for calculating KDo=wavenumber*delta
O_recon=fltarr(Ntimes, Nrings)
KDo=fltarr(Ntimes, Nrings)
for tm=0,Ntimes-1 do begin $
  O_recon(tm,0)=O(tm,rings) & $
  for p=1,Nrings-1 do begin & $
    diff=O_recon(tm,p)-O_recon(tm,p-1) & $
    if (diff gt !pi) then O_recon(tm,p)=O_recon(tm,p:*)-twopi & $
    if (diff lt -!pi) then O_recon(tm,p)=O_recon(tm,p:*)+twopi & $
  endfor & $
  KDo(tm,0)=deriv(ar, O_recon(tm,*))*delta & $
endfor

;first density wavelength
p=where(w_recon(Ntimes-1,*) lt (w_recon(Ntimes-1,0)-twopi))
p=p(0)
lambda_w=ar(p)-ar(0)

;first bending wavelength
p=where(O_recon(Ntimes-1,*) gt (O_recon(Ntimes-1,0)+twopi))
p=p(0)
lambda_o=ar(p)-ar(0)

;expected KD=wavenumber*delta
D=0.87
epsilon=1.0
muc=(5.25*!pi*J2)*(Rcenter*delta)^2
x=aj(rings)-aj(rings(0))
KD_exp=( epsilon + muc/mud*(1.0+x/delta))/(!pi*D)
KD0=KD_exp(0)

;expected wave amplitude
amplitude=mp(0)*(epsilon + muc/mud)/(2.0*!pi*D*mud*delta)

;satellite's e&i decay rate = (de/dt)/e/n
C=0.32
edot=deriv(t,e(*,planet0))/e(*,planet0)/twopi
idot=deriv(t,i(*,planet0))/i(*,planet0)/twopi
edot_exp=-(C/twopi)*mp(0)/(delta^2)

;print stats
da=ar(rings(1))-ar(rings(0))
print,'mp,  mr  = ',mp(0), mj(rings(0)) & $
print,'mud, muc = ',mud, muc & $
print,'Delta, da, h/da= ', delta, da, h_over_da & $
print,'h/Delta = ', h_over_da*da/Delta & $
print,'J2, KD0_exp = ', J2, KD0 & $
print,'C, D = ', C, D & $
print,'lambda_w, lambda_o = ',lambda_w, lambda_o & $
print,'amplitude, edot_exp = ', amplitude, edot_exp & $
print,'runtime (hours) = ', runtime_hours

;plot e and w evolution
skip=1
a_adjust=0.9995
ei_adjust=1.1
amin=min(aj)*a_adjust
amax=max(aj)/a_adjust
ei_max=ei_adjust*max(e(*,rings))/e(0,planet0)
window,xsize=650,ysize=850,retain=2
!p.multi=[0,1,4]
for tm=0,Ntimes-1,skip do begin $
  time='t='+strtrim(string(t(tm)),2)+' orbits' & $
  plot,aj(rings),e(tm,rings)/e(0,planet0), $
    yrange=[0,ei_max], charsize=2,xtitle='semimajor axis a/a!ls!n', $
    title=time,ytitle='eccentricity e/e!ls!n',psym=3, $
    xrange=[amin,amax],xstyle=1,ystyle=1 & $
  oplot,ar,ar*0+amplitude,linestyle=2 & $
  if (rings_in(0) ne -1) then oplot,aj(rings_in),e(tm,rings_in)/e(0,planet0),psym=10 & $
  if (rings_out(0) ne -1) then oplot,aj(rings_out),e(tm,rings_out)/e(0,planet0),psym=10 & $
  plots,aj(planets),e(tm,planets),psym=8 & $
  plot,aj(rings),w_rel(tm,rings)/!pi,yrange=[-1,1], $
    psym=1,charsize=2,xrange=[amin,amax],xstyle=1,$
    xtitle='semimajor axis a/a!ls!n', $
    ytitle='longitude !7x!3/!7p!3' & $
  if (rings_in(0) ne -1) then oplot,aj(rings_in),w_rel(tm,rings_in)/!pi,color=128 & $
  if (rings_out(0) ne -1) then oplot,aj(rings_out),w_rel(tm,rings_out)/!pi,color=128 & $
  plots,aj(planet0),w_rel(tm,planet0)/!pi,psym=8 & $
  oplot,aj,aj*0+w_rel(tm,planet0)/!pi,color=128 & $
  plot,aj(rings),KDw(tm,rings),yrange=[-1,5], $
    psym=3,charsize=2,xrange=[amin,amax],xstyle=1, ystyle=1,$
    xtitle='semimajor axis a/a!ls!n',ytitle='wavenumber |k|a!7D!3' & $
  oplot,aj(rings),KD_exp,linestyle=2 & $
  plot,t,idot,xtitle='time t    (orbits)', ytitle='damping rate (de/dt)/(en)',charsize=2, $
    yrange=max(abs(edot))*[-1.2, 0.1],xstyle=1 & $
  oplot,t,t*0+edot_exp,color=128 & $
  wait,0.05 & $
endfor
!p.multi=0
