;plot_system.pro

pflag=0

restore,'results.dat'
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

;longitude relative to planet0
twopi=2.0*!pi
O_rel=fltarr(Ntimes, Ntotal)
for p=0,Ntotal-1 do begin $
  Orel=O(*,p)-O(*,planet0) & $
  j=where(Orel gt !pi) & $
  if (j(0) ne -1) then Orel(j)=Orel(j)-2.0*!pi & $
  j=where(Orel lt -!pi) & $
  if (j(0) ne -1) then Orel(j)=Orel(j)+2.0*!pi & $
  O_rel(0,p)=Orel & $
endfor

;reconstructed longitudes O_recon,
;for calculating KD=wavenumber*delta
O_recon=fltarr(Ntimes, Nrings)
KD=fltarr(Ntimes, Nrings)
for tm=0,Ntimes-1 do begin $
  O_recon(tm,0)=O(tm,rings) & $
  for p=1,Nrings-1 do begin & $
    diff=O_recon(tm,p)-O_recon(tm,p-1) & $
    if (diff gt !pi) then O_recon(tm,p)=O_recon(tm,p:*)-twopi & $
    if (diff lt -!pi) then O_recon(tm,p)=O_recon(tm,p:*)+twopi & $
  endfor & $
  KD(tm,0)=deriv(ar, O_recon(tm,*))*delta & $
endfor

;first wavelength
p=where(O_recon(Ntimes-1,*) gt (O_recon(Ntimes-1,0)+twopi))
p=p(0)
lambda0=ar(p)-ar(0)

;expected KD=wavenumber*delta
D=0.87
gamma=(5.25*!pi*J2/mud)*(Rcenter*delta)^2
x=aj(rings)-aj(rings(0))
KD_exp=( 1.0 + gamma*(1.0+x/delta))/(!pi*D)
KD0=KD_exp(0)

;expected wave amplitude
i_exp=0.5*(mp(0)/mud)*(KD0/delta)

;satellite's inclination decay rate beta=(di/dt)/i/n
C=0.32
beta=deriv(t,i(*,planet0))/i(*,planet0)/twopi
beta_exp=-(C/twopi)*mp(0)/(delta^2)

;print stats
print,'mp, mr, mud = ',mp(0), mj(rings(0)), mud & $
print,'delta, h, delta/h = ', delta,  hr0, delta/hr0 & $
print,'gamma = ', gamma & $
print,'C, D = ', C, D & $
print,'mp/2/mud/delta = ', 0.5*mp(0)/mud/delta & $
print,'lambda0/delta = ',lambda0/delta

;plot e and w evolution
skip=1
a_adjust=0.9995
e_adjust=0.8
amin=min(aj)*a_adjust
amax=max(aj)/a_adjust
emin=0.0
emax=e_adjust*max(i(*,rings))/i(0,planet0)
window,xsize=650,ysize=850,retain=2
!p.multi=[0,1,4]
for tm=0,Ntimes-1,skip do begin $
  time='t='+strtrim(string(t(tm)),2)+' orbits' & $
  plot,aj(rings),i(tm,rings)/i(tm,planet0), $
    yrange=[emin,emax], charsize=2,xtitle='semimajor axis a', $
    title=time,ytitle='inclination i/i!ls!n',psym=3, $
    xrange=[amin,amax],xstyle=1,ystyle=1 & $
  oplot,ar,ar*0+i_exp,linestyle=2 & $
  if (rings_in(0) ne -1) then oplot,aj(rings_in),i(tm,rings_in)/i(tm,planet0),psym=10 & $
  if (rings_out(0) ne -1) then oplot,aj(rings_out),i(tm,rings_out)/i(tm,planet0),psym=10 & $
  plots,aj(planets),i(tm,planets),psym=8 & $
  plot,aj(rings),O_rel(tm,rings)/!pi,yrange=[-1,1], $
    psym=1,charsize=2,xrange=[amin,amax],xstyle=1,$
    xtitle='semimajor axis a', $
    ytitle='longitude !7x!3/!7p!3' & $
  if (rings_in(0) ne -1) then oplot,aj(rings_in),O_rel(tm,rings_in)/!pi,color=128 & $
  if (rings_out(0) ne -1) then oplot,aj(rings_out),O_rel(tm,rings_out)/!pi,color=128 & $
  plots,aj(planet0),O_rel(tm,planet0)/!pi,psym=8 & $
  oplot,aj,aj*0+O_rel(tm,planet0)/!pi,color=128 & $
  plot,aj(rings),KD(tm,rings),yrange=[0,5], $
    psym=0,charsize=2,xrange=[amin,amax],xstyle=1, ystyle=1,$
    xtitle='semimajor axis a',ytitle='wavenumber K!7D!3' & $
  oplot,aj(rings),KD_exp,linestyle=2 & $
  plot,t,beta,xtitle='time t    (orbits)', ytitle='e-damping rate !7b!3',charsize=2, $
    yrange=max(abs(beta))*[-1.2, 0.1],xstyle=1 & $
  oplot,t,t*0+beta_exp,color=128 & $
  wait,0.05 & $
endfor
!p.multi=0
