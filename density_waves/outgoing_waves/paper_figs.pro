;plot_system.pro

;set flags
pflag=0

;read model output
restore,'h_over_da=1/results.dat'
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
C=0.33
edot=deriv(t,e(*,planet0))/e(*,planet0)/twopi
idot=deriv(t,i(*,planet0))/i(*,planet0)/twopi
edot_exp=-(C/twopi)*mp(0)/(delta^2)

;alternate e-damping rate = (de/dt)/e/n
efactor2=(epsilon + muc/mud)/(2.0*!pi*D*delta)
edot2=-mp(0)*efactor2^2

;print stats
da=ar(rings(1))-ar(rings(0))
print,'mp,  mr  = ',mp(0), mj(rings(0)) & $
print,'mud, muc = ',mud, muc & $
print,'Delta, da, h/da= ', delta, da, h_over_da & $
print,'h/Delta = ', h_over_da*da/Delta & $
print,'J2, KD0_exp = ', J2, KD0 & $
print,'C, D = ', C, D & $
print,'lambda_w = ',lambda_w & $
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
  plot,t,edot,xtitle='time t    (orbits)', ytitle='damping rate (de/dt)/(en)',charsize=2, $
    yrange=max(abs(edot))*[-1.2, 0.1],xstyle=1 & $
  oplot,t,t*0+edot_exp,color=128 & $
  wait,0.05 & $
endfor
!p.multi=0

;plot eccentricities at selected times
times  =[    20,     40,      60,      80,    100 ]
clrs   =[   216,    254,      48,     31,     0 ]
xt     =[1.0013,  1.0042,  1.00685,  1.01045, 1.01355 ]
yt     =[ 0.003,  0.0046,  0.0063, 0.0080, 0.0105 ]
ei_max=0.015
amin=0.999
amax=1.016
loadct,39
if (pflag eq 0) then window,xsize=550,ysize=550,retain=2
if (pflag gt 0) then setplot
thck=3*pflag+1
rings2=rings(0:Nrings-2)
for t_idx=0,n_elements(times)-1 do begin $
  tm=times(t_idx) & $
  if (t_idx eq 0) then $
    plot,aj(rings),e(tm,rings)/e(tm,planet0), $
      yrange=[0,ei_max], charsize=1.5,xtitle='semimajor axis a/a!ls!n',nodata=1,  $
      ytitle='eccentricity e/e!ls!n',psym=3,xrange=[amin,amax],xstyle=1,ystyle=1, $
      thick=thck,xthick=thck, ythick=thck, charthick=thck, $
      color=0,background=255 & $
  oplot,aj(rings2),e(tm,rings2)/e(tm,planet0),color=clrs(t_idx), $
    thick=thck & $
  xyouts, xt(t_idx), yt(t_idx), $
    't = '+strtrim(string(round(t(tm)/1000.0)),2), color=clrs(t_idx), $
  charthick=thck, charsize=1.1 & $
endfor
oplot,ar,ar*0+amplitude,linestyle=2,color=0,thick=thck
plots,1,0,psym=8,symsize=1.5, color=0
loadct,0
if (pflag gt 0) then output_plot,'e.ps'


;plot w and k
if (pflag eq 0) then window,xsize=550,ysize=750,retain=2
if (pflag gt 0) then setplot,large=1
!p.multi=[0,1,2]
tm=max(times)
amax=1.0155
;longitude
plot,aj(rings),w_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.5,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n', $
  ytitle='longitudes !7x-x!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
xyouts, 0.99615,  0.16, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
xyouts, 0.99615, -0.03, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
for j=1,Nrings-2 do begin $
  dw=w_rel(tm,rings(j-1:j))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  oplot,ar(j-1:j),dw,thick=thck & $
  dw=w_rel(tm,rings(j:j+1))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  oplot,ar(j:j+1),dw,thick=thck & $
endfor
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),w_rel(tm,planet0),psym=8,symsize=1.5
;wavenumber
plot,aj(rings),KDw(tm,rings),yrange=[-0,4],symsize=0.5, $
  psym=8,charsize=1.5,xrange=[amin,amax],xstyle=1, ystyle=1,$
  xtitle='semimajor axis a/a!ls!n',ytitle='wavenumber |k|a!7D!3',$
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
oplot,aj(rings),KD_exp,linestyle=0, thick=thck, color=128 & $
!p.multi=0
if (pflag gt 0) then output_plot,'wk.ps'


;plot w(a)
if (pflag eq 0) then window,xsize=550,ysize=550,retain=2
if (pflag gt 0) then setplot
tm=max(times)
amax=1.0155
thck=3*pflag+1
plot,aj(rings),w_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.5,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n', $
  ytitle='longitudes !7x-x!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
xyouts, 0.99615,  0.16, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
xyouts, 0.99615, -0.03, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
for j=1,Nrings-2 do begin $
  dw=w_rel(tm,rings(j-1:j))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  oplot,ar(j-1:j),dw,thick=thck & $
  dw=w_rel(tm,rings(j:j+1))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  oplot,ar(j:j+1),dw,thick=thck & $
endfor
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),w_rel(tm,planet0),psym=8,symsize=1.5
!p.multi=0
if (pflag gt 0) then output_plot,'w.ps'

;initial wavenumber
print,'expected initial wavenumber (model) = ',KD_exp(rings(0))

;e-damping rate
if (pflag eq 0) then window,xsize=550,ysize=550,retain=2
if (pflag gt 0) then setplot
tm=where(t lt max(t(times)))
tm=[0,tm+1]
plot,t(tm)/1000.0,edot(tm)/1.0e-7,xtitle='time t    (10!u3!n orbits)', $
  ytitle='damping rate de!lS!n/dt    (10!u-7!ne!lS!nn!lS!n)', $
  charsize=1.5,yrange=[-2.5, 0.5],xstyle=1,ystyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
;oplot,t/1000,t*0+edot_exp/1.0e-7,color=128, thick=thck
oplot,t/1000,t*0+edot2/1.0e-7,color=128, thick=thck, linestyle=0
if (pflag gt 0) then output_plot,'edot.ps'

;L-conservation
Lee=fltarr(Ntimes)
nj=sqrt((1.0+mj)/aj^3)
factor=0.5*mj*nj*(aj^2)
for tm=0, Ntimes-1 do begin $
  l1=factor*reform(e(tm,*))^2 & $
  Lee(tm)=total(l1) & $
endfor
dLe=abs(Lee/Lee(0)-1.0)
plot_io,t,dLe,yrange=[1.0e-7, 1.0e-3],charsize=1.3,xtitle='time t    (orbits)', $
  ytitle='!7D!3L!le!n/L!le!n'
tm=where(t lt max(t(times)))
print,'maximum dLe = ',max(dLe(tm))
oplot,t,t*0+max(dLe(tm)),linestyle=2


;plot e(a) and w(a) for reflected wave
amin=0.9995
amax=1.0155
emax=0.06
tm=350
if (pflag eq 0) then window,xs=500,ys=800,retain=2
thck=2*pflag+1
!p.multi=[0,1,2]
;e(a)
if (pflag gt 0) then setplot,large=1
time='time t=2.1x10!u5!n orbits'
plot,aj(rings),e(tm,rings)/e(tm,planet0), $
  yrange=[0.0,emax], charsize=1.3,xtitle='semimajor axis a/a!ls!n', $
  title=time,ytitle='eccentricity e/e!ls!n',psym=0, $
  xrange=[amin,amax],xstyle=1,ystyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
oplot,ar,ar*0+amplitude,linestyle=2
plots,aj(planets),e(tm,planets),psym=8,symsize=1.5
;omega(a)
plot,aj(rings),w_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.3,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n', $
  ytitle='longitudes !7x-x!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
for j=1,Ntotal-2 do begin $
  dw=w_rel(tm,j-1:j)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  if ((j-1 ne planet0) and (j ne planet0)) then $
    oplot,aj(j-1:j),dw,thick=thck,color=128 & $
  dw=w_rel(tm,j:j+1)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  if ((j ne planet0) and (j+1 ne planet0))then $
    oplot,aj(j:j+1),dw,thick=thck,color=128 & $
endfor
oplot,aj(rings),w_rel(tm,rings)/!pi, psym=8, symsize=0.5 & $
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),w_rel(tm,planet0),psym=8,symsize=1.5
xyouts, 0.99717,  0.128, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
xyouts, 0.99717, -0.03, '!14_!3', orientat=90.0, charsize=1.0, charthick=thck
if (pflag gt 0) then output_plot,'ew.ps'
