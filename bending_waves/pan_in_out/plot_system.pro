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
  KD(tm,rings_out)=deriv(ar(rings_out), O_recon(tm,rings_out))*delta & $
endfor

;first wavelength
p=where(O_recon(Ntimes-1,rings_out) gt (O_recon(Ntimes-1,rings_out(0))+twopi))
p=p(0)
lambda0=aj(rings_out(p))-aj(rings_out(0))

;expected KD=wavenumber*delta
epsilon=2.0                                          ;epsilon=1(2) for one(two)-sided disk
D=0.87
Jc=0.06063*mud/(Rcenter*delta)^2
x=aj(rings)-aj(rings_out(0))
KD_exp= epsilon + (J2/Jc)*(1.0+x/delta)
KD_exp=KD_exp/(!pi*D)
KD0=KD_exp(rings_out(0))

;expected wave amplitude
fudge_factor=0.25
i_exp=fudge_factor*(mp(0)/mud)*(epsilon + J2/Jc)/(twopi*D*delta)

;satellite's inclination decay rate beta=(di/dt)/i/n
C=0.40
beta=deriv(t,i(*,planet0))/i(*,planet0)/twopi
beta_exp=-fudge_factor*(C/twopi)*mp(0)/(delta^2)
Is_exp=exp(beta_exp*twopi*t)

;print stats
print,'mp, mr, mud = ',mp(0), mj(rings(0)), mud & $
print,'delta, h, delta/h = ', delta,  hr0, delta/hr0 & $
print,'C, D = ', C, D & $
print,'fudge_factor = ',fudge_factor & $
print,'mp/2/mud/delta = ', 0.5*mp(0)/mud/delta & $
print,'KD0, lambda0 = ', KD0, lambda0 & $
print,'J2, Jc, J2/Jc = ', J2, Jc, J2/Jc

;plot e and w evolution
skip=1
a_adjust=0.9995
e_adjust=0.5
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
  plot,aj(rings_out),KD(tm,rings_out),yrange=[0,15], $
    psym=0,charsize=2,xrange=[amin,amax],xstyle=1, ystyle=1,$
    xtitle='semimajor axis a',ytitle='wavenumber K!7D!3' & $
  oplot,aj(rings_out),KD_exp(rings_out),linestyle=2 & $
  plot,t,i(*,planet0)/i(0,planet0), $
    xtitle='time t    (orbits)', ytitle='i!ls!n',charsize=2, xstyle=1 & $
  oplot,t,Is_exp,color=128 & $
  wait,0.05 & $
endfor
!p.multi=0

stop

;generate figures for paper
pflag=0
amin=0.9944
amax=1.0205
emax=0.05
tm=63
if (pflag eq 0) then window,xs=500,ys=800,retain=2
thck=2*pflag+1
!p.multi=[0,1,2]

;I(a) at early time
if (pflag gt 0) then setplot,large=1
time='time t=7.5x10!u4!n orbits'
plot,aj(rings),i(tm,rings)/i(tm,planet0), $
  yrange=[emin,emax], charsize=1.3,xtitle='semimajor axis a/a!ls!n', $
  title=time,ytitle='inclination i/i!ls!n',psym=3, $
  xrange=[amin,amax],xstyle=1,ystyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
oplot,ar,ar*0+i_exp,linestyle=2
if (rings_in(0) ne -1) then oplot,aj(rings_in),i(tm,rings_in)/i(tm,planet0),thick=thck
if (rings_out(0) ne -1) then oplot,aj(rings_out),i(tm,rings_out)/i(tm,planet0),thick=thck
plots,aj(planets),i(tm,planets),psym=8,symsize=1.3

;Omega(a) at early time
plot,aj(rings),O_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.3,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n',ytitle='longitudes !7X-X!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
for j=1,Ntotal-2 do begin $
  dw=O_rel(tm,j-1:j)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  if ((j-1 ne planet0) and (j ne planet0)) then oplot,aj(j-1:j),dw,thick=thck & $
  dw=O_rel(tm,j:j+1)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  if ((j ne planet0) and (j+1 ne planet0))then oplot,aj(j:j+1),dw,thick=thck & $
endfor
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),O_rel(tm,planet0),psym=8,symsize=1.3
if (pflag gt 0) then output_plot,'ipan1.ps'

;I(a) at later time
tm=125
if (pflag gt 0) then setplot,large=1
time='time t=1.5x10!u5!n orbits'
plot,aj(rings),i(tm,rings)/i(tm,planet0), $
  yrange=[emin,emax], charsize=1.3,xtitle='semimajor axis a/a!ls!n', $
  title=time,ytitle='inclination i/i!ls!n',psym=3, $
  xrange=[amin,amax],xstyle=1,ystyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
;oplot,ar,ar*0+i_exp,linestyle=2
if (rings_in(0) ne -1) then oplot,aj(rings_in),i(tm,rings_in)/i(tm,planet0),thick=thck
if (rings_out(0) ne -1) then oplot,aj(rings_out),i(tm,rings_out)/i(tm,planet0),thick=thck
plots,aj(planets),i(tm,planets),psym=8,symsize=1.3

;Omega(a) at later time
plot,aj(rings),O_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.3,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n',ytitle='longitudes !7X-X!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
for j=1,Ntotal-2 do begin $
  dw=O_rel(tm,j-1:j)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  if ((j-1 ne planet0) and (j ne planet0)) then $
    oplot,aj(j-1:j),dw,thick=thck & $
  dw=O_rel(tm,j:j+1)/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  if ((j ne planet0) and (j+1 ne planet0))then oplot,aj(j:j+1),dw,thick=thck & $
endfor
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),O_rel(tm,planet0),psym=8,symsize=1.3
!p.multi=0
if (pflag gt 0) then output_plot,'ipan2.ps'


;plot i vs Omega at later time
window,xs=550,ys=550,retain=2
plot_io,O_rel(tm,*)/!pi,i(tm,*),psym=3,xtitle='longitudes !7X-X!3!ls!n    (!7p!3)', $
  ytitle='sin(I)',charsize=1.3
plots,O_rel(tm,planet0)/!pi,i(tm,planet0),psym=8

;satellite's inclination
pflag=0
thck=2*pflag+1
plot,t,i(*,planet0)/i(0,planet0),xrange=[0,1.6e5], $
  xtitle='time t    (orbits)', ytitle='satellite inclination i!ls!n/i!ls!n(0)',charsize=1.3, xstyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
oplot,t,Is_exp,color=128,thick=thck

;L-conservation
Li=fltarr(Ntimes)
nj=sqrt((1.0+mj)/aj^3)
factor=0.5*mj*nj*(aj^2)
for tm=0, Ntimes-1 do begin $
  l1=factor*reform(i(tm,*))^2 & $
  Li(tm)=total(l1) & $
endfor
dLi=abs(Li/Li(0)-1.0)
plot_io,t,dLi,yrange=[1.0e-7, 1.0e-4],charsize=1.3,xtitle='time t    (orbits)', $
  ytitle='!7D!3L!li!n/L!li!n'
;tm=where(t lt 7.5e4)
tm=where(t lt 1.5e5)
print,'maximum dLi = ',max(dLi(tm))
