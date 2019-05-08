;paper_figs.pro
;  generate figures for bending wave paper
pflag=0

;read model output
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
pp=where(O_recon(Ntimes-1,*) gt (O_recon(Ntimes-1,0)+twopi))
pp=pp(0)
lambda0=ar(pp)-ar(0)

;expected KD=wavenumber*delta
D=0.87
epsilon=1.0                                                            ;for one-sided disk
Jc=4.0*mud/(Rcenter*delta)^2/(21.0*!pi)
x=aj(rings)-aj(rings(0))
KD_exp=( epsilon+(J2/Jc)*(1.0+x/delta) )/(!pi*D)
KD0=KD_exp(0)

;expected wave amplitude I/Is
i_exp=0.5*(mp(0)/mud)*(epsilon+J2/Jc)/(!pi*D*delta)

;satellite's inclination decay rate idot=(di/dt)/i/n
C=0.32
idot=deriv(t,i(*,planet0))/i(*,planet0)/twopi
idot_exp=-(C/twopi)*mp(0)/(delta^2)

;print stats
da=ar(rings(1))-ar(rings(0))
print,'mp, mr, mud = ',mp(0), mj(rings(0)), mud & $
print,'Delta, da, h/da= ', delta, da, h_over_da & $
print,'J2, Jc, J2/Jc = ', J2, Jc, J2/Jc & $
print,'C, D = ', C, D & $
print,'mp/2/mud/delta = ', 0.5*mp(0)/mud/delta & $
print,'lambda0, delta lambda0/delta= ',lambda0, delta, lambda0/delta & $
print,'KD0(formula), i_exp, idot = ', KD0, i_exp, idot_exp

;plot inclinations at selected times
times=[    13,      50,    100,    150,    187]
clrs   =[  216,    254,      48,      31,        0]
xt     =[1.0011, 1.0053, 1.01067, 1.0156, 1.019]
imax=0.013
amin=0.999
amax=1.0225
loadct,39
if (pflag eq 0) then window,xsize=550,ysize=550,retain=2
if (pflag gt 0) then setplot
thck=3*pflag+1
rings2=rings(0:Nrings-2)
for t_idx=0,n_elements(times)-1 do begin $
  tm=times(t_idx) & $
  if (t_idx eq 0) then $
    plot,aj(rings),i(tm,rings)/i(tm,planet0), $
      yrange=[0,imax], charsize=1.5,xtitle='semimajor axis a/a!ls!n',nodata=1,  $
      title=time,ytitle='inclination I/I!ls!n',psym=3,xrange=[amin,amax],xstyle=1,ystyle=1, $
      thick=thck,xthick=thck, ythick=thck, charthick=thck, $
      color=0,background=255 & $
  oplot,aj(rings2),i(tm,rings2)/i(tm,planet0),color=clrs(t_idx), thick=thck & $
  xyouts, xt(t_idx), 0.0077, 't = '+strtrim(string(round(t(tm)/1000.0)),2), color=clrs(t_idx), $
    charthick=thck, charsize=1.1 & $
endfor
oplot,ar,ar*0+i_exp,linestyle=2,color=0,thick=thck
plots,1,0,psym=8,symsize=1.5, color=0
loadct,0
if (pflag gt 0) then output_plot,'i.ps'

;plot w and k
if (pflag eq 0) then window,xsize=550,ysize=750,retain=2
if (pflag gt 0) then setplot,large=1
!p.multi=[0,1,2]
tm=max(times)
amax=1.0201
;node
plot,aj(rings),O_rel(tm,rings)/!pi,yrange=[-1,1], $
  psym=3,charsize=1.5,xrange=[amin,amax],xstyle=1,$
  xtitle='semimajor axis a/a!ls!n',ytitle='longitudes !7X-X!3!ls!n    (!7p!3)', $
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
for j=1,Nrings-2 do begin $
  dw=O_rel(tm,rings(j-1:j))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(1)=1.0 & $
  if (ddw gt 1.0 ) then dw(1)=-1.0 & $
  oplot,ar(j-1:j),dw,thick=thck & $
  dw=O_rel(tm,rings(j:j+1))/!pi & $
  ddw=dw(1)-dw(0) & $
  if (ddw lt -1.0) then dw(0)=-1.0 & $
  if (ddw gt 1.0 ) then dw(0)=1.0 & $
  oplot,ar(j:j+1),dw,thick=thck & $
endfor
oplot,[amin, amax],[0,0],color=128,thick=thck
plots,aj(planet0),O_rel(tm,planet0),psym=8,symsize=1.5
;wavenumber
plot,aj(rings),KD(tm,rings),yrange=[0,4], $
  psym=0,charsize=1.5,xrange=[amin,amax],xstyle=1, ystyle=1,$
  xtitle='semimajor axis a/a!ls!n',ytitle='wavenumber |k|a!7D!3',$
  thick=thck,xthick=thck, ythick=thck, charthick=thck & $
oplot,aj(rings),KD_exp,linestyle=2, thick=thck+1, color=128 & $
!p.multi=0
if (pflag gt 0) then output_plot,'wk.ps'

;initial wavenumber
print,'initial wavenumber (model) = ',avg(KD(tm,rings(0)))

;i-damping rate
if (pflag eq 0) then window,xsize=550,ysize=550,retain=2
if (pflag gt 0) then setplot
plot,t,idot/1.0e-7,xtitle='time t    (orbits)', $
  ytitle='damping rate dI!ls!n/dt    (10!u-7!nI!ls!nn!ls!n)', $
  charsize=1.5,yrange=[-2.5, 0.5],xstyle=1,ystyle=1, $
  thick=thck,xthick=thck, ythick=thck, charthick=thck
oplot,t,t*0+idot_exp/1.0e-7,color=128, thick=thck
if (pflag gt 0) then output_plot,'idot.ps'

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
tm=where(t lt max(t(times)))
print,'maximum dLi = ',max(dLi(tm))
oplot,t,t*0+max(dLi(tm)),linestyle=2
