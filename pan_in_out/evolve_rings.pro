pro evolve_rings, t, Ejk, Ijk, g_i, f_i, beta_i, gamma_i, e, i, w, O

;Calculate the orbital elements of the bodies having nonzero masses (eg,
;planets and/or rings). The inputs are t=vector array of times (in years) at which
;each rings' orbit elements (e, i, w, & O) are calculated. Other inputs are
;Ejk, Ijk, g_i, f_i, beta_i, & gamma_i, which are computed by eigen_rings.pro.
;The output orbit elements (e, i, w, & O) are 2D arrays, where the first index
;is the ring index, and the second is the time index.

print,'Evolving rings over time'

;number of time samplings
Ntimes=n_elements(t)

;number of planets
Ntotal=n_elements(beta_i)

;compute the rings' eccentricity and inclination versus time
h=fltarr(Ntimes,Ntotal)
k=fltarr(Ntimes,Ntotal)
p=fltarr(Ntimes,Ntotal)
q=fltarr(Ntimes,Ntotal)
for jj=0,Ntotal-1 do begin $
  for kk=0,Ntotal-1 do begin $
    h(0,jj)=h(*,jj)+Ejk(jj,kk)*sin(g_i(kk)*t+ beta_i(kk)) & $
    k(0,jj)=k(*,jj)+Ejk(jj,kk)*cos(g_i(kk)*t+ beta_i(kk)) & $
    p(0,jj)=p(*,jj)+Ijk(jj,kk)*sin(f_i(kk)*t+gamma_i(kk)) & $
    q(0,jj)=q(*,jj)+Ijk(jj,kk)*cos(f_i(kk)*t+gamma_i(kk)) & $
  endfor & $
endfor
e=sqrt(h^2+k^2)		;eccentricity
w=atan(h,k)			;longitude of pericenter
i=sqrt(p^2+q^2)		;inclination in radians
O=atan(p,q)			;longitude of asc. node in radians
j=where(w lt 0.0)
if (j(0) ne -1) then w(j)=w(j)+2.0*!pi
j=where(O lt 0.0)
if (j(0) ne -1) then O(j)=O(j)+2.0*!pi

return
end
