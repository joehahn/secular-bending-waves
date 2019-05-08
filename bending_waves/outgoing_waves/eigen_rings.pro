pro eigen_rings, Mcenter, mj, aj, ej, ij, Oj, wj, hj, $
  g_i, f_i, Ejk, Ijk, beta_i, gamma_i, J2=J2, Rcenter=Rcenter, $
  Ajk=Ajk, Bjk=Bjk

;eigen_rings.pro
;version 4.0
;December 3, 2004
;by Joe Hahn
;Saint Mary's University
;jhahn@ap.smu.ca
;
;This proceedure calculates the eigenvalues and eigenvectors that are used
;in the solution for secular evolution that is given
;in "Methods of Celestial Mechanics" by Brouwer and Clemence,
;1961, chapter 16, and also in "Solar System Dynamics" by Murray
;and Dermott, 1999, Chapter 7. The notation used here is similar
;to that used by Murray and Dermott, and this code is also described in
;Hahn, 2003, Ap.J, 595, 531. If you find these codes useful, 
;or find any bugs, or have any suggestions, then please send me an email.
;
;The inputs are: Mcenter = mass of central body in solar units (a scalar, usually unity),
;mj = masses of each orbiting rings (a vector, in solar units),
;aj = rings' semimajor axes (vector), ej = rings' eccentricity, ij = inclinations, 
;Oj = longitude of ascending node, wj = longitude of periapse, 
;hj = ring's half-width (in AU). The optional inputs are
;J2 = the primary's zonal harmonic (default = 0.0), and 
;Rcenter=radius of central body in AU. All masses are in 
;solar units, lengths are in AU, times are in years, and angles are in radians.
;
;The outputs are g_i & f_i = system's eccentricity & inclination 
;eigenfrequencies (radians/year), while Ejk & Ijk are 2D arrays containing
;the e and i eigenvectors, and the vectors beta_i and gamma_i are angles
;(in radians) pertaining to initial conditions---see Murray and Dermott for
;details. The optional outputs are the Ajk and Bjk arrays, in radians/year.
;
;Changes since version 1.0:
;The rings now have a vertical half-thickness per radius hj. This
;has the effect of softening their gravitational potentials,
;which results in the softening of the laplace coefficients over the
;dimensionless scale hj/aj.
;
;Also, it is not necessary to distinguish between an interior
;and exterior perturber. This is effected by setting the
;alpha_bar of Murray and Dermott to unity.

;Changes since version 2.0:
;The effects of the J2 term due to planetary oblateness has been added.
;Interactions between close rings having alpha=1+x where |x|<<1 are more
;accurately calculated and do not suffer roundoff errors. In particular,
;improved values for the f, g, and b1_3half functions are obtained when
;x and h are both small.

;Changes since version 3.0:
;The old eigen_rings.pro, which was an IDL batch script, is converted into
;an IDL proceedure. Also, the code gives reliable results when there is only
;one planet in the system.

print,'Computing system eigenvalues & eigenvectors'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Formulae used in Chapter 7 of "Solar System Dynamics".;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;planets' mean motions in units of radians/year
nj=2.0*!pi*sqrt(Mcenter+mj)/aj^1.5

;J2 and Rcenter default to zero
if (keyword_set(J2) eq 0) then J2=0.0
if (keyword_set(Rcenter) eq 0) then Rcenter=0.0

;Ajk arrays
Ntotal=n_elements(aj)
Ajk=dblarr(Ntotal,Ntotal)
for j=0,Ntotal-1 do begin $
  ;these are the j<>k elements
  for k=0,Ntotal-1 do begin $
    if (k ne j) then begin $
      alpha=aj(k)/aj(j) & $
      h=hj(j) & $
      hprime=hj(k) & $
      H2=0.5*(h^2+hprime^2) & $
      root_H2=sqrt(H2) & $
      x=alpha-1.d0 & $
      if ((abs(x) lt 3.0*root_H2) and (root_H2 lt 1.0e-6)) then begin $
        ;approximate form for g---more accurate for small x and h
        g= (2.0/!pi)*(2.0*H2-x^2)/((2.0*H2+x^2)^2) & $
      endif else begin $
        ;get g from Laplace coefficients
        g=-alpha*b2_3half(alpha,h,hprime)+ $
          3.0*(alpha^2)*H2*(2.0+H2)*b1_5half(alpha,h,hprime) & $
      endelse & $
      Ajk(j,k)=0.25*g*nj(j)*mj(k)/(Mcenter+mj(j)) & $
    endif & $
  endfor & $
endfor
for j=0,Ntotal-1 do begin $
  ;these are the j=k elements
  for k=0,Ntotal-1 do begin $
    if (k ne j) then begin $
      ;these are the j=k elements
      alpha=aj(k)/aj(j) & $
      h=hj(j) & $
      hprime=hj(k) & $
      H2=0.5*(h^2+hprime^2) & $
      root_H2=sqrt(H2) & $
      x=alpha-1.d0 & $
      if ((abs(x) lt 3.0*root_H2) and (root_H2 lt 1.0e-6)) then begin $
        ;approximate form for f---more accurate for small x and h
        f=-(2.0/!pi)*(2.0*H2-x^2)/((2.0*H2+x^2)^2) & $
      endif else begin $
        ;get f from Laplace coefficients
        f=alpha*b1_3half(alpha,h,hprime)- $
          3.0*(alpha^2)*H2*(2.0+H2)*b0_5half(alpha,h,hprime) & $
      endelse & $
      Ajk(j,j)=Ajk(j,j)+0.25*f*nj(j)*mj(k)/(Mcenter+mj(j)) & $
    endif & $
  endfor & $
  ;oblateness contribution
  if (J2 gt 0.0) then Ajk(j,j)=Ajk(j,j)+1.5*J2*((Rcenter/aj(j))^2)*nj(j) & $
endfor

;Bjk arrays
Bjk=fltarr(Ntotal,Ntotal)
for j=0,Ntotal-1 do begin $
  ;these are the j<>k elements
  for k=0,Ntotal-1 do begin $
    if (k ne j) then begin $
      alpha=aj(k)/aj(j) & $
      h=hj(j) & $
      hprime=hj(k) & $
      H2=0.5*(h^2+hprime^2) & $
      root_H2=sqrt(H2) & $
      x=alpha-1.d0 & $
      if ((abs(x) lt 3.0*root_H2) and (root_H2 lt 1.0e-6)) then begin $
        ;approximate form for b1_3half---more accurate for small x and h
        ab1=alpha*(2.0/!pi)/(x^2+2.0*H2) & $
      endif else begin $
        ;get b1_3half from the usual Laplace coefficient formala
        ab1=alpha*b1_3half(alpha,h,hprime) & $
      endelse & $
      Bjk(j,k)=0.25*ab1*nj(j)*mj(k)/(Mcenter+mj(j)) & $
    endif & $
  endfor & $
endfor
for j=0,Ntotal-1 do begin $
  ;these are the j=k elements
  for k=0,Ntotal-1 do begin $
    if (k ne j) then Bjk(j,j)=Bjk(j,j)-Bjk(j,k) & $
  endfor & $
  ;oblateness contribution
  if (J2 gt 0.0) then Bjk(j,j)=Bjk(j,j)-1.5*J2*((Rcenter/aj(j))^2)*nj(j) & $
endfor

;get eigenvalues & eigenvectors of the Ajk and Bjk matrices.
;the eigenvalues g_i and f_i have units of radians/year and the
;evectors are in radians. Note that IDL wants the rows & columns
;in these matrices transposed first. Note that the sign on
;the matrix of eigenvectors (which is ambigious at this stage)
;is flipped to agree with Murray & Dermott's results for
;Jupiter+Saturn (see page 282).

;solution for 1-planet system
if (Ntotal eq 1) then begin $
  g_i=0.0 & $
  f_i=0.0 & $
  Ejk=ej & $
  Ijk=ij & $
  beta_i=wj & $
  gamma_i=Oj & $
endif

;for 2+ planets
if (Ntotal gt 1) then begin $
  Ajk_trans=transpose(Ajk) & $			;transpose
  g_i=hqr(elmhes(Ajk_trans),double=1) & $	;eigenvalues & eigenvec's
  Ejk=-float(eigenvec(Ajk_trans,g_i,double=1,residual=Eres)) & $
  g_i=float(g_i) & $				;keep real part
  Bjk_trans=transpose(Bjk) & $			;transpose
  f_i=hqr(elmhes(Bjk_trans),double=1) & $	;eigenvalues & eigenvec's
  Ijk=float(eigenvec(Bjk_trans,f_i,double=1,residual=Ires)) & $
  f_i=float(f_i) & $				;keep real part

  ;the planets' p's and q's at time t=0
  hj0=ej*sin(wj) & $
  kj0=ej*cos(wj) & $
  pj0=ij*sin(Oj) & $
  qj0=ij*cos(Oj) & $

  ;solve for the S and T and beta and gamma vectors
  Ejk_trans=transpose(Ejk) & $
  S_sin_beta=cramer(Ejk_trans,hj0,double=1) & $
  S_cos_beta=cramer(Ejk_trans,kj0,double=1) & $
  beta_i=atan(S_sin_beta,S_cos_beta) & $
  S_i=sqrt(S_sin_beta^2+S_cos_beta^2) & $
  sign=S_sin_beta*sin(beta_i) & $
  j=where(sign lt 0.0) & $
  if (j(0) ne -1) then S_i(j)=-S_i(j) & $

  Ijk_trans=transpose(Ijk) & $
  T_sin_gamma=cramer(Ijk_trans,pj0,double=1) & $
  T_cos_gamma=cramer(Ijk_trans,qj0,double=1) & $
  gamma_i=atan(T_sin_gamma,T_cos_gamma) & $
  T_i=sqrt(T_sin_gamma^2+T_cos_gamma^2) & $
  sign=T_sin_gamma*sin(gamma_i) & $
  j=where(sign lt 0.0) & $
  if (j(0) ne -1) then T_i(j)=-T_i(j) & $

  ;rescale the Ejk and Ijk arrays
  for j=0,Ntotal-1 do begin $
    for k=0,Ntotal-1 do begin $
      Ejk(j,k)=Ejk(j,k)*S_i(k) & $
      Ijk(j,k)=Ijk(j,k)*T_i(k) & $
    endfor & $
  endfor & $

endif

return
end
