;ring_master.pro
;evolve the system of planets and low/zero-mass rings.

;Configure planets' and rings orbits
@inits.pro

;Compute system's eigenfrequencies and eigenvectors
.run laplace_coefficients.pro
eigen_rings, Mcenter, mj, aj, ej, ij, Oj, wj, hj, $
  g_i, f_i, Ejk, Ijk, beta_i, gamma_i, Ajk=Ajk, Bjk=Bjk, $
  J2=J2, Rcenter=Rcenter

;Compute the system's evolution over time
evolve_rings, t, Ejk, Ijk, g_i, f_i, beta_i, gamma_i, e, i, w, O

;store results
save,t,aj,mj,hj,idj,e,i,O,w,mud,delta,Ejk, Ijk, g_i, f_i, beta_i, gamma_i, $
  J2, Rcenter, mp, h_over_da, filename='results.dat'

;Show system evolution
@plot_system.pro
