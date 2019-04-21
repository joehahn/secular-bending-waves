;inits.pro
;set initial masses and coordinates of massive planets/satellites
;and the low/zero mass rings.

print,'Configuring perturbers and rings'

;time-evolution
t_start= 0.0                                                        ;start time (years)
t_stop = 3.0e5                                                   ;stop time  (years)
Ntimes = 251                                                    ;number of outputs in time-evolution
t=(t_stop-t_start)*findgen(Ntimes)/(Ntimes-1)+t_start

;mass of central body in solar units, Radius (in AU), and J2
Mcenter = 1.0
Rcenter=0.452
J2=0.016298

;number of planets/satellites
Nplanets=1

;planet/satellite initial conditions.

  ;planets' two-letter ID code
  idp=['P' ] 

  ;planet masses in solar units
  mp= [ 8.7e-12 ] 

  ;planet semimajor axes in AU
  ap= [ 1.0 ] 

  ;planet eccentricities
  ep= [ 0.0 ] 

  ;planet inclinations in *degrees*
  ip= [ 1.0e-5 ] 

  ;planet *argument* of pericenter in *degrees*
  wp= [ 0.0 ] 

  ;the planets longitude of ascending node in *degrees*
  Opp=[ 0.0 ] 

  ;convert argument of peri to longitude of peri
  for p=0, Nplanets-1 do begin $
    wp[p]=wp[p]+Opp[p] & $
    if (wp[p] gt 360.0) then wp[p]=wp[p]-360.0 & $
  endfor 

  ;convert angles to radians;
  r2d=180.0/!pi 
  for p=0, Nplanets-1 do begin $
    wp[p]=wp[p]/r2d & $
    ip[p]=ip[p]/r2d & $
    Opp[p]=Opp[p]/r2d & $
  endfor 

  ;the half-width of the planets is zero
  hp=fltarr(Nplanets) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Configure low-mass rings.

;Total number of rings in the disk
;use Nrings<250 to obtain quick results in a few minutes
Nrings=550                                     ;total number of rings
Nr_inner=50                                    ;number of rings orbiting interior to perturber
Nr_outer=Nrings-Nr_inner              ;number of rings orbiting exterior to perturber
a_min=0.99                                    ;inner ring radius (AU)
a_max=1.02                                   ;outer ring radius (AU)
delta=0.0012                                     ;gap half width (AU)
mud=5.0e-8                                     ;normalized disk mass
h_over_da=1.0                                 ;disk thickness in units of ring seperation

;set the radii (in AU) and masses (in solar units) of the rings in the inner disk,
;which orbits interior to the planet, as needed
if (Nr_inner gt 0) then begin $
  a_inner_edge=1.0-delta & $
  inner_width=a_inner_edge-a_min & $
  ar_inner=(a_inner_edge-a_min)*findgen(Nr_inner)/(Nr_inner-1)+a_min & $
  mr_inner0=2.0*mud*inner_width/Nr_inner & $                  ;individual ring mass 
  mr_inner=fltarr(Nr_inner)+mr_inner0 & $
  ;ring vertical half-width in units of semimajor axis
  ;make sure rings are overlapping some (ie, set hr0=a few*ring spacing)
  da_inner=ar_inner(1)-ar_inner(0) & $ 
  hr0_inner=h_over_da*da_inner & $
  hr_inner=hr0_inner+fltarr(Nr_inner) & $ 
endif

;set the radii (in AU) and masses (in solar units) of the rings in the outer disk,
;which orbits exterior to the planet, as needed
if (Nr_outer gt 0) then begin $
  a_outer_edge=1.0+delta & $
  outer_width=a_max-a_outer_edge & $
  ar_outer=(a_max-a_outer_edge)*findgen(Nr_outer)/(Nr_outer-1)+a_outer_edge & $
  mr_outer0=2.0*mud*outer_width/Nr_outer & $                  ;individual ring mass 
  mr_outer=fltarr(Nr_outer)+mr_outer0 & $
  ;ring vertical half-width in units of semimajor axis
  ;make sure rings are overlapping some (ie, set hr0=a few*ring spacing)
  da_outer=ar_outer(1)-ar_outer(0) & $ 
  hr0_outer=h_over_da*da_outer & $
  hr_outer=hr0_outer+fltarr(Nr_outer) & $ 
endif

;combine inner & outer disks
ar=-1
hr=-1
mr=-1
if (Nr_inner gt 0) then begin $
  ar=[ar, ar_inner] & $
  mr=[mr, mr_inner] & $
  hr=[hr, hr_inner] & $
endif
if (Nr_outer gt 0) then begin $
  ar=[ar, ar_outer] & $
  mr=[mr, mr_outer] & $
  hr=[hr, hr_outer] & $
endif & $
ar=ar[1:*]
mr=mr[1:*]
hr=hr[1:*]

;adopt uniform initial eccentricities e0 and inclinations i0 and random longitudes
e0=0.0 
i0=0.0 
er=fltarr(Nrings)+e0 
ir=fltarr(Nrings)+i0 
wr =2.0*!pi*randomu(seed,Nrings) 
Orr=2.0*!pi*randomu(seed,Nrings) 

;set ring ID
idr='r'+strarr(Nrings) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;merge planets and rings into arrays (idj, mj, aj, ei, etc)
Ntotal=Nplanets+Nrings
idj=[idp,idr]
mj=[mp,mr]
aj=[ap,ar]
ej=[ep,er]
ij=[ip,ir]
Oj=[Opp,Orr]
wj=[wp,wr]
hj=[hp,hr]

;sort according to semimajor axis
s=sort(aj)
idj=idj(s)
mj=mj(s)
aj=aj(s)
ej=ej(s)
ij=ij(s)
Oj=Oj(s)
wj=wj(s)
hj=hj(s)
