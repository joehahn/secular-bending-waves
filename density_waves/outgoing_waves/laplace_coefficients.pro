;laplace_coefficients.pro
;version 1.0
;November 6, 2002
;Joe Hahn
;hahn@lpi.usra.edu

function EllipticK,arg
;compute the complete elliptic integral of the first kind K(m).
;valid for 0<=m<1. See Abromowitz & Stegun. Fractional errors are
;smaller than 10^(-8) when compared to Maple.
m=arg^2
m1=1.0d0-m
K=1.38629436112d0     + 0.09666344259d0*m1  +0.03590092383d0*m1^2 +$
  0.03742563713d0*m1^3+ 0.01451196212d0*m1^4+                      $
[ 0.5d0               + 0.12498593597d0*m1  +0.06880248576d0*m1^2 +$
  0.03328355346d0*m1^3+ 0.00441787012d0*m1^4 ]*alog(1.0d0/m1)
return,K
end

function EllipticE,arg
;compute the complete elliptic integral of the second kind E(m).
;valid for 0<=m<1. See Abromowitz & Stegun. Fractional errors are
;smaller than 10^(-8) when compared to Maple.
m=arg^2
m1=1.0d0-m
E = 1.0d0+0.44325141463d0*m1  +0.06260601220d0*m1^2 + $
          0.04757383546d0*m1^3+0.01736506451d0*m1^4 + $
        [ 0.24998368310d0*m1  +0.09200180037d0*m1^2 + $
          0.04069697526d0*m1^3+0.00526449639d0*m1^4 ]*alog(1.0d0/m1)
return,E
end

function b0_1haf,b
;laplace coefficient b^(0)_1/2(b)
lc=b*0.0d0
j=where(b lt 1.0d0)
if (j(0) ne -1) then lc(j)=1.273239544735163d0*EllipticK(b(j))
j=where(b gt 1.0)
if (j(0) ne -1) then $
  lc(j)=1.273239544735163d0*EllipticK(1.0d0/b(j))/b(j) & $
return,lc
end

function b1_1haf,b
;laplace coefficient b^(1)_1/2(b)
lc=b*0.0d0
j=where(b lt 1.0d0)
if (j(0) ne -1) then $
  lc(j)=(1.273239544735163d0/b(j))*(EllipticK(b(j))-EllipticE(b(j)))
j=where(b gt 1.0d0)
if (j(0) ne -1) then lc(j)=(1.273239544735163d0)* $
  (EllipticK(1.0d0/b(j))-EllipticE(1.0d0/b(j)))
return,lc
end

function b1_3haf,b
;laplace coefficient b^(1)_3/2(b)
b_min=0.002
b_max=100.0
lc=b*0.0
j=where(b lt b_min)
if (j(0) ne -1) then begin
  ;use approximate form lc=3*b when b<b_min 
  bj=b(j)
  lc(j)=3.0*bj
endif
j=where((b ge b_min) and (b le b_max))
if (j(0) ne -1) then begin
  ;use exact form when b>b_min
  bj=b(j)
  lc(j)=(2.0*bj*b0_1haf(bj)-(1.0+bj^2)*b1_1haf(bj))/((1.0-bj^2))^2
endif
j=where(b gt b_max)
if (j(0) ne -1) then begin
  ;use approximate form lc=3/b^4 when b>b_max
  bj=b(j)
  lc(j)=3.0/(bj^4)
endif
return,lc
end

function b0_3haf,b
;laplace coefficient b^(0)_3/2(b)
lc=(b1_1haf(b)+(1.0+b^2)*b1_3haf(b))/(2*b)
return,lc
end

function b2_1haf,b
;laplace coefficient b^(2)_1/2(b)
lc=(1.0+b^2)*b1_1haf(b)/(2.0*b)-((1.0-b^2)^2)*b1_3haf(b)/(6.0*b)
return,lc
end

function b2_3haf,b
;laplace coefficient b^(2)_3/2(b)
b_min=0.06
b_max=10.0
lc=b*0.0
j=where(b lt b_min)
if (j(0) ne -1) then begin $
  ;use approximate form lc=3.75*b^2 when b<b_min 
  bj=b(j) & $
  lc(j)=3.75*(bj^2) & $
endif
j=where((b ge b_min) and (b le b_max))
if (j(0) ne -1) then begin $
  ;use exact form when b>b_min
  bj=b(j) & $
  lc(j)=(6.0*bj*b1_1haf(bj)-3.0*(1.0+bj^2)*b2_1haf(bj))/(1.0-bj^2)^2 & $
endif
j=where(b gt b_max)
if (j(0) ne -1) then begin $
  ;use approximate form lc=3.75/b^5 when b>b_max
  bj=b(j) & $
  lc(j)=3.75/(bj^5) & $
endif
return,lc
end

function b0_5haf,b
;laplace coefficient b^(0)_5/2(b)
b_min=0.002
b_max=1.0e3
lc=b*0.0
j=where(b lt b_min)
if (j(0) ne -1) then begin
  ;use approximate form lc=2 when b<b_min 
  bj=b(j)
  lc(j)=2.0
endif
j=where((b ge b_min) and (b le b_max))
if (j(0) ne -1) then begin
  ;use exact form when b>b_min
  bj=b(j)
  lc(j)=(2.0*bj*b1_3haf(bj)/3.0 $
    +(1.0+bj^2)*b0_3haf(bj))/((1.0-bj^2)^2)
endif
j=where(b gt b_max)
if (j(0) ne -1) then begin
  ;use approximate form lc=2/b^5 when b<b_max
  bj=b(j)
  lc(j)=2.0/(bj^5)
endif
return,lc
end

function b1_5haf,b
;laplace coefficient b^(1)_5/2(b)
b_min=1.0e-3
b_max=1.0e3
lc=b*0.0
j=where(b lt b_min)
if (j(0) ne -1) then begin
  ;use approximate form lc=5b when b<b_min 
  bj=b(j)
  lc(j)=5.0*bj
endif
j=where((b ge b_min) and (b le b_max))
if (j(0) ne -1) then begin
  ;use exact form when b>b_min
  bj=b(j)
  lc(j)=(2.0*bj*b0_3haf(bj) $
    +(1.0+bj^2)*b1_3haf(bj)/3.0)/((1.0-bj^2)^2)
endif
j=where(b gt b_max)
if (j(0) ne -1) then begin
  ;use approximate form lc=5/b^6 when b<b_max
  bj=b(j)
  lc(j)=5.0/(bj^6)
endif
return,lc
end

function integrand,x,params
;this function returns the integrand that is to be integrated by
;the function qsimp2(). Inputs are x=dependent variable,
;and parameters=array of parameters whose use is defined here.
theta=x
beta=params(0)
m=params(1)
r=params(2)
h=params(3)
hp=params(4)
H2=0.5d0*(h^2+hp^2)
den=(1.d0+beta^2)*(1.d0+H2)-2.d0*beta*cos(theta)
i=0.6366197723675814d0*cos(m*theta)/(den^r)
return,i
end

function trapzd,lower,upper,params,soln,n
;This function is called by qsimp2() which uses the trapezoidal
;rule to integrate the function integrand(). Inputs are
;lower=lower integration limit, upper=upper integration limit,
;params=array of parameters to be passed onto integrand(),
;and n=n^th call to this function. Input & output
;soln=n^th solution.
if (n eq 1) then begin
  soln=0.5d0*(upper-lower)* $
    (integrand(lower,params)+integrand(upper,params))
endif else begin
  it=(2l)^(n-2)
  tnm=it
  del=(upper-lower)/tnm
  x=lower+0.5d0*del
  sum=0.0d0
  for j=1l,it do begin
    sum=sum+integrand(x,params)
    x=x+del
  endfor
  soln=0.5d0*(soln+(upper-lower)*sum/tnm)
endelse
return,soln
end

function qsimp2,lower,upper,params,error=error,jmax=jmax
;Integrate a function using Simpson's trapezoidal rule. Inputs
;are: lower=lower integration limit, upper=upper integration
;limit, and params=array of parameters that is passes on to
;the integrand() function, which is the desired function that
;is to be integrated. The optional inputs are error=fractional
;error allowed in the solution, and jmax->maximum number of
;iterations=2^(jmax-1). This is implemented using Numerical
;Recipes subroutines qsimp and trapdz.
if (keyword_set(error) eq 0) then error=1.0d-6
if (keyword_set(jmax)  eq 0) then jmax=20
ost=-1.0d30
os =-1.0d30
j=1
flag=0
st=0.d0
while (flag eq 0) do begin
  st=trapzd(lower,upper,params,st,j)
  s=(4.d0*st-ost)/3.d0
  if (j gt 5) then $
    if ( (abs(s-os) lt error*abs(os)) or $
      ((s eq 0.d0) and (os eq 0.d0)) ) then flag=1
  if (j eq jmax) then begin
    print,'QSIMP2: too many iterations---increase jmax?'
    flag=1
  endif
  os=s
  ost=st
  j=j+1
endwhile
return,s
end

function Lap,b,m,r,h,hp,error=error
;Numerically compute the softened laplace coefficient
;~b^(m)_r(b,h) using qsimp2() which uses Simpson's trapezoidal
;rule. The softening length is h and the calculation is done to
;a fractional precision=error.
if (keyword_set(qd) eq 0) then begin
  if (keyword_set(error) eq 0) then error=1.0d-6
  lc=b*0.0d0
  params=double([0.0,m,r,h,hp])
  Nb=n_elements(b)
  pi=3.141592653589793d0
  for k=0,Nb-1 do begin
    params(0)=b(k)
    lc(k)=qsimp2(0.d0,pi,params,error=error)
  endfor
endif
return,lc
end

function b0_1half,a,h,hp
;softened laplace coefficient ~b^(0)_1/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/pi)/sqrt(2.d0*a*(g+1.d0))
b=EllipticK(X)
b=b*f
return,b
end

function b1_1half,a,h,hp
;softened laplace coefficient ~b^(1)_1/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/pi)/sqrt(2.d0*a*(g+1.d0))
b=g*EllipticK(X)-(g+1.d0)*EllipticE(X)
b=b*f
amin=0.004d0
amax=200.d0
bad=where((a lt amin) or (a gt amax))
if (bad(0) ne -1) then begin $
  bb=a/((2.d0*a*g)^1.5d0) & $
  b(bad)=bb(bad) & $
endif
return,b
end

function b0_3half,a,h,hp
;softened laplace coefficient ~b^(0)_3/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/pi)/sqrt(2.d0*a*(g+1.d0))/(2*a*(g-1.d0))
b=EllipticE(X)
b=b*f
return,b
end

function b1_3half,a,h,hp
;softened laplace coefficient ~b^(1)_3/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/pi)/sqrt(2.d0*a*(g+1.d0))/(2*a*(g-1.d0))
b=-(g-1.d0)*EllipticK(X)+g*EllipticE(X)
b=b*f
amin=0.003d0
amax=400.d0
bad=where((a lt amin) or (a gt amax))
if (bad(0) ne -1) then begin $
  bb=3.d0*a/((2.d0*a*g)^2.5d0) & $
  b(bad)=bb(bad) & $
endif
return,b
end

function b2_3half,a,h,hp
;softened laplace coefficient ~b^(2)_3/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/pi)/sqrt(2.d0*a*(g+1.d0))/(2*a*(g-1.d0))
b=-4.d0*g*(g-1.d0)*EllipticK(X)+(4.d0*(g^2)-3.d0)*EllipticE(X)
b=b*f
amin=0.025d0
amax=50.d0
bad=where((a lt amin) or (a gt amax))
if (bad(0) ne -1) then begin $
  bb=3.75d0*(a^2)/((2.d0*a*g)^3.5d0) & $
  b(bad)=bb(bad) & $
endif
return,b
end

function b0_5half,a,h,hp
;softened laplace coefficient ~b^(0)_5/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/(3.d0*pi))/((2.d0*a)^2.5d0)/((g+1.d0)^1.5d0)/(g-1.d0)^2
b=-(g-1.d0)*EllipticK(X)+4.d0*g*EllipticE(X)
b=b*f
return,b
end

function b1_5half,a,h,hp
;softened laplace coefficient ~b^(1)_3/2(a,h,hp)
pi=3.141592653589793d0
H2=0.5d0*(h^2+hp^2)
g=0.5d0*(1.d0+H2)*(1.d0/a+a)
X=sqrt(2.d0/(g+1.d0))
f=(4.d0/(3.d0*pi))/((2.d0*a)^2.5d0)/((g+1.d0)^1.5d0)/(g-1.d0)^2
b=-g*(g-1.d0)*EllipticK(X)+(g^2+3.d0)*EllipticE(X)
b=b*f
amin=0.0015d0
amax=800.d0
bad=where((a lt amin) or (a gt amax))
if (bad(0) ne -1) then begin $
  bb=5.d0*a/(2.d0*a*g)^3.5d0 & $
  b(bad)=bb(bad) & $
endif
return,b
end
