c#################################################################
      subroutine flat1(h,ro,vp,vs,h1,ro1,vp1,vs1,n,kind)
c     subroutine applies earth-flattening corrections by
c     scaling the model.for reference see Biswas 1972.
c     kind=1 for Love wave calculations
c     kind=2 for Rayleigh wave calculations
c     h -is thickness
c-input h,ro,vp,vs,n,kind; output h1,ro1,vp1,vs1.
C     TRANSFORMATION
C-----R(i)-radius of i-th boundary (free surfase is ignored)
C-----h(i) by R0*(ln(R0/R(i))-ln(R0/R(i-1));R(0)=R0
C-----vs,vp mult by factor "dif"=R0*(1/R(i+1)-1/R(i))/ln(R(i)/R(i+1))
C-----ro  mult by factor (R(i)**pwr-R(i+1)**pwr)/(ln(R(i)/R(i+1))*pwr*R0**pwr)
C-----for Rayleigh pwr=2.275 
C-----for Love     pwr=5.000  
      implicit none
      real*8 h(2),ro(2),vp(2),vs(2),hh(1000)
      real*8 h1(2),ro1(2),vp1(2),vs1(2)
c --
      integer*4 i,ii,n,kind,nm
      real*8    dif,a,fact,difr,z0,fltd,ht,hs,pwr,z1
      data a/6371.0d0/
      do i=1,n
      hh(i)=h(i)
      enddo
      pwr=2.2750d0
C     pwr=1.2750
      if(kind.eq.1) pwr=5.0d0
      nm=n-1
      hs=0.0d0
C-----transfer thickness in radius------
      do 5 i=1,n
      ht=hs
      hs=hs+hh(i)
    5 hh(i)=a-ht
C-----now hh(i) are radii  (starting from the first boundary)
c -layers' scaling
      do 10 i=1,nm
      ii=i+1
      fltd=dlog(hh(i)/hh(ii))
      dif=(1.0d0/hh(ii)-1.0d0/hh(i))*a/fltd
      difr=hh(i)**pwr-hh(ii)**pwr
      ro1(i)=ro(i)*difr/(fltd*(a**pwr)*pwr)
      vp1(i)=vp(i)*dif
   10 vs1(i)=vs(i)*dif
c half space scaling
      fact=a/hh(n)
      vp1(n)=vp(n)*fact
      vs1(n)=vs(n)*fact
      ro1(n)=ro(n)*(1.0d0/fact)**pwr
      z0=0.0d0
C-----new thicknesses----------------------------
      do 15 i=2,n
      z1=a*dlog(a/hh(i))
      h1(i-1)=z1-z0
   15 z0=z1
      h1(n)=0.0d0
      return
      end
