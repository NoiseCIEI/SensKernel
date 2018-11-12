czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c       Program SURF_PERTURB
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c-----------------------------------------------------------------------
c       This program will accept one liquid layer at the surface,
c       in wich case ellipticity of Rayleigh wave is that at the
c       top of solid array. Love wave calculations ignore liquid layer.
c       The root determination section is one of interval halving
c       once a zero crossing has been found.
c       To follow modes the program initially finds the phase velocities
c       of the 'mode' number of modes for period t1, then it finds the
c       phase velocities for period t1+dt using previous results to stay
c       on the same mode. It is preferable 'dc' positive since the
c       program will follow a mode up to its cutoff value. For 'dc'
c       negative the program cannot pick up any modes.
c       It optionally  accounts for attenuation and produces group velocity
c       derivatives by layer parameters;
c       It includes sphericity corrections by default
c       Love waves and Rayleigh waves are treated separately
c-----------------------------------------------------------------------
c         INPUT ARGUMENTS:
c-------infile outfile wavetype kmin kmax tmin tmax dt [-a] [-f] 
c                         [-d par_type  nl d_par] [-e[1/2]] [-n]
c    a. NECESSARY ARGUMENTS
c-------infile : input parameters file name----------------------------
c-------outfile: output file name----------------------------
c-------wavetype: R/L
c-------kmin,kmax - radial orders of modes to be computed
c-------tmin, tmax, tstep -  the range of periods and increment are defined
c        b. OPTIONAL  ARGUMENTS
c-------[-a] if present, attenuation is included; 
c-------[-f] if present, flattenning is included; 
c-------[-d par_type nl d_par] partial derivatives of phase and group velocities
c-------by parameter "par_type" (vp/vs/rho/thick) of the nl-th layer would be 
c-------found using +/- dx absolute perturbation.
c-------[-e[1/2/n]] if -e is present, eigenfunctions will be found; if 1 or 2 
c-------are present after -e 1st/1st and 2nd derivatives by depth will be found.
c-------if n is present after e, eigenfunctions are normalized by their maxima 
c         INPUT FILE: 
c-------file "infile" with model data and parameters of calculations; 
c-------see an example in the README file.
c         OUTPUT FILES: 
c-------grv.file : output group velocity xy-file;filename=outfile//".R.grv" or
c-------outfile//".L.grv" depending on wavetype.
c-----------------------------------------------------------------------
c-------phv.file : output phase velocity xy-file;filename=outfile//".R.phv" or
c-------outfile//".L.phv" depending on wavetype.
c-----------------------------------------------------------------------
c-------att.file : output apparent Q file;filename=outfile//".R.att" or
c-------outfile//".L.att" depending on wavetype. This file will appear only 
c-------with [-a] option!
c-----------------------------------------------------------------------
c-------eigen.file: output normalized eigenfunction as xy-file;
c-------filename=outfile//".R_VER.eig" and outfile//".R_HOR.eig" or 
c-------outfile//".L.eig". All periods are represented in one file:
c-------each period sequence starts from "T=",period value and number of depth
c-------This file will appear only with [-e] option!
c-----------------------------------------------------------------------
c-------outfile: general output file containing eigenfunctions, integrals,
c-------lagrangian, phase velocity partial derivatives, eigenfunctions 
c-------derivatives by depth, etc; to be specified by the option in the input 
c-------file (see README).
c-----------------------------------------------------------------------
        implicit none
        include "version.h"
        integer*4 nsize,nper,kind
        parameter (nsize=1000,nper=1000)
c ---
        real*8     a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  a,b,rho,d,qs
c ---
        real*8     c(nper,10),t(nper),ratio(nper,10)
        common/o/  c,t,ratio
c ---
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        real*8    t1,dt,c1,dc,cinit
        integer*4 lstop,iq,istru
        common/newly/ t1,dt,c1,dc,cinit,lstop,iq,istru
c ---
        real*8    depth(nsize),amprfi(nsize),ampz(nsize),strz(nsize),
     *            strrfi(nsize)
        integer*4 mmm
        common/rar/ depth,amprfi,ampz,strz,strrfi,mmm
c ---
        real*8 dcda(nsize),dcdb(nsize),dcdr(nsize),dwx(nsize),
     *         g1(nsize),g2(nsize)
        common/rar1/ dcda,dcdb,dcdr,dwx,g1,g2
c ---
        real*8   ccc,cvar,ugr,wvno,rat,arle
        common/rco/ ccc,cvar,ugr,wvno,rat,arle
c ---
        real*8  sumi0,sumi1,sumi2,sumi3,flagr
        common/rco1/ sumi0,sumi1,sumi2,sumi3,flagr
c ---
        integer*4 nderiv,ndpth
        real*8    dpth(nsize),dspz(nsize),dsprfi(nsize),drvz(5,nsize),
     *            drvrfi(5,nsize)
        common/deriv/ nderiv,ndpth,dpth,dspz,dsprfi,drvz,drvrfi
c ----
        logical KEY_DERIV,KEY_RHO,KEY_TH,KEY_VP,KEY_VS,KEY_ATTEN,KEY_FLAT
        logical KEY_EIGEN,KEY_EIG_NORM,KEY_EIGEN_DER1,KEY_EIGEN_DER2,KEY_CINIT
        common/log/ KEY_ATTEN,KEY_FLAT,KEY_DERIV,KEY_EIGEN,KEY_EIG_NORM,
     *              KEY_EIGEN_DER1,KEY_EIGEN_DER2,KEY_VP,KEY_VS,KEY_TH,KEY_RHO
        common/i_nit/KEY_CINIT
c ---
        real*8 a_ref(nsize),b_ref(nsize),rho_ref(nsize),d_ref(nsize),qs_ref(nsize)
        common/ref/ a_ref,b_ref,rho_ref,d_ref,qs_ref
c-------------------------------------------------------------------------
        integer*4 ilay,i,j,ncount,ndiv_store,imax(10)
        integer*4 number,numbel,mmax_store,nlay_deriv
        real*8  pi,t_base,b_corr,dx,qq,t1_store
c-------------------------------------------------------------------------
c       fact - number of wavelenghts below first layer where phase
c       velocity is less than shear velocity that we may consider
c       to be effective halfspace.
        data pi/3.1415927d0/,t_base/1.0d0/
        pi=datan(1.0d0)*4.0d0
        ncount = 0
        dx=0.0d0
        nlay_deriv=0
        write(*,'(20(1H=),1x,a,1x,20(1H=))') VERSION
        call INIT(dx,nlay_deriv,kind)
           ndiv_store=ndiv
           t1_store=t1
           mmax_store=mmax
c----------------------------------------------------------------------
c       kmax - number of periods; t1 - initial starting period;
c       dt - period increment; c1 - initial phase velocity guess;
c       dc - phase velocity increment; istru<1 - depth dependent
c       parameters will be calculated in midpoints of sublayers,
c       else program STRUT will be called.
c----------------------------------------------------------------------
c       nderiv - highest order of derivative of eigenfunction;
c       ndpth - number of depthes for which eigenfunction and its
c       derivatives are desired.
c----------------------------------------------------------------------
c       dpth(j) - array of depths for calculation of eigenfunction
c       and its derivatives.
c----------------------------------------------------------------------
C---------computations start-------------------------------------------
C##########################################################################
1     ncount=ncount+1
      ndiv=ndiv_store
      mmax=mmax_store
      t1=t1_store
      do i=1,nsize
        depth(i)=0.0d0
      enddo
      do i=1,nper
        do j=1,10    
          c(i,j)=0.0d0
          ratio(i,j)=0.0d0
        enddo
      enddo
c----------------------------------------------------------------------
c       ndiv - number of subdivisions for one layer;
c       mode - number of modes for which dispersion curves are desired;
c !!!       iedit = 0 - only dispersion, iedit = 1 - dispersion, amplitudes
c       of eigenfunctions and derivatives of phase velocity,
c----------------------------------------------------------------------
C----------------for perturbations--------------------------------S
        if(KEY_DERIV)                        then
               if(ncount.eq.2) then
        if(KEY_TH) d_ref(nlay_deriv)=d_ref(nlay_deriv)-dx
        if(KEY_VP) a_ref(nlay_deriv)=a_ref(nlay_deriv)-dx
        if(KEY_VS) b_ref(nlay_deriv)=b_ref(nlay_deriv)-dx
        if(KEY_RHO) rho_ref(nlay_deriv)=rho_ref(nlay_deriv)-dx
        KEY_EIGEN=.false.
                               endif
                if(ncount.eq.3) then
        if(KEY_TH) d_ref(nlay_deriv)=d_ref(nlay_deriv)+2*dx
        if(KEY_VP) a_ref(nlay_deriv)=a_ref(nlay_deriv)+2*dx
        if(KEY_VS) b_ref(nlay_deriv)=b_ref(nlay_deriv)+2*dx
        if(KEY_RHO) rho_ref(nlay_deriv)=rho_ref(nlay_deriv)+2*dx
                               endif
                                              endif
         if(KEY_FLAT.and..not.KEY_ATTEN)then
C--------flattening without attenuation-------------------------------
         call flat1(d_ref,rho_ref,a_ref,b_ref,d,rho,a,b,mmax,kind)
                                        endif
         if(.not.KEY_FLAT.and..not.KEY_ATTEN)then
C---------no flattening and no  attennuation-----------------
              do i=1,mmax
                d(i)=d_ref(i)
                a(i)=a_ref(i)
                b(i)=b_ref(i)
                rho(i)=rho_ref(i)
              enddo
                              endif
c----------------------------------------------------------------------
        nmax=mmax
        number=0
        numbel=0
        lstop=0
        idrop=0
        iq=0
        imax(1)=0
        ilay=1
        if(b_ref(1).eq.0.0) ilay=2
        b_corr=0.0d0
        if(KEY_ATTEN) b_corr=qs_ref(ilay)*dlog(t_base/t1)/pi
        qq=b_ref(ilay)
        if(kind.eq.2)qq=0.9d0*qq
        if(KEY_CINIT)qq=cinit
        c1=qq*(1.0d0+b_corr)
c       if(ilay.eq.2.and.kind.eq.2)c1=1.5
        write(*,*) 'START c1=',c1
c-----------------------start of main calculations------------
c       if(KEY_ATTEN.and.kind.eq.1) write(*,*) '   PERIOD (s)    ALPHA_L (1/1000km)           Q_L '
c       if(KEY_ATTEN.and.kind.eq.2) write(*,*) '   PERIOD (s)    ALPHA_R (1/1000km)           Q_R '
        call calcul(ncount,dx,imax,number,numbel,kind,t_base,*1)
              write(*,*) 'ALL IS DONE'
                              END
