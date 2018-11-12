         subroutine init(dx,nlay_deriv,kind)
C---------to initialize the program--------------------------
        implicit none
        integer*4 nsize,nper,nlay_deriv,kind
        real*8    dx
c ---
        parameter (nsize=1000,nper=1000)
c ---
        real*8     a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  a,b,rho,d,qs
c --- 
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        real*8     c(nper,10),t(nper),ratio(nper,10)
        common/o/  c,t,ratio
c ---
        real*8    t1,dt,c1,dc,cinit
        integer*4 lstop,iq,istru
        common/newly/ t1,dt,c1,dc,cinit,lstop,iq,istru
c ---
        integer*4 nderiv,ndpth
        real*8    dpth(nsize),dspz(nsize),dsprfi(nsize),drvz(5,nsize),
     *            drvrfi(5,nsize)
        common/deriv/ nderiv,ndpth,dpth,dspz,dsprfi,drvz,drvrfi
c ---
        logical KEY_DERIV,KEY_RHO,KEY_TH,KEY_VP,KEY_VS,KEY_ATTEN,KEY_FLAT
        logical KEY_EIGEN,KEY_EIG_NORM,KEY_EIGEN_DER1,KEY_EIGEN_DER2,KEY_CINIT
        common/log/ KEY_ATTEN,KEY_FLAT,KEY_DERIV,KEY_EIGEN,KEY_EIG_NORM,
     *              KEY_EIGEN_DER1,KEY_EIGEN_DER2,KEY_VP,KEY_VS,KEY_TH,KEY_RHO
        common/i_nit/KEY_CINIT
c ---
        real*8 a_ref(nsize),b_ref(nsize),rho_ref(nsize),d_ref(nsize),qs_ref(nsize),dpth_ref(1000)
        common/ref/ a_ref,b_ref,rho_ref,d_ref,qs_ref,dpth_ref
C-------------------------------------------------------------------------
        integer*4 i,i_d,init_c,k_min,k_max,ln,le,ls,lt
        integer*4 ncount,lns,lqq,narg
        real*8    dpth0,ddepth,r0,tmax,perturb
        character*10 symbol,type
        character*1  symb
        integer*4    lnblnk
        character*80 infile,outfile
        character*80 der_phvfile1,der_phvfile2,der_grvfile1,der_grvfile2
        character*80 phvfile,grvfile
        character*80 eigenfileV,eigenfileH
        character*80 attfile
        data ndiv/5/,fact/4.0/,istru/0/,dc/0.01/,iedit/2/,nderiv/1/,perturb/1.0/
        data ndpth/1000/,dpth0/0.0/,ddepth/1.0/,r0/6371.0/
c       fact - number of wavelenghts below first layer where phase
c       velocity is less than shear velocity that we may consider
c       to be effective halfspace.
C--------------INITIATION------------------------------------------S
        KEY_FLAT=      .FALSE.
        KEY_EIGEN=     .FALSE.
        KEY_EIG_NORM=  .FALSE.
        KEY_EIGEN_DER1=.FALSE.
        KEY_EIGEN_DER2=.FALSE.
        KEY_DERIV=     .FALSE.
        KEY_RHO=       .FALSE.
        KEY_VP=        .FALSE.
        KEY_VS=        .FALSE.
        KEY_TH=        .FALSE.
        KEY_ATTEN=     .FALSE.
        KEY_CINIT=     .FALSE.
        ncount=0
        i_d=0
c----------------------------------------------------------------------
        write(*,*) 'surface wave forward calculations for earthquake seismology problems:'
        write(*,*) 'options for flattenning, attenuation, partial derivatives,'
        write(*,*) 'eigenfunctions and their derivatives by depth'
C-------interpreting arguments--------------------------------S
        narg=iargc()
        if (narg.lt.8.or.narg.gt.17)then
        write(*,*) 'USAGE: SURF_PERTURB infile  outfile R/L kmin kmax tmin tmax tstep'
        write(*,*) 'f] [-a] [-s depth_step] [-d par_type nl dx] [-e[1/2/n]] [-c  cinit]
     *  [-p perturb]'
        STOP
        endif
        kind=0
        call GETARG(1,infile)
        call GETARG(2,outfile)
        ln=lnblnk(outfile)
        open(1,file=infile,STATUS='OLD')
        call GETARG(3,symb)
        open(2,file=outfile(1:ln)//'.'//symb)
        if(symb.eq.'L') kind=1
        if(symb.eq.'R') kind=2
        PRint*,symb, kind
        if(symb.ne.'R'.and.symb.ne.'L') STOP' WRONG WAVE TYPE'
          call GETARG(4,symbol)
        read(symbol,*) k_min
          call GETARG(5,symbol)
        read(symbol,*) k_max
          call GETARG(6,symbol)
        read(symbol,*) t1    
          call GETARG(7,symbol)
        read(symbol,*) tmax
          call GETARG(8,symbol)
        read(symbol,*) dt
        if(narg.gt.8)  then
        PRint*,'dt=',dt
                do i=9,narg
        call GETARG(i,symbol)
        if (symbol(1:2).eq.'-a')KEY_ATTEN=.true.
        if (symbol(1:2).eq.'-f')KEY_FLAT=.true.
        if (symbol(1:2).eq.'-c')then
          KEY_CINIT=.true.
          init_c=i+1
                                 endif
        if (symbol(1:2).eq.'-s') then
          call GETARG(i+1,symbol)
          read(symbol,*)ddepth
                                 endif
        if (symbol(1:2).eq.'-d') then
          KEY_DERIV=.true.
          i_d=i
                                 endif  
        if (symbol(1:2).eq.'-p') then
          call GETARG(i+1,symbol)
           read(symbol,*) perturb
                                 end if

        if (symbol(1:2).eq.'-e')then
          KEY_EIGEN=.true.
          le=lnblnk(symbol)
          if(le.eq.2) go to 1111
           if(symbol(3:3).eq.'1'.or.symbol(3:3).eq.'2') KEY_EIGEN_DER1=.true.
           if(symbol(3:3).eq.'2') KEY_EIGEN_DER2=.true.
           if(symbol(3:3).eq.'n') KEY_EIG_NORM=.true.
1111                           endif        
                             enddo
                  endif
C-------interpreting arguments--------------------------------E
C--------DEFINING and OPENING of OUTPUT FILES------------------------S
C--------for velocities---------------------------------S
           grvfile(1:ln+6)=outfile(1:ln)//'.'//symb//'.grv'
           phvfile(1:ln+6)=outfile(1:ln)//'.'//symb//'.phv'
           open(21,file=phvfile(1:ln+6))
           open(22,file=grvfile(1:ln+6))
C--------for velocities---------------------------------E
C--------for initial phase velocity---------------------S
         if(KEY_CINIT) then
         call GETARG(init_c,symbol)
         read(symbol,*) cinit
         write(*,*) KEY_CINIT,'starting phase velocity=',cinit
                       endif
C--------for initial phase velocity---------------------E
C--------for partial derivatives found by numerical technique-S
           if(KEY_DERIV) then
           call GETARG(i_d+1,type)
           lt=lnblnk(type)
           if(type(1:lt).eq.'vs'.or.type(1:lt).eq.'VS')KEY_VS=.true.
           if(type(1:lt).eq.'vp'.or.type(1:lt).eq.'VP')KEY_VP=.true.
           if(type(1:lt).eq.'rho'.or.type(1:lt).eq.'RHO')KEY_RHO=.true.
           if(type(1:lt).eq.'th'.or.type(1:lt).eq.'TH')KEY_TH=.true.
           call GETARG(i_d+3,symbol)
           read(symbol,*) dx     
           call GETARG(i_d+2,symbol)
           ls=lnblnk(symbol)
         read(symbol,*) nlay_deriv
         lqq=ln+lt+ls+12 
         der_phvfile1(1:lqq)=outfile(1:ln)//'_'//type(1:lt)//
     *                       '_'//symbol(1:ls)//'.'//symb//'.phv_1st'
         der_phvfile2(1:lqq)=outfile(1:ln)//'_'//type(1:lt)//
     *                       '_'//symbol(1:ls)//'.'//symb//'.phv_2nd'
         der_grvfile1(1:lqq)=outfile(1:ln)//'_'//type(1:lt)//
     *                       '_'//symbol(1:ls)//'.'//symb//'.grv_1st'
         der_grvfile2(1:lqq)=outfile(1:ln)//'_'//type(1:lt)//
     *                       '_'//symbol(1:ls)//'.'//symb//'.grv_2nd'
        open(23,file=der_phvfile1(1:lqq))
        open(24,file=der_phvfile2(1:lqq))
        open(25,file=der_grvfile1(1:lqq))
        open(26,file=der_grvfile2(1:lqq))
                         endif
C---------to make depths------------------------------S
           dpth(1)=dpth0
           dpth_ref(1)=dpth(1)
           do i=2,ndpth
             dpth(i)=dpth(i-1)+ddepth
             dpth_ref(i)=dpth(i)
           enddo
           if(KEY_FLAT) then
           do i=1,ndpth
             dpth_ref(i)=dpth(i)
             dpth(i)=R0*dlog(R0/(R0-dpth(i)))
           enddo
                        endif
C---------to read depths------------------------------E
C--------for partial derivatives found by numerical technique--E
C--------for eigenfunctions------------------------------------S
         if(KEY_EIGEN) then
         eigenfileH=outfile(1:ln)//'.'//symb//'_HOR.eig'
         lns=lnblnk(eigenfileH)
         open(27,file=eigenfileH(1:lns))
         if(symb.eq.'R') then
         eigenfileV=outfile(1:ln)//'.'//symb//'_VER.eig'
         lns=lnblnk(eigenfileV)
         open(28,file=eigenfileV(1:lns))
                      endif
                      endif
C--------for eigenfunctions------------------------------------E
C--------for attenuation--------------------------------------S
         if(KEY_ATTEN) then
        attfile=outfile(1:ln)//'.'//symb//'.att'
        lns=lnblnk(attfile)
        open(29,file=attfile(1:lns))
                      endif
C--------for attenuation---------------------------------------E
C--------DEFINING and OPENING of OUTPUT FILES------------------------E
        if(KEY_EIGEN) iedit=2
        if(KEY_EIGEN_DER2) nderiv=2
C-------reading of the input model into "ref_files"-----S
c----------------------------------------------------------------------
c       d - thickness of layer in kilometers,
c       a - compressional wave velocity in km/sec,
c       b - transvers wave velocity in km/sec,
c       rho - density in gm/sm**3.
c       qs - Q factor for S-wave
c----------------------------------------------------------------------
         PRint*,'NSIZE=',nsize
        do i=1,nsize
        if(KEY_ATTEN) read(1,*,end=99)d_ref(i),a_ref(i),b_ref(i),rho_ref(i),qs_ref(i)
        if(.NOT.KEY_ATTEN) read(1,*,end=99)d_ref(i),a_ref(i),b_ref(i),rho_ref(i)
        end do
99      mmax=i-1
        PRint*,'init_mmax=',mmax
        if (mmax.lt.2) STOP' LESS THAN 2  LAYERS'
        d(mmax)=0.0d0
        if(KEY_ATTEN) then
        do i=1,mmax
        if(qs_ref(i).eq.0.0) qs_ref(i)=10000.
        qs_ref(i)=1./qs_ref(i)
        enddo
                     endif
C-------reading of the input model into "ref_files"-----E
        mode=k_max-k_min+1
c-------to define an array of periods--------------------------------S
        kmax=int((tmax-t1)/dt)+1
        PRint*,'kmax=',kmax
        if(kmax.gt.nper) then
        write(*,*) 'Too many periods=',kmax,'>',nper
        kmax=nper
                         endif
        write(*,*) 'NLAYERS=', mmax,' NUMB. OF PERIODS=',kmax
        t(1)=t1
                                do i=2,kmax
        t(i)=t(i-1)+dt
                        enddo
CNEW-----------------------------------------------S
              do i=1,kmax
        t(i)=t(i)*perturb
              enddo
CNEW-----------------------------------------------E

c-------to define an array of periods--------------------------------E
        write(*,*) '    YOUR CHOICE:'
        write(*,*) 'DDEPTH=',ddepth, ' PERTURB=',perturb
        IF(symb.EQ.'R') write(*,*) 'WAVE TYPE = RAYLEIGH   '
        IF(symb.EQ.'L') write(*,*) 'WAVE TYPE = LOVE       '
        IF(KEY_FLAT) write(*,*) 'FLATTENNING          '
        IF(KEY_ATTEN) write(*,*) 'ANELASTICITY         '
        IF(KEY_EIGEN) write(*,*) 'EIGEN FUNCTIONS      '
        IF(KEY_DERIV) then
        write(*,*) 'PARTIAL DERIVATIVES for the ',nlay_deriv,'th LAYER'
        IF(KEY_VP) write(*,*) 'by VP with perturbation=',dx,'km/s'
        IF(KEY_VS) write(*,*) 'by VS with perturbation=',dx,'km/s'
        IF(KEY_RHO) write(*,*) 'by RHO with perturbation=',dx,'g/cm**3'
        IF(KEY_TH) write(*,*) 'by THICKNESS with perturbation=',dx,'km'
                     endif
        return
        end
