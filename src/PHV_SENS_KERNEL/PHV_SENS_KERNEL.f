             include "version.h"
             character*80 cross_sect,surf_inp_file,outfile,bred
             character*63 bred1
             character*1 wavetype,sym_mode
             character*3 sym_per 
             real*4 vel_p(400),vel_s(400),rho(400),border(400)
             real*4 depth(3000),a(3000),b(3000),rh(3000),thick(3000)
             real*4 v1(3000),dv1(3000),v2(3000),dv2(3000),v3(3000),dv3(3000)
             real*4 dcdb(3000),dcda(3000),dcdrh(3000)
             real*4  period(10,100)
             integer nper(10)
             data marg/4/,eps/0.0001/,perturb/1./
C------------------------INIT-----------------------S
             PRINT*,'to extract R/L phase velocity partial derivatives by Vs,Vp,rho from output of SURF_PERTURB'
             narg=iargc()
             if(narg.ne.marg)then
             PRINT*,'USAGE: PHV_SENS_KERNEL cross_sect surf_inp_file  wavetype(R/L) outfile'
             STOP
                            endif
             write(*,'(20(1H=),1x,a,1x,20(1H=))') VERSION
             call GETARG(1,cross_sect)
             call GETARG(2,surf_inp_file)
             lsi=lnblnk(surf_inp_file)
             lsymb=0
             if(surf_inp_file(lsi:lsi).eq.'+'.or.surf_inp_file(lsi:lsi).eq.'-')lsymb=1
             call GETARG(3,wavetype)
             if(wavetype.ne.'R'.and.wavetype.ne.'L')STOP'WRONG WAVETYPE'
             call GETARG(4,outfile)
             lso=lnblnk(outfile)
      PRINT*,'Wavetype=',Wavetype,' perturb=',perturb
             open(1,file=cross_sect,status='OLD')
      mode_end=0
C------------------------INIT-----------------------E
C------cross_section reading ----S
             border(1)=0.0
             do i=1,100000
             read(1,*,end=99) thick(i),vel_p(i),vel_s(i),rho(i)
             border(i+1)=border(i)+thick(i)
             enddo
99           nd=i-1
             PRINT *,'n_depth_points=',nd
C------cross_section reading ----E
C-----to determine requested modes and periods from central .phv file--------S
             if(lsymb.eq.0)then
             open(22,file=surf_inp_file(1:lsi)//'.'//wavetype//'.phv',status='OLD')
                           else
             open(22,file=surf_inp_file(1:lsi-1)//'.'//wavetype//'.phv',status='OLD')
                              endif
444         DO ip=1,1000
      read(22,'(a)',end=99999)bred
      if(bred(5:10).eq.'    ')then
      read(22,'(1X)',end=99999)
      go to 555
                               endif
      read(bred,*)period(mode_end+1,ip)
      enddo
555   mode_end=mode_end+1
      nper(mode_end)=ip-1
      go to 444
99999 close(22)  
C-----to determine requested modes and periods from .phv file--------E
C-------------------------Loop for modes---------------------------------------S
      dO imode=1,mode_end
      mode_req=imode-1
      write(sym_mode,'(I1)')mode_req
C            open(22,file=surf_inp_file(1:lsi)//'.'//wavetype//'.phv',status='OLD')
      if(wavetype.eq.'R')open(2,file=surf_inp_file(1:lsi)//'.R',status='OLD')
      if(wavetype.eq.'L')open(2,file=surf_inp_file(1:lsi)//'.L',status='OLD')
C-------------------------Loop for periods---------------------------------------S
      kk=3000
      DO ip=1,nper(imode)
      iper=int(period(imode,ip)/perturb+0.5)
      period_req=period(imode,ip)
      sym_per='   '
      if(iper.lt.10)write(sym_per(1:1),'(i1)')iper
      if(iper.ge.10.and.iper.lt.100)write(sym_per(1:2),'(i2)')iper
      if(iper.ge.100)write(sym_per(1:3),'(i3)')iper
      lp=lnblnk(sym_per)
      if(wavetype.eq.'R')open(8,file=outfile(1:lso)//'.R_'//sym_mode//'_'//sym_per(1:lp))
      if(wavetype.eq.'L')open(8,file=outfile(1:lso)//'.L_'//sym_mode//'_'//sym_per(1:lp))
             Do j=1,1000000
C------------search for mode-----------------------------S
             IF(wavetype.eq.'R')THEN
31      if(ip.ne.1)go to 30
             read(2,'(a)',end=888)bred             
             if(bred(33:45).ne.'Rayleigh mode')go to 31
             read(bred(48:48),*) mode
             if(mode.lt.mode_req+1)go to 31
             if(mode.gt.mode_req+1)go to 888
C------------search for mode-----------------------------E
C------------the mode is found----------------------->
C------------search for period---------------------------S
30           read(2,'(a)')bred
             if(bred(3:3).ne.'@')go to 30
             read(2,'(a)')bred1
             read(bred1,'(4(E14.7,2X))')per_q,c_r,u_r,wvn_r
             if(lsymb.eq.0.and.abs(per_q-period_req).gt.eps)go to 30
C------------search for period---------------------------E
C------------period is found------------------------>
             read(2,'(F15.4)')sumI0_r
C------------loop for depths-------S
                 do i=1,3000
C------------reading V1----S
             read(2,'(a)',end=777)bred
             if(bred(3:3).eq.'$')go to 34
             read(bred,'(3(E14.7,2X))')depth(i),v2(i),dv2(i)
                 enddo
777          STOP'Unexpected end met'
C------------reading V1----E
C------------reading V2----S
34           idepth=i-1
             do i=1,idepth
             read(2,'(a)')bred
             read(bred,'(3(E14.7,2X))')depth(i),v1(i),dv1(i)
             enddo
             d_step=depth(2)-depth(1)
             go to 35
C------------reading V2----E
C-----------------Love waves------------------------S
             ELSE
C------------search for mode-----------------------------S
             if(ip.eq.1.and.imode.eq.1)then
             read(2,'(a)',end=888)bred             
             read(2,'(a)')bred
                                           endif
             read(2,'(4(E14.7,2X))')per,c_l,u_l,wvn_l
             read(2,'(F12.4)')sumI0_l
C------------loop for depths-------S
                 do i=1,kk 
C------------reading V3----S
             read(2,'(a)',end=33)bred
             if(bred(3:3).eq.'@')go to 33
             if(bred(35:38).eq.'Love')then
             read(2,'(a)') bred
             go to 33
                                      endif
             read(bred,'(3(E14.7,2X))',end=33)depth(i),v3(i),dv3(i)
                 enddo
C------------reading V3----E
33           idepth=i-1    
             d_step=depth(2)-depth(1)
             go to 35
                STOP'Unexpected end2 met'
C------------finding parameters-S
             ENDIF
             EndDo
35           continue 
                do l=1,idepth
                do jj=1,nd-1
             if(depth(l).ge.border(jj).and.depth(l).lt.border(jj+1))then
             a(l)=vel_p(jj)
             b(l)=vel_s(jj)
             rh(l)=rho(jj)
             elseif(depth(l).ge.border(nd-1))then
             a(l)=vel_p(nd)
             b(l)=vel_s(nd)
             rh(l)=rho(nd)
                                    endif
                enddo
                enddo
             go to 36
C-----------finding parameters-E
36          continue 
C------------writing output------S
             if(wavetype.eq.'R')then
C-----------------------Rayleigh waves--------------S
             do k=1,idepth
            dcdb(k)=b(k)*rh(k)/u_r/sumI0_r*((v1(k)+dv2(k)/wvn_r)**2 +4./wvn_r*dv1(k)*v2(k))
            dcda(k)=a(k)*rh(k)/u_r/sumI0_r*((v2(k)-dv1(k)/wvn_r)**2)
            dcdrh(k)=(0.5/rh(k)*(dcda(k)*a(k)+dcdb(k)*b(k))-c_r**2/2./u_r/sumI0_r*(v1(k)**2+v2(k)**2))
      if(k.eq.1)write(8,'(F10.4,6(1X,E14.5),I2)') depth(k),dcdb(k)*b(k)/c_r,dcda(k)
     #*a(k)/c_r,dcdrh(k)*rh(k)/c_r,per_q,c_r,u_r,mode_req
      if(k.ne.1)write(8,'(F10.4,3(1X,E14.5))') depth(k),dcdb(k)*b(k)/c_r,dcda(k)*a(k)/c_r,dcdrh(k)*rh(k)/c_r
             enddo
C-----------------------Rayleigh waves--------------E
             else
C---------------------------Love waves--------------S
             do k=1,idepth
             q1=v3(k)**2
             q2=(1./wvn_l*dv3(k))**2
             dcdb(k)=b(k)*rh(k)/u_l/sumI0_l*(q1+q2)
             dcdrh(k)=0.5/rh(k)*dcdb(k)*b(k)-0.5*c_l**2/u_l/sumI0_l*v3(k)**2 
             if(k.eq.1)write(8,'(F10.4,5(1X,E14.5),I2)') depth(k),dcdb(k)*b(k)/c_l,dcdrh(k)*rh(k)/c_l,per,c_l,u_l,mode_req
             if(k.ne.1)write(8,'(F10.4,2(1X,E14.5))') depth(k),dcdb(k)*b(k)/c_l,dcdrh(k)*rh(k)/c_l
             enddo
C---------------------------Love waves--------------E
             endif
             close(8)
        ENDDO
C-------------------------Loop for periods---------------------------------------E
             enddO
C-------------------------Loop for modes---------------------------------------E
888          PRINT*,'ALL IS DONE FOR ',mode_end,' MODES'
             END
