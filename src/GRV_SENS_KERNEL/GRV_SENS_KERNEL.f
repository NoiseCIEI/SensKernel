             include "version.h"
             character*80 inname,bred
             character*3 symb
             character*1 symo
             character*1 wavetype,mode
             integer nper(10)
             real*4 period(10,100)
             data marg/3/,pi2/6.2831853/
C--domega= log(T2/T0)---difference in periods between central and
C--------- leftside periods/central period
C---to run this program one should make:
C---three calculations of Rayleigh/Love waves with shifted periods
C---three calculations of phase sensitivity kernels with shifted periods
C---leftside period files with +
C---rightside periodfiles  with -
C---central period file no signs (no + or -)
C------------------------INIT-----------------------S
             narg=iargc()
             print*,'to find R/L normalized group velocity derivatives by Vs,Vp,rho'
             if(narg.ne.marg.and.narg.ne.marg+1)STOP'USAGE: GRV_SENS_KERNEL inname mode_end wavetype'
             write(*,'(20(1H=),1x,a,1x,20(1H=))') VERSION
             call GETARG(1,inname)
             lin=lnblnk(inname)
             call GETARG(2,symo)
             read(symo,'(I1)')mode_end
             call GETARG(3,wavetype)
             mode=' '
             mode_curr=0
C------------------------INIT-----------------------E
C-----to determine requested modes and periods from .phv file--------S
             open(22,file=inname(1:lin)//'.'//wavetype//'.phv',status='OLD')
444         DO ip=1,1000
      read(22,'(a)',end=99999)bred
      if(bred(5:10).eq.'    ')then
      read(22,'(1X)',end=99999)
      go to 555
                              endif
      read(bred,*)period(mode_curr+1,ip)
      enddo
555   mode_curr=mode_curr+1
      nper(mode_curr)=ip-1
      go to 444
99999 close(22)
C-----to determine requested modes and periods from .phv file--------E
C------------Loop for modes---------------------------S
                   dO imode=1,mode_end+1
             mode_req=imode-1
             write(mode,'(I1)')mode_req
C-------------period    loop------S
             DO IP=1,nper(imode)
             period_req=period(imode,ip)
             symb='   '
             iper=int(period_req)
             if(iper.lt.10)write(symb(1:1),'(I1)')iper
             if(iper.ge.10.and.iper.lt.100)write(symb(1:2),'(I2)')iper
             if(iper.ge.100)write(symb(1:3),'(I3)')iper
             lsy=lnblnk(symb)
             open(1,file=inname(1:lin)//'-.phv.'//wavetype//'_'//mode//'_'//symb(1:lsy),status='OLD')
             open(2,file=inname(1:lin)//'.phv.'//wavetype//'_'//mode//'_'//symb(1:lsy),status='OLD')
             open(3,file=inname(1:lin)//'+.phv.'//wavetype//'_'//mode//'_'//symb(1:lsy),status='OLD')
             open(4,file=inname(1:lin)//'.grv.'//wavetype//'_'//mode//'_'//symb(1:lsy))
C------------loop for depth-Freche derivatives-----S
C-----reading phase derivatives for given mode and period----S
             do i=1,2500
             if(wavetype.eq.'R')THEN
             if(i.eq.1)then
             read(2,'(F10.4,6E15.5,I2)',end=99) depth,dcdb,dcda,dcdrho,T0,c,u,mod
             read(1,'(F10.4,6E15.5,I2)',end=99) depth,dcdb_min,dcda_min,dcdrho_min,T1,c_m,u_m,mod
             read(3,'(F10.4,6E15.5,I2)',end=99) depth,dcdb_plu,dcda_plu,dcdrho_plu,T2,c_b,u_b,mod
                        else
             read(2,'(F10.4,3E15.5)',end=99) depth,dcdb,dcda,dcdrho
             read(1,'(F10.4,3E15.5)',end=99) depth,dcdb_min,dcda_min,dcdrho_min
             read(3,'(F10.4,3E15.5)',end=99) depth,dcdb_plu,dcda_plu,dcdrho_plu
                        endif
                              ELSE
             if(i.eq.1)then
             read(2,'(F10.4,5E15.5,I2)',end=99) depth,dcdb,dcdrho,T0,c,u,mod
             read(1,'(F10.4,5E15.5,I2)',end=99) depth,dcdb_min,dcdrho_min,T1,c_m,u_m,mod
             read(3,'(F10.4,5E15.5,I2)',end=99) depth,dcdb_plu,dcdrho_plu,T2,c_b,u_b,mod
                       else
             read(2,'(F10.4,2E15.5)',end=99) depth,dcdb,dcdrho
             read(1,'(F10.4,2E15.5)',end=99) depth,dcdb_min,dcdrho_min
             read(3,'(F10.4,2E15.5)',end=99) depth,dcdb_plu,dcdrho_plu
                      endif
                             ENDIF
             dcdb=dcdb*c
             dcdrho=dcdrho*c
             dcdb_plu=dcdb_plu*c_b
             dcdrho_plu=dcdrho_plu*c_b
             dcdb_min=dcdb_min*c_m
             dcdrho_min=dcdrho_min*c_m
                if(wavetype.eq.'R')then
             dcda=dcda*c
             dcda_plu=dcda_plu*c_b
             dcda_min=dcda_min*c_m
                                    endif
C-----reading phase derivatives for given mode and period----E
             
             if(i.eq.1)then
             domega=log(T2/T0)
             c0=(c_m+c_b)/2.
             u0=(u_m+u_b)/2.
             u_c=u/c
                        endif
C-----------calculating group velocity derivatives------S
      dudb=(u_c/2.*(2.-u_c)*(dcdb_plu+dcdb_min) -0.5*u_c**2*(dcdb_plu-dcdb_min)/domega)/u
      dudrho=(u_c/2.*(2.-u_c)*(dcdrho_plu+dcdrho_min)+0.5*u_c**2*(dcdrho_plu-dcdrho_min)/domega)/u
      if(wavetype.eq.'R') duda=(u_c/2.*(2.-u_c)*(dcda_plu+dcda_min)-0.5*u_c**2*(dcda_plu-dcda_min)/domega)/u
C-----------calculating group velocity derivatives------E
C-----------writing output-----S
      if(i.eq.1.and.wavetype.eq.'R')write(4,'(F10.4,6E15.5,I2)')depth,dudb,duda,dudrho,T0,c0,u0,mode_req
      if(i.eq.1.and.wavetype.eq.'L')write(4,'(F10.4,5E15.5,I2)')depth,dudb,dudrho,T0,c0,u0,mode_req
      if(i.gt.1.and.wavetype.eq.'R')write(4,'(F10.4,3E15.5)')depth,dudb,duda,dudrho
      if(i.gt.1.and.wavetype.eq.'L')write(4,'(F10.4,2E15.5)')depth,dudb,dudrho
C-----------writing output-----S
             enddo
99           continue
C------------loop for depth-Freche derivatives-----E
             close(1)
             close(2)
             close(3)
             close(4)
             ENDDO
C------------Loop for periods---------------------------E
             enddO
C------------Loop for modes---------------------------E
             PRint*,'ALL IS DONE FOR ',IP-1,' PERIODS and ',mode_end+1,' modes'
             STOP
             END
