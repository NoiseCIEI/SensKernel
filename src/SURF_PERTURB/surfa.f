czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        subroutine NEVILL (t,c1,c2,del1,del2,ifunc,cc)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        implicit none
        integer*4 nsize,ifunc
        real*8    t,c1,c2,del1,del2,cc
        parameter (nsize=1000)
c ---
        real*8     a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  a,b,rho,d,qs
c ---
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        real*8    t1,dt,cc1,dc,cinit
        integer*4 lstop,iqq,istru
        common/newly/t1,dt,cc1,dc,cinit,lstop,iqq,istru
        real*8     x(20),y(20)
c--------------------------------------------------------------------
        integer*4  ic,j,m,nev,kk
        real*8     accur1,accur2,c3,del3,s1,s13,ss2,ss1,s2,s32,dltar
        accur1=0.1d-5
        accur2=0.1d-7
        m=0
        ic=0
        c3=(c1+c2)/2.0d0
        del3=dltar(c3,t,ifunc)
        nev=1
1310    ic=ic+1
             if(ic.lt.5000)           goto 777
c------------------------------------------------------
c         TOO  MANY  CYCLES        ????? STOP ??????
c------------------------------------------------------
c                                                       STOP
        lstop=lstop+10
        write(*,*) 'LSTOP=', lstop
                   goto 2000
c------------------------------------------------------
c
777               continue
cMB             if(c1-c3) 1320,1320,1330
cMB1320         if(c2-c3) 1344,1344,1000
cMB1330         if(c2-c3) 1000,1344,1344
             if(c1.gt.c3) goto 1330
             if(c2.le.c3) goto 1344
             if(c2.gt.c3) goto 1000
1330         if(c2.ge.c3) goto 1344
1000    s13=del1-del3
        s32=del3-del2
cMB          if(sign(1.d0,del3)*sign(1.d0,del1)) 1441,1441,1443
             if(sign(1.d0,del3)*sign(1.d0,del1).gt.0.0d0) goto 1443
        c2=c3
        del2=del3
                  goto 1444
1443    c1=c3
        del1=del3
1444              continue
cMB             if(dabs(c1-c2)-accur1) 20,20,22
             if((dabs(c1-c2)-accur1).le.0.0d0) goto 20
cxx 22                continue
             if(sign(1.d0,s13).ne.sign(1.d0,s32)) nev=0
        ss1=dabs(del1)
        s1=0.1d0*ss1
        ss2=dabs(del2)
        s2=0.1d0*ss2
             if(s1.gt.ss2.or.s2.gt.ss1)           goto 1344
             if(nev.eq.0)           goto 1344
             if(nev.eq.2)           goto 1350
        x(1)=c1
        y(1)=del1
        x(2)=c2
        y(2)=del2
        m=1
                  goto 1355
1344    c3=(c1+c2)/2.0d0
        del3=dltar(c3,t,ifunc)
        nev=1
        m=1
                  goto 1310
1350    x(m+1)=c3
        y(m+1)=del3
1355              do 1360 kk=1,m
        j=m-kk+1
             if(dabs(y(m+1)-y(j)).le.accur2)          goto 1344
        x(j)=(-y(j)*x(j+1)+y(m+1)*x(j))/(y(m+1)-y(j))
1360              continue
cxx 21                continue
        c3=x(1)
        del3=dltar(c3,t,ifunc)
        nev=2
        m=m+1
             if(m.gt.10) m=10
                  goto 1310
20                continue
        cc=c3
2000              continue
                                                        RETURN
                                                        END
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        function DLTAR (cc,tt,kk)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        implicit none
        integer*4 nsize,nper,kk
        real*8    cc,tt
        parameter (nsize=1000,nper=1000)
c ---
        real*8    a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/ a,b,rho,d,qs
c ---
        real*8    c(nper,10),t(nper),ratio(nper,10)
        common/o/ c,t,ratio
c ---
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        integer*4 ii
        real*8    dlt_old,dmax,sum,dltar,dltar1,dltar4
c-----------------------------------------------------------------------
cMB          if(idrop) 899,899,905
             if(idrop.gt.0) goto 905
cxx 899               continue
        dmax=fact*cc*tt
        mmax=nmax
        sum=0
                  do 900 ii=1,nmax
cMb          if(cc-b(ii)) 901,900,900
             if(cc-b(ii).ge.0.0d0) goto 900
        sum=sum+d(ii)
cMB          if(sum-dmax) 900,900,902
             if(sum-dmax.le.0.0d0) goto 900
        mmax=ii
             goto 904
900               continue
904     idrop=1
             if(mmax.lt.2) mmax=2
905               continue
             goto (1,2,3,4,5),kk
c---------------------------------------------
c       LOVE  WAVE  PERIOD  equation
c------------------
1       dltar=dltar1(cc,tt,1)
                                                        RETURN
c---------------------------------------------
c       RAYLEIGH  WAVE  PERIOD  equation
c------------------
2       dltar=dltar4(cc,tt,1)
        dlt_old=dltar
                                                        RETURN
c---------------------------------------------
c       RAYLEIGH  WAVE  ellipticity
c------------------
3       dltar=dltar4(cc,tt,2)
                                                        RETURN
c---------------------------------------------
c       RAYLEIGH  WAVE  amplitude  response  component
c------------------
4       dltar=dltar4(cc,tt,3)
                                                        RETURN
c---------------------------------------------
c       LOVE  WAVE  amplitude  response  component
c------------------
5       dltar=dltar1(cc,tt,2)
                                                        RETURN
                                                        END
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        function DLTAR1 (c,t,mup)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c-----------------------------------------------------------------------
c       Haskell-Thompson Love wave formulation from halfspace to surface
c-----------------------------------------------------------------------
        implicit none
        integer*4 nsize,mup
        real*8   c,t
        parameter (nsize=1000)
c ---
        real*8     a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  a,b,rho,d,qs
c ---
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        integer*4  k,m,mmm1
        real*8    ett,wvno,covb,h,exqm,eut,cosq,exqp,tt,sinq,q
        real*8    rb,z,ut,y,dltar1
c ---
        ett=0.0d0
        wvno=6.2831853d0/(c*t)
        covb=c/b(mmax)
        h=rho(mmax)*b(mmax)*b(mmax)
        rb=dsqrt(dabs(covb**2-1.0d0))
        ut=1.0d0
        tt=h*rb
        mmm1=mmax-1
                                do 1340 k=1,mmm1
        m=mmax-k
                        if(b(m).eq.0.0d0)                 go to 1340
        covb=c/b(m)
        rb=dsqrt(dabs(covb**2-1.0d0))
        h=rho(m)*b(m)*b(m)
        q=-wvno*d(m)*rb
                        if(rb.lt.0.1d-20)               go to 1221
cMB                     if(c-b(m)) 1209,1221,1231
                        if(c-b(m).lt.0.0d0) goto 1209
                        if(c-b(m).eq.0.0d0) goto 1221
        sinq=dsin(q)
        y=sinq/rb
        z=rb*sinq
        cosq=dcos(q)
                                                        go to 1242
1221    y=-wvno*d(m)
        z=0
        cosq=1
                                                        go to 1242
1209    exqp=dexp(q)
        exqm=1.0d0/exqp
        y=(exqp-exqm)/(2.0d0*rb)
        z=-rb*rb*y
        cosq=(exqp+exqm)/2.0d0
1242    eut=cosq*ut-y*tt/h
        ett=h*z*ut+cosq*tt
        ut=eut
        tt=ett
1340                                    continue
                                                        go to (1,2),mup
1       dltar1=-ett
                                                                RETURN
2       dltar1=ut
                                                                RETURN
                                                                END
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        function DLTAR4 (c,t,mup)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        implicit none
        integer*4 nsize,mup
        real*8    c,t
        parameter (nsize=1000)
c ---
        real*8     p(nsize),s(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  p,s,rho,d,qs
c ---
        integer*4 nmax,mmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/fact,nmax,mmax,kmax,idrop,iedit,ndiv,mode
c ---
        integer*4 jump,m
        real*8    a11,a12,a13,a14,a15,a21,a22,a23,a24,a31,a32,a33
        real*8    a41,a42,a51,accur,accurs,arga,argb,b1,b2,b3,b4
        real*8    b5,bb1,bb2,bb3,bb4,bb5,cc,ccm,cosp,cosq,csq,dltar4
        real*8    g,g1,g1s,gg1,gm,gra,gs,pm,pp,ppp,qm,r12,ra,rad,rb
        real*8    rba,rhoc,rhocs,rhp,rr,rs1,rs2,rs3,rs4,rsinp,rsinq
        real*8    sinpr,sinqr,ss,sss,suu,wvno
c-------------------------------------------------------------------
        ra=0.0d0
        rb=0.0d0
        g=0.0d0
        g1=0.0d0
        r12=0.0d0
        accur=1.d-8
        accurs=1.d-8
        wvno=6.28318531/(c*t)
        csq=c*c
        jump=1
             if(mup.gt.1) jump=2
1000    b1=0.0d0
        b2=0.0d0
        b3=0.0d0
        b4=0.0d0
        b5=0.0d0
cMB          if(jump-2) 1,2,3
             if(jump.gt.2) goto 3
             if(jump.eq.2) goto 2
        b1=1.0d0
             goto 4
2       b2=1.0d0
             goto 4
3       b3=1.0d0
             goto 4
4                 continue
                  do 50 m=1,mmax
        arga=1.0d0-csq/p(m)**2
        ra=dsqrt(dabs(arga))
             if(arga.gt.0.0d0) ra=-ra
             if(dabs(s(m)).gt.accurs)      goto 101
c---------------------------------------------
c       LIQUID  SURFACE  LAYER
c-------------------------------
        pm=wvno*ra*d(m)
             if(mup.gt.1)      goto 50
        rhoc=rho(m)*csq
             if(dabs(ra).lt.accur)      goto 312
cMB             if(ra) 313,312,314
             if(ra.gt.0.0d0) goto 314
             if(ra.lt.0.0d0) goto 313
312     sinpr=wvno*d(m)
        rsinp=0.0d0
        cosp=1.0d0
             goto 315
313     sinpr=(dexp(pm)-dexp(-pm))/(2.0d0*ra)
        rsinp=-ra*ra*sinpr
        cosp=0.5d0*(dexp(pm)+dexp(-pm))
             goto 315
314     sinpr=dsin(pm)/ra
        rsinp=ra*dsin(pm)
        cosp=dcos(pm)
315               continue
        a11=cosp
        a21=rhoc*sinpr
        a31=0.0d0
        a41=0.0d0
        a51=0.0d0
        a12=0.0d0
        a22=0.0d0
        a32=0.0d0
        a42=0.0d0
        a13=0.0d0
        a23=0.0d0
        a33=0.0d0
        a14=0.0d0
        a24=0.0d0
        a15=0.0d0
             goto 1001
c-------------------------------------------------
101     argb=1.0d0-csq/s(m)**2
        rb=dsqrt(dabs(argb))
             if(argb.gt.0.0d0) rb=-rb
        g=2.0d0*s(m)**2/csq
        g1=g-1.0d0
cMB          if(mmax-m) 40,52,40
             if(mmax.eq.m) goto 52
        rhoc=rho(m)*csq
        pm=wvno*ra*d(m)
        qm=wvno*rb*d(m)
cMB          if(ra) 213,212,214
             if(ra.gt.0.0d0) goto 214
             if(ra.lt.0.0d0) goto 213
        rsinp=0.0
        sinpr=wvno*d(m)
        cosp=1.0d0
             goto 215
213     rsinp=-ra*0.5d0*(dexp(pm)-dexp(-pm))
        sinpr=-rsinp/(ra**2)
        cosp=0.5d0*(dexp(pm)+dexp(-pm))
             goto 215
214     rsinp=ra*dsin(pm)
        sinpr=rsinp/(ra**2)
        cosp=dcos(pm)
215                 continue
             if(dabs(rb).lt.accur)      goto 218
             if(rb.gt.0.0d0)      goto 217
        rsinq=-rb*0.5d0*(dexp(qm)-dexp(-qm))
        sinqr=-rsinq/(rb**2)
        cosq=0.5d0*(dexp(qm)+dexp(-qm))
             goto 219
217     rsinq=rb*dsin(qm)
        sinqr=rsinq/(rb**2)
        cosq=dcos(qm)
             goto 219
218     rsinq=0.0d0
        sinqr=wvno*d(m)
        cosq=1.0d0
             goto 219
219     rr=rsinp*rsinq
        ss=sinpr*sinqr
        cc=cosp*cosq
        rs1=rsinp*cosq
        rs2=sinqr*cosp
        rs3=sinpr*cosq
        rs4=rsinq*cosp
        gm=2.0d0*g-1.0d0
        gs=g*g
        g1s=g1*g1
        ccm=1.0d0-cc
        gg1=g*g1
        rhocs=rhoc*rhoc
        suu=gs*rr+g1s*ss
        a11=2.0d0*gs-gm
        a11=a11*cc-suu-2.0d0*gg1
        a12=-(rs1+rs2)/rhoc
        a13=gm*ccm+g1*ss+g*rr
        a13=-2.0d0*a13/rhoc
        a14=(rs3+rs4)/rhoc
        a15=2.0d0*ccm+rr+ss
        a15=a15/rhocs
        a21=rhoc*(g1s*rs3+gs*rs4)
        a22=cc
        a23=2.0d0*(g*rs4+g1*rs3)
        a24=sinpr*rsinq
        a31=rhoc*(gg1*gm*ccm+g1s*g1*ss+gs*g*rr)
        a32=g1*rs2+g*rs1
        a33=1.0d0+2.0d0*(2.0d0*gg1*ccm+suu)
        a41=-rhoc*(g1s*rs2+gs*rs1)
        a42=rsinp*sinqr
        a51=rhocs*(2.0d0*gs*g1s*ccm+gs*gs*rr+g1s*g1s*ss)
1001              continue
c--------------------------------------------------
c         MATRIX  MULTIPLICATION
c               effect  of  immaginary  elements  included
c---------------------------
        bb1=a11*b1+a12*b2+a13*b3+a14*b4+a15*b5
        bb2=a21*b1+a22*b2+a23*b3+a24*b4-a14*b5
        bb3=a31*b1+a32*b2+a33*b3-0.5*a23*b4+0.5*a13*b5
        bb4=a41*b1+a42*b2-2.*a32*b3+a22*b4-a12*b5
        bb5=a51*b1-a41*b2+2.*a31*b3-a21*b4+a11*b5
        b1=bb1
        b2=bb2
        b3=bb3
        b4=bb4
        b5=bb5
50                continue
c--------------------------------------------------
c         the  FOLLOWING  EXPRESSION  is  VALID   rb=0
c------------------------------
52                continue
        pp=p(m)
        sss=s(m)**2
        ppp=pp**2
        rhp=rho(m)*pp
        gra=g*ra
        g1s=g1*g1
        rba=rb-1.0d0/ra
        a11=-2.0d0*rb*sss/ppp+csq*g1s/ppp/gra
        a12=rhp*pp
        a13=-rb/a12+g1/a12/gra
        a14=rb/a12/gra
        a15=rba/rhp/rhp/csq/g
        a12=-1.0d0/g/a12
        bb1=a11*b1+a12*b2+2.*a13*b3+a14*b4+a15*b5
             goto (501,502,503),mup
501     dltar4=-bb1
                                                        RETURN
502          if(jump.eq.2) r12=bb1
        jump=jump+1
             if(jump.eq.3)      goto 1000
        dltar4=0.5d0*bb1/r12
                                                        RETURN
503     dltar4=dabs(bb1)
             if(dabs(s(1)).gt.accurs)                    RETURN
        ra=c/p(1)
        rad=wvno*d(1)*dsqrt(dabs(ra*ra-1.0d0))
        dltar4=dabs(bb1*dcos(rad))
                                                        RETURN
                                                        END
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        subroutine LEIGEN(num,ixa,jxa)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c----------------------------------------------------------------------
c       num = 0 - the first call, num > 0 - following calls;
c       ixa - index of period; jxa - index of mode;
c       jmax < 0 or =0 - end of program, jmax > 0 - number of layers;
c       d - layer thickness; a - Pvel; b - Svel; rho - density;
c       depth - depth of the midpoint of layer; xmu,xlamb - Lame const;
c       t - period; c - phase vel; cvar, ugr - phase and group vel 
c       obtained from integrals; wvnum - wave number.
c----------------------------------------------------------------------
        implicit none
        integer*4 nsize,nper,num,ixa,jxa
        parameter (nsize=1000,nper=1000)
c ---
        real*8     a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/  a,b,rho,d,qs
c ---
        real*8    cc(nper,10),tt(nper),rat(nper,10)
        common/o/ cc,tt,rat
c ---
        integer*4 nmax,jmax,kmax,idrop,iedit,ndiv,mode
        real*8    any
        common/c/ any,nmax,jmax,kmax,idrop,iedit,ndiv,mode
c ---
        real*8    dept1(nsize),amp(nsize),ampuz(nsize),stresz(nsize),
     *            stress(nsize)
        integer*4 mmax
        common/rar/ dept1,amp,ampuz,stresz,stress,mmax
c ---
        real*8 dcda(nsize),dcdb(nsize),dcdr(nsize),dwx(nsize),
     *         g(nsize),g2(nsize)
        common/rar1/ dcda,dcdb,dcdr,dwx,g,g2
c ---
        real*8    c,cvar,ugr,wvno,ratio,ale
        common/rco/ c,cvar,ugr,wvno,ratio,ale
c ---
        real*8    sumi0,sumi1,sumi2,sumi3,flagr
        common/rco1/ sumi0,sumi1,sumi2,sumi3,flagr
c ---
        integer*4 nderiv,ndpth
        real*8    dpth(nsize),dspz(nsize),dsp(nsize),drvz(5,nsize),drv(5,nsize)
        common/deriv/ nderiv,ndpth,dpth,dspz,dsp,drvz,drv
c ---
        logical KEY_ATTEN,KEY_FLAT,KEY_DERIV
        common/log/ KEY_ATTEN,KEY_FLAT,KEY_DERIV
c ---
        real*8    dmm(5),smm(5)
        real*8    depth(nsize)
        real*8    xmu(nsize)
        real*8    pr(5*nsize)
c ---
        integer*4 i,ii,jj,mm1,k,kk,maxu,j,ivre,l,ldiv,m,mmm1,nderv
        real*8    const_lim,const_lim1,xxmin,base,ampl,cosq,covb,div
        real*8    dldk,dldm,dmax,dm,dldr,ett,dz,eut,expq,exqm,exqp
        real*8    h,rab2,omegsq,q,z,rab1,rab3,ut,sm,sinq,rb,stres,t
        real*8    sum,tq,wvnum,ut0,uttq,wvar,wvnosq,xkk,y
c ---
        data const_lim/1.d+10/,const_lim1/1.d+5/
        xxmin=1.0d-20
        maxu=0
c
c
                        if(num.gt.0.and..not.KEY_ATTEN)     go to 503
        mmax=jmax
        nmax=mmax
cMB                     if(mmax) 777,777,778
                        if(mmax.le.0) goto 777
cxx 778                                     continue
        base=0.0d0
        mm1=mmax-1
        ivre=(nsize-1)/mm1
                        if(ndiv.gt.ivre) ndiv=ivre
        div=float(ndiv)
                        if(ndiv.le.1)                   go to 20003
        jj=1
                        if(b(1).le.0.1d-10) jj=2
                                do 10001 j=jj,mm1
        ldiv=(j-jj)*ndiv
                                do 10001 i=1,ndiv
        pr(ldiv+i)=d(j)/div
        pr(ldiv+i+nsize)=b(j)
        pr(ldiv+i+2*nsize)=rho(j)
        if(KEY_ATTEN)pr(ldiv+i+3*nsize)=qs(j)
10001                                   continue
        mmax=(mm1-jj+1)*ndiv+jj
        d(mmax)=0.0d0
        a(mmax)=a(nmax)
        b(mmax)=b(nmax)
        rho(mmax)=rho(nmax)
        if(KEY_ATTEN)qs(mmax)=qs(nmax)
        nmax=mmax
        mm1=mmax-1
                                do 10002 j=jj,mm1
        d(j)=pr(j-jj+1)
        b(j)=pr(j-jj+nsize+1)
        rho(j)=pr(j-jj+2*nsize+1)
        if(KEY_ATTEN)qs(j)=pr(j-jj+3*nsize+1)
10002                                   continue
20003                                   continue
c
c
c
c
c
                                do 22 i=1,mmax
        base=base+d(i)
        depth(i)=base-d(i)*0.5d0
        xmu(i)=rho(i)*b(i)*b(i)
cMB                     if(i-mmax)22,24,24
                        if(i.gt.mmax) goto 24
22                                      continue
24      depth(mmax)=base-d(mmax)
c
c
503                                     continue
                                do 500 j=1,nsize
        amp(j)=0.0d0
        stress(j)=0.0d0
        dcdr(j)=0.0d0
        dcdb(j)=0.0d0
500     g(j)=0.0d0
        t=tt(ixa)
        c=cc(ixa,jxa)
        mmax=nmax
c-----------------------------------------------------------------------
c       Layer dropping procedure; any is the same as fact.
c-----------------------------------------------------------------------
                        if(any.le.0.0d0) any=7.0d0
        dmax=any*c*t
        sum=0.0d0
                                do 900 ii=1,mmax
        maxu=ii
cMB                     if(c-b(ii)) 901,900,900
                        if(c.ge.b(ii)) goto 900
        sum=sum+d(ii)
                        if(ii.eq.mmax)                  go to 902
                        if(sum.le.dmax)                 go to 900
cMB                     if(b(ii+1)-b(ii)) 902,900,90009
                        if(b(ii+1)-b(ii).lt.0.0d0) goto 902
                        if(b(ii+1)-b(ii).gt.0.0d0) goto 902
900                                     continue
        maxu=maxu+1
902     mmax=maxu
        wvno=6.2831853d0/(c*t)
        wvnum=wvno
        wvnosq=wvno*wvno
        omegsq=(6.2831853d0/t)**2
C-------the halfspace is handled-------------------------
        ut0=1.0d0
7777    ut=ut0
        covb=c/b(mmax)
        h=rho(mmax)*b(mmax)*b(mmax)
        rb=wvno*dsqrt(abs(covb**2-1.0d0))
        tq=-h*rb*ut0
        amp(mmax)=ut
        stress(mmax)=tq
        g(mmax)=tq/xmu(mmax)
cMB                     if(rb) 1001,1000,1001
                        if(rb.gt.0.0d0) goto 1001
        dm=1.0e25
        sm=0.0d0
                                                        go to 1002
1001    dm=0.5d0/rb
        sm=0.5d0*rb
1002                                    continue
        dldr=omegsq*dm
        dldm=-(wvnosq*dm+sm)
        dcdb(mmax)=2.0d0*rho(mmax)*b(mmax)*c*dldm/wvno
        dcdr(mmax)=(c/wvno)*(dldr+b(mmax)*b(mmax)*dldm)
        mmm1=mmax-1
        sumi0=rho(mmax)*dm
        sumi1=h*dm
        sumi2=h*sm
C---------------loop for knots-----------S
                                do 1340 k=1,mmm1
         IF(abs(ut).gt.const_lim) THEN
         ut0=ut0/const_lim1
         goto 7777 
                                  ENDIF
        m=mmax-k
                        if(b(m).eq.0.0)                 go to 1340
        covb=c/b(m)
        rb=wvno*dsqrt(abs(covb**2-1.0d0))
        h=rho(m)*b(m)*b(m)
        dz=d(m)/4.0d0
        dmm(1)=ut*ut
        smm(1)=(tq/h)**2
                                do 1339 kk=2,5
        xkk=kk-1
        q=rb*dz*xkk
cMB                     if(c-b(m)) 1207,1221,1231
                        if(c-b(m).lt.0.0d0) goto 1207
                        if(c-b(m).eq.0.0d0) goto 1221
        sinq=dsin(q)
        y=sinq/rb
        z=-rb*sinq
        cosq=dcos(q)
                                                        go to 1242
1221    y=dz*xkk
        z=0.0d0
        cosq=1.0d0
                                                        go to 1242
1207    exqp=dexp(q)
        exqm=1.0d0/exqp
        y=(exqp-exqm)/(2.0d0*rb)
        z=rb*rb*y
        cosq=(exqp+exqm)/2.0d0
1242    eut=cosq*ut-y*tq/h
        ett=-h*z*ut+cosq*tq
        dmm(kk)=eut*eut
        smm(kk)=(ett*ett)/(h*h)
cMB                     if(kk-3) 1339,1301,1339
                        if(kk.ne.3) goto 1339
        amp(m)=eut
        stress(m)=ett
        g(m)=ett/xmu(m)
1339                                    continue
        ut=eut
        tq=ett
        dm=(dz/22.5d0)*(7.0d0*(dmm(1)+dmm(5))+32.0d0*(dmm(2)+dmm(4))+12.0d0*dmm(3))
        sm=(dz/22.5d0)*(7.0d0*(smm(1)+smm(5))+32.0d0*(smm(2)+smm(4))+12.0d0*smm(3))
        dldm=-(wvnosq*dm+sm)
        dldr=omegsq*dm
        dcdb(m)=2.0d0*rho(m)*b(m)*c*dldm/wvno
        dcdr(m)=(c/wvno)*(dldr+b(m)*b(m)*dldm)
        sumi0=sumi0+rho(m)*dm
        sumi1=sumi1+h*dm
        sumi2=sumi2+h*sm
1340                                    continue
C---------------loop for knots-----------E
        amp(nsize)=1.0d0
        stress(nsize)=tq/ut
        dcdb(nsize)=0.0d0
        dcdr(nsize)=0.0d0
        depth(nsize)=0.0d0
                        if(b(1).eq.0.0d0) depth(nsize)=d(1)
        dldk=-2.0d0*wvno*sumi1
                                do 1500 l=1,mmax
                        if(b(l).eq.0.0d0)                 go to 1500
        amp(l)=amp(l)/ut
        stress(l)=stress(l)/ut
        dcdb(l)=dcdb(l)/dldk
        dcdr(l)=dcdr(l)/dldk
        g(l)=stress(l)/xmu(l)
c--------------------------------------------------------
c       Exclusion of low amplitudes from final output
c--------------------------------------------------------
cMB                     if(b(l)-b(mmax)) 1506,1507,1507
                        if(b(l)-b(mmax).lt.0.0d0) goto 1506
                                        continue
cMB                     if(abs(amp(l))-xxmin) 1501,1502,1502
                        if(abs(amp(l))-xxmin.gt.0.0d0) goto 1502
        amp(l)=0.0d0
        stress(l)=0.0d0
        dcdb(l)=0.0d0
        dcdr(l)=0.0d0
        g(l)=0.0d0
                                                        go to 1500
1502                                    continue
1506                                    continue
1500                                    continue
        sumi0=sumi0/(ut*ut)
        sumi1=sumi1/(ut*ut)
        sumi2=sumi2/(ut*ut)
        wvar=(omegsq*sumi0-sumi2)/sumi1
        cvar=dsqrt(omegsq/wvar)
        wvar=dsqrt(wvar)
        ugr=sumi1/(c*sumi0)
        flagr=omegsq*sumi0-wvnosq*sumi1-sumi2
        ale=1.0d0/(2.0d0*c*ugr*sumi0)/dsqrt(6.28318d0)*1.d-15
                        if(b(1).le.0.0d0)                 go to 5556
                                do 5555 l=1,mmax
        i=mmax-l+1
        j=i+1
        dept1(j)=depth(i)
        amp(j)=amp(i)
        stress(j)=stress(i)
        dcdb(j)=dcdb(i)
        dcdr(j)=dcdr(i)
5555    g(j)=g(i)
        mmax=mmax+1
                                                        go to 501
5556                            do 5557 l=1,mmax
5557    dept1(l)=depth(l)
501     dept1(1)=0.0d0
                        if(b(1).le.0.0) dept1(1)=d(1)
        amp(1)=1.0d0
        stress(1)=0.0d0
        dcdb(1)=0.0d0
        dcdr(1)=0.0d0
        g(1)=0.0d0
c----------------------------------------------------------------------
c       Calculation of eigenfunction and its derivatives
c----------------------------------------------------------------------
                        if(iedit.ne.2)                  go to 800
                                do 799 i=1,5
                                do 799 j=1,nsize
        dsp(j)=0.0d0
        drv(i,j)=0.0d0
799                                     continue
                                do 801 l=1,ndpth
                        if(dpth(l).lt.dept1(1))         go to 801
c----------------------------------------------------------------------
c       If dpth(l) is above the free surface, the output will be zero
c----------------------------------------------------------------------
                                do 802 ii=2,mmax
        rab1=dept1(ii)-dpth(l)
                        if(abs(rab1).lt.1.0d-5) rab1=0.0d0
                        if(rab1.le.0.0d0)                  go to 802
        i=ii-1
        j=ii
                        if(b(1).le.0.0d0) i=ii
        rab2=d(i)/2.0d0
        rab1=rab1-rab2
                        if(abs(rab1).lt.1.0d-5) rab1=0.0d0
                        if(rab1.le.0.0d0)                  go to 803
        j=j-1
        i=i-1
                                                        go to 803
802                                     continue
        j=mmax
        i=mmax-1
                        if(b(1).le.0.0d0) i=j
803     dz=dpth(l)-dept1(j)
        ampl=amp(j)
        stres=stress(j)
        covb=c/b(i)
        rb=wvno*dsqrt(abs(covb**2-1.0d0))
cMB                        if(c-b(i)) 810,811,812
                        if(c-b(i).gt.0.0d0) goto 812
                        if(c-b(i).eq.0.0d0) goto 811
        nderv=nderiv+1
                                do 804 ii=1,nderv
        k=ii-1
        rab1=rb**k
        q=rb*dz
        expq=dexp(q)
        uttq=(ampl+stres/(xmu(i)*rb))/ampl
        rab2=rab1*expq*uttq*ampl/2.0d0
        uttq=abs(uttq)
                        if(uttq.lt.1.0d-5) rab2=0.0d0
        rab3=(-1)**k*rab1*(ampl/2.0d0-stres/(2.0d0*xmu(i)*rb))/expq
                        if(k.ne.0)                      go to 805
        dsp(l)=rab2+rab3
                                                        go to 804
805     drv(k,l)=rab2+rab3
804                                     continue
                                                        go to 801
811     dsp(l)=ampl+stres*dz/xmu(i)
                        if(nderiv.lt.1)                 go to 801
        drv(1,l)=stres/xmu(i)
                                                        go to 801
812     nderv=nderiv+1
                                do 806 ii=1,nderv
        k=ii-1
        rab1=rb**k
        q=rb*dz
        q=q+k*3.1415926d0/2.0d0
        cosq=dcos(q)
        sinq=dsin(q)
        rab2=rab1*ampl*cosq
        rab3=stres*rab1*sinq/(xmu(i)*rb)
                        if(k.ne.0)                      go to 807
        dsp(l)=rab2+rab3
                                                        go to 806
807     drv(k,l)=rab2+rab3
806                                     continue
801                                     continue
800                                                             RETURN
777                                     continue
                                                                STOP
                                                                END
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        subroutine REIGEN (num,ixa,jxa)
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        implicit none
        integer*4 nsize,nper,num,ixa,jxa
        parameter (nsize=1000,nper=1000)
c ---
        real*8 aur1,auz1,
     1  atz1,atr1,aur2,auz2,atz2,atr2,aa,bb,daur1,dauz1,datz1,datr1,
     2  daur2,dauz2,datz2,datr2,eaur1,eauz1,eatz1,eatr1,eaur2,eauz2,
     3  eatz2,eatr2,saur1,sauz1,satz1,satr1,saur2,sauz2,satz2,satr2,
     4  xnorm,dmmr,dmmz,smmz,smmr,drsz,dzsr,dldr,dldm,dldl,dldk
        complex*16 cra,cr1,cr2,cr3
c ---
        real*8    a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/ a,b,rho,d,qs
c ---
        real*8    cc(nper,10),tt(nper),rat(nper,10)
        common/o/ cc,tt,rat
c ---
        integer*4 nmax,jmax,kmax,idrop,iedit,ndiv,mode
        real*8    fact
        common/c/   fact,nmax,jmax,kmax,idrop,iedit,ndiv,mode 
c ---
        real*8    te1,dete,ce1,dece,cinit
        integer*4 lstop,iqq,istru
        common/newly/ te1,dete,ce1,dece,cinit,lstop,iqq,istru
c ---
        real*8     dept1(nsize),ampur(nsize),ampuz(nsize),stresz(nsize),
     *             stresr(nsize)
        integer*4   mmax
        common/rar/ dept1,ampur,ampuz,stresz,stresr,mmax
c ---
        real*8    dcda(nsize),dcdb(nsize),dcdr(nsize),dwx(nsize),
     *            g1(nsize),g2(nsize)
        common/rar1/ dcda,dcdb,dcdr,dwx,g1,g2
c ---
        real*8     c,cvar,ugr,wvno,ratio,are
        common/rco/ c,cvar,ugr,wvno,ratio,are
c ---
        real*8      sumi0,sumi1,sumi2,sumi3,flagr
        common/rco1/ sumi0,sumi1,sumi2,sumi3,flagr
c ---
        integer*4 nderiv,ndpth
        real*8    dpth(nsize),dspz(nsize),dspr(nsize),drvz(5,nsize),
     *            drvr(5,nsize)
        common/deriv/ nderiv,ndpth,dpth,dspz,dspr,drvz,drvr
c ---
        logical KEY_ATTEN,KEY_FLAT,KEY_DERIV
        common/log/ KEY_ATTEN,KEY_FLAT,KEY_DERIV
c ---
        real*8   wwt(4),wt(4)
        real*8 yy1(nsize,5),yy2(nsize,5),yy3(nsize,5),yy4(nsize,5)
        real*8 yz1(nsize,5),yz2(nsize,5),yz3(nsize,5),yz4(nsize,5)
        real*8 dmrsmz(5),dmzsmr(5),depth(nsize)
        real*8 dmr(5),dmz(5),smr(5),smz(5)
        real*8 pr(5*nsize),xlamb(nsize),xmu(nsize)
c ---
        integer*4 i,ii,ivre,iter,j,jj,k,kk,l,ll,ldiv,m,mm,mm1,mmm1,max,n,nderv
        real*8    omega,a1,base,a4,a21,a13,a12,a2,a31,a24,a3,a34
        real*8    amplz,a43,a42,amplr,ap,atz,atr,auz,aur,brkt,bp
        real*8    c1,c2,c3,cosq,c4,cos2rm,covb,cova,cosra,ddz,del
        real*8    det,div,dmax,durdz,duzdz,dz,gam,fac1,expq,fac3
        real*8    fac2,fac4,gamm1,h,rab4,omegsq,rab3,q,rab1,rab2
        real*8    ra,rb,sin2ra,sinq,strsr,strsz,sum,t,tzz,wvar,uttq
        real*8    pi2,wvnum,wvnosq,xdiv,xtest,xxmin
        data wwt/0.0d0,0.5d0,0.5d0,1.0d0/ 
cMB    *       wt/0.1666667,0.3333333,0.3333333,0.1666667/
        wt(1)=1.0d0/6.0d0
        wt(2)=1.0d0/3.0d0
        wt(3)=wt(2)
        wt(4)=wt(1)
        pi2=datan(1.0d0)*8.0d0
        xnorm=0.0d0
        max=0
        tzz=0.0d0
        xdiv=0.0d0
        aur=0.0d0
        auz=0.0d0
        atz=0.0d0
        atr=0.0d0
c----------------------------------------------------------------------
                  do 99999 i=1,nsize
                  do 99999 j=1,5
        yy1(i,j)=1.d0
        yy2(i,j)=1.d0
        yy3(i,j)=1.d0
        yy4(i,j)=1.d0
        yz1(i,j)=2.d0
        yz2(i,j)=1.d0
        yz3(i,j)=1.d0
        yz4(i,j)=1.d0
99999             continue
        xxmin=1.0d-15
             if(num.gt.0.and..not.KEY_ATTEN)      goto 503
        mmax=jmax
        nmax=mmax
             if(fact.le.0.0d0) fact=4.
             if(num.eq.-1)      goto 88888
cMB          if(mmax) 777,777,778
             if(mmax.le.0) goto 777
c
cxx  778             continue
             if(istru.lt.1)     goto 20005
        call strut
             goto 20004
20005             continue
        base=0.0d0
        mm1=mmax-1
        ivre=(nsize-1)/mm1
             if(ndiv.gt.ivre) ndiv=ivre
        div=dble(ndiv)
        jj=1
             if(ndiv.le.1)      goto 20003
             if(b(1).le.0.1d-10) jj=2
                  do 10001 j=jj,mm1
        ldiv=(j-jj)*ndiv
                  do 10001 i=1,ndiv
        pr(ldiv+i)=d(j)/div
        pr(ldiv+nsize+i)=a(j)
        pr(ldiv+2*nsize+i)=b(j)
        pr(ldiv+3*nsize+i)=rho(j)
        if(KEY_ATTEN) pr(ldiv+4*nsize+i)=qs(j)
10001             continue
        mmax=(mm1-jj+1)*ndiv+jj
        d(mmax)=0.0d0
        a(mmax)=a(nmax)
        b(mmax)=b(nmax)
        rho(mmax)=rho(nmax)
        if(KEY_ATTEN)qs(mmax)=qs(nmax)
        nmax=mmax
        mm1=mmax-1
                  do 10002 j=jj,mm1
        d(j)=pr(j-jj+1)
        a(j)=pr(j-jj+nsize+1)
        b(j)=pr(j-jj+2*nsize+1)
        rho(j)=pr(j-jj+3*nsize+1)
        if(KEY_ATTEN) qs(j)=pr(j-jj+4*nsize+1)
10002             continue
20003             continue
c
c
20004             continue
        mmax=nmax
88888   base=0.0d0
                  do 22 i=1,mmax
        base=base+d(i)
        depth(i)=base-d(i)*0.5d0
        xmu(i)=rho(i)*b(i)*b(i)
        xlamb(i)=rho(i)*(a(i)*a(i)-2.0d0*b(i)*b(i))
cMB          if(i-mmax) 22,24,24
             if(i-mmax.ge.0) goto 24
22                continue
24      depth(mmax)=base-d(mmax)
c
c
503               continue
                  do 500 j=1,nmax
        dcda(j)=0.0d0
        dcdb(j)=0.0d0
        dcdr(j)=0.0d0
        ampur(j)=0.0d0
        ampuz(j)=0.0d0
        stresr(j)=0.0d0
        stresz(j)=0.0d0
        g1(j)=0.0d0
        g2(j)=0.0d0
500               continue
        t=tt(ixa)
        c=cc(ixa,jxa)
        ratio=rat(ixa,jxa)
        mmax=nmax
        dmax=fact*t*c
        sum=0.0d0
                  do 900 ii=1,mmax
        max=ii
cMB          if(c-b(ii)) 901,900,900
             if(c-b(ii).ge.0.0d0) goto 900
        sum=sum+d(ii)
             if(ii.eq.mmax)      goto 902
             if(sum.le.dmax)      goto 900
cMB          if(a(ii+1)-a(ii)) 902,903,90009
             if(a(ii+1)-a(ii).gt.0.0d0) goto 90009
             if(a(ii+1)-a(ii).lt.0.0d0) goto 902
cMB 903          if(b(ii+1)-b(ii)) 902,900,90009
             if(b(ii+1)-b(ii).gt.0.0d0) goto 90009
             if(b(ii+1)-b(ii).lt.0.0d0) goto 902
900               continue
90009   max=max+1
902     mmax=max
        sumi0=0.0d0
        sumi1=0.0d0
        sumi2=0.0d0
        sumi3=0.0d0
cMB     wvno=6.2831853072/(c*t)
        wvno=pi2/(c*t)
        wvnosq=wvno*wvno
        wvnum=wvno
cMB     omega=6.2831853072/t
        omega=pi2/t
        omegsq=omega*omega
             if(b(1).gt.0.0d0)      goto 4000
        ra=c/a(1)
        cr1=ra*ra-1.0d0
cMB     cra=wvno*cdsqrt(cr1)
        cra=cdsqrt(cr1)*wvno
                        if(cdabs(cra).le.1.0d-35)        go to 4001
                                                        go to 4002
4001    sumi0=rho(1)*d(1)
        sumi1=0.0d0
        sumi2=0.0d0
        sumi3=0.0d0
        tzz=0.0d0
             goto 4000
4002    cr1=cdsin(2.0d0*cra*d(1))
        cr2=4.0d0*cra
        cr3=cr1/cr2
        sin2ra=dreal(cr3)
        cr1=cdcos(cra*d(1))
        cosra=dreal(cr1)
        cos2rm=1.0d0/(cosra*cosra)
        fac1=(0.5d0*d(1)+sin2ra)*cos2rm
        fac3=wvno*(0.5d0*d(1)-sin2ra)*cos2rm
        cr1=cra*cra
        rab1=dreal(cr1)
        fac2=wvno*fac3/(rab1)
        fac4=rab1*fac3/wvno
        sumi0=rho(1)*(fac1+fac2)
        sumi1=xlamb(1)*fac2
        sumi2=xlamb(1)*fac3
        sumi3=xlamb(1)*fac4
        cr1=cdsin(cra*d(1))
        cr2=cr1/cra
        rab1=dreal(cr2)
        tzz=-rho(1)*omegsq*rab1/cosra
4000              continue
        cova=c/a(mmax)
        covb=c/b(mmax)
        gam=2.0d0/covb**2
        gamm1=gam-1.0d0
        ra=wvno*sqrt(dabs(cova**2-1.0d0))
        rb=wvno*sqrt(dabs(covb**2-1.0d0))
        det=wvnosq-ra*rb
C       write(*,*) num,' DET ',det,ra,rb
        h=rho(mmax)*omegsq
        brkt=-gamm1*wvno+gam*ra*rb/wvno
        iter=0
        aur1=1.0d0
        auz1=0.0d0
        atz1=-h*brkt/det
        atr1=-h*ra/det
        mmm1=mmax-1
                  do 1346 mm=1,mmm1
        m=mmax-mm
             if(b(m).le.0.0d0)      goto 1346
        xdiv=1.0d0
        ddz=-d(m)/(4.0d0*xdiv)
        a12=1.0d0/(xlamb(m)+2.0d0*xmu(m))
        a13=wvno*xlamb(m)*a12
        a21=-omegsq*rho(m)
        a24=wvno
        a31=-wvno
        a34=1.0d0/xmu(m)
        a42=-a13
        a43=a21+4.0d0*wvnosq*xmu(m)*(xlamb(m)+xmu(m))*a12
        yy3(m,5)=aur1
        yy1(m,5)=auz1
        yy2(m,5)=atz1
        yy4(m,5)=atr1
                  do 1338 kk=2,5
        k=6-kk
        eaur1=aur1
        eauz1=auz1
        eatz1=atz1
        eatr1=atr1
        daur1=0.0d0
        dauz1=0.0d0
        datz1=0.0d0
        datr1=0.0d0
                  do 1401 ll=1,4
        saur1=aur1+wwt(ll)*ddz*daur1
        sauz1=auz1+wwt(ll)*ddz*dauz1
        satz1=atz1+wwt(ll)*ddz*datz1
        satr1=atr1+wwt(ll)*ddz*datr1
        daur1=a31*sauz1+a34*satr1
        dauz1=a12*satz1+a13*saur1
        datz1=a21*sauz1+a24*satr1
        datr1=a42*satz1+a43*saur1
        eaur1=eaur1+wt(ll)*ddz*daur1
        eauz1=eauz1+wt(ll)*ddz*dauz1
        eatz1=eatz1+wt(ll)*ddz*datz1
        eatr1=eatr1+wt(ll)*ddz*datr1
1401              continue
        aur1=eaur1
        auz1=eauz1
        atz1=eatz1
        atr1=eatr1
cxx 1400              continue
        yy1(m,k)=auz1
        yy2(m,k)=atz1
        yy3(m,k)=aur1
        yy4(m,k)=atr1
1338              continue
1346              continue
             if(b(1).gt.0.0d0)      goto 1347
        yy1(1,1)=yy1(2,1)
        yy2(1,1)=yy2(2,1)
        yy3(1,1)=yy3(2,1)
        yy4(1,1)=yy4(2,1)
1347              continue
4003    aur2=0.0d0
        auz2=1.0d0
        atz2=-h*rb/det
        atr2=-h*brkt/det
             if(iter.eq.0)      goto 4004
        aur1=1.0d0
        auz1=0.0d0
        atz1=-h*brkt/det
        atr1=-h*ra/det
        aur2=aur2+xnorm*aur1
        auz2=auz2+xnorm*auz1
        atz2=atz2+xnorm*atz1
        atr2=atr2+xnorm*atr1
4004              do 2346 mm=1,mmm1
        m=mmax-mm
             if(b(m).le.0.0d0)      goto 2346
        ddz=-d(m)/(4.0d0*xdiv)
        a12=1.0d0/(xlamb(m)+2.0d0*xmu(m))
        a13=wvno*xlamb(m)*a12
        a21=-omegsq*rho(m)
        a24=wvno
        a31=-wvno
        a34=1.0d0/xmu(m)
        a42=-a13
        a43=a21+4.0d0*wvnosq*xmu(m)*(xlamb(m)+xmu(m))*a12
        yz3(m,5)=aur2
        yz1(m,5)=auz2
        yz2(m,5)=atz2
        yz4(m,5)=atr2
                do 2338 kk=2,5
        k=6-kk
        eaur2=aur2
        eauz2=auz2
        eatz2=atz2
        eatr2=atr2
        daur2=0.0d0
        dauz2=0.0d0
        datz2=0.0d0
        datr2=0.0d0
                  do 2401 ll=1,4
        saur2=aur2+wwt(ll)*ddz*daur2
        sauz2=auz2+wwt(ll)*ddz*dauz2
        satz2=atz2+wwt(ll)*ddz*datz2
        satr2=atr2+wwt(ll)*ddz*datr2
        daur2=a31*sauz2+a34*satr2
        dauz2=a12*satz2+a13*saur2
        datz2=a21*sauz2+a24*satr2
        datr2=a42*satz2+a43*saur2
        eaur2=eaur2+wt(ll)*ddz*daur2
        eauz2=eauz2+wt(ll)*ddz*dauz2
        eatz2=eatz2+wt(ll)*ddz*datz2
        eatr2=eatr2+wt(ll)*ddz*datr2
2401              continue
        aur2=eaur2
        auz2=eauz2
        atz2=eatz2
        atr2=eatr2
cxx 2400              continue
        yz1(m,k)=auz2
        yz2(m,k)=atz2
        yz3(m,k)=aur2
        yz4(m,k)=atr2
2338              continue
2346              continue
             if(b(1).gt.0.0d0)      goto 2347
        yz1(1,1)=yz1(2,1)
        yz2(1,1)=yz2(2,1)
        yz3(1,1)=yz3(2,1)
        yz4(1,1)=yz4(2,1)
2347              continue
        aa=yz3(1,1)-ratio*yz1(1,1)
        bb=ratio*yy1(1,1)-yy3(1,1)
             if(dabs(bb).lt.1.0d-10) bb=dsign(1.0d-10,bb)
        xnorm=aa/bb
        bb=xnorm*yy1(1,1)+yz1(1,1)
             if(dabs(bb).lt.1.0d-10) bb=dsign(1.0d-10,bb)
        ampur(nsize)=(xnorm*yy3(1,1)+yz3(1,1))/bb
        ampuz(nsize)=(xnorm*yy1(1,1)+yz1(1,1))/bb
        stresz(nsize)=(xnorm*yy2(1,1)+yz2(1,1))/bb
        stresr(nsize)=(xnorm*yy4(1,1)+yz4(1,1))/bb
        iter=iter+1
             if(iter.gt.1)      goto 2011
        xtest=dabs(ampur(nsize)/ratio-1.0d0)
             if(xtest.ge.0.00001d0)      goto 4003
2011    dcdb(nsize-1)=iter
        dcda(nsize)=0.0d0
        dcdb(nsize)=0.0d0
        dcdr(nsize)=0.0d0
        ampur(nsize-1)=yy3(1,1)
        ampur(nsize-2)=yz3(1,1)
        ampuz(nsize-1)=yy1(1,1)
        ampuz(nsize-2)=yz1(1,1)
        stresz(nsize-1)=yy2(1,1)
        stresz(nsize-2)=yz2(1,1)
        stresr(nsize-1)=yy4(1,1)
        stresr(nsize-2)=yz4(1,1)
        dcdb(nsize-2)=bb
        dcda(nsize-2)=0.0d0
        dcda(nsize-1)=0.0d0
        dcdr(nsize-2)=0.0d0
        dcdr(nsize-1)=0.0d0
                  do 7000 m=1,mmax
             if(b(m).le.0.0d0)      goto 7000
cMB          if(m-mmax) 7001,77777,77777
             if(m.ge.mmax) goto 77777
        dz=d(m)/4.
                  do 1339 kk=1,5
        aur=(xnorm*yy3(m,kk)+yz3(m,kk))/bb
        auz=(xnorm*yy1(m,kk)+yz1(m,kk))/bb
        atz=(xnorm*yy2(m,kk)+yz2(m,kk))/bb
        atr=(xnorm*yy4(m,kk)+yz4(m,kk))/bb
        durdz=atr/xmu(m)-wvno*auz
        duzdz=(atz+wvno*xlamb(m)*aur)/(xlamb(m)+2.0d0*xmu(m))
        dmr(kk)=aur*aur
        dmz(kk)=auz*auz
        smr(kk)=durdz*durdz
        smz(kk)=duzdz*duzdz
        dmrsmz(kk)=aur*duzdz
        dmzsmr(kk)=auz*durdz
cMB          if(kk-3) 1339,1301,1339
             if(kk.ne.3) goto 1339
        ampur(m)=aur
        ampuz(m)=auz
        stresz(m)=atz
        stresr(m)=atr
        g1(m)=(atz+aur*xlamb(m)*wvno)/(xlamb(m)+2.0d0*xmu(m))
        g2(m)=atr/xmu(m)
1339              continue
        dmmr=(dz/22.5d0)*(7.0d0*(dmr(1)+dmr(5))+32.0d0*(dmr(2)+
     *        dmr(4))+12.0d0*dmr(3))
        dmmz=(dz/22.5d0)*(7.0d0*(dmz(1)+dmz(5))+32.0d0*(dmz(2)+
     *        dmz(4))+12.0d0*dmz(3))
        smmz=(dz/22.5d0)*(7.0d0*(smz(1)+smz(5))+32.0d0*(smz(2)+
     *        smz(4))+12.0d0*smz(3))
        smmr=(dz/22.5d0)*(7.0d0*(smr(1)+smr(5))+32.0d0*(smr(2)+
     *        smr(4))+12.0d0*smr(3))
        drsz=(dz/22.5d0)*(7.0d0*(dmrsmz(1)+dmrsmz(5))+32.0d0*(dmrsmz(2)
     *  +dmrsmz(4))+12.0d0*dmrsmz(3))
        dzsr=(dz/22.5d0)*(7.0d0*(dmzsmr(1)+dmzsmr(5))+32.0d0*(dmzsmr(2)
     *  +dmzsmr(4))+12.0d0*dmzsmr(3))
        sumi0=sumi0+rho(m)*(dmmr+dmmz)
        sumi1=(xlamb(m)+2.0d0*xmu(m))*dmmr+xmu(m)*dmmz+sumi1
        sumi2=xmu(m)*dzsr-xlamb(m)*drsz+sumi2
        sumi3=(xlamb(m)+2.0d0*xmu(m))*smmz+xmu(m)*smmr+sumi3
        dldl=-wvnosq*dmmr+2.0d0*wvno*drsz-smmz
        dldm=-wvnosq*(2.0d0*dmmr+dmmz)-2.0d0*wvno*dzsr-(2.0d0*smmz+smmr)
        dldr=omegsq*(dmmr+dmmz)
        dcdb(m)=2.0d0*rho(m)*b(m)*c*(dldm-2.0d0*dldl)/wvno
        dcda(m)=2.0d0*rho(m)*a(m)*c*dldl/wvno
        dcdr(m)=(c/wvno)*(dldr+xlamb(m)*dldl/rho(m)+xmu(m)*dldm/rho(m))
cMB          if(dabs(auz)+dabs(aur)-xxmin) 7002,7002,7000
             if((dabs(auz)+dabs(aur)-xxmin).le.0.0d0) goto 7002
7000              continue
77777             continue
             if((b(1).gt.0.1d-10).or.m.ne.2)      goto 7002
        aur=ratio
        auz=1.0d0
        atr=0.0d0
        atz=tzz
7002    ampur(m)=aur
        ampuz(m)=auz
        stresr(m)=atr
        stresz(m)=atz
        g1(m)=(atz+aur*xlamb(m)*wvno)/(xlamb(m)+2.0d0*xmu(m))
        g2(m)=atr/xmu(m)
        ap=-rho(m)*(wvno*aur+rb*auz)/det
        bp=-rho(m)*(-ra*aur/wvno-auz)/det
        a1=-wvno*ap/rho(m)
        a2=-wvno*rb*bp/rho(m)
        a3=ra*ap/rho(m)
        a4=wvnosq*bp/rho(m)
cMB          if(rb) 7005,7006,7005
             if(rb.eq.0.0d0) goto 7006
        dmmr=a1*a1/(2.0d0*ra)+2.0d0*a1*a2/(ra+rb)+a2*a2/(2.0d0*rb)
        dmmz=a3*a3/(2.0d0*ra)+2.0d0*a3*a4/(ra+rb)+a4*a4/(2.0d0*rb)
        smmz=ra*a3*a3/2.0d0+2.0d0*ra*rb*a3*a4/(ra+rb)+rb*a4*a4/2.0d0
        smmr=ra*a1*a1/2.0d0+2.0d0*ra*rb*a1*a2/(ra+rb)+rb*a2*a2/2.0d0
        drsz=-a1*a3/2.0d0-(a1*a4*rb+a2*a3*ra)/(ra+rb)-a2*a4/2.0d0
        dzsr=-a1*a3/2.0d0-(a1*a4*ra+a2*a3*rb)/(ra+rb)-a2*a4/2.0d0
             goto 7010
7006    ugr=b(m)
        flagr=0.00d0
        sumi0=rho(m)*10.0d0**(25)
        sumi1=xmu(m)*10.0d0**(25)
        sumi2=0.0d0
        sumi3=0.0d0
        are=0.0d0
        dcdb(m)=-2.0d0*wvno*10.0d0**(25)
             goto 531
7010              continue
        sumi0=sumi0+rho(m)*(dmmr+dmmz)
        sumi1=(xlamb(m)+2.0d0*xmu(m))*dmmr+xmu(m)*dmmz+sumi1
        sumi2=xmu(m)*dzsr-xlamb(m)*drsz+sumi2
        sumi3=(xlamb(m)+2.0d0*xmu(m))*smmz+xmu(m)*smmr+sumi3
        dldr=omegsq*(dmmr+dmmz)
        dldm=-wvnosq*(2.0d0*dmmr+dmmz)-2.0d0*wvno*dzsr-(2.0d0*smmz+smmr)
        dldl=-wvnosq*dmmr+2.0d0*wvno*drsz-smmz
        dcda(m)=2.0d0*rho(m)*a(m)*c*dldl/wvno
        dcdb(m)=2.0d0*rho(m)*b(m)*c*(dldm-2.0d0*dldl)/wvno
        dcdr(m)=(c/wvno)*(dldr+xlamb(m)*dldl/rho(m)+xmu(m)*dldm/rho(m))
cxx 7011              continue
        ugr=(wvno*sumi1+sumi2)/(omega*sumi0)
        flagr=omegsq*sumi0-wvnosq*sumi1-2.0d0*wvno*sumi2-sumi3
        wvar=(-sumi2+sqrt(dabs(sumi2**2-sumi1*(sumi3-
     *  omegsq*sumi0))))/sumi1
        cvar=omega/wvar
        are=1.0d0/(2.0d0*c*ugr*sumi0)/dsqrt(pi2)*1.0d-15
531               continue
cxx 501               continue
c
c
c
c
c
c
        n=1
             if(b(1).le.0.0d0) n=2
                  do 5010 m=n,mmax
        dldk=-2.0d0*(wvnum*sumi1+sumi2)
        dcdr(m)=dcdr(m)/dldk
        dcda(m)=dcda(m)/dldk
        dcdb(m)=dcdb(m)/dldk
        dwx(m)=(dcda(m)*4.0d0/3.0d0*b(m)/a(m)+dcdb(m))*b(m)
5010                continue
                        if(b(1).le.0.0d0)         go to 5556
        do 5555 m=1,mmax
        i=mmax-m+1
        j=i+1
        dept1(j)=depth(i)
        ampur(j)=ampur(i)
        ampuz(j)=ampuz(i)
        stresz(j)=stresz(i)
        stresr(j)=stresr(i)
        dcda(j)=dcda(i)
        dcdb(j)=dcdb(i)
        dcdr(j)=dcdr(i)
        dwx(j)=dwx(i)
        g1(j)=g1(i)
        g2(j)=g2(i)
5555                continue
        mmax=mmax+1
                                                go to 5558
5556                            do 5557 m=1,mmax
5557    dept1(m)=depth(m)
5558    dept1(1)=0.0d0
        if(b(1).le.0.0d0) dept1(1)=d(1)
        ampur(1)=ratio
        ampuz(1)=1.0d0
        stresz(1)=0.0d0
        if(b(1).le.0.0d0) stresz(1)=tzz
        stresr(1)=0.0d0
        dcda(1)=0.0d0
        dcdb(1)=0.0d0
        dcdr(1)=0.0d0
        dwx(1)=0.0d0
        g1(1)=0.0d0
        g2(1)=0.0d0
c----------------------------------------------------------------------
c       Calculation of eigen function and its derivatives
c----------------------------------------------------------------------
                        if(iedit.ne.2)                  go to 800
                                do 799 i=1,5
                                do 799 j=1,nsize
        dspz(j)=0.0d0
        dspr(j)=0.0d0
        drvz(i,j)=0.0d0
        drvr(i,j)=0.0d0
799                                     continue
                                do 801 l=1,ndpth
                        if(dpth(l).lt.dept1(1))         go to 801
c----------------------------------------------------------------------
c       If dpth(l) is above the free surface, the output will be zero
c----------------------------------------------------------------------
                                do 802 ii=2,mmax
        rab1=dept1(ii)-dpth(l)
                        if(dabs(rab1).lt.1.0d-5) rab1=0.0d0
                        if(rab1.le.0.0d0)                  go to 802
        i=ii-1
        j=ii
                        if(b(1).le.0.0d0) i=ii
        rab2=d(i)/2.0d0
        rab1=rab1-rab2
                        if(dabs(rab1).lt.1.0d-5) rab1=0.0d0
                        if(rab1.le.0.0d0)                  go to 803
        j=j-1
        i=i-1
                                                        go to 803
802                                     continue
        j=mmax
        i=mmax-1
                        if(b(1).le.0.0d0) i=j
803     dz=dpth(l)-dept1(j)
        amplz=ampuz(j)
        amplr=ampur(j)
        strsz=stresz(j)
        strsr=stresr(j)
        cova=c/a(i)
        covb=c/b(i)
        rb=wvno*sqrt(dabs(covb**2-1.0d0))
        ra=wvno*sqrt(dabs(cova**2-1.0d0))
        del=2.0d0*b(i)*b(i)-c*c
c----------------------------------------------------------------
c       Calculation of terms, connected with a-velocity
c----------------------------------------------------------------
        rab1=(c-a(i))/c
                        if(dabs(rab1).lt.1.0d-5) rab1=0.0d0
cMB                     if(rab1) 810,811,812
                        if(rab1.ge.0.0d0) goto 812
                        if(rab1.eq.0.0d0) goto 811
        nderv=nderiv+1
        q=ra*dz
        expq=dexp(q)
        rab2=2.0d0*xmu(i)*wvno*ra*amplr
        rab3=del*rho(i)*wvno*wvno*amplz
        rab4=2.0d0*rho(i)*c*c*ra*wvno
        c1=(rab2-rab3+wvno*strsr-ra*strsz)/rab4
        uttq=dabs(c1/amplr)
                        if(uttq.lt.1.0d-5) c1=0.0d0
        c2=(rab2+rab3-wvno*strsr-ra*strsz)/rab4
                                do 804 ii=1,nderv
        k=ii-1
        rab1=ra**k
        rab2=rab1*c1*expq
        rab3=(-1)**k*rab1*c2/expq
                        if(k.ne.0)                      go to 805
        dspr(l)=rab2+rab3
        dspz(l)=(ra*rab2-ra*rab3)/wvno
                                                        go to 804
805     drvr(k,l)=rab2+rab3
        drvz(k,l)=(ra*rab2-ra*rab3)/wvno
804                                     continue
                                                        go to 700
811     rab2=(2.0d0*xmu(i)*wvno*amplr-strsz)/(rho(i)*c*c*wvno)
        rab3=(del*rho(i)*wvno*amplz-strsr)/(rho(i)*c*c)
        dspr(l)=rab2-rab3*dz
        dspz(l)=-rab3/wvno
                        if(nderiv.lt.1)                 go to 700
        drvr(1,l)=-rab3
        drvz(1,l)=0.0d0
                                                        go to 700
812     nderv=nderiv+1
        c1=(2.0d0*xmu(i)*wvno*amplr-strsz)/(rho(i)*c*c*wvno)
        c2=(del*rho(i)*wvno*amplz-strsr)/(rho(i)*c*c*wvno)
                                do 806 ii=1,nderv
        k=ii-1
        rab1=ra**k
cMB     q=ra*dz+k*3.1415926/2.
        q=ra*dz+k*pi2/4.0d0
        cosq=dcos(q)
        sinq=dsin(q)
                        if(k.ne.0)                      go to 807
        dspr(l)=c1*cosq-wvno*c2*sinq/ra
        dspz(l)=-ra*c1*sinq/wvno-c2*cosq
                                                        go to 806
807     drvr(k,l)=(c1*cosq-wvno*c2*sinq/ra)*rab1
        drvz(k,l)=-(ra*c1*sinq/wvno+c2*cosq)*rab1
806                                     continue
c----------------------------------------------------------------------
c       Calculation of terms, connected with b-velocity
c----------------------------------------------------------------------
700     rab1=(c-b(i))/c
                        if(dabs(rab1).lt.1.0d-5) rab1=0.0d0
cMB                     if(rab1) 710,711,712
                        if(rab1.gt.0.0d0) goto 712
                        if(rab1.eq.0.0d0) goto 711
        nderv=nderiv+1
        q=rb*dz
        expq=dexp(q)
        rab2=2.0d0*xmu(i)*wvno*rb*amplz
        rab3=del*rho(i)*wvno*wvno*amplr
        rab4=2.0d0*rho(i)*c*c*wvno*wvno
        c3=(rab2-rab3+wvno*strsz-rb*strsr)/rab4
        uttq=dabs(c3/amplz)
                        if(uttq.lt.1.0d-5) c3=0.0d0
        c4=(-rab2-rab3+wvno*strsz+rb*strsr)/rab4
                                do 704 ii=1,nderv
        k=ii-1
        rab1=rb**k
        rab2=rab1*c3*expq
        rab3=(-1)**k*rab1*c4/expq
                        if(k.ne.0)                      go to 705
        dspr(l)=dspr(l)+rab2+rab3
        dspz(l)=dspz(l)+(rab2-rab3)*wvno/rb
                                                        go to 704
705     drvr(k,l)=drvr(k,l)+rab2+rab3
        drvz(k,l)=drvz(k,l)+(rab2-rab3)*wvno/rb
704                                     continue
                                                        go to 801
711     rab2=(2.0d0*xmu(i)*wvno*amplz-strsr)/(rho(i)*c*c*wvno)
        rab3=(del*rho(i)*wvno*amplr-strsz)/(rho(i)*c*c)
        dspz(l)=dspz(l)+rab2-rab3*dz
        dspr(l)=dspr(l)-rab3/wvno
                        if(nderiv.lt.1)                 go to 801
        drvz(1,l)=drvz(1,l)-rab3
                                                        go to 801
712     nderv=nderiv+1
        c3=(2.0d0*xmu(i)*wvno*amplz-strsr)/(rho(i)*c*c*wvno)
        c4=(del*rho(i)*wvno*amplr-strsz)/(rho(i)*c*c*wvno)
                                do 706 ii=1,nderv
        k=ii-1
        rab1=rb**k
cMB     q=rb*dz+k*3.1415926/2.
        q=rb*dz+k*pi2/4.0d0
        cosq=dcos(q)
        sinq=dsin(q)
                        if(k.ne.0)                      go to 707
        dspz(l)=dspz(l)+c3*cosq-wvno*c4*sinq/rb
        dspr(l)=dspr(l)-rb*c3*sinq/wvno-c4*cosq
                                                        go to 706
707     drvz(k,l)=drvz(k,l)+(c3*cosq-wvno*c4*sinq/rb)*rab1
        drvr(k,l)=drvr(k,l)-(rb*c3*sinq/wvno+c4*cosq)*rab1
706                                     continue
801                                     continue
800                                                             RETURN
777                    continue
                                                                STOP
                                                                END
c///////////////////////////////////////////////////////////////////
c-------------------------------------------------------------------
        subroutine STRUT
c-------------------------------------------------------------------
        implicit none
        integer*4 nsize
        parameter (nsize=1000)
c ---
        integer*4 nst,idum
        real*8    accur,ds(nsize)
        common/st/ nst,idum,accur,ds
c ---
        real*8    a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/d/ a,b,rho,d,qs
c ---
        integer*4 nmax,q1,q2,q3,q4,q5,q6
        real*8    q7
        common/c/ q7,nmax,q1,q2,q3,q4,q5,q6
c ---
        real*8    r(nsize),rr(nsize)
        integer*4 i,j,k,n,l,m,mmax
        real*8    h,f
        h=0.0
c-------------------------------------------------------------------
c
c    ---- NST GT 0 ---- ---- IF B(1)=0   NMAX GT 0 --------
c
c-------------------------------------------------------------------
        d(nmax)=0.
             if(b(1).ge.0.1d-10)      goto 100
                                      do 1  i=1,nst
             if(ds(i).ge.d(1))      goto 100
        ds(i)=d(1)
1                                     continue
100                         continue
        mmax=nmax-1
             if(mmax.eq.0)      goto 101
        r(1)=d(1)
             if(mmax.eq.1)      goto 101
                                      do 2  i=2,mmax
        j=i-1
        r(i)=r(j)+d(i)
2                                     continue
101                         continue
        rr(1)=0.
        l=1
        m=1
        n=1
                                      do 3  i=1,98
             if(n.gt.nst)      goto 103
             if(m.gt.mmax)      goto 102
        h=ds(n)-r(m)
             if(h.ge.accur)      goto 103
102                         continue
        n=n+1
             if(l.gt.1)      goto 104
             if(abs(ds(n-1)).le.accur)      goto 3
             goto 106
104                         continue
        f=ds(n-1)-abs(rr(l))
             if(abs(f).le.accur)      goto 3
106                         continue
             if(m.gt.mmax)      goto 107
             if(abs(h).le.accur)      goto 103
107                         continue
        l=l+1
        rr(l)=ds(n-1)
             goto 3
103                         continue
             if(m.gt.mmax)      goto 105
        l=l+1
        rr(l)=-r(m)
        m=m+1
3                                     continue
105                         continue
        k=nmax
        rr(l+1)=rr(l)
                                      do 4  i=1,l
        j=l-i+1
        a(j)=a(k)
        b(j)=b(k)
        rho(j)=rho(k)
             if(rr(j).ge.0.0)      goto 108
        rr(j)=-rr(j)
        k=k-1
108                         continue
        d(j)=rr(j+1)-rr(j)
        
4                                     continue
        nmax=l
                                                        RETURN
                                                        END



