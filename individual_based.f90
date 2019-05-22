  PROGRAM chaos ! AGENT-BASED LOGISTIC MODEL
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: X(:,:),X0(:),DR(:),DRP(:,:),XAV(:) 
  DOUBLE PRECISION, ALLOCATABLE :: DAV(:),H(:),BR(:)
  LOGICAL, ALLOCATABLE :: LV(:) 
  INTEGER, ALLOCATABLE :: IDM(:) 
  INTEGER N,K,K1,K2,K3,I,IMAX,I0,IL,IT,IW,NH,plotcount,DMC,DMTOT,KC
  DOUBLE PRECISION T,DT,TMAX,TMEAS,RAN3,EPS,MUT,DRT,S,W,TVF,RR,R,DX,XTMFE
  DOUBLE PRECISION DTM,TV,aa,comp,CC,DAVT,CCA,LHMAX, DIST,SA,SB,CB, XCA
  DOUBLE PRECISION ATM,XTM,BRT
  LOGICAL GAUSS
  N=5 ! SPACE DIMENSIONALITY
  IMAX=30000 ! MAX POPULATION
  MUT=1.D-2 ! MUTATION AMPLITUDE
  EPS=MUT ! SPREAD AROUND THE INITIAL POSITION 
  I0=300 ! I0=K0=TARGET NUMBER OF INDIVIDUALS 
  NH=50 ! NUMBER OF POINTS IN COORDINATE HISTOGRAM
  LHMAX=2.D0 ! MAXIMUM OF HISTOGRAM FOR COORDINATE
  TMAX=80000.D0 ! END TIME
  DTM=10.D0 ! MEASURE TIME
  SA=0.5D0   ! WIDTH OF THE GAUSSIAN PART OF COMP. KERN.
  SB=0.15d0!0.16D0   ! WIDTH OF THE GAUSSIAN OF THE BIRTH RATE
  CB=0.75d0!0.9d0 ! BIRTH COEFFICIENT
  XCA=1.0D0
  GAUSS= .TRUE. ! IF TRUE, THE G. PART WITH SA IS IN COMP. KERNEL 
  ALLOCATE(X(N,IMAX),DR(IMAX),BR(IMAX),DRP(IMAX,IMAX),LV(IMAX))
  ALLOCATE(X0(N),XAV(N),DAV(N),H(-NH:NH),IDM(N))
  open(unit=12,file='rand.dat') ! READING THE RANDOM SEED
  read(12,*)iw
  close(12)


  LV=.FALSE.
  IT=INT(I0*0.1) ! 0.6 IS IMPIRICAL
  DO I=1,IT ! SEEDING THE INITIAL CONDITIONS AROUND X0
     DO K=1,N
        X(K,I)=(2*ran3(iw)-1.d0)*EPS
     END DO
     LV(I)=.TRUE.
  END DO
  DR=0.D0 ! DEATH RATES
  DRP=0.D0
  DRT=0.D0
  BRT=0.D0
  DO I=1,IT ! DEATH RATE OF I
     DO K=1,IT ! DEATH RATE OF I CAUSED BY K
  DRP(I,K)=comp(x,i,k,i0,n,imax,GAUSS,XCA)
  DR(I)=DR(I)+DRP(I,K)
     END DO
     DRT=DRT+DR(I)
     atm=1.d0
     do k=2,n  ! BIRTH RATE
  xtm=x(k,i)   
  atm=atm*(dexp(-xtm*xtm/(2.d0*sb*sb))*(1.d0-cb) + cb)   
  end do
  BR(I)=atm
  BRT=BRT+ATM
  END DO

  IL=IT ! TOTAL NUMBER OF INDIVIDUALS=NUMBER OF LIVE INDIVIDUALS
  H=0.D0
  T=0.D0
  TMEAS=0.d0
  plotcount=0
!  OPEN(UNIT=15,FILE='movie_cyl.dat')
!  OPEN(UNIT=10,FILE='movie.dat')
  OPEN(UNIT=20,FILE='dimensions_ind_015_075_300_1.dat')
  OPEN(UNIT=30,FILE='xy_av1.dat')
  OPEN(UNIT=40,FILE='dispersion1.dat')
  OPEN(UNIT=90, FILE='hist1.dat') ! 

  DO WHILE (T.LE.TMAX.AND.IT.LT.IMAX) !!!! BEGINNING OF UPDATE !!!!!
  TV=DRT+BRT ! TOTAL BIRTH+DEATH RATE
  R=ran3(iw)
  IF (R.LE.BRT/TV) THEN ! BIRTH EVENT
        RR=ran3(iw)
        S=0.D0
        K=0
        DO WHILE (S.LE.RR) ! CHOOSING THE FATHER NUMBER K
          K=K+1
          IF (LV(K)) S=S+BR(K)/BRT
       END DO
!!$     
!!$4    K=int(ran3(iw)*IT)+1 ! CHOOSING THE FATHER
!!$     IF (.NOT.LV(K)) GO TO 4
     I=1
     DO WHILE (LV(I)) ! FINDING A SPARE SPOT
     I=I+1
     END DO
     IF(I.GT.IT) IT=I ! ADDING A NEWBORN TO THE RIGHT END OF THE LIST
     IL=IL+1 ! INCREMENTING THE NUMBER OF LIVE INDIVIDUALS
     LV(I)=.TRUE. ! ASSINGING A SEAT
     DO K1=1,N
        X(K1,I)=X(K1,K) + MUT*(2*ran3(iw)-1.D0)
     END DO

     DO K=1,IT ! INCREMENT IN DEATH RATE OF K CAUSED BY NEW I
        IF(LV(K)) THEN
   DRP(K,I)=comp(x,k,i,i0,n,imax,GAUSS,XCA)
   DR(K)=DR(K)+DRP(K,I)
   DRT=DRT+DRP(K,I)
   END IF
   END DO

   DR(I)=0.D0
     DO K=1,IT ! DEATH RATE OF NEW I CAUSED BY K
        IF(LV(K)) THEN
  DRP(I,K)=comp(x,i,k,i0,n,imax,GAUSS,XCA)
  DR(I)=DR(I)+DRP(I,K)
 END IF
 END DO
 DRT=DRT+DR(I)
  atm=1.d0
  do k=2,n  ! BIRTH RATE
  xtm=x(k,i)   
  atm=atm*(dexp(-xtm*xtm/(2.d0*sb*sb))*(1.d0-cb) + cb)   
  end do
  BR(I)=atm
  BRT=BRT+BR(I)
  
    ELSE   IF (R.GT.IL/TV) THEN ! ! DEATH EVENT
        RR=ran3(iw)
        S=0.D0
        K=0
        DO WHILE (S.LE.RR) ! CHOOSING VICTIM'S NUMBER K
          K=K+1
          IF (LV(K)) S=S+DR(K)/DRT
       END DO
        IL=IL-1
        LV(K)=.FALSE.
        DRT=DRT-DR(K)
        BRT=BRT-BR(K)
        DR(K)=0.D0
        BR(K)=0.D0
        DO I=1,IT ! REDUCING OTHER DEATH RATES
           IF(LV(I)) THEN
           DR(I)=DR(I)-DRP(I,K)
           DRT=DRT-DRP(I,K)
           DRP(I,K)=0.D0
           DRP(K,I)=0.D0
           END IF
        END DO
     END IF

6     DT=10/TV*ran3(iw) ! ADVANCING TIME
      IF (DEXP(-DT*TV).LE.ran3(iw)) GO TO 6
      T=T+DT

      IF(T.GE.TMEAS) THEN ! MEASUREMENT
         DRT=0.D0 ! RECALCULATING THE DEATH RATE
         BRT=0.D0
         DO K=1,IT
            IF(LV(K)) THEN
               DRT=DRT+DR(K)
               BRT=BRT+BR(K)
            END IF   
         END DO
!    WRITE(10,*) (plotcount*100,K1=1,2)
!   WRITE(15,*) (plotcount*100,K1=1,2)
     plotcount=plotcount+1

     TMEAS=TMEAS+DTM

     DMTOT=0 ! MEASURING DIMENSION
     IDM=0
        DO KC=1,IT
            IF(LV(KC)) THEN
               DMC=N
               DO K1=2,N
                  IF (DABS(X(K1,KC)).LE.SB*2.D0) then 
                     DMC=DMC-1
                  END IF
               END DO
                  IDM(DMC)=IDM(DMC)+1              
                  DMTOT=DMTOT+DMC
            END IF
         END DO
         WRITE(*,*)'AVERAGE DIMENSION',1.D0*DMTOT/IL
         DO K1=1,N
            WRITE(*,*) K1,IDM(K1)
         END DO   
         XAV=0.D0 ! AVERAGE COORDINATE
         DAV=0.D0 ! DISPERSION OF COORDINATE
         DO K=1,IT
            IF(LV(K)) THEN
               S=0
               DO K1=2,N
               S=S+X(K1,K)**2   
               END DO
!               WRITE(15,*) X(1,K),DSQRT(S)
!               WRITE(10,*) (X(K1,K),K1=1,2)
               DO K1=1,N
                  XAV(K1)=XAV(K1)+X(K1,K)/IL
                  DAV(K1)=DAV(K1)+X(K1,K)*X(K1,K)/IL
                   K2=NINT(NH*X(K1,K)/LHMAX) ! HISTOGRAM OF DENSITY
                  IF(iabs(K2).LE.NH) H(K2)=H(K2)+1
               END DO
             END IF
         END DO
         DAVT=0.D0
         DIST=0.D0
          DO K1=1,N
            DIST=DIST+(X0(K1)-XAV(K1))**2
            DAV(K1)=DAV(K1)-XAV(K1)*XAV(K1)
            DAVT=DAVT+DAV(K1) ! TOTAL DISPERSION
            CC=CC-(XAV(K1)**4)/4
         END DO
         CCA=DEXP(CC) ! AV. CARRYING CAPACITY         
!         WRITE(15,*) '&'
         WRITE(20,*) t,(IDM(K1),K1=1,N),1.D0*DMTOT/IL
         WRITE(30,*)XAV(1),XAV(2), DSQRT(DAV(1)),DSQRT(DAV(2)) 
         WRITE(40,*) T,  DSQRT(DAVT), IL
         WRITE(*,*)  T,  DSQRT(DAVT), IL
      END IF

   END DO

!   CLOSE(15)
!   CLOSE(10)
   CLOSE(20)
   CLOSE(30)
   CLOSE(40)

  DO K1=-NH,NH
     WRITE(90,*)K1*LHMAX/NH,H(K1)
  END DO
  CLOSE(90)



  open(unit=12, file='rand.dat')
  write(12,*)nint(-ran3(iw)*1.d+6)
  close(12)         
  open(unit=35, file='final_points.dat') !writing last snapshot
  write(35,*)t
  DO K=1,IT
  IF(LV(K)) WRITE(35,*) (X(K1,K),K1=1,N)
  END DO
  close (35)
  write(*,*) 'survivors, it, il', it, il
END PROGRAM CHAOS

  function comp(x,i,k,i0,n,imax,GAUSS,XCA) ! compeeting effect of
  integer i,k,i0,n,imax ! k on i with ith carrying capacity
  integer k1 ! that is, death rate of i from k
  double precision x(n,imax), s, w, comp, sk, sa , XCA
  logical GAUSS
  sa=.5d0 ! width of the competition kernel
   comp=0.d0
   if(i.ne.k) then
        W=0.D0
        S=0.D0
  do k1=1,n
     dx=x(k1,k)-x(k1,i)
     if(GAUSS) w=w-dx*dx/(2*sa*sa) ! ADDING A GAUSSIAN COMPETITION KERNEL
     S=S-((X(K1,I)-XCA)**4)/4
  end do
  comp=DEXP(W-S)/I0
    end if
  return
  end function comp

!==========================================================================
!  ran3 from Numerical Recipes, 2nd edition
!--------------------------------------------------------------------------
      function ran3(idum)
!==========================================================================

!         implicit real*4(m)

      integer idum
      integer mbig,seed,mz
!      real mbig,seed,mz
      double precision ran3,fac

!         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
!      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
 
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
!      real  mj,mk,ma(55)

      save iff,inext,inextp,ma

      data iff /0/

1      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if(ran3.le.0.or.ran3.ge.1) goto 1
!!      if(ran3.le.0.or.ran3.ge.1) Then
      !! write(6,*)'RAN3 failed: ', ran3
      ! stop
!!      endif
      return
      end function ran3



