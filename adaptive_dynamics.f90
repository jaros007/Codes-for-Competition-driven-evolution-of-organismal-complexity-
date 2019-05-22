  PROGRAM chaos !Multicluster AD for chaotic attractors
  IMPLICIT NONE ! general model with arbitrary coefficients
  DOUBLE PRECISION, ALLOCATABLE :: X(:),DX(:),XN(:),X0(:),DV(:),VLT(:)
  DOUBLE PRECISION, ALLOCATABLE :: CC(:,:),ME(:,:),CR(:),CN(:),DC(:)
  DOUBLE PRECISION, ALLOCATABLE :: B(:,:),DME(:),X2M(:),XC(:),BRT(:)   
  DOUBLE PRECISION, ALLOCATABLE :: B1(:,:),DSEED(:),SI(:),AD(:,:),VD(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: PC(:),PA(:),DP(:),PN(:),X1M(:)
  INTEGER N,K,L,IW,K1,K2,K3,MC,K4,IC,KC,NH,plotcount,nctot,iw0,NM,km, DIN
  INTEGER KC1,KC2,JE,NE,J1, J2, JL,MCIN,NR,NRTOT,per,glass,conv
  INTEGER NDK,MX,DMC,DMTOT
  DOUBLE PRECISION T, DT, TFIN, TMEAS, DTM, DIF,ran3,TR,VAV,WW,pt,pm
  DOUBLE PRECISION SA,SB,COMP,DTP,TP,S7,V,VMIN, POPMIN,CB,POPBIRTH
  DOUBLE PRECISION DISP,SIG,C1,C2,DIST,dw,lw,d1,VLTM,WUW,TAV,TAVV
  DOUBLE PRECISION GS,XP,GW,mu,TSEED,DTSEED,DCC,DCMIN,DCFIN,VGLASS
  DOUBLE PRECISION VFIN,DTMERGE,TMERGE,AA,LSEED, CBIN,NCFAR,KCFAR,PDIM
  DOUBLE PRECISION PW,DMAV,DMA,DMFIN,coeffbirth,XCA,PDECR,BRATE,LMERGE
  LOGICAL GAUSS 
  LOGICAL, ALLOCATABLE :: LIVE(:),NEAR(:)
  INTEGER, ALLOCATABLE :: CLDM(:),DML(:),DDD(:) 
  GS(XP,GW)=DEXP(-XP*XP/(GW*GW*2.d0)) ! GAUSS
  N=5    !SPACE DIMENSIONALITY
!  DIN=1 ! DIMENSIONALITY OF INITIAL CLUSTERS
!  PDIM=0.5d0! PROBABILITY TO INCREASE DIMENSION ON EVERY SPLITTING
!  PDECR=0.0d0! PROBABILITY TO DECREASE DIMENSION ON EVERY SPLITTING
  TFIN=6000.d0
  DTSEED=1.D0 ! TIME TO SEED A NEW CLUSTER
  DTMERGE=DTSEED ! TIME TO MERGE CLUSTERS
  DTM=1.d0  ! MEASURE TIME
  TAV=1.d0  ! AVERAGING TIME FOR THE DISPLACEMENT
  DT=1.d-2  ! POSITION INTEGRATION TIME STEP
  DTP=1.d-2 ! POPULATION INTEGRATION TIME STEP
  MC=300   ! MAXIMUM NUMBER OF CLUSTERS
  MCIN=1   ! INITIAL NUMBER OF CLUSTERS
  NE=10.D0/DTP ! NUMBER OF TIME STEPS IN ECOLOGY INTEGRATIOn  
  GW=1.D0   ! WIDTH OF A GAUSSIAN DISTRIBUTION FOR COEFFICIENTS
  SA=0.5D0   ! WIDTH OF THE GAUSSIAN PART OF COMP. KERN.
  SB=0.17D0   ! WIDTH OF THE GAUSSIAN OF THE BIRTH RATE
  coeffbirth=1.d0!0.9d0 ! POWER
  XCA=0.9D0
  GAUSS=.true. ! IF TRUE, THE G. PART WITH SA IS IN COMP. KERNEL
  CB= 1.D-2 ! DISTANCE BETWEEN THE OFFSPRING AND THE ANSESTOR 
  DCMIN=CB ! MINIMAL DISTANCE BELOW WHICH THE CLUSTERS ARE MERGED
  DCFIN=0.1D0 ! FINAL VIZIBLE MINIMAL DISTANCE BETWEEN CLUSTERS 
  VMIN=1.D-4 ! MINIMUM VELOCITY UPON WHICH ALL STOPS
  VGLASS=1.D-4 ! A THRESHOLD BETWEEN STEADY AND GLASSY BEHAVIOUR.
  POPMIN=1.D-8! THRESHOLD FOR THE POPULATION
  POPBIRTH=2*POPMIN ! RARE MUTANT POPULATION
  mu=0.5D0 !0.0001d0!fraction of mutant pop size
  NRTOT=1 ! TOTAL NUMBER OF REPETITIONS
  MX=NINT(TFIN/DTM)+1 ! NUMBER OF MEASUREMENT POINTS TO ALLOCATE MEM
  ALLOCATE(X(N*MC),DX(N*MC),XN(N*MC),X0(N),DV(N*MC),VLT(N*MC))
  ALLOCATE(B(N,N),CC(N,MC),LIVE(MC),NEAR(MC),DSEED(N))
  ALLOCATE(B1(N,N),SI(N),CLDM(N),AD(N,MX),VD(N,MX),DML(N))
  ALLOCATE(ME(MC,MC),PC(MC),PA(MC),DP(MC),PN(MC),BRT(MC))
  ALLOCATE(CR(N*MC),CN(N*MC),DC(N*MC),X1M(N),X2M(N),XC(N),DDD(N))
  open(unit=12,file='rand.dat')
  read(12,*)iw
!  iw=-536456!3546345
  close(12)
  per=0
  conv=0
  glass=0
  nctot=0
  ncfar=0
  vfin=0.D0
  AD=0.D0
  VD=0.D0

  
  DO NR=1,NRTOT

!  open(unit=10,file='coeff2.dat') ! coefficients from adpt. dyn.
  do k1=1,n
     do k2=1,n
!        read (10,*) b1(k1,k2)
43     XP=(RAN3(IW)*2-1.D0)*4*GW ! GAUSSIAN ASSIGNMENT
       IF(GS(XP,GW).LE.RAN3(IW)) GO TO 43
     B1(K1,K2)=XP/dsqrt(1.d0*N) ! THIS IS INSTEAD OF READING IT FROM coeff*.dat FI
!        write (10,*) b1(k1,k2)
       end do
       XP=(RAN3(IW)*2-1.D0)*XCA
       XC(K1)=XP ! CENTER OF THE CARRYING CAPACITY
  end do
  b1=0.d0
  si=sa!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xc=XCA
  !  close(10)
77 format(10f8.3)
  X0=0.D0 ! IF WE DON'T READ THE INITIAL POSITION FROM THE FILE
  X=0.d0
  CC=0.D0
  LIVE=.FALSE.
  PA=0.D0
  plotcount=0
  
  DO KC=1,MCIN ! SEEDING CLUSTERS IN THE FIRST DIMENSION
     LIVE(KC)=.TRUE.
     PA(KC)=0.5D0
     DO K1=1,N ! ON
   CC(K1,KC)=X0(K1)+(2.d0*ran3(iw)-1.d0)*CB ! RATHER RANDOM CLUSTER PLACING
    X((kc-1)*n+k1)=CC(K1,KC)  ! converting it into a single vector
     END DO
  END DO
  T=0.D0
!  TMEAS=50*dtm
  TMEAS=dtm
  TAVV=TAV
  TSEED=0.d0
  TMERGE=DTMERGE
  V=VMIN*2 ! TO PREVENT IT FROM STOPPING
  JL=MCIN ! NUMBER OF LIVE CLUSTERS
  LMERGE=DCMIN
  OPEN(UNIT=50, FILE='velocity_np.dat')!, position='append') ! writing trajectory
  OPEN(UNIT=53, FILE='cluster_dim_np.dat')!, position='append') ! writing cluster dims
  OPEN(UNIT=15, FILE='movie_np.dat') ! writing a movie

  DO WHILE (T.LE.TFIN)!.AND.V.GE.VMIN) 
    VLT=0.D0
  DO WHILE (T.LE.TAVV)  
     VAV=0.D0
     DMAV=0.D0
     SIG=0.D0
     dist=0.d0
     NM=0 ! NUMBER OF FAR CLUSTERS
     KM=0
   DO WHILE (T.LE.TMEAS)! LOOP ON A MEASUREMENT
     IF(T.GE.TMERGE) THEN ! MERGING CLOSE CLUSTERS
        TMERGE=T+DTMERGE
        KM=KM+1
     NEAR=.FALSE.
     DO KC=1,MC  ! REMOVING CLUSTERS WHICH ARE TOO CLOSE
        IF(live(kc)) THEN ! removing coinciding clusters
           NEAR(KC)=.TRUE.
              DO KC2=1,KC-1
               IF(LIVE(KC2)) THEN
                  DCC=0.D0
                  NDK=0
                    DO J1=1,N
                       DCC=DCC+(CC(J1,KC)-CC(J1,KC2))**2
                    END DO
                  DCC=DSQRT(DCC)
                  IF (DCC.LE.DCFIN.AND.NEAR(KC2)) THEN 
             NEAR(KC)=.FALSE. ! LABELING A CLOSE CLUSTER
          END IF
                  IF (DCC.LE.LMERGE) THEN !LIQUIDATING CLUSTER KC
        DO J1=1,N ! PRESERVING CENTRE OF MASS POSITION
           CC(J1,KC2)=(CC(J1,KC)*PA(KC)+CC(J1,KC2)*PA(KC2))/(PA(KC)+PA(KC2)) 
          X((KC2-1)*n+J1)=CC(J1,KC2)
        END DO
                     LIVE(KC)=.FALSE.
                     PA(KC2)=PA(KC2)+PA(KC)
                     PA(KC)=0.D0                    
                  END IF
               END IF
               END DO
            END IF
         END DO

         DO KC=1,MC ! SUMMING UP THE NUMBER OF WELL-SEPARATED CLUSTERS
            IF(NEAR(KC))NM=NM+1
         END DO
     END IF

     IF(JL.LT.MC.AND.T.GE.TSEED) THEN ! SEEDING A NEW CLUSTER
        JL=JL+1
        TSEED=T+DTSEED
25        KC1=INT(RAN3(IW)*MC)+1
        IF(.NOT.LIVE(KC1)) GO TO 25 ! FINDING THE FATHER
        KC2=1
        DO WHILE (LIVE(KC2)) ! PLACE FOR OFFSPRING
           KC2=KC2+1
        END DO
        DSEED=0.D0
        LSEED=0.D0
        DO J1=1,N ! MAKING THE DISTANCE EXACTLY CB
  3     XP=(RAN3(IW)*2-1.D0)*4*CB ! GAUSSIAN ASSIGNMENT
       IF(GS(XP,CB).LE.RAN3(IW)) GO TO 3
        DSEED(J1)=XP       
       !           DSEED(J1)=2.d0*RAN3(IW)-1.d0
           LSEED=LSEED+DSEED(J1)*DSEED(J1)
        END DO
        LSEED=DSQRT(LSEED)
!         LSEED=DSQRT(LSEED)/CB
       DO J1=1,N
          CC(J1,KC2)=CC(J1,KC1)+(1.d0-mu)*DSEED(J1) ! preserving CM
!          CC(J1,KC2)=CC(J1,KC1)+(1.d0-mu)*DSEED(J1)/LSEED ! preserving CM
          CC(J1,KC1)=CC(J1,KC1)-mu*DSEED(J1) !
!          CC(J1,KC1)=CC(J1,KC1)-mu*DSEED(J1)/LSEED !
          X((KC2-1)*N+J1)=CC(J1,KC2)
          X((KC1-1)*N+J1)=CC(J1,KC1)
      END DO
      PA(KC2)=PA(KC1)*mu ! SPLITTING THE ANCESTOR
      PA(KC1)=PA(KC1)*(1.d0-mu) 
      LIVE(KC2)=.TRUE.
   END IF
   LMERGE=LSEED

   
   ME=0.D0 ! POPULATION DYNAMICS
   BRT=0.D0 ! CLUSTER'S BIRTH RATE 
     DO KC1=1,MC  ! COMPUTING THE TABLE OF DEATH RATES FOR CLUSTER KC1
        IF(LIVE(KC1)) THEN
           BRT(KC1)=brate(kc1,cc,mc,n,sb,coeffbirth)
        DO KC2=1,MC ! affected by the cluster kc2
           IF(LIVE(KC2)) THEN
    ME(KC1,KC2)=comp(kc1,kc2,cc,mc,n,b1,si,gauss,xc) !TABLE
           END IF
        END DO
        END IF
     END DO
     PC=PA

     DO JE=1,NE
      CALL DERIVS_POP (PC,DP,MC,ME,LIVE,BRT) 
1     CALL RK_POP (PC,DP,MC,DTP,PN,ME,LIVE,BRT)
     PA=PA+PN
     PC=PN
     END DO
     PA=PA/(NE+1)
     S7=0.D0
     JL=MC
     PM=POPMIN ! MAXIMUM POPULATION
     PT=0.D0 ! TOTAL POPULATION
     DO KC=1,MC
        S7=S7+DP(KC)*DP(KC)
        IF(PA(KC).LT.POPMIN) THEN 
           JL=JL-1
           pa(kc)=0.d0
           LIVE(KC)=.FALSE.
        END IF
        PT=PT+PA(KC)
        IF (PA(KC).GE.PM) PM=PA(KC)
    END DO
     S7=DSQRT(S7)
     CALL DERIVS_MOV(x,dx,pa,n,mc,b1,si,gauss,LIVE,PM,SB,coeffbirth,XC)
     CALL RK_MOV(x,dx,n,dt,xn,mc,b1,pa,si,gauss,LIVE,PM,SB,coeffbirth,XC)
!     T=T+DT/PM
     T=T+DT
     DV=(XN-X)/DT 
     X=XN
     V=0.D0
     DMA=0.D0
     X1M=0.D0
     X2M=0.D0
     DISP=0.D0
   do kc1=1,mc
      WW=0.D0
      PW=0.D0
     do j1=1,n
    cc(j1,kc1)=x((kc1-1)*n+j1) ! reconstructing the coordinate array
    WW=WW+DV((kc1-1)*n+j1)**2
    PW=PW+1.D0
    VLT((kc1-1)*n+j1)=VLT((kc1-1)*n+j1)+DV((kc1-1)*n+j1)*PA(KC1)/PT*DT
    X1M(J1)=X1M(J1)+x((kc1-1)*n+j1)*PA(KC1)/PT
    X2M(J1)=X2M(J1)+x((kc1-1)*n+j1)**2*PA(KC1)/PT
     end do
     WW=DSQRT(WW)*PA(KC1)/PT
     PW=PW*PA(KC1)/PT
     V=V+WW
     DMA=DMA+PW
   end do
   
   DO J1=1,N
   DISP=DISP+X2M(J1)-X1M(J1)**2
   END DO
   VAV=VAV+V*DT
   DMAV=DMAV+DMA*DT
   IF(DISP.GE.0.D0) SIG=SIG+DSQRT(DISP)*DT
   dw=0.d0
   do kc1=1,mc
      d1=0.d0
    do kc2=1,mc
       if(live(kc1).and.live(kc2)) then
          lw=0.d0
        do j1=1,n
    lw=lw+(x((kc1-1)*n+j1)-x((kc2-1)*n+j1))**2
        end do
        lw=dsqrt(lw)
       end if
       d1=d1+lw*pa(kc2)/pt
    end do
    dw=dw+d1*pa(kc1)/pt
   end do
      DIST=DIST+DW*dt/2 ! since each distance is counted twice
   END DO   ! END CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN MEASUREMENTS
   
    KCFAR=1.D0*NM/KM ! averaged over tmeas
    VAV=VAV/DTM
    DMAV=DMAV/DTM
    SIG=SIG/DTM
    DIST=DIST/DTM
    TMEAS=TMEAS+DTM

    WRITE(15,*) (plotcount*100,K1=1,4)
     plotcount=plotcount+1

     CLDM=0
          WRITE(*,*)

        DDD=0
        DO KC=1,MC
        IF(LIVE(KC)) THEN
           DO K1=2,N
              IF (DABS(CC(K1,KC)).GE.SB*2.D0) THEN
                 DDD(K1)=DDD(K1)+1
              END IF   
           END DO
        END IF
        END DO


        JL=0 ! COUNTING LIVE CLUSTERS AGAIN
        DMTOT=0
        DO KC=1,MC
            IF(LIVE(KC)) THEN
!!! this needs to be changed           CLDM(BR(KC))=CLDM(BR(KC))+1
               DMC=N
               DML=1
               DO K1=2,N
                  IF (DABS(CC(K1,KC)).LE.SB*2.D0) then !.AND.DDD(K1).LE.4) THEN
                     DMC=DMC-1
                     DML(K1)=0
                  END IF
               END DO
                  WRITE(*,*)(DML(K1),K1=1,N)
                  WRITE(*,*)(CC(K1,KC),K1=1,N)
                  CLDM(DMC)=CLDM(DMC)+1              
                  JL=JL+1
                  DMTOT=DMTOT+DMC
            WRITE(15,*)(CC(J1,KC),J1=1,2),PA(KC),T
            END IF
         END DO

         DO K1=1,N
            AD(K1,PLOTCOUNT)=AD(K1,PLOTCOUNT)+CLDM(K1)
            VD(K1,PLOTCOUNT)=VD(K1,PLOTCOUNT)+CLDM(K1)**2
         END DO
         
 
     WRITE(50,*) T,JL,DMTOT*1.d0/JL,VAV
     WRITE(*,*) T,JL,DMTOT*1.d0/JL,V
     DO J1=1,N
     IF(CLDM(J1).GE.1) WRITE(*,*) J1,CLDM(J1)
     END DO
     WRITE(53,*) T,(CLDM(J1),J1=1,N)

  END DO
  TAVV=TAVV+TAV
    VLT=VLT/TAV ! COMPUTING COARSE-GRAINED DISPLACEMENT
      VLTM=0.D0
     DO KC=1,MC
      WUW=0.D0
     do j1=1,n
    WUW=WUW+VLT((kc-1)*n+j1)**2
     end do
     WUW=DSQRT(WUW)
     VLTM=VLTM+WUW
     END DO
!     WRITE(*,*)T-TAV/2,V,VAV,VLTM
  END DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of the run
  WRITE(*,*) 'FINAL T,V', T, VAV, JL
  write(50,*)'&'
  write(50,*)'&', 'coeffbirth  ', coeffbirth, ' Sb ', sb, 'offset', xc(1)
  write(53,*)'&', 'coeffbirth  ', coeffbirth, ' Sb ', sb, 'offset', xc(1)
  CLOSE(50)  
  CLOSE(53)  
  CLOSE(15)
    DO KC=1,MC  ! FINALLY REMOVING CLUSTERS WHICH ARE TOO CLOSE
        IF(live(kc)) THEN ! removing coinciding clusters
              DO KC2=1,KC-1
               IF(LIVE(KC2)) THEN
                  DCC=0.D0
                    DO J1=1,N
                       DCC=DCC+(CC(J1,KC)-CC(J1,KC2))**2
                    END DO
                  DCC=DSQRT(DCC)
                  IF (DCC.LE.DCFIN) THEN !LIQUIDATING CLUSTER KC
                     LIVE(KC)=.FALSE.
                     PA(KC2)=PA(KC2)+PA(KC)
                     PA(KC)=0.D0
                  END IF
                  END IF
               END DO
            END IF
         END DO
  OPEN(UNIT=20,FILE='last_snapshot.dat')
  KC=0
      DO KC1=1,MC
         if(live(kc1)) THEN
      KC=KC+1
      WRITE(20,*)(CC(J1,KC1),J1=1,N),PA(KC1)
         end if
     END DO
  CLOSE(20)
  write(*,*) kc
  do kc1=1,mc
  IF(LIVE(KC1)) write(*,*) pa(kc1)
  end do
  nctot=nctot+jl ! clusters that are alive and have not converged
  ncfar=ncfar+kcfar ! clusters that are further apart
  vfin=vfin+vav ! arithmetic average of total per capita speed
  dmfin=dmfin+dmav
  IF(V.GE.VGLASS) THEN 
     PER=PER+1
 ELSE IF (V.LT.VGLASS.AND.V.GE.VMIN) THEN
      GLASS=GLASS+1
 ELSE IF (V.LT.VMIN) THEN
      CONV=CONV+1
  END IF
  write(*,*)'STATS', PER, GLASS, CONV, vfin/NR, dmfin/NR, 1.*nctot/NR, 1.*ncfar/NR 
END DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of cycle on repetitions
!!$OPEN(UNIT=25,FILE='stat_neq.dat')!, position='append')
!!$  write(25,*)N,mc,sa,PER, GLASS, CONV, vfin/NR, dmfin/NR, 1.*nctot/NR, 1.*ncfar/NR 
!!$  close(25)


  OPEN(UNIT=27,FILE='clusters_stat.dat')!, position='append')
  OPEN(UNIT=29,FILE='clusters_stat_simple.dat')!, position='append')
  AD=AD/NRTOT
  VD=VD/NRTOT
  DO J1=1,PLOTCOUNT
     DO K1=1,N
        VD(K1,J1)=VD(K1,J1)-AD(K1,J1)**2
        IF(VD(K1,J1).GT.0.D0) THEN
           VD(K1,J1)=DSQRT(VD(K1,J1))
        ELSE IF (VD(K1,J1).LE.0.D0) THEN
           VD(K1,J1)=0.D0
        END  IF
      END DO
  write(27,*) J1*DTM,(AD(K1,J1),VD(K1,J1),K1=1,N) 
  write(29,*) J1*DTM,(AD(K1,J1),K1=1,N) 
   END DO
  close(27)
  close(29)

  

  open(unit=12, file='rand.dat')
  write(12,*)nint(-ran3(iw)*1.d+6)
  close(12)   
END PROGRAM CHAOS

!---------------------------------------------------------------------
  function comp(c1,c2,cc,mc,n,b1,si,gauss,xc) ! competing effect
  integer n ! of k on i with ith carrying capacity
  integer j1,j2,j3,c1,c2,mc,ngauss,br(mc)
  double precision cc(n,mc),s,w,comp,sa,dy
  double precision b1(n,n),si(n),xc(n)
  logical GAUSS
        W=0.D0
        S=0.D0
   do j1=1,n ! cluster c1 is affected by c2
      dy=cc(j1,c2)-cc(j1,c1)
     do j2=1,n
        W=W-b1(j1,j2)*(cc(j1,c1)-cc(j1,c2))*cc(j2,c2)
        end do
     if(GAUSS) W=W-DY*DY/(2*si(j1)*si(j1)) ! ADDING GAUSSIAN PART TO A PARTICULAR. D
     S=S-((cc(j1,c1)-xc(j1))**4)/4 ! it's present independent of dm
    end do
  comp=DEXP(W-S)
  return
  end function comp
!-----------------------------------------------------------------------
  function brate(kc,cc,mc,n,sb,cb)
    implicit none
  integer mc,n,kc,k  
  double precision cc(n,mc),sb,cb,a,x,brate
  a=1.d0
  do k=2,n
  x=cc(k,kc)   
  a=a*(dexp(-x*x/(2.d0*sb*sb))*(1.d0-cb) + cb)   
  end do
  brate=a
  end function brate  
  !--------------------------------------------------------------------
  SUBROUTINE RK_POP(y,dydx,n,h,yout,me,live,brt)
! Y=Y_i, DYDX= Y'_i, n-dim,  h -step, YOUT=Y_i+1, derivs calculates derivs
! FOR AUTONOMOUS SYSTEM WITHOUT TIME DEPENDENCE
      INTEGER n,NMAX
      DOUBLE PRECISION h,dydx(n),y(n),yout(n),me(n,n),brt(n)
      logical live(n)
      EXTERNAL derivs_pop
      PARAMETER (NMAX=400)
      INTEGER i
      DOUBLE PRECISION h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5d0
      h6=h/6.d0
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs_pop(yt,dyt,n,me,live,brt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs_pop(yt,dym,n,me,live,brt)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs_pop(yt,dyt,n,me,live,brt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))

14    continue
      return
      END SUBROUTINE RK_POP

!---------------------------------------------------------------------
       SUBROUTINE derivs_pop(y,dy,n,me,live,brt)
! Y=Y_i, DY= Y'_i, n-dim
      INTEGER n,k1,k2
      DOUBLE PRECISION y(n),dy(n),me(n,n),brt(n)
      LOGICAL LIVE(N)
      do k1=1,n
         dy(k1)=0.d0
         IF(LIVE(K1))THEN
         do k2=1,n
            IF(LIVE(K2)) THEN
            dy(k1)=dy(k1) - me(k1,k2)*y(k2)
            END IF
         end do ! changeable below, birth rate now set as 1/dim
            dy(k1)=(dy(k1)+brt(k1))*y(k1)  
          end if
            end do
      RETURN
       END SUBROUTINE DERIVS_POP

!--------------------------------------------------------------------
      SUBROUTINE RK_MOV(y,dydx,n,h,yout,mc,b1,pa,si,gauss,LIVE,PM,SB,CB,XC)
! Y=Y_i, DYDX= Y'_i, n-dim,  h -step, YOUT=Y_i+1, derivs calculates derivs
! FOR AUTONOMOUS SYSTEM WITHOUT TIME DEPENDENCE
      INTEGER n,NMAX
      DOUBLE PRECISION dydx(n*mc),y(n*mc),yout(n*mc),pa(mc),si(n)
      double precision b1(n,n),XC(N),sb,cb
      LOGICAL GAUSS, LIVE(MC)
      EXTERNAL derivs_mov
      PARAMETER (NMAX=2000)
      INTEGER i
      DOUBLE PRECISION H, h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX),PM,SA
      hh=h*0.5d0
      h6=h/6.d0
      do 11 i=1,n*mc
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs_mov(yt,dyt,pa,n,mc,b1,si,gauss,LIVE,PM,SB,CB,XC)
      do 12 i=1,n*mc
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs_mov(yt,dym,pa,n,mc,b1,si,gauss,LIVE,PM,SB,CB,XC)
      do 13 i=1,n*mc
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs_mov(yt,dyt,pa,n,mc,b1,si,gauss,LIVE,PM,SB,CB,XC)
      do 14 i=1,n*mc
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))

14    continue
      return
      END SUBROUTINE RK_MOV


!-------------------------------------------------------------------- 
  SUBROUTINE DERIVS_MOV(y,dy,pa,n,mc,b1,si,gauss,LIVE,PM,SB,CB,XC)
       implicit none
  integer j1,j2,j3,c1,c2,mc,n,nmax,mcmax
  double precision y(n*mc),dy(n*mc),pa(mc)
  double precision b1(n,n),si(n),XC(N)
  parameter(nmax=10, mcmax=400) ! largest dimension and number of clusters
  double precision me(mcmax,mcmax),cc(nmax,mcmax),sg(nmax,mcmax,mcmax)
  double precision s,w,comp,ddy,aa,bb,v,pm,sb,cb,ww
  logical GAUSS,LIVE(MC)
  do c1=1,mc
     do j1=1,n
        cc(j1,c1)=y((c1-1)*n+j1) ! reconstructing the coordinate array
     end do
  end do
  me=0.d0
  do c1=1,mc
     if(live(c1)) then
     do c2=1,mc
        if(live(c2)) then
       W=0.D0
       S=0.D0
   do j1=1,n ! cluster c1 is affected by c2
     SG(J1,C1,C2)=0.D0
     DDY=cc(j1,c2) - cc(j1,c1)
     do j2=1,n
        W=W+b1(j1,j2)*DDY*cc(j2,c2)
        bb=b1(j1,j2)*cc(j2,c2)
        SG(J1,C1,C2)=SG(J1,C1,C2)+bb
        end do
     if(GAUSS) THEN ! ADDING GAUSSIAN PART TO CERTAIN D
         W=W-DDY*DDY/(2*SI(J1)*SI(J1)) ! ADDING GAUSSIAN PART
         SG(J1,C1,C2)=SG(J1,C1,C2)-DDY/SI(J1)/SI(J1) ! SELECTION GRAD.
     end if
    S=S-((cc(J1,c1)-XC(J1))**4)/4 ! adding carrying capacity
    SG(J1,C1,C2)=SG(J1,C1,C2)-(cc(j1,c1)-XC(J1))**3    
end do
  ME(C1,C2)=DEXP(W-S) ! C2 AFFECTS C1
  end if
end do
end if
end do

  do c1=1,mc
     do j1=1,n
        ww=cc(j1,c1)
        v=0.d0
        if(j1.ge.2) v= -(1.D0-CB)*ww/(sb*sb)*dexp(-ww*ww/(2*sb*sb)) !from birth rate
      do c2=1,mc
         v=v+me(c1,c2)*PA(C2)*SG(J1,C1,C2)
      end do
!      v=v*PA(C1)/PM ! mutation rate prop. to population
      v=v*PA(C1) ! mutation rate prop. to population
      dy((c1-1)*n+j1)=v ! assinging the derivative
   end do
  end do
  RETURN
     END SUBROUTINE DERIVS_MOV
         
!---------------------------------------------------------------------
      function ran3(idum)
      integer idum
      integer mbig,seed,mz
      double precision ran3,fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
       integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
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
      return
      end function ran3

!---------------------------------------------------------------------
