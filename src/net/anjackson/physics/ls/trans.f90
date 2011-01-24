! --> trans.f90
! Hard sphere HCPFCC lattice-switch code:
! Indirect addressing, system size can be any N = (n*6)**3.
! NVT or NPT.
!
! v1.29:(13th May 1999) Andrew N Jackson, original code from GJA & NBW.
! ---------------------------------------------------------------------------------------
! Code Switches (via cpp):
! ---------------------------------------------------------------------------------------
#define _false_          0       /* false flag                                         */
#define _true_           1       /* true flag                                          */
! ---------------------------------------------------------------------------------------
! Simulation Options:
#define SystemSize       2       /* System size, NM = (SystemSize*6)**3.              */
#define MoveType         2       /* Particle move type;                               */
                                 /*  1 = Top-Hat, 2 = Random Walk, 3 = Gaussian,      */
                                 /*  4 = RW + cubic cut-off, 5 = RW + spherical cutoff*/
#define ConstPres        _false_ /* Which ensemble; 0 = NVT, 1 = NPT.                 */
#define FixedAspectRatio _false_ /* Fix the aspect ratio of the system (NPT only).    */
#define LattSwitch       _true_  /* Enable the lattice-switch moves                   */
#define RandomSwitch     _false_ /* Makes the L-S numbering random                    */
#define RandPlaneSwitch  _false_ /* Makes the L-S plane-mapping random                */
#define LockInGateway    _false_ /* Lock the system in the M=0 macrostate             */
#define UseRndAtm        _true_  /* Random or Typewriter atom selection               */
#define LockedPos        _false_ /* Lock the spheres in place                         */
#define CheckForDeath    _false_ /* Check if any spheres have move to a different site*/
#define StayInCoMFrame   _false_ /* Keep in the centre-of-mass frame                  */
#define NNsCompare       _false_ /* Compare 1st and 2nd nearest neighbour simulation  */
#define Check1NNsOnly    _false_ /* Only check 1st (not 2nd) NNs for overlaps         */
#define UseDodeca        _false_ /* Uses hard-dodecahedra instead of spheres          */
#define PolyDisperse     _true_ /* Simulate a polydisperse system.                   */

! Transition probability matrix calc:
#define TransProb        _false_ /* Macrostate transition probability matrix analysis.*/
#define TP_sml_max       60      /* Small system goes out to  M=+/-60                 */
#define TP_med_max       350     /* Medium system goes out to M=+/-350                */
#define TP_lrg_max       1500    /* Large system goes out to  M=+/-1500               */

! Simulation parameters:
#define InitPhase        2       /* Initial phase:  1 = hcp, 2 = fcc                  */
#define MeltCheckTime	 500     /* Melt checking period (in MCS)                     */
#define ConfSaveTime	 5000    /* Configuration saving period (in MCS)              */

! Output Options:
#define ClosestDist      _false_ /* Collect data for closest approach distribution    */
#define NNSepDist        _false_ /* Collect data for NN seperation distribution       */
#define LogCoMVec        _false_ /* Log the center of mass vector                     */
#define LogPlaneDispl    _false_ /* Log the plane disps and the squeeze factor        */
#define CalcDodecVecH    _false_ /* Bin the dodecahedra perp/para pressure data       */
#define OutputXYZConf    _true_  /* Output a file suitable for XMol in SaveConfig     */
#define OutputRotM       _false_ /* Output the rotated-u order parameter              */
#define SpewLocalDist    _false_ /* Dumps the displacements to a file for plotting    */
#define CalcNNDist       _false_ /* Calcs P(u.R) for NNs in both lattices              */
#define CalcNNChange     _false_ /* Calcs # of 1st NNs which change due to the switch  */

! ---------------------------------------------------------------------------------------
! Notes:
!
! - Contigious simulation of polydisperse systems may require saving the radii.(?)
! - Optimize RNG by creating 3*NM RDNS at once.
!
! ---------------------------------------------------------------------------------------
! Global variable storage for the simulation code:
! ---------------------------------------------------------------------------------------
! Define all the system variables in a place where all routines can see them: 
! Added by Anj to make life easier.
  MODULE SystemDef
    IMPLICIT NONE
    INTEGER smallNM, mediumNM,largeNM
    PARAMETER (smallNM=216, mediumNM=1728, largeNM=5832)
    INTEGER NM,NumWHTS
! Set system size here according to SystemSize:
#if   SystemSize == 1
    PARAMETER (NM=smallNM, NumWHTS=200)
#elif SystemSize == 2
    PARAMETER (NM=mediumNM, NumWHTS=1000)
#elif SystemSize == 3
    PARAMETER (NM=largeNM, NumWHTS=10000)
#endif
! Local overlap order parameter:
    INTEGER true,false
    PARAMETER (true=1)
    PARAMETER (false=0)
    INTEGER NOVER(NM,19),IOVER(NM)
! Neighbour lists:
    INTEGER NoMoreNeighbours
    PARAMETER (NoMoreNeighbours=99999)
    INTEGER N0LIST(NM,7,2),NLIST(NM,19,2),INDL(NM,19,2)
! Comparing 1st and 2nd NN simulations:
    INTEGER NLInfo(NM,19,2),N0LInfo(NM,7,2)
    INTEGER OldIJover,Old2NNs,New2NNs,Current2NNs,Count2NNs,Count2NNsMCS
! MCMC weights:
    REAL*8 WHTS(NumWHTS)
! Lattice-sites:
    REAL*8 XLAT(NM,2),YLAT(NM,2),ZLAT(NM,2)
! Displacements:
    REAL*8 DDX(NM),DDY(NM),DDZ(NM)
! Trail displacements
    REAL*8 EPSX,EPSY,EPSZ
! Box dimensions:
    REAL*8 AL0,BL0,CL0,AL,BL,CL
! Miscellaneous:
    INTEGER NCOUNT,ILAT,ICLAT
    INTEGER ISTART,ISEED,NEQBM,NSAMP,NMOVES
    REAL*8 DIAMETER,DIAMETER2,STEP,cond,cutoff,cutoff2,halfcutoff
    REAL*8 Pres,Volstep
    INTEGER VolTried,VolAcc1,VolAcc2
    LOGICAL overlapped
    INTEGER IZERO,OVOLD,OVNEW,NOLAP
    REAL*8 IJDr2,Uperp2,Upara
    REAL*8 density,pi,maxnns
!Anj: Variables for collecting the `closest approach' distribution:
    REAL*8 CloseDr2
    INTEGER CloseI,CloseJ
!Anj: Variable to hold probabilities for gaussian position selection scheme:
    REAL*8 GPo,GPn,Rold,Rnew,GPratio
!Anj: Switches to allow sections of the code to be turned on and off:
    LOGICAL TransProbFin
!Anj: Melt checking flag:
    LOGICAL deathbymeltdown
!Anj: Look-up table for hard dodecahedra:
    REAL*8 lnnx(NM,19,2),lnny(NM,19,2),lnnz(NM,19,2)
!Anj: Histogram handles for HS pressure:
    INTEGER upara_hh,uperp_hh
!Anj: Look-up table for each plane of which plane (below=-1 or above=+1) changes:
    INTEGER PlaneChanges(SystemSize*6),WhichPlane(NM),WhichXRow(NM),WhichYRow(NM)
!Anj: The copy-space for the U rotations:
    REAL*8 Dcp(NM,3)
!Anj: NN dist (u.R) histogram collection handles:
    INTEGER nnd_hists(10)
!Anj: Variables for the polydispersity calulation:
    REAL*8 polyd,radii(NM)
  END MODULE SystemDef


! -----------------------------------------------------------------------------
! The Core Code:
! -----------------------------------------------------------------------------
  PROGRAM fcchcp
    USE SystemDef
    IMPLICIT NONE
    REAL*8 delx,dely,delz
    REAL RVEC(100)
    REAL*8 Zoutfac
    INTEGER I,J,IMOVE,NACT,NAC,NFLIP,NREJ,JX,JCOUNT,NCHK,NJX,IJover
    INTEGER ILOOP,ios,maxM,outM
    LOGICAL acceptmove,overlaps2NN
!Anj: Gaussian position function type declaration:
    REAL*8 GaussPos,GaussProb
!Anj: Analysis time step flags:
    INTEGER MCTStep, NTStep
!Anj: Trial variables for closest approach code:
    REAL*8 TryCDr2
    INTEGER TryCI,TryCJ
    LOGICAL TryC
!Anj: Neighbour vector storage:
    REAL*8 DUvec(18,2)
!Anj: Proper new-vector creation:
    REAL*8 costheta,sintheta,phi,rmagn,udx,udy,udz,udr
! For histograms:
    INTEGER Init_Histogram

!Anj: TransProb Analysis finished flag:
    TransProbFin = .FALSE.

! Start in one or other phase. Note ilat=1 is HCP, ilat=2 is FCC
! ICLAT is reference to the conjugate lattice (FCC->HCP or HCP->FCC)
    ILAT  = InitPhase
    ICLAT = 3 - ILAT

! Read the defining parameters from the standard input:
#if UseDodeca == _false_
 WRITE(*,FMT='(A,I4)') ' Hard-Sphere Lattice-Switch Monte Carlo code v1.29: NM=',NM
#else
 WRITE(*,FMT='(A,I4)') ' Hard-Dodecahedra Lattice-Switch Monte Carlo code v1.29: NM=',NM
#endif
#if UseDodeca == _true_ && Check1NNsOnly == _false_
 WRITE(*,*) 'WARNING: Dodecahedra implemented with 1st AND 2ND nearest neighbours!'
#endif
#if UseDodeca == _true_ && MoveType != 2
 WRITE(*,*) 'WARNING: Dodecahedra implemented without proper random walk!'
#endif
#if ConstPres == _true_
    WRITE(*,*)'Enter pressure,volstep,density,nmoves,step,nsamp,neqbm,istart,iseed'
    READ(*,*) Pres,Volstep,density,nmoves,step,nsamp,neqbm,istart,iseed
#elif MoveType == 4 || MoveType == 5
    WRITE(*,*)' Enter density, nmoves, step, cut-off, nsamp, neqbm, istart, iseed'
    READ(*,*) density,nmoves,step,cutoff,nsamp,neqbm,istart,iseed
#elif PolyDisperse == _true_
    WRITE(*,*)' Enter density, polydispersity, nmoves, step, nsamp, neqbm, istart, iseed'
    READ(*,*) density,polyd,nmoves,step,nsamp,neqbm,istart,iseed
#else
    WRITE(*,*)' Enter density, nmoves, step, nsamp, neqbm, istart, iseed'
    READ(*,*) density,nmoves,step,nsamp,neqbm,istart,iseed
#endif

! Initialise some variables and counters
    pi = 2.0d0*ASIN(1.0d0)
    IZERO=NumWHTS/2
    NAC = 0
    NCOUNT = 0
    NFLIP = 0
    NREJ = 0
    NACT = 0
    VolTried = 0; VolAcc1 = 0; VolAcc2 = 0
    DIAMETER = density**(1.0d0/3.0d0)
    DIAMETER2 = DIAMETER*DIAMETER
    maxM = (numWHTS/2) - 1
    Zoutfac = DBLE(NM)*density/SQRT(2.0d0)
#if MoveType == 4 || MoveType == 5
    cutoff2 = cutoff*cutoff/4.0d0
    halfcutoff = cutoff/2.0d0
#endif
#if NNsCompare == _true_
    Current2NNS = 0
    Count2NNs = 0
    Count2NNsMCS = 0
#endif

! Initialise RNG
    CALL RMARIN(iseed)

! Set up perfect crystal lattices for HCP and FCC, and the lattice-switch:
    CALL MakeLattices
#if CalcNNChange == _true_
    CALL CountNNChanges    
#endif

#if RandomSwitch == _true_
! Randomize the lattice-switch particle mapping:
    CALL MakeLSTotallyRandom
#endif
#if RandPlaneSwitch == _true_
! Randomize the lattice-switch plane mapping:
    CALL MakeLSPlanesRandom
#endif
#if CalcNNChange == _true_
    CALL CountNNChanges    
#endif

! Make neighbour lists for both HCP and FCC
    CALL MAKELIST

#if TransProb == _true_
! Initialise transition probability measuring code:
    CALL TPInitialise
#endif

#if LogPlaneDispl == _true_
    CALL InitPlaneSep
#endif

#if SpewLocalDist == _true_
          CALL InitDisplOutput
#endif

#if CalcNNDist == _true_
          CALL InitNNDistCalc
#endif

#if CalcDodecVecH == _true_
          upara_hh = Init_Histogram(-1.0d0,4.0d0,100)
          uperp_hh = Init_Histogram(-1.0d0,4.0d0,100)
#endif

! Initialise the values of the displacements and the overlap array:
    DO I = 1,NM
       DDX(I) = 0.0
       DDY(I) = 0.0
       DDZ(I) = 0.0
       DO J=1,maxnns
          NOVER(I,J)=false
       END DO
    END DO

! Read in old configuration if instructed to do so:
    IF(ISTART.EQ.1) THEN
       WRITE(6,*) 'READING IN OLD CONFIG'
       CALL LoadConfig
    ELSE
       WRITE(6,*) 'COLD START'
    ENDIF

#if PolyDisperse == _true_
! Create the radii array for the polydisperse simulation:
    CALL Initialise_Polydisperse_Radii
#endif

#if ClosestDist == _true_
! Initialise closest pair distribution code:
    CALL FindClosestPair
#endif

! Write header to the output stream:
#if ConstPres == _true_
    WRITE(*,*) 'Pressure=',Pres,' Volstep=',Volstep
#endif
    WRITE(*,*)'particle diameter=',diameter,' nmoves=',nmoves
#if MoveType == 4 || MoveType == 5
    WRITE(*,*)'step=',step,'cut-off=',cutoff,' nsamp=',nsamp,' neqbm=',neqbm
#else
    WRITE(*,*)'step size=',step,' nsamp=',nsamp,' neqbm=',neqbm
#endif
    WRITE(*,*)'istart=',istart,' iseed=',iseed
#if ConstPres == _true_
    WRITE(*,FMT='(A,F12.8,A,F12.8,A,F12.8)') ' NPT: AL0xBL0xCL0 = ',AL0,'x',BL0,'x',CL0
    WRITE(*,FMT='(A,F12.8,A,F12.8,A,F12.8)') ' NPT: ALxBLxCL = ',AL,'x',BL,'x',CL
    WRITE(*,FMT='(A,F12.6)') ' NPT: V = ',AL*BL*CL
#else
    WRITE(*,FMT='(A,F12.8,A,F12.8,A,F12.8)') ' NVT: ALxBLxCL = ',AL,'x',BL,'x',CL
#endif
#if NNsCompare == _true_
    WRITE(*,*) 'NNC: Current2NNs = ',Current2NNs
#endif

! Read in weights, if file::'weights' doesn't exist then
! initialize the weights to 1.0:
    OPEN (UNIT=55,FILE='weights',STATUS='OLD',IOSTAT=ios)
    IF(ios.EQ.0) THEN
       DO I=1,NumWHTS
          READ(55,*,IOSTAT = ios) WHTS(I)
       END DO
       CLOSE(UNIT=55)
    ELSE
       DO I=1,NumWHTS
          WHTS(I) = 1.0
       END DO
    END IF
! Write the weights to the standard output :
    DO I=1,NumWHTS
       WRITE(*,*) 'W',WHTS(I)
    END DO

! Main loop
    DO IMOVE=1,NMOVES

!Anj: Periodic Configuration Save:
       IF(MOD(IMOVE,ConfSaveTime).EQ.0) CALL SaveConfig

#if CheckForDeath == _true_
!Anj: Periodic Melt check:
       IF(MOD(IMOVE,MeltCheckTime).EQ.0) THEN
          deathbymeltdown = .FALSE.
          CALL MeltCheck(deathbymeltdown,IMOVE)
          IF (deathbymeltdown) GOTO 999
       END IF
#endif

! Sample the data every IMOVE = N * NSAMP loops:
       IF((IMOVE.GT.NEQBM.AND.MOD(IMOVE,NSAMP).EQ.0)) THEN
          WRITE(*,FMT='(A,I12,I8,I3)') ' A ',IMOVE,NCOUNT,ILAT
          outM = (ILAT*2-3)*NCOUNT
#if ConstPres == _true_
          WRITE(*,FMT='(A,I12,I8,F11.5,4F10.6)') ' V ',IMOVE,outM,AL*BL*CL,AL,BL,CL,Zoutfac/(AL*BL*CL)
#endif
#if ClosestDist == _true_
          WRITE(*,FMT='(A,I11,I8,2I5,F12.7,E12.5)') ' CP ',IMOVE,outM,CloseI,CloseJ,SQRT(CloseDr2),SQRT(CloseDr2)-DIAMETER
#endif
#if NNSepDist == _true_
          CALL MeasureNNSeperations(IMOVE)
#endif
#if LogCoMVec == _true_
!Anj: Log CoM vector:
          CALL LogCoM
#endif
#if LogPlaneDispl == _true_
!Anj: Log Plane displacement and seperation:
          CALL LogPlaneSep
#endif
#if NNsCompare == _true_
          WRITE(*,*) 'NN',IMOVE,outm,Current2NNs
#endif
#if OutputRotM == _true_
          CALL TestRotatedM(IMOVE)
#endif
#if SpewLocalDist == _true_
          CALL OutputDisplacements
#endif
#if CalcNNDist == _true_
          CALL CalcNNDistributions
#endif
       ENDIF

#if NNsCompare == _true_
       IF (Current2NNs > 0) Count2NNsMCS = Count2NNsMCS + 1
#endif

!  Loop over all atoms
       DO ILOOP = 1,NM

#if LattSwitch == _true_
! If NCOUNT=0 permit change of basis:
          IF(NCOUNT==0) THEN
             ICLAT=ILAT
             ILAT=-ILAT+3
             IF(IMOVE.GT.NEQBM) NFLIP = NFLIP+1
          ENDIF
#endif

#if ConstPres == _true_ && FixedAspectRatio == _false_
! Anj: If NPT, change the volume (freely) once every sweep (on average):
          CALL RanMar(RVEC,1)
          IF ((INT(NM*RVEC(1))+1).EQ.1) CALL NPTChangeVolumeFree
#endif
#if ConstPres == _true_ && FixedAspectRatio == _true_
! Anj: If NPT, change the volume (not the aspect) once every sweep (on average):
          CALL RanMar(RVEC,1)
          IF ((INT(NM*RVEC(1))+1).EQ.1) CALL NPTChangeVolumeFixed
#endif

#if UseRndAtm == _true_
! Anj: Random atom selection:
          CALL RANMAR(RVEC,1)
          I = INT(NM*RVEC(1)) + 1
#else
! Anj: Typewriter atom selection:
          I = ILOOP
#endif

! Periodically dump the local vectors::
#if CalcDodecVecH == _true_ 
             IF (IMOVE>NEQBM .AND. MOD(IMOVE,NSAMP)==0) THEN
                DO JCOUNT=1,12
                   J = NLIST(I,JCOUNT,ILAT)
                   CALL Interaction(I,J,JCOUNT,ILAT,.FALSE.,IJover)
!                   IF ( Upara < -1.0d0 + 0.1d0 .AND. J<I) WRITE(*,*) 'DV ',Upara,Uperp2
                   IF (I<J) THEN
                      CALL Add_To_Histogram(upara_hh,Upara)
                      CALL Add_This_To_Histogram(uperp_hh,Upara,Uperp2)
                   END IF
                END DO
                IF (I==NM) THEN
                   CALL Output_Histogram(upara_hh,'upara.h.dat')
                   CALL Output_Histogram(uperp_hh,'uperp.h.dat')
                END IF
             END IF
#endif

! NOTE: STEP is absolute displacement, but
! ESPX and DDX are fractional changes (relative to box size)
#if     MoveType == 1
!Anj: Top Hat move selection technique:
          CALL RANMAR(RVEC,3)
          EPSX = (DBLE(RVEC(1))-0.5d0)*STEP/AL
          EPSY = (DBLE(RVEC(2))-0.5d0)*STEP/BL
          EPSZ = (DBLE(RVEC(3))-0.5d0)*STEP/CL
#elif MoveType == 2
!Anj: Random walk selection technique:
#if UseDodeca == _false_
          CALL RANMAR(RVEC,3)
          EPSX = DDX(I)+((DBLE(RVEC(1))-0.5d0)*STEP)/AL
          EPSY = DDY(I)+((DBLE(RVEC(2))-0.5d0)*STEP)/BL
          EPSZ = DDZ(I)+((DBLE(RVEC(3))-0.5d0)*STEP)/CL
#else
          CALL RANMAR(RVEC,3)
          EPSX = DDX(I)+((DBLE(RVEC(1))-0.5d0)*STEP)
          EPSY = DDY(I)+((DBLE(RVEC(2))-0.5d0)*STEP)
          EPSZ = DDZ(I)+((DBLE(RVEC(3))-0.5d0)*STEP)
!-----------------------------------------------
! Two different algorithms to generate moves in a sphere:
!          udr = STEP*STEP
!          DO WHILE (udr > (STEP*STEP)/4.0d0)
!             CALL RANMAR(RVEC,3)
!             udx = ((DBLE(RVEC(1))-0.5d0)*STEP)
!             udy = ((DBLE(RVEC(2))-0.5d0)*STEP)
!             udz = ((DBLE(RVEC(3))-0.5d0)*STEP)
!             udr = udx*udx+udy*udy+udz*udz
!          ENDDO
!          EPSX = DDX(I)+udx
!          EPSY = DDY(I)+udy
!          EPSZ = DDZ(I)+udz
!          IF (ILOOP==1) WRITE(69,*) 'C ',udx,udy,udz
!-----------------------------------------------
!          CALL RANMAR(RVEC,3)
!          costheta = -1.0d0 + 2.0d0*DBLE(RVEC(1))
!          sintheta = SIN(ACOS(costheta))
!          phi =   2*pi*DBLE(RVEC(2))
!          rmagn = STEP*(DBLE(RVEC(3))**(1.0d0/3.0d0))
!          EPSX = DDX(I)+rmagn*sintheta*COS(phi)
!          EPSY = DDY(I)+rmagn*sintheta*SIN(phi)
!          EPSZ = DDZ(I)+rmagn*costheta
!          IF (ILOOP==1) WRITE(69,*) 'C ',EPSX-DDX(I),EPSY-DDY(I),EPSZ-DDZ(I)
#endif
#elif MoveType == 3
!Anj: Gaussian distribution:
          EPSX = GaussPos(STEP)/AL
          EPSY = GaussPos(STEP)/BL
          EPSZ = GaussPos(STEP)/CL
          Rold = DDX(I)*DDX(I)*AL*AL + DDY(I)*DDY(I)*BL*BL + DDZ(I)*DDZ(I)*CL*CL
          Rnew = EPSX*EPSX*AL*AL + EPSY*EPSY*BL*BL + EPSZ*EPSZ*CL*CL
          GPratio = -0.5d0*(Rold - Rnew)/(STEP*STEP)
#elif MoveType == 4
! Anj: Random Walk with a cuboid cut-off:
          CALL RANMAR(RVEC,3)
          EPSX = DDX(I)*AL+((DBLE(RVEC(1))-0.5d0)*STEP)
          EPSY = DDY(I)*BL+((DBLE(RVEC(2))-0.5d0)*STEP)
          EPSZ = DDZ(I)*CL-((DBLE(RVEC(3))-0.5d0)*STEP)
     IF (ABS(EPSX)<halfcutoff .AND. ABS(EPSY)<halfcutoff .AND. ABS(EPSZ)<halfcutoff) THEN
             EPSX = EPSX/AL
             EPSY = EPSY/BL
             EPSZ = EPSZ/CL
          ELSE
             EPSX = DDX(I)
             EPSY = DDY(I)
             EPSZ = DDZ(I)
             NAC = NAC -1
          END IF
#elif MoveType == 5
! Anj: Random Walk with a spherical cut-off:
          CALL RANMAR(RVEC,3)
          EPSX = DDX(I)+((DBLE(RVEC(1))-0.5d0)*STEP)/AL
          EPSY = DDY(I)+((DBLE(RVEC(2))-0.5d0)*STEP)/BL
          EPSZ = DDZ(I)+((DBLE(RVEC(3))-0.5d0)*STEP)/CL
          Rnew = EPSX*EPSX*AL*AL + EPSY*EPSY*BL*BL + EPSZ*EPSZ*CL*CL
          IF (Rnew>cutoff2) THEN
             EPSX = DDX(I)
             EPSY = DDY(I)
             EPSZ = DDZ(I)
             NAC = NAC -1
          END IF
#endif

#if LockedPos == _true_
! Anj: If locking the particles positions, then undo the effect of the MoveType code:
          EPSX = DDX(I)
          EPSY = DDY(I)
          EPSZ = DDZ(I)
#endif

!  Begin loop over neighbours: Checking for overlaps:
#if NNsCompare == _true_
          Old2NNs = 0
          New2NNs = 0
#endif
#if ClosestDist == _true_
          TryC = .FALSE.
#endif
          overlapped = .FALSE.
          JCOUNT = 1
          J = NLIST(I,JCOUNT,ILAT)
          DO WHILE (J.NE.NoMoreNeighbours .AND. .NOT. overlapped)
#if NNsCompare == _false_
             CALL Interaction(I,J,JCOUNT,ILAT,.TRUE.,IJover)
             IF (IJover.EQ.1) overlapped = .TRUE.
#else
             CALL Interaction(I,J,JCOUNT,ILAT,.FALSE.,OldIJover)
             IF (OldIJover == 1 .AND. NLInfo(I,JCOUNT,ILAT)==2) THEN
                Old2NNs = Old2NNs + 1
             END IF

             CALL Interaction(I,J,JCOUNT,ILAT,.TRUE.,IJover)
             IF (IJover == 1 .AND. NLInfo(I,JCOUNT,ILAT)==1) overlapped = .TRUE.
             IF (IJover == 1 .AND. NLInfo(I,JCOUNT,ILAT)==2) New2NNs = New2NNs + 1
#endif
#if ClosestDist == _true_
             IF (IJDr2 < CloseDr2) THEN
                TryC = .TRUE.
                TryCDr2 = IJDr2
                TryCI = I
                TryCJ = J
             END IF
#endif
             JCOUNT=JCOUNT + 1
             J = NLIST(I,JCOUNT,ILAT)
          ENDDO

#if NNsCompare == _true_
! Update 1NN and 2NNcounters:
          IF (overlapped == .FALSE.) Current2NNs = Current2NNs + (New2NNs-Old2NNs)
          IF (Current2NNs > 0) Count2NNs = Count2NNs + 1
          IF (Current2NNs < 0) THEN
             WRITE(*,*) 'ACK!: No. of overlaps with 2NNs is negative!'
             STOP
          ENDIF
#endif

#if TransProb == _false_
! IF there is no reason to calculate dM for invalid moves, jump to the next sphere:
          IF (overlapped) THEN
             NREJ = NREJ + 1
             CYCLE
          END IF
#endif

! Calculate number of overlaps in other basis
          DO JX = 1,maxnns
             IOVER(JX)=NOVER(I,JX)
          ENDDO
          CALL OVERLAP(I)
          OVOLD=NCOUNT
          OVNEW=NCOUNT+NOLAP
!Anj: Avoid buggering the bounds of the weights array:
          IF (OVOLD>maxM) OVOLD = maxM
          IF (OVNEW>maxM) OVNEW = maxM
!Anj: Set the sign right according to the phase:
          IF(ILAT.EQ.1) THEN
             OVOLD=-OVOLD
             OVNEW=-OVNEW
          ENDIF

! Calc accept/reject probability with preweighting function
! This expression assumes we read in the log of the weights.
#if MoveType == 3
          cond = (WHTS(IZERO+OVOLD)-WHTS(IZERO+OVNEW))+GPratio
#else
          cond = WHTS(IZERO+OVOLD)-WHTS(IZERO+OVNEW)
#endif
          IF (cond > 0.0d0) cond = 0.0d0
          cond = exp(cond)

#if LockInGateway == _true_
! Lock the simulation in the gateway region:
          IF (OVOLD.EQ.0 .AND. OVNEW.NE.0) cond = 0.0d0
#endif

!Anj: Actual accept/reject result placed into LOGICAL variable 'acceptmove':
          acceptmove = .FALSE.
          CALL RANMAR(RVEC,1)
          IF (RVEC(1).LT.cond) acceptmove=.TRUE.

#if TransProb == _true_
! Call the code which measures the transition probabilities using a barrier:
          CALL TPAccControl(OVOLD,OVNEW,IMOVE,acceptmove)
#endif

! IF move is good, and the weights allow, accept it:
          IF(.NOT.overlapped .AND. acceptmove) THEN
! Fix the center of mass; Calc the change in the displacement:
#if StayInCoMFrame == _true_
             delx = EPSX - DDX(I)
             dely = EPSY - DDY(I)
             delz = EPSZ - DDZ(I)
#endif
!  Update displacements etcetera:
             NCOUNT = NCOUNT+NOLAP
             NAC = NAC + 1
             IF (OVOLD.NE.OVNEW) NACT=NACT+1
             DDX(I) = EPSX
             DDY(I) = EPSY
             DDZ(I) = EPSZ
! Fix the center of mass:
#if StayInCoMFrame == _true_
             CALL FixCoMCompletely(delx,dely,delz)
#endif
! Update overlap matrix:
             DO JX = 1,maxnns
                NJX = NLIST(I,JX,ICLAT)
                NOVER(I,JX)  =  IOVER(JX)
                NOVER(NJX,indl(I,jx,iclat))  =  IOVER(JX)
             ENDDO
#if ClosestDist == _true_
! Keep tabs of the closest pair:
          IF ( I == CloseI .OR. I == CloseJ ) THEN
             Call FindClosestPair
          ELSE IF (TryC) THEN
             CloseDr2 = TryCDr2
             CloseI = TryCI
             CloseJ = TryCJ
          END IF
#endif
          ELSE
! Increase #rejected counter:
             NREJ=NREJ+1
          ENDIF

! End of atom loop...
       END DO

! End of total moves loop...
     END DO


! Simulation finished.
! Perform post-run consistency check:
 999 CALL ConsistencyCheck
       
! Save final displacements and overlap matrix
     CALL SaveConfig

! Output move statistics:
#if LockedPos == _false_
     IF (nmoves.GT.0) THEN
        WRITE(6,*) 'Move acceptance rate= ',NAC/DBLE(NMOVES*NM)
        WRITE(6,*) 'Move to new M acceptance rate = ',NACT/DBLE(NMOVES*NM)
     END IF
     WRITE(6,*) 'Total number of basis transformations= ',NFLIP
#endif

! Anj: Volume move acceptance:
#if ConstPres == _true_
     IF (VolTried.GT.0) THEN
 WRITE(*,*) 'Volumes moves leading to no olaps [%] =',100.0*REAL(VolAcc1)/REAL(VolTried)
 WRITE(*,*) 'Total vol moves acc. inc. weights (%) =',100.0*REAL(VolAcc2)/REAL(VolTried)
     END IF
#endif

! Anj: 1st v 2nd NN calculations
#if NNsCompare == _true_
     WRITE(*,*) '--NNs Comparison---------------------------------------------'
     WRITE(*,*) 'Total number of attempted moves = ',NM,'x',NMOVES
     WRITE(*,*) ' - Number of which broke 2NN constraint = ',Count2NNs
     WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
     WRITE(*,*) 'Total number of sweeps = ',NMOVES
     WRITE(*,*) ' - Number of which broke 2NN constraint = ',Count2NNsMCS
     WRITE(*,*) '--NNs Comparison---------------------------------------------'
#endif

! Density checking:
     WRITE(*,*) '--Density Checks---------------------------------------------'
     density = (NM*4.0d0*pi*((diameter/2.0d0)**3)/3.0d0)/(AL*BL*CL)
     WRITE(*,*) 'absolute density = ',density
     WRITE(*,*) 'close-packed density(1) = ',(SQRT(2.0d0)*pi/6.0d0)
     WRITE(*,*) 'close-packed density(2) = ',(NM*4.0d0*pi*(0.5d0**3)/3.0d0)/(AL*BL*CL)
     WRITE(*,*) 'rel. density(1) = ',density/(SQRT(2.0d0)*pi/6.0d0)
     WRITE(*,*) 'rel. density(2) = ',NM*diameter**3/(SQRT(2.0d0)*AL*BL*CL)
     WRITE(*,*) '--Density Checks---------------------------------------------'

  END PROGRAM fcchcp


! -----------------------------------------------------------------------------
! The simulation subroutines:
! -----------------------------------------------------------------------------

! Initialises the particle radii list for polydisperse simulation:
! ----------------------------------------------------------------
  SUBROUTINE Initialise_Polydisperse_Radii
    Use SystemDef
    IMPLICIT NONE
    INTEGER i,jc,pd_M,n_tries
    LOGICAL try_new_set
    REAL RVEC(NM)
! REAL*8 GaussPos
! Gaussian distriubution: radii(i) = (GaussPos(polyd*DIAMETER) + DIAMETER)/2.0d0

! Choose maximum polydispersity if polyd=0.0:
    IF (polyd==-1.0d0) polyd = (1.0d0 - DIAMETER)/DIAMETER
    
! Choose a new set of radii at random,
    CALL RANMAR(RVEC,NM)
    DO i = 1,NM
       radii(i) = (DBLE(RVEC(i)-0.5)*polyd*DIAMETER + DIAMETER)/2.0d0
    END DO

! Re-select radii until there are no overlaps:
    CALL CountOverAtV(AL,BL,CL,ILAT,pd_M,.FALSE.)
    IF ( pd_M.NE.0 ) THEN
       WRITE(*,*) 'erk: polydispersity is too large and has caused overlaps!'
       STOP
    END IF

! Output the radii information:
    WRITE(*,*) 'H Simulating a polydisperse system, polyd = ',polyd
    OPEN(UNIT=88,FILE='polyd.radii',STATUS='UNKNOWN')
    DO i = 1,NM
       WRITE(88,*) i,radii(i),2.0d0*radii(i)
    END DO
    CLOSE(UNIT=88)

  END SUBROUTINE Initialise_Polydisperse_Radii

! Randomizes the mapping of the lattice-switch:
! Swaps NM*swaps random pairs of sphere positions in lattice 2.
! -------------------------------------------------------------
  SUBROUTINE MakeLSTotallyRandom
    Use SystemDef
    IMPLICIT NONE
    INTEGER swi,swj,swaps,nswi
    PARAMETER(swaps = 10)
    REAL*8 lx,ly,lz
    REAL rands(2)

! Loop over pairs:
    DO nswi = 1, NM*swaps
! Choose a pair at random:
       CALL RANMAR(rands,2)
       swi = NM*rands(1) + 1
       swj = NM*rands(2) + 1
! Swap the data between swi and swj:
                lx = XLAT(swi,2);          ly = YLAT(swi,2);          lz = ZLAT(swi,2)
       XLAT(swi,2) = XLAT(swj,2); YLAT(swi,2) = YLAT(swj,2); ZLAT(swi,2) = ZLAT(swj,2)
       XLAT(swj,2) = lx;          YLAT(swj,2) = ly;          ZLAT(swj,2) = lz
! On to the next one...
    END DO

! Warn user that we are using the random mapping:
    WRITE(*,*) 'H Using random LS mapping, which should mess up the PlaneSep calc.'

! Write out the lattices:
    CALL OutputLattices

  END SUBROUTINE MakeLSTotallyRandom

! Randomizes the mapping of the lattice-switch:
! Swaps NM*swaps random pairs of planes in lattice 2.
! -------------------------------------------------------------
  SUBROUTINE MakeLSPlanesRandom
    Use SystemDef
    IMPLICIT NONE
    INTEGER pswi,pswj,swi,swj,swaps,nswi,insw
    PARAMETER(swaps = 5)
    REAL*8 lx,ly,lz
    REAL rands(2)

! Loop over pairs:
!    DO nswi = 1, NM*swaps
    DO nswi = 1,INT(NM**(1.0/3.0)/2.0)
! Choose a pair of planes at random:
!       CALL RANMAR(rands,2)
!       pswi = (NM**(1.0/3.0))*rands(1) + 1
!       pswj = (NM**(1.0/3.0))*rands(2) + 1
       pswi = nswi; pswj = INT(NM**(1.0/3.0)) - pswi + 1
       WRITE(*,*) pswi,pswj
! Swap the planes over:
       DO insw = 1,(NM**(2.0/3.0))
          swi = (pswi-1.0)*(NM**(2.0/3.0)) + insw
          swj = (pswj-1.0)*(NM**(2.0/3.0)) + insw
                lx = XLAT(swi,2);          ly = YLAT(swi,2);          lz = ZLAT(swi,2)
       XLAT(swi,2) = XLAT(swj,2); YLAT(swi,2) = YLAT(swj,2); ZLAT(swi,2) = ZLAT(swj,2)
       XLAT(swj,2) = lx;          YLAT(swj,2) = ly;          ZLAT(swj,2) = lz

       END DO
! On to the next one...
    END DO

! Warn user that we are using the random mapping:
    WRITE(*,*) 'H Using random LS plane mapping, which should mess up the PlaneSep calc.'

! Write out the lattices:
    CALL OutputLattices

  END SUBROUTINE MakeLSPlanesRandom

! Counts the number of 1stNNs that change during the switch:
! ----------------------------------------------------------
  SUBROUTINE CountNNChanges
    Use SystemDef
    IMPLICIT NONE
    INTEGER i,jc,j,dirn_ok,total_nn, total_ok
    REAL*8 dx(2),dy(2),dz(2),prec
    PARAMETER (prec = 1e-8)

! Ensure the NN lists are up to date:
    CALL MAKELIST

    total_nn = 0; total_ok = 0
    DO i = 1,NM
       DO jc = 1,18
          j = NLIST(i,jc,ILAT)
          IF (j>i .AND. NLInfo(i,jc,ILAT)==1) THEN
             dx(ILAT)  = XLAT(j,ILAT)   - XLAT(i,ILAT)
             dx(ICLAT) = XLAT(j,ICLAT) - XLAT(i,ICLAT)
             dy(ILAT)  = YLAT(j,ILAT)   - YLAT(i,ILAT)
             dy(ICLAT) = YLAT(j,ICLAT) - YLAT(i,ICLAT)
             dz(ILAT)  = ZLAT(j,ILAT)  - ZLAT(i,ILAT)
             dz(ICLAT) = ZLAT(j,ICLAT) - ZLAT(i,ICLAT)

             dx(ILAT)   = dx(ILAT)  - AINT(dx(ILAT))
             dx(ICLAT)  = dx(ICLAT) - AINT(dx(ICLAT))
             dy(ILAT)   = dy(ILAT)  - AINT(dy(ILAT))
             dy(ICLAT)  = dy(ICLAT) - AINT(dy(ICLAT))
             dz(ILAT)   = dz(ILAT)  - AINT(dz(ILAT))
             dz(ICLAT)  = dz(ICLAT) - AINT(dz(ICLAT))
             dx(ILAT)   = dx(ILAT)  - AINT(2.0*dx(ILAT))
             dx(ICLAT)  = dx(ICLAT) - AINT(2.0*dx(ICLAT))
             dy(ILAT)   = dy(ILAT)  - AINT(2.0*dy(ILAT))
             dy(ICLAT)  = dy(ICLAT) - AINT(2.0*dy(ICLAT))
             dz(ILAT)   = dz(ILAT)  - AINT(2.0*dz(ILAT))
             dz(ICLAT)  = dz(ICLAT) - AINT(2.0*dz(ICLAT))

             dirn_ok = 0
             IF (INT(dx(ILAT)/prec) == INT(dx(ICLAT)/prec)) dirn_ok = dirn_ok + 1
             IF (INT(dy(ILAT)/prec) == INT(dy(ICLAT)/prec)) dirn_ok = dirn_ok + 1
             IF (INT(dz(ILAT)/prec) == INT(dz(ICLAT)/prec)) dirn_ok = dirn_ok + 1
             total_nn = total_nn + 1
             IF (dirn_ok==3) total_ok = total_ok + 1
          END IF
       END DO
    END DO

    WRITE(*,*) 'H OK NNs (tot,ok,not ok) ',total_nn,total_ok,total_nn - total_ok

  END SUBROUTINE CountNNChanges

! This should fix the center of mass.
! Built to be applied every time we move a sphere.
! ------------------------------------------------
  SUBROUTINE FixComCompletely(deltax,deltay,deltaz)
    Use SystemDef
    IMPLICIT NONE
    INTEGER js
    REAL*8 deltax,deltay,deltaz,dxx,dyy,dzz

    dxx = deltax/DBLE(NM)
    dyy = deltay/DBLE(NM)
    dzz = deltaz/DBLE(NM)

    DO js = 1,NM
       DDX(js) = DDX(js) - dxx
       DDY(js) = DDY(js) - dyy
       DDZ(js) = DDZ(js) - dzz
    END DO

  END SUBROUTINE FixComCompletely

! Code to implment the MC stepping of the volume changing constant pressure
! simulation. Aspect ratio is free to change.
! -------------------------------------------------------------------------
  SUBROUTINE NPTChangeVolumeFree
    USE SystemDef
    IMPLICIT NONE
    REAL RVEC(4)
    REAL*8 Anew,Bnew,Cnew,Vold,Vnew
    INTEGER dirn,Mold,Mnew,Iolaps,Wold,Wnew,I,J,Jc

! Initialise random numbers, record initial volume/MCMC weight:
    CALL RANMAR(RVEC,4)
    Vold = AL*BL*CL
    Wold = IZERO - (2*ILAT-3)*NCOUNT

! Random walk the volume, allowing the cell shape to vary:
    Anew = AL; Bnew = BL; Cnew = CL
    Anew = Anew + Volstep*(RVEC(1)-0.5d0)
    Bnew = Bnew + Volstep*(RVEC(2)-0.5d0)
    Cnew = Cnew + Volstep*(RVEC(3)-0.5d0)
    Vnew = Anew*Bnew*Cnew

! Definite expansion:
    IF (Anew > AL .AND. Bnew > BL .AND. Cnew > Cl) THEN
       Iolaps = 0
! Possible contraction:
    ELSE
       CALL CountOverAtV(Anew,Bnew,Cnew,ILAT,Iolaps,.FALSE.)
    END IF

    VolTried = VolTried + 1
! If no overlaps have been created:
    IF (Iolaps.EQ.0) THEN
       VolAcc1 = VolAcc1 + 1

! Generate acceptance condition:
       CALL CountOverAtV(Anew,Bnew,Cnew,ICLAT,Mnew,.FALSE.)
       Wnew = IZERO - (2*ILAT-3)*Mnew
       cond = WHTS(Wold)-WHTS(Wnew) - Pres*(Vnew-Vold) + REAL(NM)*LOG(Vnew/Vold)
       IF (cond > 0.0d0) cond = 0.0d0

! Accept?
       IF (RVEC(4).LT.EXP(cond)) THEN
          VolAcc2 = VolAcc2 + 1
          AL = Anew; BL = Bnew; CL = Cnew
          NCOUNT = Mnew
! Update the local overlap array:
          CALL CountOverAtV(AL,BL,CL,ICLAT,Mnew,.TRUE.)
       END IF

    END IF

    RETURN
  END SUBROUTINE NPTChangeVolumeFree


! Code to implment the MC stepping of the volume changing constant pressure
! simulation. Aspect ratio is FIXED.
! -------------------------------------------------------------------------
  SUBROUTINE NPTChangeVolumeFixed
    USE SystemDef
    IMPLICIT NONE
    REAL RVEC(4)
    REAL*8 Anew,Bnew,Cnew,Vold,Vnew,Vchange
    INTEGER dirn,Mold,Mnew,Iolaps,Wold,Wnew,I,J,Jc

! Initialise random numbers, record initial volume/MCMC weight:
    CALL RANMAR(RVEC,2)
    Vold = AL*BL*CL
    Wold = IZERO - (2*ILAT-3)*NCOUNT

! Random walk the volume, allowing the cell shape to vary:
    Anew = AL*(1.0d0+Volstep*(RVEC(1)-0.5d0))
    Bnew = BL*(1.0d0+Volstep*(RVEC(1)-0.5d0))
    Cnew = CL*(1.0d0+Volstep*(RVEC(1)-0.5d0))
    Vnew = Anew*Bnew*Cnew

! Expansion:
    IF (Vnew > Vold) THEN
       Iolaps = 0
! Contraction:
    ELSE
       CALL CountOverAtV(Anew,Bnew,Cnew,ILAT,Iolaps,.FALSE.)
    END IF

    VolTried = VolTried + 1
! If no overlaps have been created:
    IF (Iolaps.EQ.0) THEN
       VolAcc1 = VolAcc1 + 1

! Generate acceptance condition:
       CALL CountOverAtV(Anew,Bnew,Cnew,ICLAT,Mnew,.FALSE.)
       Wnew = IZERO - (2*ILAT-3)*Mnew
       cond = WHTS(Wold)-WHTS(Wnew) - Pres*(Vnew-Vold) + REAL(NM)*LOG(Vnew/Vold)

! Accept?
       IF (RVEC(2).LT.EXP(cond)) THEN
          VolAcc2 = VolAcc2 + 1
          AL = Anew; BL = Bnew; CL = Cnew
          NCOUNT = Mnew
! Update the local overlap array:
          CALL CountOverAtV(AL,BL,CL,ICLAT,Mnew,.TRUE.)
       END IF

    END IF

    RETURN
  END SUBROUTINE NPTChangeVolumeFixed


! Code to count the total number of overlap pairs for the current set of 
! displacements and the given trial volume:
! ----------------------------------------------------------------------
  SUBROUTINE CountOverAtV(An,Bn,Cn,latt,olaps,UpdateNOVER)
    USE SystemDef
    IMPLICIT NONE
    REAL*8 An,Bn,Cn,Ao,Bo,Co
    INTEGER latt,olaps,vI,vJ,vJc,vIJo
    LOGICAL UpdateNOVER

! Change to trial volume, storing the current state:
    Ao = AL; Bo = BL; Co = CL
    AL = An; BL = Bn; CL = Cn

! Loop over all particles and all nearest neighbours to find overlaps:
    olaps = 0
    DO vI=1,NM
       vJc = 1
       vJ = NLIST(vI,vJc,latt)
       DO WHILE (vJ.NE.NoMoreNeighbours)
          CALL Interaction(vI,vJ,vJc,latt,.FALSE.,vIJo)
          olaps = olaps + vIJo
          IF (UpdateNOVER) NOVER(vI,vJc) = vIJo
          vJc=vJc + 1
          vJ = NLIST(vI,vJc,latt)
       END DO
    END DO
! If the number of overlapping spheres is not even, DIE HORRIBLY!
    IF (MOD(olaps,2).NE.0) THEN
       WRITE(*,*) 'Oh, bugger: Odd number of overlapping spheres in CountOverAtV.'
       WRITE(*,*) ' erk-',An,Bn,Cn,latt,olaps,UpdateNOVER
       STOP
    ENDIF
! We want the number of particle pairs, so divide by two:
    olaps = olaps/2

! Restore the original volume:
    AL = Ao; BL = Bo; CL = Co

    RETURN
  END SUBROUTINE CountOverAtV


! Calculates the interaction energy between two particles, in this case
! it counts the number of pairs of overlaps, ie it returns 0 or 1 in olap.
! It also places the actual seperation^2 in SystemDef:IJDr2
! ------------------------------------------------------------------------
  SUBROUTINE Interaction(I,J,Jc,latt,trial,olap)
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,Jc,olap,latt
    LOGICAL trial
    REAL*8 DX,DY,DZ,DX2,DY2,DZ2,IX,IY,IZ
    REAL*8 Uij(3),nij(3),nijmag,test_sep
    
! Deal with trial moves or the cuurent positions:
    IF (trial) THEN 
       IX = EPSX;  IY = EPSY;  IZ = EPSZ
    ELSE
       IX = DDX(I); IY = DDY(I); IZ = DDZ(I)
    ENDIF

#if UseDodeca == _true_
! Create displacement difference vector:
    Uij(1) = (DDX(J) - IX)
    Uij(2) = (DDY(J) - IY)
    Uij(3) = (DDZ(J) - IZ)

    Upara = Uij(1)*lnnx(I,Jc,latt) + Uij(2)*lnny(I,JC,latt) + Uij(3)*lnnz(I,Jc,latt)
! This is the parallel component of U, the perpendicular bit is:
#if CalcDodecVecH == _true_
! Uperp**2 = U**2 - Upara**2
    Uperp2 = (Uij(1)*Uij(1)+Uij(2)*Uij(2)+Uij(3)*Uij(3)) - Upara*Upara
#endif

    IF(Upara > -1.0d0) THEN
       olap = 0
    ELSE
       olap = 1
    END IF
    IJDr2 = Upara

#else
! Calculate position difference vector:
    DX =  XLAT(I,latt) + IX - XLAT(J,latt) - DDX(J)
    DX = DX - AINT(DX)
    DX = DX - AINT(DX+DX)

    DY =  YLAT(I,latt) + IY - YLAT(J,latt) - DDY(J)
    DY = DY - AINT(DY)
    DY = DY - AINT(DY+DY)

    DZ =  ZLAT(I,latt) + IZ - ZLAT(J,latt) - DDZ(J)
    DZ = DZ - AINT(DZ)
    DZ = DZ - AINT(DZ+DZ)

    DX2 = DX*DX*AL*AL
    DY2 = DY*DY*BL*BL
    DZ2 = DZ*DZ*CL*CL
    IJDr2 = DX2+DY2+DZ2

#if PolyDisperse == _true_
    test_sep = (radii(I)+radii(J))
    test_sep = test_sep*test_sep
!    WRITE(*,*) 'H Simulating a polydisperse system, C',test_sep
#else
    test_sep = DIAMETER2
#endif

! Return overlap cost of this configuration:
    IF(IJDr2.LT.test_sep) THEN
       olap = 1
    ELSE
       olap = 0
    END IF
#endif
    RETURN
  END SUBROUTINE Interaction


! Do a consistency check on the number of conjugate overlaps
! ----------------------------------------------------------
  SUBROUTINE MeltCheck(meltdown,time)
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,JCOUNT,meltnum,time
    REAL*8 ilix,iliy,iliz,ilisep2,iljsep2,iljx,iljy,iljz
    LOGICAL melted,meltdown,KillTheCode
    PARAMETER ( KillTheCode = .FALSE.)

    meltnum = 0

    DO I = 1,NM
       melted = .FALSE.
       ilix = DDX(I)
       iliy = DDY(I)
       iliz = DDZ(I)
       ilisep2 = ilix*ilix*AL*AL + iliy*iliy*BL*BL + iliz*iliz*CL*CL
       JCOUNT = 1
       J = NLIST(I,JCOUNT,ILAT)
       DO WHILE (J.NE.NoMoreNeighbours)

          iljx = XLAT(I,ILAT) + DDX(I) - XLAT(J,ILAT)
          iljx = iljx - AINT(iljx)
          iljx = iljx - AINT(iljx+iljx)

          iljy = YLAT(I,ILAT) + DDY(I) - YLAT(J,ILAT)
          iljy = iljy - AINT(iljy)
          iljy = iljy - AINT(iljy+iljy)

          iljz = ZLAT(I,ILAT) + DDZ(I) - ZLAT(J,ILAT)
          iljz = iljz - AINT(iljz)
          iljz = iljz - AINT(iljz+iljz)

          iljsep2 = iljx*iljx*AL*AL + iljy*iljy*BL*BL + iljz*iljz*CL*CL
          IF (iljsep2.LT.ilisep2) THEN
             IF (.NOT.melted) THEN
                meltnum=meltnum+1
                melted=.TRUE.
             END IF
           WRITE(62,*) 't=',time,': particle',I,' is nearer the site',J,' than it`s own!'
          END IF

       JCOUNT = JCOUNT + 1
       J = NLIST(I,JCOUNT,ILAT)
       END DO
    END DO

    IF (meltnum.GT.0 .AND. KillTheCode) THEN
       WRITE(*,*)
      WRITE(*,*) 'Fatal error in MeltCheck: ',meltnum,' particles have left their sites.'
       WRITE(*,*) '  i.e. Sample appears to have melted.'
       WRITE(*,*)
       meltdown = .TRUE.
    ELSE
       meltdown = .FALSE.
    END IF
     
  END SUBROUTINE MeltCheck

! Do a consistency check on the number of conjugate overlaps
! ----------------------------------------------------------
  SUBROUTINE ConsistencyCheck
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,JCOUNT,NCHK
    
    WRITE(*,*) 'Doing consistency checks'
    CALL CountOverAtV(AL,BL,CL,ICLAT,NCHK,.FALSE.)
    IF(NCHK.NE.NCOUNT) THEN
       WRITE(*,*) 'CONSISTENCY ERROR IN NCOUNT',NCHK,NCOUNT
       STOP
    ELSE
       WRITE(*,*) 'OK'
    END IF
     
  END SUBROUTINE ConsistencyCheck


! This subroutine counts the change in the number of overlaps if we switch basis.
! SWITCHING OVER TO N0LIST DOESN'T WORK ANYMORE.
! -------------------------------------------------------------------------------
  SUBROUTINE OVERLAP(I)
    USE SystemDef
    IMPLICIT NONE
    INTEGER J,I,JCOUNT,IJover

    NOLAP = 0
    JCOUNT = 1
    J = NLIST(I,JCOUNT,ICLAT)
!   J = N0LIST(I,JCOUNT,ICLAT)
    DO WHILE (J.NE.NoMoreNeighbours)
       IF(IOVER(JCOUNT)==true) THEN
          IOVER(JCOUNT) = false
          NOLAP = NOLAP-1
       ENDIF
       CALL Interaction(I,J,JCOUNT,ICLAT,.TRUE.,IJover)
       IF (IJover.EQ.1) THEN
          NOLAP = NOLAP +1
          IOVER(JCOUNT) = true
       ENDIF
    JCOUNT=JCOUNT + 1
    J = NLIST(I,JCOUNT,ICLAT)
!   J = N0LIST(I,JCOUNT,ICLAT)
    ENDDO

    RETURN
  END SUBROUTINE OVERLAP


! Load in the old configuration:
! ------------------------------
  SUBROUTINE LoadConfig
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,IT1,IT2,ILEN

    OPEN (unit=25,file='conf_in',FORM='FORMATTED',STATUS='UNKNOWN')

    READ(25,*) NCOUNT,ILAT,iseed
    ICLAT = 3 - ILAT

    DO I=1,NM
#if PolyDisperse == _true_
       READ(25,*)DDX(I),DDY(I),DDZ(I),radii(I)
#else
       READ(25,*)DDX(I),DDY(I),DDZ(I)
#endif
    ENDDO
    READ(25,*) ILEN
    DO I=1,ILEN
       READ(25,*) IT1,IT2
!Anj: Bug corrention to make the read-in data match the output data:
       DO J = 1,maxnns
          IF (NLIST(IT1,J,ICLAT).EQ.IT2) THEN
             NOVER(IT1,J) = true
          ENDIF
       ENDDO
    ENDDO

#if ConstPres == _true_
    READ(25,*) AL,BL,CL
#endif
#if NNsCompare == _true_
    READ(25,*) Current2NNs
#endif
    CLOSE(UNIT=25)

! Perform consistency check on the loaded configuration:
    CALL ConsistencyCheck

    RETURN
  END SUBROUTINE LoadConfig


! Save the new configuration:
! ---------------------------
  SUBROUTINE SaveConfig
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,ILEN,K
!Anj: New RNG seed generator function:
    INTEGER RMnewseed

    OPEN (unit=26,file='conf_out',FORM='FORMATTED',STATUS='UNKNOWN')

!Anj: header: no. overlaps in other phase, current lattice, RNG seed:
    iseed = RMnewseed()
    WRITE(26,*) NCOUNT,ILAT,iseed

!Anj: positions (and radii if polydisperse):
    DO I =1,NM
#if PolyDisperse == _true_
       WRITE(26,*)DDX(I),DDY(I),DDZ(I),radii(I)
#else
       WRITE(26,*)DDX(I),DDY(I),DDZ(I)
#endif
    ENDDO

!Anj: Overlap matrix element !=0:
    ILEN=0
    DO I=1,NM
       DO J=1,maxnns
          IF(NOVER(I,J)==true) ILEN=ILEN+1
       ENDDO
    ENDDO
    WRITE(26,*)ILEN
    DO I=1,NM
       DO J=1,maxnns
          IF(NOVER(I,J)==true) WRITE(26,*) I,NLIST(I,J,ICLAT)
       ENDDO
    ENDDO

!Anj: Extras:
#if ConstPres == _true_
    WRITE(26,*) AL,BL,CL
#endif
#if NNsCompare == _true_
    WRITE(26,*) Current2NNs
#endif
    CLOSE(UNIT=26)
       
! Write out in a format suitable for visualisation.
#if OutputXYZConf == _true_
       OPEN (unit=27,file='conf.xyz',FORM='FORMATTED',STATUS='UNKNOWN')
       WRITE(27,*)NM+8
       WRITE(27,*)
       DO I =1,NM
#if UseDodeca == _true_
         WRITE(27,*) 'C ',XLAT(I,ILAT)*AL+DDX(i)*(1.0d0-DIAMETER),YLAT(I,ILAT)*BL+DDY(i)*(1.0d0-DIAMETER),ZLAT(I,ILAT)*CL+DDZ(I)*(1.0d0-DIAMETER)
#else
         WRITE(27,*) 'C ',AL*(XLAT(I,ILAT)+DDX(i)),BL*(YLAT(I,ILAT)+DDY(i)),CL*(ZLAT(I,ILAT)+DDZ(I))
#endif
       ENDDO
! Write out the corners of the periodic cell:
       DO I = -1,1,2
          DO J = -1,1,2
             DO K = -1,1,2
                WRITE(27,*) 'H ',0.5d0*I*AL,0.5d0*J*BL,0.5d0*K*CL
             END DO
          END DO
       END DO
       CLOSE(UNIT=27)
#endif

    RETURN
  END SUBROUTINE SaveConfig


! MakeLattices: (Anj,17/Mar/99)
!    Construct the lattices:    X/Y/ZLAT(N,1/2)
!    Input is NM, 
!        which must corresponds to a cube with a linear dimension divisible by 6.
!    The fcc (2) is created from the hcp (1) by shifting pairs of planes.
!    These lattices are then scaled and shifted onto a unit cube centred about the origin
!    The initial system dimensions are also defined here (AL,BL,CL).
!----------------------------------------------------------------------------------------
  SUBROUTINE MakeLattices
    USE SystemDef
    IMPLICIT NONE
    INTEGER ilatt,ix,iy,iz,in
    INTEGER Nx,Ny,Nz
    REAL*8 xs,ys,zs,t
    REAL*8 xo,yo,zo,xoo,yoo
    REAL*8 fccshift(2),hcpshift(2)

! Define N-dimensions from overall system size:
    Nx = NINT(DBLE(NM)**(1.0d0/3.0d0))
    Ny = Nx
    Nz = Nx
    IF (MOD(Nz,6).NE.0) THEN
       WRITE(*,*) 'MakeLattices: No. of x-y planes must be n*6 for this lattice builder!'
       STOP
    END IF

! Define seperation in the x direction of the origin of two planes A and B
    t = SQRT(3.0d0)/3.0d0
! Define distances between those atoms which share the same plane:
    xs = SQRT(3.0d0)/2.0d0
    ys = 1.0d0
    zs = SQRT(2.0d0/3.0d0)

! Initialize particle counter:
    in = 1
! Loop in z, over the x-y planes:
    DO iz = 1,Nz
! Build look up table of which planes changes (above=+1, below=-1):
       PlaneChanges(iz) = 1 - 2*MOD(iz,2)
! Define origin of hcp planes:
       IF (MOD(iz,2)==1) xo = 0.0d0        ! A
       IF (MOD(iz,2)==0) xo = +t           ! B
       yo = 0.0d0
       zo = DBLE(iz-1)*zs
! Shift these to get fcc from hcp:                                      ! hcp -> fcc
       IF (MOD(iz,6)==1 .OR. MOD(iz,6)==2) THEN
          hcpshift(1) = +0.0d0*t; hcpshift(2) = +0.0d0
          fccshift(1) = +0.0d0*t; fccshift(2) = +0.0d0
       ELSE IF (MOD(iz,6)==3 .OR. MOD(iz,6)==4) THEN
          hcpshift(1) = +0.0d0*t; hcpshift(2) = +0.0d0
          fccshift(1) = +2.0d0*t; fccshift(2) = +0.0d0
       ELSE IF (MOD(iz,6)==5 .OR. MOD(iz,6)==0) THEN
          hcpshift(1) = +0.0d0*t; hcpshift(2) = +0.0d0
          fccshift(1) = +1.0d0*t; fccshift(2) = +0.0d0
       ELSE
          WRITE(*,*) 'MakeLattices:  Ooer, plane belongs to no pair!'
       END IF
! Loop within the x-y planes:
       DO ix = 1,Nx
          xoo = xo + DBLE(ix-1)*xs
          yoo = yo + DBLE(MOD(ix-1,2))*ys/2.0d0
          DO iy=1,Ny
            XLAT(in,1) = xoo + hcpshift(1)
            YLAT(in,1) = yoo + DBLE(iy-1)*ys + hcpshift(2)
            ZLAT(in,1) = zo
            XLAT(in,2) = xoo + fccshift(1)
            YLAT(in,2) = yoo + DBLE(iy-1)*ys + fccshift(2)
            ZLAT(in,2) = zo
            WhichPlane(in) = iz
            WhichXRow(in) = ix
            WhichYRow(in) = 2*(iy-1) + MOD(ix-1,2) + 1
!            WRITE(*,*) in,WhichPlane(in),WhichXRow(in),WhichYRow(in)
            in = in + 1
          END DO
       END DO
    END DO

! Calculate the lattice size:
    AL0 = DBLE(Nx)*xs
    BL0 = DBLE(Ny)*ys
    CL0 = DBLE(Nz)*zs
    AL = AL0; BL = BL0; CL = CL0

! Map lattice into the simulation cube, by scaling and by periodic boundary conds:
    DO ilatt = 1,2
       DO in = 1,Nx*Ny*Nz
          XLAT(in,ilatt) = (XLAT(in,ilatt)/AL0)-0.5d0
          XLAT(in,ilatt) = XLAT(in,ilatt)-AINT(XLAT(in,ilatt))
          XLAT(in,ilatt) = XLAT(in,ilatt)-AINT(XLAT(in,ilatt)+XLAT(in,ilatt))

          YLAT(in,ilatt) = (YLAT(in,ilatt)/BL0)-0.5d0
          YLAT(in,ilatt) = YLAT(in,ilatt)-AINT(YLAT(in,ilatt))
          YLAT(in,ilatt) = YLAT(in,ilatt)-AINT(YLAT(in,ilatt)+YLAT(in,ilatt))

          ZLAT(in,ilatt) = (ZLAT(in,ilatt)/CL0)-0.5d0
          ZLAT(in,ilatt) = ZLAT(in,ilatt)-AINT(ZLAT(in,ilatt))
          ZLAT(in,ilatt) = ZLAT(in,ilatt)-AINT(ZLAT(in,ilatt)+ZLAT(in,ilatt))
       END DO
    END DO

! Write out the lattices:
    CALL OutputLattices

    RETURN
  END SUBROUTINE MakeLattices

!   Dump the lattices to the files fcclatt and hcplatt
!-----------------------------------------------------
  SUBROUTINE OutputLattices
    USE SystemDef
    IMPLICIT NONE
    INTEGER ilatt,in

! Open output file to store the lattices in:
    OPEN (unit=31,file='hcplatt',FORM='FORMATTED',STATUS='UNKNOWN')
    OPEN (unit=32,file='fcclatt',FORM='FORMATTED',STATUS='UNKNOWN')

! Map lattice into the simulation cube, by scaling and by periodic boundary conds:
    DO ilatt = 1,2
       DO in = 1,NM
          WRITE(30+ilatt,*) XLAT(in,ilatt),YLAT(in,ilatt),ZLAT(in,ilatt)
       END DO
    END DO

! Close lattice files:
    CLOSE(UNIT=31)
    CLOSE(UNIT=32)

  END SUBROUTINE OutputLattices

! This subroutine set up the (constant!) neighbour list for all atoms.
! 15th/March/1999: Altered to store 1st/2nd NN information (in NLInfo).
! ---------------------------------------------------------------------
  SUBROUTINE MAKELIST
    USE SystemDef
    IMPLICIT NONE
    INTEGER I,J,K,NF,NL,NL0,NH,inlatt,inclatt,nnorder
    REAL*8 DX0,DX,DY,DZ,DX2,DY2,DZ2,DSEP2,dsep

! Open the files to store the NNlists in:
    OPEN (unit=50,file='N0list',FORM='FORMATTED',STATUS='UNKNOWN')
    OPEN (unit=51,file='NHlist',FORM='FORMATTED',STATUS='UNKNOWN')
    OPEN (unit=52,file='NFlist',FORM='FORMATTED',STATUS='UNKNOWN')

! Initialising the NN lists:
    DO I=1,NM
       DO J=1,19
          NLIST(I,J,1)=NoMoreNeighbours
          NLIST(I,J,2)=NoMoreNeighbours
       ENDDO
       DO J=1,7
          N0LIST(I,J,1)=NoMoreNeighbours
          N0LIST(I,J,2)=NoMoreNeighbours
       ENDDO
    ENDDO
    maxnns = 0

!Anj: Loop over lattices:
    DO inlatt = 1,2
       inclatt = 3 - inlatt
! Loop over atoms:
       DO I = 1,nm
          NL=0
          NL0=0
! Loop over atoms over atoms:
          DO J = 1,nm
! Do not neighbour an atom to itself:
             IF (I.NE.J) THEN
! Calculate seperation:
                DX =  XLAT(I,inlatt)-XLAT(J,inlatt)
                DX = DX - AINT(DX)
                DX = DX - AINT(DX+DX)
                DX2 = DX*DX*AL*AL

                DY =  YLAT(I,inlatt)-YLAT(J,inlatt)
                DY = DY - AINT(DY)
                DY = DY - AINT(DY+DY)
                DY2 = DY*DY*BL*BL

                DZ =  ZLAT(I,inlatt)-ZLAT(J,inlatt)
                DZ = DZ - AINT(DZ)
                DZ = DZ - AINT(DZ+DZ)
                DZ2 = DZ*DZ*CL*CL

                DSEP2 = DX2+DY2+DZ2

! Check that perfect lattice doesn't overlap
                IF(DSEP2.LT.0.999999d0) THEN
               IF (inlatt==1) WRITE(6,*) 'HCP neighbours too close in MAKELIST',I,J,DSEP2
               IF (inlatt==2) WRITE(6,*) 'FCC neighbours too close in MAKELIST',I,J,DSEP2
               STOP
                ENDIF

! Calculate if it's a first or second nearest neighbour:
                nnorder = 0
                IF(DSEP2.LT.1.2d0) nnorder = 1
                IF(DSEP2.GE.1.2d0 .AND. DSEP2.LT.2.1d0) nnorder = 2
! Add to neighbour lists:
#if Check1NNsOnly == _true_
                IF (nnorder==1) THEN
#else
                IF (nnorder==1 .OR. nnorder==2) THEN
#endif
                   NL = NL+1
                   NLIST(I,NL,inlatt) = J
                   NLInfo(I,NL,inlatt) = nnorder
                   WRITE(50+inlatt,*) I,NL,J,nnorder
! Store lattice neighbour normal vectors:
                   dsep = SQRT(DSEP2)
                   lnnx(I,NL,inlatt) = -DX*AL/dsep
                   lnny(I,NL,inlatt) = -DY*BL/dsep
                   lnnz(I,NL,inlatt) = -DZ*CL/dsep
!                   IF (inlatt==1) WRITE(67,FMT='(A,3F14.7)') 'C ',lnnx(I,NL,inlatt),lnny(I,NL,inlatt),lnnz(I,NL,inlatt)
!                   IF (inlatt==2) WRITE(68,FMT='(A,3F14.7)') 'C ',lnnx(I,NL,inlatt),lnny(I,NL,inlatt),lnnz(I,NL,inlatt)
! Also pick out the neighbours on planes above and below.
                   DX0 =  XLAT(I,inclatt)-XLAT(J,inclatt)
                   DX0 = DX0 - AINT(DX0+DX0)
                   IF(ABS(DX-DX0).GT.0.01) THEN 
                      NL0 = NL0+1
                      N0LIST(I,NL0,inlatt) = J
                      N0LInfo(I,NL0,inlatt) = nnorder
                      WRITE(50,*) I,NL0,J,nnorder
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
          IF ( NL > maxnns) maxnns = NL
       ENDDO
    ENDDO

    CLOSE(UNIT=50)
    CLOSE(UNIT=51)
    CLOSE(UNIT=52)

! Construct cross-referencing arrays:
    DO i=1,nm
       DO j=1,maxnns
          NF = nlist(i,j,1)
          NH = nlist(i,j,2)
          DO k=1,maxnns
             IF(nlist(nf,k,1).EQ.i) indl(i,j,1)=k   
             IF(nlist(nh,k,2).EQ.i) indl(i,j,2)=k   
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE MAKELIST


! Transition probability routines and data storage:
!----------------------------------------------------------------------------
  MODULE TPData
    IMPLICIT NONE
    INTEGER IMLimit
    INTEGER, PARAMETER :: Itpwidth = 8
    INTEGER, ALLOCATABLE :: TPM(:,:),TPMvisits(:)
    INTEGER :: Itptimer,Itpcycle,Ibarrier,oldtime
    LOGICAL, PARAMETER :: TwoSweep = .FALSE.
  END MODULE TPData

! Initialisation:
!----------------
  SUBROUTINE TPInitialise
    USE SystemDef
    USE TPData
    IMPLICIT NONE
    INTEGER Itpa,Itpb

    IF (NM == smallNM)  IMLimit = TP_sml_max
    IF (NM == mediumNM) IMLimit = TP_med_max
    IF (NM == largeNM)  IMLimit = TP_lrg_max

    WRITE(*,*) 'TPM: Simulation must be at least ',NEQBM+(IMLimit+1)*NSAMP,' MCS long.'

    ALLOCATE(TPM(0:IMLimit,-Itpwidth:Itpwidth), TPMvisits(0:IMLimit))

    oldtime=0; Itptimer = 0; Itpcycle = 0; Ibarrier = 1

    DO Itpa = 0,IMLimit
       TPMvisits(Itpa) = 0
       DO Itpb = -Itpwidth,Itpwidth
          TPM(Itpa,Itpb) = 0
       END DO
    END DO
    
    RETURN
  END SUBROUTINE TPInitialise

! Acceptance control and data logging:
!-------------------------------------
  SUBROUTINE TPAccControl(Mold,Mnew,time,accept)
    USE SystemDef
    USE TPData
    IMPLICIT NONE
    LOGICAL accept, lockmove
    INTEGER Mold,Mnew,time

! If this is already over, return:
    IF (TransProbFin) RETURN

! Default to not locking the particle:
    lockmove = .FALSE.

! If the move _would_ be accepted normally:
    IF (.NOT.overlapped .AND. accept) THEN

! If the barrier is moving to greater |M|:
       IF (Itpcycle==0) THEN
! Calc the current barrier position
          IF (time-NEQBM>oldtime) THEN
             oldtime = time-NEQBM
             Itptimer = Itptimer + 1
          END IF
          IF (Itptimer>NSAMP) THEN
             Ibarrier = Ibarrier + 1
             Itptimer = 0
! Periodically output the TPM matrix, to see how things are going:
             CALL TPMOutput
          END IF
! Decide whether to accept:
          IF (ABS(Mnew)>Ibarrier) lockmove = .TRUE.
! If we have reached the end, move to the second cycle:
          IF (Ibarrier>IMlimit) THEN
             Itptimer = 0
             Itpcycle = 1
             RETURN
          END IF

! Else, the barrier is moving to lower |M|:
       ELSE
! Calc the current barrier position
          IF (time-NEQBM>oldtime) THEN
             oldtime = time-NEQBM
             Itptimer = Itptimer + 1
          END IF
          IF (Itptimer>NSAMP .AND. ABS(Mold)<=Ibarrier-1) THEN
             Ibarrier = Ibarrier - 1
             Itptimer = 0
          END IF
! Decide whether to accept:
          IF (ABS(Mnew)>Ibarrier .AND. ABS(Mold)<=Ibarrier) lockmove = .TRUE.
! If we have reached the end, output the results:
          IF (Ibarrier==1 .OR. .NOT.TwoSweep) THEN
             CALL TPMOutput
             TransProbFin = .TRUE.
             RETURN
          END IF
       END IF
! If not going to be accepted, then record the 'no-move':
    ELSE
       Mnew = Mold
    END IF

! If equilibrated, record the move in the TP matrix:
    IF (time>NEQBM .AND. ABS(Mnew-Mold)<=Itpwidth) THEN
       IF (Mold>0) THEN
          TPM(Mold,Mnew-Mold) = TPM(Mold,Mnew-Mold) + 1
          TPMvisits(Mold) = TPMvisits(Mold) + 1
       END IF
       IF (Mold==0) THEN
          TPM(0,ABS(Mnew)) = TPM(0,ABS(Mnew)) + 1
          TPMvisits(0) = TPMvisits(0) + 1
       END IF
       IF (Mold<0) THEN
          TPM(-Mold,-(Mnew-Mold)) = TPM(-Mold,-(Mnew-Mold)) + 1
          TPMvisits(-Mold) = TPMvisits(-Mold) + 1
       END IF
    END IF

! Lock the move out, if that's what we want to do:
    IF (lockmove) accept = .FALSE.

    RETURN
  END SUBROUTINE TPAccControl

! TPM output routine:
!--------------------
  SUBROUTINE TPMOutput
    USE SystemDef
    USE TPdata
    IMPLICIT NONE
    INTEGER Itpa,Itpb

! First, output the transition probability matrix to 'TPM.dat':
    OPEN(UNIT=50, FILE='TPM.dat', STATUS='UNKNOWN')
    WRITE(50,*) 0,IMLimit,Itpwidth
    DO Itpa = 0,IMLimit
       WRITE(50,*)
       DO Itpb = -Itpwidth,Itpwidth
          IF (TPMvisits(Itpa).NE.0) THEN
             WRITE(50,*) DBLE(TPM(Itpa,Itpb))/DBLE(TPMvisits(Itpa))
          ELSE
             WRITE(50,*) 0.0d0
          ENDIF
       END DO
    END DO
    CLOSE(UNIT=50)

    RETURN
  END SUBROUTINE TPMOutput


!----------------------------------------------------------------------------
! Analysis routines and data storage:
!----------------------------------------------------------------------------
  MODULE AnalysisData
    IMPLICIT NONE
!Anj: General analysis parameters:
    INTEGER M,INACvM
!Anj: Plane definition data:
    REAL*8 Pdisplz(12)
    INTEGER dvMSet(12,144),IdvMSet,IdvMNum
  END MODULE AnalysisData

!-------------------------------------------------------------------
!  Specific analysis routines called from the general analysis code:
!-------------------------------------------------------------------

! Calcs data from u.R for NNs which change during L-S:
! -------------------------------------------------------
  SUBROUTINE InitNNDistCalc
    Use SystemDef
    IMPLICIT NONE
    INTEGER Init_Histogram,hnum

    nnd_hists(1) = Init_Histogram(-4.0d0, 4.0d0, 60)
    nnd_hists(2) = Init_Histogram(-4.0d0, 4.0d0, 60)
    nnd_hists(3) = Init_Histogram(-4.0d0, 4.0d0, 60)
    nnd_hists(4) = Init_Histogram(-0.5d0, 24.5d0, 25)
    nnd_hists(5) = Init_Histogram(-1.0d0, 25.0d0, 13)

  END SUBROUTINE InitNNDistCalc
! ------------------------------
  SUBROUTINE CalcNNDistributions
    Use SystemDef
    IMPLICIT NONE
    INTEGER i,jc,j,ill,il,plch,ipl,jpl
    INTEGER ipy,ipx
    REAL*8 udotr,sfactor
    REAL*8 correl,cdx,cdy,cdr

! Output is scaled by the appropriate factor so that units are sphere-seperation.
    sfactor = 1.0d0/(1.0d0-0.7778**(1.0d0/3.0d0))


! Loop over all atoms:
    DO i = 1,NM
! Loop over both lattices
       DO ill = 1,2
          IF (ill==1) il = ILAT
          IF (ill==2) il = ICLAT
          ipl = WhichPlane(i)
          plch = PlaneChanges(ipl)
! Loop over neighbours:
          jc = 1
          j = NLIST(i,jc,il)
          DO WHILE (j.NE.NoMoreNeighbours)
! Does this neighbour change during LS?:
             jpl = WhichPlane(j)
             IF ( ipl.NE.jpl .AND. (j-i)*plch > 0 .AND. NLInfo(i,jc,il) == 1 ) THEN
! If so, calculate u.R and spew it:
        udotr = DDX(i)*AL*lnnx(i,jc,il) + DDY(i)*BL*lnny(i,jc,il) + DDZ(i)*CL*lnnz(i,jc,il)
!                IF (udotr > 0.0d0) WRITE(47+ill,FMT='(E14.5)') udotr
!                IF (udotr > 0.0d0) CALL Add_To_Histogram(nnd_hists(ill),udotr)
                CALL Add_To_Histogram(nnd_hists(ill),udotr*sfactor)
             END IF
             jc = jc + 1
             j = NLIST(i,jc,il)
          END DO

! Also, collect histogram of u_z:
! HS:
!          CALL Add_To_Histogram(nnd_hists(3),sfactor*CL*DDZ(i)*PlaneChanges(WhichPlane(i)))
! HD:
          CALL Add_To_Histogram(nnd_hists(3),DDZ(i)*PlaneChanges(WhichPlane(i)))

       END DO

! Also, collect the in-plane autocorrelation function for planes:
       DO j = 1,NM
          IF (WhichPlane(i)==WhichPlane(j) .AND. WhichXRow(i)==WhichXRow(j)) THEN
             DO ipy = -1,1
                cdy = YLAT(j,ILAT)+DBLE(ipy) - YLAT(i,ILAT)
! HS:
!                correl = DDZ(j)*DDZ(i)*CL*CL*sfactor*sfactor
! HD:
                correl = DDZ(j)*DDZ(i)
                CALL Add_This_To_Histogram(nnd_hists(4),ABS(cdy)*BL,correl)
             END DO
          END IF
          IF (WhichPlane(i)==WhichPlane(j) .AND. WhichYRow(i)==WhichYRow(j)) THEN
             DO ipx = -1,1
                cdx = (XLAT(j,ILAT)+DBLE(ipx) - XLAT(i,ILAT))/(SQRT(3.0d0)/2.0d0)
! HS:
!                correl = DDZ(j)*DDZ(i)*CL*CL*sfactor*sfactor
! HD:
                correl = DDZ(j)*DDZ(i)
                CALL Add_This_To_Histogram(nnd_hists(5),ABS(cdx)*AL,correl)
             END DO
          END IF
       END DO
       
    END DO

    CALL Output_Histogram(nnd_hists(1),'nnd.this.dat')
    CALL Output_Histogram(nnd_hists(2),'nnd.other.dat')
    CALL Output_Histogram(nnd_hists(3),'u_z.dat')
    CALL Output_Histogram(nnd_hists(4),'u_z0.u_zY.dat')
    CALL Output_Histogram(nnd_hists(5),'u_z0.u_zX.dat')

  END SUBROUTINE CalcNNDistributions

! Routine to test what happens is we use a rotated M:
! NB Count no. overlaps between and within layers.
!
! ---------------------------------------------------
  SUBROUTINE TestRotatedM(itime)
    Use SystemDef
    IMPLICIT NONE
    INTEGER i,ipl,rM,itime,tot_shifted
    REAL*8 phi,xyzrad,lo_phi,hi_phi,shift_frac,max_rot
    PARAMETER (lo_phi = 0.0d0*3.141582654d0/180.0d0, hi_phi = 60.0*3.141592654d0/180.0d0)
    PARAMETER (max_rot = 60.0d0*3.141582654d0/180.0d0)

! Copy off the current displacement array, and rotate the displacement array:
    tot_shifted = 0
    DO i = 1,NM
! Copy:
       Dcp(i,1) = DDX(i)*AL
       Dcp(i,2) = DDY(i)*BL
       Dcp(i,3) = DDZ(i)*CL
! Rotate:
       ipl = WhichPlane(i)
! Only change displacements towards the changing plane:
       IF (SIGN(Dcp(i,3),DBLE(PlaneChanges(ipl)))==Dcp(i,3)) THEN
          xyzrad = SQRT(Dcp(i,1)*Dcp(i,1)+Dcp(i,2)*Dcp(i,2)+Dcp(i,3)*Dcp(i,3))
          IF (xyzrad.NE.0.0d0) THEN
             phi = ABS(ASIN(Dcp(i,3)/xyzrad))
          ELSE
             phi = 0.0d0
          END IF
! Determine fraction of total rotation based on phi:
          shift_frac = 0.0d0
          IF (phi > lo_phi .AND. phi < hi_phi) shift_frac = (phi-lo_phi)/(hi_phi-lo_phi)
          IF (phi >= hi_phi) shift_frac = 1.0d0

! Actually rotate the displ. in the x-y plane:
          DDX(i) = (Dcp(i,1)*COS(max_rot*shift_frac) - Dcp(i,2)*SIN(max_rot*shift_frac))/AL
          DDY(i) = (Dcp(i,1)*SIN(max_rot*shift_frac) + Dcp(i,2)*COS(max_rot*shift_frac))/BL
          DDZ(i) = Dcp(i,3)/CL
          IF (shift_frac .NE. 0.0d0) THEN
             tot_shifted = tot_shifted + 1
! Output to check the calculation:
             WRITE(*,FMT='(A,2I5,8E14.4)') 'R',i,ipl,Dcp(i,1),Dcp(i,2),Dcp(i,3),phi,max_rot*shift_frac,DDX(i)*AL,DDY(i)*BL,DDZ(i)*CL
          END IF

       END IF
    END DO

! Calculate the number of overlaps:
    CALL CountOverAtV(AL,BL,CL,ICLAT,rM,.FALSE.)

! Output both the M's
    WRITE(*,*) 't-M-rM',itime,NCOUNT*(2*ILAT-3),rM*(2*ILAT-3),tot_shifted

! Copy original displacements back into the array:
    DO i = 1,NM
       DDX(i) = Dcp(i,1)/AL
       DDY(i) = Dcp(i,2)/BL
       DDZ(i) = Dcp(i,3)/CL
    END DO

  END SUBROUTINE TestRotatedM

! This routine takes the current particle arrangement and finds the 
! seperations of all nearest neighbours.
! -----------------------------------------------------------------
  SUBROUTINE MeasureNNSeperations(imove)
    Use SystemDef
    IMPLICIT NONE
    INTEGER imove,outM,I,Jc,J,IJo

    outM = (ILAT*2-3)*NCOUNT
! Loop over spheres:
    DO I = 1,NM
! Loop over neighbours:
       Jc = 1
       J = NLIST(I,Jc,ILAT)
       DO WHILE (J.NE.NoMoreNeighbours)

          IF (J > I .AND. NLInfo(I,Jc,ILAT) == 1) THEN
             CALL Interaction(I,J,Jc,ILAT,.FALSE.,IJo)
             WRITE(*,FMT='(A,I11,I8,2I5,F12.7,E12.5)') ' CP ',imove,outM,I,J,SQRT(IJDr2),SQRT(IJDr2)-DIAMETER
          END IF

          Jc = Jc + 1
          J = NLIST(I,Jc,ILAT)
       ENDDO
    ENDDO

  END SUBROUTINE MeasureNNSeperations

! This routine takes the current particle arrangement and finds the 
! seperation and identity of the closest pair of particles.
! -----------------------------------------------------------------
  SUBROUTINE FindClosestPair
    Use SystemDef
    IMPLICIT NONE
    INTEGER CfindI,CfindJ,CfindJc,CfindIJover
    REAL*8 CfindDr2

    CloseDr2 = 0.0d0
! Loop over spheres:
    DO CfindI = 1,NM
! Loop over neighbours:
       CfindJc = 1
       CfindJ = NLIST(CfindI,CfindJc,ILAT)
       DO WHILE (CfindJ.NE.NoMoreNeighbours)

          CALL Interaction(CfindI,CfindJ,CfindJc,ILAT,.FALSE.,CfindIJover)
          IF ((CfindI == 1 .AND. CfindJc == 1) .OR. IJDr2<CloseDr2 ) THEN
             CloseDr2 = IJDr2
             CloseI = CfindI
             CloseJ = CfindJ
          END IF

          CfindJc = CfindJc + 1
          CfindJ = NLIST(CfindI,CfindJc,ILAT)
       ENDDO
    ENDDO
    
  END SUBROUTINE FindClosestPair

! Output local particle displacements to a file:
!-----------------------------------------------
  SUBROUTINE InitDisplOutput
    USE SystemDef
    IMPLICIT NONE

    OPEN(UNIT=46, FILE='displ.out',STATUS='UNKNOWN')

  END SUBROUTINE InitDisplOutput
!-------------------------------
  SUBROUTINE OutputDisplacements
    USE SystemDef
    IMPLICIT NONE
    INTEGER i,j,jc,olps,alatt
    REAL*8 dx,dy,dz,dr
    
    alatt=ILAT
    DO i=1,NM
! Which planes changes (above=+1, below=-1):
       IF ( PlaneChanges(WhichPlane(i)) == -1 ) THEN
! Count up local overlaps:
          olps = 0
          jc = 1
          j = NLIST(i,jc,alatt)
          DO WHILE (j.NE.NoMoreNeighbours)
             IF (NLInfo(i,jc,alatt) == 1 .AND. NOVER(i,jc)==true) olps = olps + 1
             jc = jc + 1
             j = NLIST(i,jc,alatt)
          ENDDO

          jc = 1
          j = NLIST(i,jc,alatt)
          DO WHILE (j.NE.NoMoreNeighbours)
             IF (NLInfo(i,jc,alatt) == 1 .AND. olps>0) THEN
!             IF (NOVER(i,jc)==true) THEN
             dx = (DDX(j)+XLAT(j,alatt)-DDX(i)-XLAT(i,alatt))
             dx = dx - AINT(dx)
             dx = dx - AINT(dx+dx)
             dy = (DDY(j)+YLAT(j,alatt)-DDY(i)-YLAT(i,alatt))
             dy = dy - AINT(dy)
             dy = dy - AINT(dy+dy)
             dz = (DDZ(j)+ZLAT(j,alatt)-DDZ(i)-ZLAT(i,alatt))
             dz = dz - AINT(dz)
             dz = dz - AINT(dz+dz)
!r             dr = SQRT(dx*dx*AL*AL + dy*dy*BL*BL + dz*dz*CL*CL)
!r             WRITE(46,FMT='(3E12.4)') dx*AL,dy*BL,dr
             WRITE(46,FMT='(3E12.4)') dx*AL,dy*BL,dz*CL
!             WRITE(46,FMT='(3E12.4)') DDX(i)*AL,DDY(i)*BL,DDZ(i)*CL
             END IF
             jc = jc + 1
             j = NLIST(i,jc,alatt)
          ENDDO
!          WRITE(46,FMT='(3E12.4)') DDX(i)*AL,DDY(i)*BL,DDZ(i)*CL
!          IF ( olps == 0 ) WRITE(46,FMT='(3E12.4)') DDX(i)*AL,DDY(i)*BL,DDZ(i)*CL
!          IF ( olps > 0 ) WRITE(46,FMT='(3E12.4)') DDX(i)*AL,DDY(i)*BL,DDZ(i)*CL
       END IF
    END DO

  END SUBROUTINE OutputDisplacements


!  Routine to log the Centre of Mass behaviour:
!--------------------------------------------------
  SUBROUTINE LogCoM
    USE SystemDef
    USE AnalysisData
    IMPLICIT NONE
    REAL*8 CoMx,CoMy,CoMz,CoMr,cx,cy,cz,cr
    INTEGER Im,Ibuc,Iatm,Mflag
    LOGICAL Out

!Anj: Set-up M:
    INACvM = 2*ILAT-3
    M = INACvM*NCOUNT


! Initialise the CoM vector:
    CoMx = 0.0d0; CoMy = 0.0d0; CoMz = 0.0d0; CoMr = 0.0d0

! Calc the CoM as the average vector sum of the set {r}:
    DO Iatm = 1,NM
       cx = DDX(Iatm)*AL
       cy = DDY(Iatm)*BL
       cz = DDZ(Iatm)*CL
       CoMx = CoMx + cx
       CoMy = CoMy + cy
       CoMz = CoMz + cz
    END DO
    CoMx = CoMx/DBLE(NM)
    CoMy = CoMy/DBLE(NM)
    CoMz = CoMz/DBLE(NM)
    CoMr = CoMr/DBLE(NM)
    CoMr = SQRT(CoMx*CoMx + CoMy*CoMy + CoMz*CoMz)

    WRITE(*,FMT='(A,I5,4e14.6)') ' C ',M,CoMx,CoMy,CoMz,CoMr

  END SUBROUTINE LogCoM

!  Define the set of particles from which to examine plane separation:
!----------------------------------------------------------------------
  SUBROUTINE InitPlaneSep
    USE SystemDef
    USE AnalysisData
    IMPLICIT NONE
    CHARACTER*32 filename
    INTEGER iset,inum,im,idir,indx,isc,Ibuc

! Particles are numbered such that 1-N**(2/3) are in the first plane, etc..
! Loop over planes
    IdvMset = SystemSize*6
    IdvMNum = IdvMSet*IdvMSet
    DO iset = 1,IdvMSet
       DO isc = 1,IdvMNum
          dvMset(iset,isc) = (iset-1)*IdvMNum + isc
!          WRITE(72,*) iset,isc,dvMSet(iset,isc),ZLAT(dvMSet(iset,isc),1)
       END DO
    END DO

    RETURN
  END SUBROUTINE InitPlaneSep


! Output the current displacements etc for each of the
! sets of particles defined in PlaneSepInit.
!-----------------------------------------------------
  SUBROUTINE LogPlaneSep
    USE SystemDef
    USE AnalysisData
    IMPLICIT NONE
    REAL*8 dispz,difactor,diifactor
    INTEGER inum,iset,iiset,iiiset,I,IM,Ibuc,Mflag,normfac

!Anj: Set-up M:
    INACvM = 2*ILAT-3
    M = INACvM*NCOUNT

!  Loop over the set of particles and find the average z displacement for each plane:
    DO iset = 1,IdvMSet
       dispz = 0.0
       DO inum = 1,IdvMNum
          I = dvMSet(iset,inum)
          dispz = dispz + DDZ(I)*CL
       END DO
       Pdisplz(iset) = dispz/DBLE(IdvMNum)
! Store results for the planes in an array:
!       WRITE(*,FMT='(A,I5,E14.6)') 'Z ',M,Pdisplz(iset)

    END DO

! Calculate the d factor:
    difactor = 0.0d0; diifactor = 0.0d0
    normfac = 0
    DO iset = 1,IdvMSet,2
       iiset = iset + 1
       IF (iiset>IdvMSet) iiset = 1
       iiiset = iset - 1
       IF (iiiset<1) iiiset = IdvMSet
       difactor = difactor + (Pdisplz(iiset) - Pdisplz(iset))
       diifactor = diifactor + (Pdisplz(iset) - Pdisplz(iiiset))
       normfac = normfac + 1
    ENDDO
    difactor = difactor/DBLE(normfac)
    diifactor = diifactor/DBLE(normfac)
    WRITE(*,FMT='(A,I5,2E14.6)') 'D ',M,difactor,diifactor

    RETURN
  END SUBROUTINE LogPlaneSep


!----------------------------------------------------------------------------
! Library routines:
!----------------------------------------------------------------------------
! Generic Histogram taking routines, 
! supports up to ten histograms, each of up to 200 buckets:
! ---------------------------------------------------------
  MODULE Histo_Mod
    INTEGER HistMAX,BinMAX
    PARAMETER(HistMAX = 10, BinMAX = 200)
    INTEGER UnitHist,AveHist
    PARAMETER(UnitHist = 1, AveHist = 2)
    INTEGER MAXMoments
    PARAMETER(MAXMoments = 4)
    INTEGER nhistos
    DATA nhistos/0/
    REAL*8 hlob(HistMAX),hhib(HistMAX)
    INTEGER hnb(HistMAX)
    REAL*8 hdata(HistMAX,BinMAX)
    INTEGER hndata(HistMAX),hnbdata(HistMAX,BinMAX)
    REAL*8 hmoments(HistMAX,MAXMoments)
    INTEGER HistType(HistMAX)
  END MODULE Histo_Mod
! ---------------------------------
  FUNCTION Init_Histogram(lo,hi,n)
    Use Histo_Mod
    INTEGER Init_Histogram,n,ib,im
    REAL*8 lo,hi

! Store as next histogram:
    nhistos = nhistos + 1
    IF (nhistos > HistMax) THEN
       WRITE(*,*) 'Init_Histogram: Exceeded max num. of histograms,',HistMAX
       nhistos = nhistos - 1
       Init_Histogram = -1
       RETURN
    END IF

! Store the bounds:
    hlob(nhistos) = lo
    hhib(nhistos) = hi
    IF ( n <= BinMAX .AND. n>0 ) THEN
       hnb(nhistos) = n
    ELSE
       IF ( n > BinMax ) WRITE(*,*) 'Init_Histogram: Too many buckets! Using ',BinMAX
       IF ( n <= 0 )     WRITE(*,*) 'Init_Histogram: <= zero buckets! Using ',BinMAX
       hnb(nhistos) = BinMAX
    END IF

! Init histogram to zero:
    hndata(nhistos) = 0
    DO ib = 1,hnb(nhistos)
       hdata(nhistos,ib) = 0
       hnbdata(nhistos,ib) = 0
    END DO
    DO im = 1,MAXMoments
       hmoments(nhistos,im) = 0.0
    END DO

! Return histgram handle:
    Init_Histogram = nhistos

  END FUNCTION Init_Histogram
! ---------------------------------
  SUBROUTINE Add_To_Histogram(hh,hx)
    Use Histo_Mod
    INTEGER hh,bin,im
    REAL*8 hx
    
! Ignore bad histogram handles
    IF (hh < 1) RETURN

! Assign type to histogram:
    HistType(hh) = UnitHist

! Calculate which bin this number belongs in:
    bin = INT(DBLE(hnb(hh))*(hx - hlob(hh))/(hhib(hh) - hlob(hh))) + 1

! Add it to the histogram if it fits inside it:
    IF (bin > 0 .AND. bin <=BinMAX) THEN
       hndata(hh) = hndata(hh) + 1
       hdata(hh,bin) = hdata(hh,bin) + 1
    END IF

! Add up the moments also:
    DO im = 1,MAXMoments
       hmoments(hh,im) = hmoments(hh,im) + hx**im
    END DO

  END SUBROUTINE Add_To_Histogram
! ---------------------------------
  SUBROUTINE Add_This_To_Histogram(hh,hx,hadd)
    Use Histo_Mod
    INTEGER hh,bin,im
    REAL*8 hx,hadd
    
! Ignore bad histogram handles
    IF (hh < 1) RETURN

! Assign type to histogram:
    HistType(hh) = AveHist

! Calculate which bin this number belongs in:
    bin = INT(DBLE(hnb(hh))*(hx - hlob(hh))/(hhib(hh) - hlob(hh))) + 1

! Add it to the histogram if it fits inside it:
    IF (bin > 0 .AND. bin <=BinMAX) THEN
       hndata(hh) = hndata(hh) + 1
       hdata(hh,bin) = hdata(hh,bin) + hadd
       hnbdata(hh,bin) = hnbdata(hh,bin) + 1
    END IF

! Add up the moments also:
    DO im = 1,MAXMoments
       hmoments(hh,im) = hmoments(hh,im) + hadd*hx**im
    END DO

  END SUBROUTINE Add_This_To_Histogram
! ---------------------------------
  SUBROUTINE Output_Histogram(hh,hname)
    Use Histo_Mod
    INTEGER hh,ib,im
    REAL*8 hx,hy,normal
    CHARACTER*(*) hname
    
! Ignore bad histogram handles
    IF (hh < 1) RETURN

! Calculate normalization:
    normal = 0.0d0
    DO ib = 1,hnb(hh)
       IF ( HistType(hh) == UnitHist ) THEN
          normal = normal + hdata(hh,ib)*(hhib(hh)-hlob(hh))/DBLE(hnb(hh))
       END IF
    END DO

! Open the file:
    OPEN(UNIT=73, FILE=hname, STATUS='UNKNOWN')

! Output the moments:
    IF (hndata(hh).NE.0) THEN
    DO im = 1,MAXMoments
       WRITE(73,FMT="('# moment ',I2,' = ',E16.5)") im,hmoments(hh,im)/DBLE(hndata(hh))
    END DO
    WRITE(73,FMT="('# std       = ',E16.5)") SQRT(ABS(hmoments(hh,2)/DBLE(hndata(hh))-(hmoments(hh,1)/DBLE(hndata(hh)))**2))
    END IF

! Output the histogram:
    DO ib = 1,hnb(hh)
       hx = hlob(hh) + (DBLE(ib)-0.5d0)*(hhib(hh)-hlob(hh))/DBLE(hnb(hh))

       IF ( normal > 0.0 .AND. HistType(hh) == UnitHist ) THEN
          hy = hdata(hh,ib)/normal
       ELSE IF ( HistType(hh) == AveHist .AND. hnbdata(hh,ib) > 0 ) THEN
          hy = (hdata(hh,ib)/DBLE(hnbdata(hh,ib)))
       ELSE
          hy = 0
       END IF
       WRITE(73,*) hx,hy
    END DO

! Close the file:
    CLOSE(UNIT=73)

  END SUBROUTINE Output_Histogram
 
! Takes a string and a number and creates a suitable filename 
! string such as 'EchoAndBounce667.dat'
! -----------------------------------------------------------
  SUBROUTINE NameAndNumber(head,num,newstr)
    IMPLICIT NONE
    CHARACTER*(*) head,newstr
    CHARACTER(4) numstr
    INTEGER num,lenh

    newstr = ' '
    lenh = LEN(TRIM(head))
    newstr(1:lenh)=head(1:lenh)
    WRITE(numstr,FMT='(I4)') num+1000
    newstr(lenh+1:lenh+3)=numstr(2:4)
    newstr(lenh+4:lenh+7)='.dat'

  END SUBROUTINE

! Code to generate a gaussian random deviate of a given standard deviation:
! std    = half-half width of the gaussian you want.
!
! Basically from Press et al (p.280), with some modification.
! -------------------------------------------------------------------------
  FUNCTION GaussPos(std)
    USE SystemDef
    IMPLICIT NONE
    REAL*8 GaussPos,std,val,amp
    PARAMETER (amp=0.3989422804d0)
    REAL RVEC(2),ran1,ran2
    LOGICAL loop
    INTEGER iset
    REAL*8 fac,gset,rsq,v1,v2
    SAVE gset,iset,ran2
    DATA iset/0/

      IF (iset.EQ.0) THEN
        loop = .TRUE.
        DO WHILE (loop)
          CALL RANMAR(RVEC,2)
          ran1 = RVEC(1)
          ran2 = RVEC(2)
          v1 = 2.0d0*ran1-1.0d0
          v2 = 2.0d0*ran2-1.0d0
          rsq = v1*v1 + v2*v2
          IF (rsq.GT.0.0d0 .AND. rsq.LT.1.0) loop = .FALSE.
        END DO
        fac = sqrt(-2.0d0*log(rsq)/rsq)
!       fac = 1.0d0
        gset = v1*fac
        val = v2*fac
        GaussPos = val*std
        iset = 1
      ELSE
        GaussPos = gset*std
        iset = 0
      END IF

    RETURN
  END FUNCTION GaussPos

! Returns the natural logarithm of the probatility of a given value from a gaussian dist:
! --------------------------------------------------------------
  FUNCTION GaussProb(std,pos)
   IMPLICIT NONE
    REAL*8 GaussProb,std,pos,amp
    PARAMETER (amp=0.3989422804d0)

!    GaussProb = log(amp/std) - 0.5d0*pos*pos/(std*std)
    GaussProb = 1.0d0

    RETURN
  END FUNCTION GaussProb

! -------------------------------------------------------------------
  INTEGER FUNCTION RMnewseed()
    IMPLICIT NONE
    REAL RV(1)
    INTEGER ns,MAX
    PARAMETER (MAX = 900000000)

    ns = 0
    DO WHILE (ns.EQ.0 .OR. ns.EQ.MAX)
       CALL RANMAR(RV,1)
       ns = RV(1)*MAX
    END DO

    RMnewseed = ns
    RETURN
  END FUNCTION RMnewseed

! -------------------------------------------------------------------
  SUBROUTINE RMARIN(IJKL)

!       INITIALIZING ROUTINE FOR RANMAR. THE INPUT VALUE SHOULD
!       BE IN THE RANGE:   0 <= IJKL <= 900 000 000
!       TO GET THE STANDARD VALUES IN THE MARSAGLIA - ZAMAN PAPER
!       (I=12, J=34, K=56, L=78) PUT  IJKL = 54217137

    INTEGER IJ,IJKL,KL,I,J,K,L,II,JJ,M,I97,J97,IM_RAN
    REAL  S,T,U,C,CD,CM
    COMMON/RASET1/U(97),C,CD,CM,I97,J97,IM_RAN(97)
    SAVE /RASET1/

    IJ = IJKL / 30082
    KL = IJKL - IJ * 30082
    I  = MOD(IJ/177,177) + 2
    J  = MOD(IJ,177) + 2
    K  = MOD(KL/169,178) + 1
    L  = MOD(KL,169)
!        WRITE(*,*) 'RANMAR INITIALIZED: ',IJKL,I,J,K,L
    DO II=1,97
       S = 0.
       T = 0.5
       DO JJ=1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF(MOD(L*M,64).GE.32) S = S + T
          T = 0.5 * T
       ENDDO
       U(II) = S
    ENDDO
    C = 362436. / 16777216.
    CD = 7654321. / 16777216.
    CM = 16777213. / 16777216.
    I97 = 97
    J97 = 33
    DO II=1,97
       IM_RAN(II) = II-1
    ENDDO
    IM_RAN(1) = 97
    RETURN
  END SUBROUTINE RMARIN


! -----------------------------------------------------------------------------
  SUBROUTINE RANMAR(RVEC,LEN)

!       UNIVERSAL RANDOM NUMBER GENERATOR PROPOSED BY MARSAGLIA
!       AND ZAMAN IN REPORT FSU-SCRI-87-50
!       GENERATES VECTOR 'RVEC' OF LENGTH 'LEN' OF PSEUDORANDOM
!       NUMBERS; THE COMMON BLOCK INCLUDES EVERYTHING NEEDED TO
!       COMPLETELY SPECIFY THE STATE OF THE GENERATOR.

    INTEGER IVEC,LEN,I97,J97,IM_RAN
    REAL  UNI,U,C,CD,CM,RVEC
    DIMENSION RVEC(*)
    COMMON/RASET1/U(97),C,CD,CM,I97,J97,IM_RAN(97)
    SAVE /RASET1/

    DO IVEC=1,LEN
       UNI = U(I97) - U(J97)
       IF(UNI.LT.0.) UNI = UNI + 1.
       U(I97) = UNI
!Shouldn't this be quicker...
       I97 = I97 - 1
       IF (I97.EQ.0) I97 = 97
       J97 = J97 - 1
       IF (J97.EQ.0) J97 = 97
!than this...
!       I97 = IM_RAN(I97)
!       J97 = IM_RAN(J97)
! Yes it is (by about 2%)!
       C = C - CD
       IF (C.LT.0.) C = C + CM
       UNI = UNI - C
       IF (UNI.LT.0.) UNI = UNI + 1.
       RVEC(IVEC) = UNI
    ENDDO
    RETURN
  END SUBROUTINE RANMAR

! -----------------------------------------------------------------------------
! <-- End of trans.f90
! -----------------------------------------------------------------------------
