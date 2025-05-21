MODULE icbdyn
   !!======================================================================
   !!                       ***  MODULE  icbdyn  ***
   !! Iceberg:  time stepping routine for iceberg tracking
   !!======================================================================
   !! History :  3.3  !  2010-01  (Martin&Adcroft)  Original code
   !!             -   !  2011-03  (Madec)  Part conversion to NEMO form
   !!             -   !                    Removal of mapping from another grid
   !!             -   !  2011-04  (Alderson)  Split into separate modules
   !!             -   !  2011-05  (Alderson)  Replace broken grounding routine with one of
   !!             -   !                       Gurvan's suggestions (just like the broken one)
   !!----------------------------------------------------------------------
   USE par_oce        ! NEMO parameters
   USE dom_oce        ! NEMO ocean domain
   USE phycst         ! NEMO physical constants
   USE in_out_manager                      ! IO parameters
   !
   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines
   USE icbdia         ! iceberg budget routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_dyn  ! routine called in icbstp.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: icbdyn.F90 15088 2021-07-06 13:03:34Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_dyn  ***
      !!
      !! ** Purpose :   iceberg evolution.
      !!
      !! ** Method  : - See Martin & Adcroft, Ocean Modelling 34, 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   !
      !
      LOGICAL  ::   ll_bounced, ll_bouncedx, ll_bouncedy, ll_silted
      REAL(wp) ::   zdrag_groundedx, zdrag_groundedy
      REAL(wp) ::   zdrag_groundedxG, zdrag_groundedyG
      REAL(wp) ::   zdrag_groundedxS, zdrag_groundedyS
      REAL(wp) ::   Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori
      REAL(wp) ::   Xdrag_ocn_stored, Xdrag_atm_stored, Xdrag_ice_stored, Xgrad_press_sc_stored, Xdrag_Cori_stored
      REAL(wp) ::   Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori
      REAL(wp) ::   Ydrag_ocn_stored, Ydrag_atm_stored, Ydrag_ice_stored, Ygrad_press_sc_stored, Ydrag_Cori_stored
      REAL(wp) ::   zuvel1 , zvvel1 , zu1, zv1, zax1, zay1, zxi1 , zyj1
      REAL(wp) ::   zuvel2 , zvvel2 , zu2, zv2, zax2, zay2, zxi2 , zyj2
      REAL(wp) ::   zuvel3 , zvvel3 , zu3, zv3, zax3, zay3, zxi3 , zyj3
      REAL(wp) ::   zuvel4 , zvvel4 , zu4, zv4, zax4, zay4, zxi4 , zyj4
      REAL(wp) ::   zuvel_n, zvvel_n, zxi_n   , zyj_n
      REAL(wp) ::   zdt, zdt_2, zdt_6, ze1, ze2, pe1, pe2
      TYPE(iceberg), POINTER ::   berg
      TYPE(point)  , POINTER ::   pt
      INTEGER  ::   ii, ii0
      INTEGER  ::   ij, ij0
      REAL(wp) ::   zT, zD, zF, zFequil, zFdisequil
      REAL(wp) :: issh_x0y0, issh_x1y1
      REAL(wp) :: ibathy_x0y0, ibathy_x0y1, ibathy_x1y0, ibathy_x1y1
      REAL(wp) :: ibathy_xx1y1, ibathy_x1yy1
      REAL(wp) :: distX1, distX2, distY1, distY2, dist
      REAL(wp) :: pe1_x0y0, pe1_x0y1, pe1_x1y0, pe1_x1y1
      REAL(wp) :: pe2_x0y0, pe2_x0y1, pe2_x1y0, pe2_x1y1
      REAL(wp) :: angleTopoMaxGradient, angleTopoTraject
      REAL(wp) :: angleTopoX, angleTopoY
      REAL(wp) :: angleTopoXlagr, angleTopoYlagr
      REAL(wp) :: angleTopoXlagrOld, angleTopoYlagrOld, angleTopolagrOld
      REAL(wp) :: gravity_termX, gravity_termY
      REAL(wp) :: TopoX1, TopoX2, TopoY1, TopoY2
      INTEGER  ::   FlagStopX, FlagStopY
      INTEGER  ::   AnyFlagStopX, AnyFlagStopY
      REAL(wp) ::   pi , pj      ! current iceberg position
      REAL(wp) ::   pi0, pj0     ! previous iceberg position
      REAL(wp) ::   pu  , pv     ! current iceberg velocities
      REAL(wp) ::   zdrag_groundedx_stored, zdrag_groundedy_stored
      REAL(wp) ::   zdrag_groundedxG_stored, zdrag_groundedyG_stored
      REAL(wp) ::   zdrag_groundedxS_stored, zdrag_groundedyS_stored      
      REAL(wp) ::   zaxe_stored, zaye_stored
      REAL(wp) ::   zaxeBeforeFG, zayeBeforeFG
      REAL(wp) ::   zaxeBeforeFG_stored, zayeBeforeFG_stored


      ! Parameter of solid-body kinetic Coulomb friction
      REAL(wp) ::   muFric 
      REAL(wp) ::   unitWeightSilt
      REAL(wp) ::   shearStrengthSilt
      REAL(wp) ::   depthSilt
      REAL(wp) ::   taperCoeffK
      REAL(wp) ::   StopThreshold
      REAL(wp) ::   SiltScourWidth
      INTEGER  ::   ibounce_method
      INTEGER  ::   berg_number

      ibounce_method = 1

      muFric       = 0.002_wp
      ! 0.00002_wp ! 0._wp !0.2_wp !0.002_wp ! 0.2_wp !0.00002_wp ! 0.002_wp ! 0.10_wp
      unitWeightSilt     = 8000._wp
      shearStrengthSilt     = 6000._wp
      depthSilt    = 8._wp
      taperCoeffK  = 200000.0_wp
      StopThreshold = 5.E-5_wp

      IF (ibounce_method/=3) THEN
       muFric       = 0._wp
       unitWeightSilt     = 0._wp
       shearStrengthSilt     = 0._wp
       depthSilt    = 0._wp
       taperCoeffK  = 200000.0_wp
       StopThreshold = 0._wp
      ENDIF


      ! Yavor Kostov: The staggered time-stepping of the legacy
      ! 4th order Runge-Kutta method below appears to be non-standard,
      ! especially when updating velocity. This may be an intended feature.
      ! It is beyond the scope of Kostov et al. (2025) to modify
      ! the staggering of the Runge-Kutta stages.
      ! As of now, all Kostov et al. (2025) updates are designed to be
      ! compatible with the legacy staggering.
      !!----------------------------------------------------------------------
      !
      ! 4th order Runge-Kutta to solve:   d/dt X = V,  d/dt V = A
      !                    with I.C.'s:   X=X1 and V=V1
      !
      !                                    ; A1=A(X1,V1)
      !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1 ; A2=A(X2,V2)
      !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2 ; A3=A(X3,V3)
      !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3 ; A4=A(X4,V4)
      !
      !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
      !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
      
      ! time steps
      zdt   = berg_dt
      zdt_2 = zdt * 0.5_wp
      zdt_6 = zdt / 6._wp

      berg => first_berg                    ! start from the first berg
      !
      DO WHILE ( ASSOCIATED(berg) )          !==  loop over all bergs  ==!
         !
         pt => berg%current_point
         berg_number = berg%number(1)
         SiltScourWidth = 90._wp ! 0.1_wp*berg%current_point%width

         ll_bounced = .FALSE.
         ll_bouncedx = .FALSE.
         ll_bouncedy = .FALSE.
         ll_silted   = .FALSE.
         zdrag_groundedx = 0._wp
         zdrag_groundedy = 0._wp
         zdrag_groundedxG = 0._wp
         zdrag_groundedyG = 0._wp
         zdrag_groundedxS = 0._wp
         zdrag_groundedyS = 0._wp
         angleTopoX=0._wp
         gravity_termX=0._wp
         angleTopoY=0._wp
         gravity_termY=0._wp
         FlagStopX=0
         FlagStopY=0
         AnyFlagStopX=0
         AnyFlagStopY=0
         zaxeBeforeFG=0._wp
         zayeBeforeFG=0._wp
         zaxeBeforeFG_stored=0._wp
         zayeBeforeFG_stored=0._wp
         Xdrag_ocn_stored=0._wp
         Xdrag_atm_stored=0._wp
         Xdrag_ice_stored=0._wp
         Xgrad_press_sc_stored=0._wp
         Xdrag_Cori_stored=0._wp
         Ydrag_ocn_stored=0._wp
         Ydrag_atm_stored=0._wp
         Ydrag_ice_stored=0._wp
         Ygrad_press_sc_stored=0._wp
         Ydrag_Cori_stored=0._wp
 
         angleTopoXlagrOld=pt%angleTopoXlagr
         angleTopoYlagrOld=pt%angleTopoYlagr
         angleTopolagrOld =pt%angleTopolagr

         ! STEP 1 !
         ! ====== !
         zxi1 = pt%xi   ;   zuvel1 = pt%uvel     !**   X1 in (i,j)  ;  V1 in m/s
         zyj1 = pt%yj   ;   zvvel1 = pt%vvel


         ii0 = INT(zxi1+0.5_wp) + (nn_hls-1)  ;  ij0 = INT(zyj1+0.5_wp) + (nn_hls-1)
         ii0 = mi1( ii0 )
         ij0 = mj1( ij0 )
         zT = berg%current_point%thickness               ! total thickness
         zD = rho_berg_1_oce * zT                        ! draught (keel depth)
         zF = zT - zD                                    ! freeboard
         zFequil = zF                                    ! unless the iceberg is
                                                      ! grounded, the freeboard
                                                      ! is consistent with
                                                      ! an equilibrium of
                                                      ! vertical forces
         zFdisequil = 0._wp
         zdrag_groundedx_stored = 0._wp
         zdrag_groundedy_stored = 0._wp
         zdrag_groundedxG_stored =  0._wp
         zdrag_groundedyG_stored = 0._wp
         zdrag_groundedxS_stored =  0._wp
         zdrag_groundedyS_stored = 0._wp         
         zaxe_stored = 0._wp
         zaye_stored = 0._wp

     
      pi0=zxi1
      pj0=zyj1
      pu= zuvel1
      pv= zvvel1
     
      CALL icb_utl_interp( pi0, pj0, pe1=pe1, pe2=pe2 )       
      pi = pi0 + (pu/pe1)*zdt_2 !The subsequent call to the grounding algorithm...
      pj = pj0 + (pv/pe2)*zdt_2 !needs a first estimate of the projected location


        IF (ibounce_method .EQ. 3) THEN
                 CALL icb_ground( kt, berg, pi, pi0, zuvel1, zuvel1,   &
            &                   pj, pj0, zvvel1, zvvel1, ll_bounced, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted,&
            &                   zdrag_groundedx, zdrag_groundedy, &
            &                   zdrag_groundedxS, zdrag_groundedyS, &
            &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
            &                   ibounce_method, muFric, unitWeightSilt, shearStrengthSilt,depthSilt,SiltScourWidth)
        ENDIF


         !                                         !**   A1 = A(X1,V1)
         CALL icb_accel( kt, berg , zxi1, ze1, zuvel1, zuvel1, zax1,     &
            &                   zyj1, ze2, zvvel1, zvvel1, zay1, zdt_2, 0.5_wp,&
            &                   ll_bouncedx, ll_bouncedy, ll_silted,&
            &            zdrag_groundedxS,zdrag_groundedyS,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth,muFric,&
            &            zdrag_groundedx,zdrag_groundedy,FlagStopX,FlagStopY, &
            & zdrag_groundedxG,zdrag_groundedyG,zFdisequil,zaxeBeforeFG, zayeBeforeFG,&
            &            Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori,&
            &            Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori)

         !
         zu1 = zuvel1 / ze1                           !**   V1 in d(i,j)/dt
         zv1 = zvvel1 / ze2

         ! At each RK4 stage, we update the weighted sum of the acceleration components, accordingly
         zdrag_groundedx_stored=zdrag_groundedx_stored + zdrag_groundedx/6._wp
         zdrag_groundedy_stored=zdrag_groundedy_stored + zdrag_groundedy/6._wp
         zdrag_groundedxG_stored=zdrag_groundedxG_stored + zdrag_groundedxG/6._wp
         zdrag_groundedyG_stored=zdrag_groundedyG_stored  + zdrag_groundedyG/6._wp
         zdrag_groundedxS_stored=zdrag_groundedxS_stored + zdrag_groundedxS/6._wp
         zdrag_groundedyS_stored=zdrag_groundedyS_stored  + zdrag_groundedyS/6._wp         
         zaxe_stored = zaxe_stored + zax1 / 6._wp
         zaye_stored = zaye_stored + zay1 / 6._wp
         zaxeBeforeFG_stored = zaxeBeforeFG_stored + zaxeBeforeFG/ 6._wp
         zayeBeforeFG_stored = zayeBeforeFG_stored + zayeBeforeFG/ 6._wp
         Xdrag_ocn_stored = Xdrag_ocn_stored + Xdrag_ocn/ 6._wp
         Xdrag_atm_stored = Xdrag_atm_stored + Xdrag_atm/ 6._wp
         Xdrag_ice_stored = Xdrag_ice_stored + Xdrag_ice/ 6._wp
         Xgrad_press_sc_stored = Xgrad_press_sc_stored + Xgrad_press_sc/ 6._wp
         Xdrag_Cori_stored = Xdrag_Cori_stored + Xdrag_Cori/ 6._wp
         Ydrag_ocn_stored = Ydrag_ocn_stored + Ydrag_ocn/ 6._wp
         Ydrag_atm_stored = Ydrag_atm_stored + Ydrag_atm/ 6._wp
         Ydrag_ice_stored = Ydrag_ice_stored + Ydrag_ice/ 6._wp
         Ygrad_press_sc_stored = Ygrad_press_sc_stored + Ygrad_press_sc/ 6._wp
         Ydrag_Cori_stored = Ydrag_Cori_stored + Ydrag_Cori/ 6._wp


         ! STEP 2 !
         ! ====== !
         !                                         !**   X2 = X1+dt/2*V1   ;   V2 = V1+dt/2*A1
         ! position using di/dt & djdt   !   V2  in m/s
         zxi2 = zxi1 + zdt_2 * zu1          ;   zuvel2 = zuvel1 + zdt_2 * zax1
         zyj2 = zyj1 + zdt_2 * zv1          ;   zvvel2 = zvvel1 + zdt_2 * zay1

         ! If any of the RK4 stages triggers a flag, it stops the iceberg motion
         IF (FlagStopX .EQ. 1) THEN
         zxi2 = zxi1                        ;   zuvel2 = 0._wp ; zu1 = 0._wp
         AnyFlagStopX = 1       
         ENDIF
         IF (FlagStopY .EQ. 1) THEN
         zyj2 = zyj1                        ;   zvvel2 = 0._wp ; zv1 = 0._wp
         AnyFlagStopY = 1
         ENDIF

         FlagStopX=0
         FlagStopY=0         
         !
        CALL icb_ground( kt, berg, zxi2, zxi1, zu1, zu1,   &
            &                   zyj2, zyj1, zv1, zv1, ll_bounced, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &                   zdrag_groundedx, zdrag_groundedy, &
            &                   zdrag_groundedxS, zdrag_groundedyS, &
            &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
            &                   ibounce_method, muFric, unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth)

         !                                         !**   A2 = A(X2,V2)
        CALL icb_accel( kt, berg , zxi2, ze1, zuvel2, zuvel1, zax2,    &
            &                   zyj2, ze2, zvvel2, zvvel1, zay2, zdt_2, 0.5_wp,&
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &            zdrag_groundedxS,zdrag_groundedyS,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth,muFric,&
            &            zdrag_groundedx,zdrag_groundedy,FlagStopX,FlagStopY, &
            &  zdrag_groundedxG,zdrag_groundedyG,zFdisequil,zaxeBeforeFG, zayeBeforeFG,&
            &            Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori,&
            &            Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori)
         !
         zu2 = zuvel2 / ze1                           !**   V2 in d(i,j)/dt
         zv2 = zvvel2 / ze2

         ! At each RK4 stage, we update the weighted sum of the acceleration components, accordingly
         zdrag_groundedx_stored=zdrag_groundedx_stored+2._wp*zdrag_groundedx/6._wp
         zdrag_groundedy_stored=zdrag_groundedy_stored+2._wp*zdrag_groundedy/6._wp
         zdrag_groundedxG_stored=zdrag_groundedxG_stored+2._wp*zdrag_groundedxG/6._wp
         zdrag_groundedyG_stored=zdrag_groundedyG_stored+2._wp*zdrag_groundedyG/6._wp
         zdrag_groundedxS_stored=zdrag_groundedxS_stored+2._wp*zdrag_groundedxS/6._wp
         zdrag_groundedyS_stored=zdrag_groundedyS_stored+2._wp*zdrag_groundedyS/6._wp         
         zaxe_stored = zaxe_stored + 2._wp*zax2 / 6._wp
         zaye_stored = zaye_stored + 2._wp*zay2 / 6._wp
         zaxeBeforeFG_stored = zaxeBeforeFG_stored + 2._wp*zaxeBeforeFG/ 6._wp
         zayeBeforeFG_stored = zayeBeforeFG_stored + 2._wp*zayeBeforeFG/ 6._wp
         Xdrag_ocn_stored = Xdrag_ocn_stored + 2._wp*Xdrag_ocn/ 6._wp
         Xdrag_atm_stored = Xdrag_atm_stored + 2._wp*Xdrag_atm/ 6._wp
         Xdrag_ice_stored = Xdrag_ice_stored + 2._wp*Xdrag_ice/ 6._wp
         Xgrad_press_sc_stored = Xgrad_press_sc_stored + 2._wp*Xgrad_press_sc/ 6._wp
         Xdrag_Cori_stored = Xdrag_Cori_stored + 2._wp*Xdrag_Cori/ 6._wp
         Ydrag_ocn_stored = Ydrag_ocn_stored + 2._wp*Ydrag_ocn/ 6._wp
         Ydrag_atm_stored = Ydrag_atm_stored + 2._wp*Ydrag_atm/ 6._wp
         Ydrag_ice_stored = Ydrag_ice_stored + 2._wp*Ydrag_ice/ 6._wp
         Ygrad_press_sc_stored = Ygrad_press_sc_stored + 2._wp*Ygrad_press_sc/ 6._wp
         Ydrag_Cori_stored = Ydrag_Cori_stored + 2._wp*Ydrag_Cori/ 6._wp

         !
         ! STEP 3 !
         ! ====== !
         !                                         !**  X3 = X1+dt/2*V2  ;   V3 = V1+dt/2*A2; A3=A(X3)
         zxi3  = zxi1  + zdt_2 * zu2   ;   zuvel3 = zuvel1 + zdt_2 * zax2
         zyj3  = zyj1  + zdt_2 * zv2   ;   zvvel3 = zvvel1 + zdt_2 * zay2
         IF (FlagStopX .EQ. 1) THEN
         zxi3 = zxi1                        ;   zuvel3 = 0._wp ; zu2 = 0._wp
         AnyFlagStopX = 1
         ENDIF
         IF (FlagStopY .EQ. 1) THEN
         zyj3 = zyj1                        ;   zvvel3 = 0._wp ; zv2 = 0._wp
         AnyFlagStopY = 1
         ENDIF

         FlagStopX=0
         FlagStopY=0
         !
         CALL icb_ground( kt, berg, zxi3, zxi1, zu1, zu2,   &
            &                   zyj3, zyj1, zv1, zv2, ll_bounced, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &                   zdrag_groundedx, zdrag_groundedy, &
            &                   zdrag_groundedxS, zdrag_groundedyS, &
            &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
            &                   ibounce_method, muFric, unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth)
         !                                         !**   A3 = A(X3,V3)
         CALL icb_accel( kt, berg , zxi3, ze1, zuvel3, zuvel1, zax3,    &
            &                   zyj3, ze2, zvvel3, zvvel1, zay3, zdt, 1._wp, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &            zdrag_groundedxS,zdrag_groundedyS,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth,muFric,&
            &            zdrag_groundedx,zdrag_groundedy,FlagStopX,FlagStopY, &
            &  zdrag_groundedxG,zdrag_groundedyG,zFdisequil,zaxeBeforeFG,zayeBeforeFG,&
            &            Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori,&
            &            Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori)
         !
         zu3 = zuvel3 / ze1                           !**   V3 in d(i,j)/dt
         zv3 = zvvel3 / ze2

         ! At each RK4 stage, we update the weighted sum of the acceleration components, accordingly
         zdrag_groundedx_stored=zdrag_groundedx_stored+2._wp*zdrag_groundedx/6._wp
         zdrag_groundedy_stored=zdrag_groundedy_stored+2._wp*zdrag_groundedy/6._wp
         zdrag_groundedxG_stored=zdrag_groundedxG_stored+2._wp*zdrag_groundedxG/6._wp
         zdrag_groundedyG_stored=zdrag_groundedyG_stored+2._wp*zdrag_groundedyG/6._wp
         zdrag_groundedxS_stored=zdrag_groundedxS_stored+2._wp*zdrag_groundedxS/6._wp
         zdrag_groundedyS_stored=zdrag_groundedyS_stored+2._wp*zdrag_groundedyS/6._wp         
         zaxe_stored = zaxe_stored + 2._wp*zax3 / 6._wp
         zaye_stored = zaye_stored + 2._wp*zay3 / 6._wp
         zaxeBeforeFG_stored = zaxeBeforeFG_stored + 2._wp*zaxeBeforeFG/ 6._wp
         zayeBeforeFG_stored = zayeBeforeFG_stored + 2._wp*zayeBeforeFG/ 6._wp
         Xdrag_ocn_stored = Xdrag_ocn_stored + 2._wp*Xdrag_ocn/ 6._wp
         Xdrag_atm_stored = Xdrag_atm_stored + 2._wp*Xdrag_atm/ 6._wp
         Xdrag_ice_stored = Xdrag_ice_stored + 2._wp*Xdrag_ice/ 6._wp
         Xgrad_press_sc_stored = Xgrad_press_sc_stored + 2._wp*Xgrad_press_sc/ 6._wp
         Xdrag_Cori_stored = Xdrag_Cori_stored + 2._wp*Xdrag_Cori/ 6._wp
         Ydrag_ocn_stored = Ydrag_ocn_stored + 2._wp*Ydrag_ocn/ 6._wp
         Ydrag_atm_stored = Ydrag_atm_stored + 2._wp*Ydrag_atm/ 6._wp
         Ydrag_ice_stored = Ydrag_ice_stored + 2._wp*Ydrag_ice/ 6._wp
         Ygrad_press_sc_stored = Ygrad_press_sc_stored + 2._wp*Ygrad_press_sc/ 6._wp
         Ydrag_Cori_stored = Ydrag_Cori_stored + 2._wp*Ydrag_Cori/ 6._wp

         ! STEP 4 !
         ! ====== !
         !                                         !**   X4 = X1+dt*V3   ;   V4 = V1+dt*A3
         zxi4 = zxi1 + zdt * zu3   ;   zuvel4 = zuvel1 + zdt * zax3
         zyj4 = zyj1 + zdt * zv3   ;   zvvel4 = zvvel1 + zdt * zay3
         IF (FlagStopX .EQ. 1) THEN
         zxi4 = zxi1                        ;   zuvel4 = 0._wp ; zu3 = 0._wp
          AnyFlagStopX = 1
         ENDIF
         IF (FlagStopY .EQ. 1) THEN
         zyj4 = zyj1                        ;   zvvel4 = 0._wp ; zv3 = 0._wp
         AnyFlagStopY = 1
         ENDIF

         FlagStopX=0
         FlagStopY=0
         CALL icb_ground( kt, berg, zxi4, zxi1, zu1, zu3,   &
            &                   zyj4, zyj1, zv1, zv3, ll_bounced, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &                   zdrag_groundedx, zdrag_groundedy, &
            &                   zdrag_groundedxS, zdrag_groundedyS, &
            &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
            &                   ibounce_method, muFric, unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth) 

         !                                         !**   A4 = A(X4,V4)
         CALL icb_accel( kt, berg , zxi4, ze1, zuvel4, zuvel1, zax4,&
            &                   zyj4, ze2, zvvel4, zvvel1, zay4, zdt, 1._wp, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted, &
            &            zdrag_groundedxS,zdrag_groundedyS,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth,muFric,&
            &            zdrag_groundedx,zdrag_groundedy,FlagStopX,FlagStopY,&
            & zdrag_groundedxG,zdrag_groundedyG,zFdisequil,zaxeBeforeFG,zayeBeforeFG,&
            &            Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori,&
            &            Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori)
         zu4 = zuvel4 / ze1                           !**   V4 in d(i,j)/dt
         zv4 = zvvel4 / ze2

         ! At each RK4 stage, we update the weighted sum of the acceleration components, accordingly
         zdrag_groundedx_stored=zdrag_groundedx_stored + zdrag_groundedx/6._wp
         zdrag_groundedy_stored=zdrag_groundedy_stored + zdrag_groundedy/6._wp
         zdrag_groundedxG_stored=zdrag_groundedxG_stored + zdrag_groundedxG/6._wp
         zdrag_groundedyG_stored=zdrag_groundedyG_stored + zdrag_groundedyG/6._wp
         zdrag_groundedxS_stored=zdrag_groundedxS_stored + zdrag_groundedxS/6._wp
         zdrag_groundedyS_stored=zdrag_groundedyS_stored + zdrag_groundedyS/6._wp         
         zaxe_stored = zaxe_stored + zax4 / 6._wp
         zaye_stored = zaye_stored + zay4 / 6._wp
         zaxeBeforeFG_stored = zaxeBeforeFG_stored + zaxeBeforeFG/ 6._wp
         zayeBeforeFG_stored = zayeBeforeFG_stored + zayeBeforeFG/ 6._wp
         Xdrag_ocn_stored = Xdrag_ocn_stored + Xdrag_ocn/ 6._wp
         Xdrag_atm_stored = Xdrag_atm_stored + Xdrag_atm/ 6._wp
         Xdrag_ice_stored = Xdrag_ice_stored + Xdrag_ice/ 6._wp
         Xgrad_press_sc_stored = Xgrad_press_sc_stored + Xgrad_press_sc/ 6._wp
         Xdrag_Cori_stored = Xdrag_Cori_stored + Xdrag_Cori/ 6._wp
         Ydrag_ocn_stored = Ydrag_ocn_stored + Ydrag_ocn/ 6._wp
         Ydrag_atm_stored = Ydrag_atm_stored + Ydrag_atm/ 6._wp
         Ydrag_ice_stored = Ydrag_ice_stored + Ydrag_ice/ 6._wp
         Ygrad_press_sc_stored = Ygrad_press_sc_stored + Ygrad_press_sc/ 6._wp
         Ydrag_Cori_stored = Ydrag_Cori_stored + Ydrag_Cori/ 6._wp

         IF (FlagStopX .EQ. 1) THEN
         zu4 = 0._wp
         AnyFlagStopX = 1
         ENDIF
         IF (FlagStopY .EQ. 1) THEN
         zv4 = 0._wp
         AnyFlagStopY = 1
         ENDIF
         FlagStopX=0
         FlagStopY=0



         ! FINAL STEP !
         ! ========== !
         !                                         !**   Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
         !                                         !**   Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
         zxi_n   = pt%xi   + zdt_6 * (  zu1  + 2.*(zu2  + zu3 ) + zu4  )
         zyj_n   = pt%yj   + zdt_6 * (  zv1  + 2.*(zv2  + zv3 ) + zv4  )
         zuvel_n = pt%uvel + zdt_6 * (  zax1 + 2.*(zax2 + zax3) + zax4 )
         zvvel_n = pt%vvel + zdt_6 * (  zay1 + 2.*(zay2 + zay3) + zay4 )




         IF (AnyFlagStopX .EQ. 1) THEN
         zuvel_n = 0._wp 
         zxi_n   = pt%xi
         ENDIF
         IF (AnyFlagStopY .EQ. 1) THEN
         zvvel_n = 0._wp
         zyj_n   = pt%yj  
         ENDIF
         IF (((AnyFlagStopX .EQ. 1).OR.(AnyFlagStopY .EQ. 1)) .AND. &
            & (SQRT(zuvel_n**2 + zvvel_n**2).LT.StopThreshold)) THEN
          zxi_n   = pt%xi
          zyj_n   = pt%yj 
          zuvel_n = 0._wp
          zvvel_n = 0._wp
         ENDIF

         CALL icb_ground( kt, berg, zxi_n, zxi1, zu1, zuvel_n,   &
            &                   zyj_n, zyj1, zv1, zvvel_n, ll_bounced, &
            &                   ll_bouncedx, ll_bouncedy, ll_silted,&
            &                   zdrag_groundedx, zdrag_groundedy, &
            &                   zdrag_groundedxS, zdrag_groundedyS, &
            &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
            &                   ibounce_method, muFric, unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth)

        ! We have to do a final calculation of the resulting freeboard disequilibrium (if any),
        ! the angle of the trajectory if sliding along topography, and the maximum topographic angle.

         CALL icb_utl_interp( zxi_n, zyj_n,     &
         &                 ibathy=ibathy_x1y1                  &   ! bathy
         &                  )
         CALL icb_utl_interp( zxi_n, zyj_n,     &
         &                 issh=issh_x1y1                  &   ! ssh
         &                  )
         pt%ibathy_x1y1 = ibathy_x1y1
         pt%issh_x1y1 = issh_x1y1

         pt%zFequil = zFequil

        ! The freeboard calculation for grounded icebergs takes into account the SSH, the background bathymetry depth,
        ! and the assumed sediment depth where iceberg keels can penetrate
        ! If the total iceberg thickness is greater than (the equilibrium freeboard) +SSH +(background bathymetry)
        ! +(sedimend thickness), then we have a freeboard out of equilibrium w.r.t. the Archimedes balance
        zFdisequil = 0._wp
        IF (((zT-depthSilt-ABS(ibathy_x1y1)-issh_x1y1).GT.zFequil).AND.(ibounce_method.EQ.3)) THEN
            zF = zT -depthSilt - ibathy_x1y1 -issh_x1y1
            zFdisequil = MAX(0._wp, zF - zFequil)
        ENDIF

         pt%zFdisequil = zFdisequil

         ! We save the acceleration components from this timestep
         pt%zdrag_groundedx_stored = zdrag_groundedx_stored  
         pt%zdrag_groundedy_stored = zdrag_groundedy_stored
         pt%zdrag_groundedxG_stored = zdrag_groundedxG_stored
         pt%zdrag_groundedyG_stored = zdrag_groundedyG_stored
         pt%zdrag_groundedxS_stored = zdrag_groundedxS_stored
         pt%zdrag_groundedyS_stored = zdrag_groundedyS_stored
         pt%zaxe_stored = zaxe_stored
         pt%zaye_stored = zaye_stored
         pt%zaxeBeforeFG_stored = zaxeBeforeFG_stored
         pt%zayeBeforeFG_stored = zayeBeforeFG_stored
         pt%Xdrag_ocn_stored = Xdrag_ocn_stored
         pt%Xdrag_atm_stored = Xdrag_atm_stored
         pt%Xdrag_ice_stored = Xdrag_ice_stored
         pt%Xgrad_press_sc_stored = Xgrad_press_sc_stored
         pt%Xdrag_Cori_stored = Xdrag_Cori_stored
         pt%Ydrag_ocn_stored = Ydrag_ocn_stored
         pt%Ydrag_atm_stored = Ydrag_atm_stored
         pt%Ydrag_ice_stored = Ydrag_ice_stored
         pt%Ygrad_press_sc_stored = Ygrad_press_sc_stored
         pt%Ydrag_Cori_stored = Ydrag_Cori_stored

         CALL icb_utl_interp( pi0, pj0,   &
         &                 ibathy=ibathy_x0y0                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi0, zyj_n,    &
         &                 ibathy=ibathy_x0y1                  &   ! bathy
         &                  )
         CALL icb_utl_interp( zxi_n, pj0,    &
         &                 ibathy=ibathy_x1y0                  &   ! bathy
         &                  )
         CALL icb_utl_interp( zxi_n, zyj_n,     &
         &                 ibathy=ibathy_x1y1                  &   ! bathy
         &                  )

         ! We calculate the vectors of the iceberg motion in the x and y directions
         ! as well as the angles of the iceberg trajectory along the topography
         ii0 = INT(pi0+0.5_wp) + (nn_hls-1)  ;  ij0 = INT(pj0+0.5_wp) + (nn_hls-1) ! initial gridpoint position (T-cell)
         ii  = INT(zxi_n +0.5_wp) + (nn_hls-1)  ;  ij  = INT(zyj_n +0.5_wp) + (nn_hls-1)
         distX1=0.5_wp*(e1t(ii,ij0)+e1t(ii0,ij0))*ABS(zxi_n-pi0)
         distX2=0.5_wp*(e1t(ii,ij)+e1t(ii0,ij))*ABS(zxi_n-pi0)
         distY1=0.5_wp*(e2t(ii0,ij)+e2t(ii0,ij0))*ABS(zyj_n-pj0)
         distY2=0.5_wp*(e2t(ii,ij)+e2t(ii,ij0))*ABS(zyj_n-pj0)
         TopoX1=ibathy_x1y0  -ibathy_x0y0
         TopoX2=ibathy_x1y1  -ibathy_x0y1
         TopoY1=ibathy_x0y1  -ibathy_x0y0
         TopoY2=ibathy_x1y1  -ibathy_x1y0

        IF ((distX1 .LT. 0.01_wp).OR.(distX2 .LT. 0.01_wp)) THEN
         angleTopoXlagr=0._wp
        ELSE
         angleTopoXlagr=atan((0.5_wp*(-TopoX1/distX1)+0.5_wp*(-TopoX2/distX2))) &
      &  /(3.14_wp/180._wp)
        END IF

        IF ((distY1 .LT. 0.01_wp).OR.(distY2 .LT. 0.01_wp)) THEN
         angleTopoYlagr=0._wp
        ELSE
         angleTopoYlagr=atan((0.5_wp*(-TopoY1/distY1)+0.5_wp*(-TopoY2/distY2))) &
      &  /(3.14_wp/180._wp)
        ENDIF

        dist = ((0.5_wp*distX1+0.5_wp*distX2)**2 + &
      & (0.5_wp*distY1+0.5_wp*distY2)**2)**0.5_wp

        IF (dist .LT. 0.01_wp) THEN
           angleTopoTraject=0._wp
        ELSE
           angleTopoTraject= atan(-(ibathy_x1y1-ibathy_x0y0)/dist) &
      &    /(3.14_wp/180._wp)
        ENDIF
        

        ! We save the topographic angles and the trajectory angles
        pt%angleTopolagr=angleTopoTraject

        pt%angleTopoXlagr=angleTopoXlagr
        pt%angleTopoYlagr=angleTopoYlagr


         pt%uvel = zuvel_n                        !** save in berg structure
         pt%vvel = zvvel_n
         pt%xi   = zxi_n
         pt%yj   = zyj_n
 

         berg => berg%next                         ! switch to the next berg
         !
      END DO                                  !==  end loop over all bergs  ==!
      !
   END SUBROUTINE icb_dyn


   SUBROUTINE icb_ground( kt, berg, pi, pi0, pu1, pu,   &
      &                         pj, pj0, pv1, pv, ld_bounced, &
      &                         ld_bouncedx, ld_bouncedy, ll_silted, &
      &                   zdrag_groundedx, zdrag_groundedy, &
      &                   zdrag_groundedxS, zdrag_groundedyS, &
      &                   zdrag_groundedxG,zdrag_groundedyG,zFdisequil, &
      &                   ibounce_method, muFric,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_ground  ***
      !!
      !! ** Purpose :   iceberg grounding.
      !!
      !! ** Method  : - adjust velocity and then put iceberg back to start position
      !!                NB three possibilities available one of which is hard-coded here
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt 
      TYPE(iceberg ), POINTER, INTENT(in   ) ::   berg             ! berg
      !
      REAL(wp), INTENT(inout) ::   pi , pj      ! current iceberg position
      REAL(wp), INTENT(in   ) ::   pi0, pj0     ! previous iceberg position
      REAL(wp), INTENT(inout) ::   pu1, pu , pv1, pv     ! current iceberg velocities
      LOGICAL , INTENT(  out) ::   ld_bounced, ll_silted   ! bounced indicator
      LOGICAL , INTENT(  out) ::   ld_bouncedx   ! bounced indicator
      LOGICAL , INTENT(  out) ::   ld_bouncedy   ! bounced indicator
      REAL(wp), INTENT(  out) ::   zdrag_groundedx, zdrag_groundedy     ! grounding drag
      REAL(wp), INTENT(  out) ::   zdrag_groundedxG, zdrag_groundedyG   ! grounding drag
      REAL(wp), INTENT(  out) ::   zdrag_groundedxS, zdrag_groundedyS   ! silt drag
      REAL(wp), INTENT(  out) ::   zFdisequil
      REAL(wp), INTENT(in   ) ::   muFric
      INTEGER,  INTENT(in   ) ::   ibounce_method
      REAL(wp), INTENT(in   ) ::   unitWeightSilt, shearStrengthSilt, depthSilt, SiltScourWidth
      REAL(wp) :: ibathy_x0y0, ibathy_x0y1, ibathy_x1y0, ibathy_x1y1
      REAL(wp) :: ibathy_x0y0pX, ibathy_x0y0mX
      REAL(wp) :: ibathy_x0y0pY, ibathy_x0y0mY
      REAL(wp) :: issh_x0y0, issh_x1y1 
      REAL(wp) :: distX1, distX2, distY1, distY2, dist
      REAL(wp) :: distXp, distYp, distXm, distYm
      REAL(wp) :: pe1_x0y0, pe1_x0y1, pe1_x1y0, pe1_x1y1
      REAL(wp) :: pe2_x0y0, pe2_x0y1, pe2_x1y0, pe2_x1y1
      REAL(wp) :: angleTopoMaxGradient, angleTopoTraject
      REAL(wp) :: angleTopoX, angleTopoY
      REAL(wp) :: angleTopoXlagr, angleTopoYlagr
      REAL(wp) :: gravity_termX, gravity_termY, zdrag_groundedS
      REAL(wp) :: gravity_termX2, gravity_termY2 
      REAL(wp) :: TopoX1, TopoX2, TopoY1, TopoY2
      REAL(wp) :: angleTopoXlagrOld, angleTopoYlagrOld, angleTopolagrOld
      TYPE(point)  , POINTER ::   pt
      INTEGER  ::   berg_number
      REAL(wp), PARAMETER ::   taperCoeffK  = 200000.0_wp

      INTEGER  ::   ii, ii0, ii0px, ii0mx
      INTEGER  ::   ij, ij0, ij0py, ij0my
      INTEGER  ::   ikb
      REAL(wp) ::   zT, zD, zF, zFequil, zDequil
      REAL(wp), DIMENSION(jpk) :: ze3t
      REAL(wp) ::   depthG

      !!----------------------------------------------------------------------
      !
      ld_bounced = .FALSE.
      ld_bouncedx = .FALSE.
      ld_bouncedy = .FALSE.
      ll_silted   = .FALSE.
      zdrag_groundedx=0._wp
      zdrag_groundedy=0._wp
      zdrag_groundedxG=0._wp
      zdrag_groundedyG=0._wp
      zdrag_groundedS=0._wp
      zdrag_groundedxS=0._wp
      zdrag_groundedyS=0._wp
      depthG = 0._wp

      pt => berg%current_point
      berg_number = berg%number(1)

      ! We read the saved topographic and trajectory angles
      angleTopoXlagrOld=pt%angleTopoXlagr
      angleTopoYlagrOld=pt%angleTopoYlagr
      angleTopolagrOld =pt%angleTopolagr

      !
      ii0 = INT(pi0+0.5_wp) + (nn_hls-1)  ;  ij0 = INT(pj0+0.5_wp) + (nn_hls-1)      ! initial gridpoint position (T-cell)
      ii  = INT(pi +0.5_wp) + (nn_hls-1)  ;  ij  = INT(pj +0.5_wp) + (nn_hls-1)      ! current     -         -
      !
      !  IF( ii == ii0  .AND.  ij == ij0  )   RETURN           ! berg remains in the same cell
      !
      ! map into current processor
      ii0 = mi1( ii0 )
      ij0 = mj1( ij0 )
      ii  = mi1( ii  )
      ij  = mj1( ij  )
      !
      ! assume icb is grounded if tmask(ii,ij,1) or tmask(ii,ij,ikb), depending of the option is not 0
      IF ( ln_M2016 .AND. ln_icb_grd ) THEN
         !
         ! draught (keel depth)
         zT = berg%current_point%thickness               ! total thickness
         zD = rho_berg_1_oce * zT                        ! draught (keel depth)
         zF = zT - zD                                    ! freeboard
         zFequil = zF                                    ! unless the iceberg is
                                                      ! grounded, the freeboard
                                                      ! is consistent with
                                                      ! an equilibrium of
                                                      ! vertical forces
         zDequil = zD
         zFdisequil = 0._wp

         !
         ! interpol needed data

         CALL icb_utl_interp( pi0, pj0,   &
         &                 issh=issh_x0y0                  &   ! ssh
         &                  )
         CALL icb_utl_interp( pi, pj,   &
         &                 issh=issh_x1y1                  &   ! ssh
         &                  )

         CALL icb_utl_interp( pi, pj, pe3t=ze3t )

         ! we calculate changes in the background bathymetry depth
         CALL icb_utl_interp( pi0, pj0,   &  
         &                 ibathy=ibathy_x0y0                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi0, pj,    &   
         &                 ibathy=ibathy_x0y1                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi, pj0,    &   
         &                 ibathy=ibathy_x1y0                  &   ! bathy
         &                  )                  
         CALL icb_utl_interp( pi, pj,     &
         &                 ibathy=ibathy_x1y1                  &   ! bathy
         &                  )
         ! we calculate the distance over which the depth varies
         distX1=0.5_wp*(e1t(ii,ij0)+e1t(ii0,ij0))*ABS(pi-pi0)
         distX2=0.5_wp*(e1t(ii,ij)+e1t(ii0,ij))*ABS(pi-pi0)
         distY1=0.5_wp*(e2t(ii0,ij)+e2t(ii0,ij0))*ABS(pj-pj0)
         distY2=0.5_wp*(e2t(ii,ij)+e2t(ii,ij0))*ABS(pj-pj0)
         TopoX1=ibathy_x1y0  -ibathy_x0y0
         TopoX2=ibathy_x1y1  -ibathy_x0y1
         TopoY1=ibathy_x0y1  -ibathy_x0y0
         TopoY2=ibathy_x1y1  -ibathy_x1y0

         
       

         CALL icb_utl_interp( pi0+0.1_wp, pj0, &
         &                 ibathy=ibathy_x0y0pX                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi0-0.1_wp, pj0, &
         &                 ibathy=ibathy_x0y0mX                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi0, pj0+0.1_wp, &
         &                 ibathy=ibathy_x0y0pY                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi0, pj0-0.1_wp, &
         &                 ibathy=ibathy_x0y0mY                  &   ! bathy
         &                  )

        ii0px= INT(pi0+0.5_wp+0.1_wp) + (nn_hls-1) ; ij0py= INT(pj0+0.5_wp+0.1_wp) +(nn_hls-1)
        ii0mx= INT(pi0+0.5_wp-0.1_wp) + (nn_hls-1) ; ij0my= INT(pj0+0.5_wp-0.1_wp) +(nn_hls-1)
         ii0px = mi1( ii0px )
         ij0py = mj1( ij0py )
         ii0mx = mi1( ii0mx )
         ij0my = mj1( ij0my )

         distXp=0.5_wp*(e1t(ii0,ij0)+e1t(ii0px,ij0))*0.1_wp
         distXm=0.5_wp*(e1t(ii0,ij0)+e1t(ii0mx,ij0))*0.1_wp
         distYp=0.5_wp*(e2t(ii0,ij0)+e2t(ii0,ij0py))*0.1_wp
         distYm=0.5_wp*(e2t(ii0,ij0)+e2t(ii0,ij0my))*0.1_wp

      ! The maximum topographic slope at the iceberg location
        angleTopoMaxGradient=atan(-(((0.5_wp*(ibathy_x0y0pX-ibathy_x0y0)/distXp+ &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mX)/distXm)**2._wp + &
      &                     (0.5_wp*(ibathy_x0y0pY-ibathy_x0y0)/distYp+ &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mY)/distYm)**2._wp)**0.5_wp)) &
      &                      /(3.14_wp/180._wp)


        dist = ((0.5_wp*distX1+0.5_wp*distX2)**2 + &
      & (0.5_wp*distY1+0.5_wp*distY2)**2)**0.5_wp

      ! The angle of the trajectory if the iceberg is sliding along topography   
        IF (dist .LT. 0.01_wp) THEN
           angleTopoTraject=0._wp
        ELSE
           angleTopoTraject= atan(-(ibathy_x1y1-ibathy_x0y0)/dist) &
      &    /(3.14_wp/180._wp)
        ENDIF

        IF ( isnan(angleTopolagrOld)) THEN
           angleTopolagrOld=angleTopoTraject
        ENDIF

         angleTopoX=atan(-(0.5_wp*(ibathy_x0y0pX-ibathy_x0y0)/distXp + &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mX)/distXm)) &
      &  /(3.14_wp/180._wp)
      ! The effect of gravity along the x-direction
         gravity_termX=SIN((3.14_wp/180._wp)*angleTopoX) &
      &               *COS((3.14_wp/180._wp)*angleTopoX)
      ! The impact of the vertical deflection of the trajectory onto horizontal velocity
        IF ((distX1 .LT. 0.01_wp).OR.(distX2 .LT. 0.01_wp)) THEN
         gravity_termX2=0._wp       
        ELSE 
         gravity_termX2=-pu*SQRT(pu**(2._wp)+pv**(2._wp))*ABS(tan(0.5_wp*(3.14_wp/180._wp)*(angleTopolagrOld+angleTopoTraject)))&
      &    *((-angleTopolagrOld+angleTopoTraject)*(3.14_wp/180._wp))&
      &    /((dist))
        END IF


         angleTopoY=atan(-(0.5_wp*(ibathy_x0y0pY-ibathy_x0y0)/distYp + &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mY)/distYm)) &
      &  /(3.14_wp/180._wp)
      ! The effect of gravity along the x-direction
         gravity_termY=SIN((3.14_wp/180._wp)*angleTopoY) &
      &               *COS((3.14_wp/180._wp)*angleTopoY)
     ! The impact of the vertical deflection of the trajectory onto horizontal velocity
        IF ((distY1 .LT. 0.01_wp).OR.(distY2 .LT. 0.01_wp)) THEN
         gravity_termY2=0._wp
        ELSE
         gravity_termY2=-pv*SQRT(pu**(2._wp)+pv**(2._wp))*ABS(tan(0.5_wp*(3.14_wp/180._wp)*(angleTopolagrOld+angleTopoTraject)))&
      &    *((-angleTopolagrOld+angleTopoTraject)*(3.14_wp/180._wp))&
      &    /((dist))
        END IF        


            IF (zDequil-ABS(ibathy_x1y1)-issh_x1y1 .GT. 0._wp) THEN
             ll_silted   = .TRUE.

      ! The resistance experienced by icebergs ploughing scours into the sediment layers:
        IF  (berg%current_point%mass .GT. 1.E-20_wp) THEN
             zdrag_groundedS = (0.5_wp*unitWeightSilt*SiltScourWidth &
            &                         *((MAX(0._wp,(MIN(depthSilt, &
            &                          zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
            &                          ))))**2._wp) &
            &                   +2._wp*shearStrengthSilt*SiltScourWidth &
            &                         *(MAX(0._wp,(MIN(depthSilt, &
            &                          zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
            &                          )))) &
            &                   +0.5_wp*SQRT(2._wp)*shearStrengthSilt &
            &                         *((MAX(0._wp,(MIN(depthSilt, &
            &                          zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
            &                          ))))**2._wp) &
            &                 )/berg%current_point%mass
        END IF

        ! Projection of the sediment resistance onto the x- and y-components of iceberg motion
             IF ((SQRT(pu**2 + pv**2).GT. 1.E-20_wp).AND.(ibounce_method==3)) THEN
            zdrag_groundedxS = &
             &  -zdrag_groundedS*(abs(pu)/SQRT(pu**2 + pv**2)) &
             & * tanh(taperCoeffK*pu)
            zdrag_groundedyS = &
             &  -zdrag_groundedS*(abs(pv)/SQRT(pu**2 + pv**2)) &
             & * tanh(taperCoeffK*pv)
             ENDIF

            ENDIF

         ! 
         !compute bottom level
         CALL icb_utl_getkb( ikb, ze3t, zD )
         !
         ! berg reach a new t-cell, but an ocean one
         ! .AND. needed in case berg hit an isf (tmask(ii,ij,1) == 0 and tmask(ii,ij,ikb) /= 0)
!         IF(  tmask(ii,ij,ikb) /= 0._wp .AND. tmask(ii,ij,1) /= 0._wp ) RETURN
!         !
!      ELSE
!         IF(  tmask(ii,ij,1)  /=   0._wp  )   RETURN           ! berg reach a new t-cell, but an ocean one
!      END IF
      IF (((zT-depthSilt-ABS(ibathy_x1y1)-issh_x1y1).LE.zFequil).AND.tmask(ii,ij,1)/=0._wp) RETURN
      
      IF((  tmask(ii,ij,1)  /=   0._wp  ).AND.ibounce_method/=3)   RETURN           ! berg reach a new t-cell, but an ocean one
      
      !
      ! From here, berg have reach land: treat grounding/bouncing
      ! -------------------------------
      ld_bounced = .TRUE.

      !! not obvious what should happen now
      !! if berg tries to enter a land box, the only location we can return it to is the start 
      !! position (pi0,pj0), since it has to be in a wet box to do any melting;
      !! first option is simply to set whole velocity to zero and move back to start point
      !! second option (suggested by gm) is only to set the velocity component in the (i,j) direction
      !! of travel to zero; at a coastal boundary this has the effect of sliding the berg along the coast

      SELECT CASE ( ibounce_method )
      CASE ( 1 )
         pi = pi0
         pj = pj0
         pu = 0._wp
         pv = 0._wp
      CASE ( 2 )
         IF( ii0 /= ii ) THEN
            pi = pi0                   ! return back to the initial position
            pu = 0._wp                 ! zeroing of velocity in the direction of the grounding
         ENDIF
         IF( ij0 /= ij ) THEN
            pj = pj0                   ! return back to the initial position
            pv = 0._wp                 ! zeroing of velocity in the direction of the grounding
         ENDIF
      CASE ( 3 )
        IF(tmask(ii,ij,1) .EQ. 0._wp) THEN
!         IF( ii0 /= ii ) THEN
            pi = pi0                   ! return back to the initial position
            pu = 0._wp                 ! zeroing of velocity in the direction of the grounding
!         ENDIF
!         IF( ij0 /= ij ) THEN
            pj = pj0                   ! return back to the initial position
            pv = 0._wp                 ! zeroing of velocity in the direction of the grounding
!         ENDIF
        ELSE
            ld_bouncedx = .TRUE.
            ld_bouncedy = .TRUE.
            ld_bounced = .TRUE.
            depthG = depthSilt+ibathy_x0y0+issh_x0y0
            zF = zT - depthG
            zFdisequil = MAX(0._wp, zF - zFequil)

         ! Solid-body Coulomb friction on icebergs whose keels reach the solid basement   
         IF (SQRT(pu**2 + pv**2).GT. 1.E-20_wp) THEN
            zdrag_groundedx = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*( muFric*(cos((3.14_wp/180._wp) &
            & * angleTopoTraject))*(cos((3.14_wp/180._wp) &
            & * angleTopoMaxGradient))*(abs(pu)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pu))
            zdrag_groundedy = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*( muFric*(cos((3.14_wp/180._wp) &
            & * angleTopoTraject))*(cos((3.14_wp/180._wp) &
            & * angleTopoMaxGradient))*(abs(pv)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pv))
         ENDIF
         
         ! Gravitational and geometric deceleration experienced by icebergs along the
         ! sloping solid basement
            zdrag_groundedxG = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*(gravity_termX) &
            &  +gravity_termX2
            zdrag_groundedyG = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*(gravity_termY) &
            &  +gravity_termY2

       END IF
      END SELECT
      END IF
      !
   END SUBROUTINE icb_ground


   SUBROUTINE icb_accel( kt, berg , pxi, pe1, puvel, puvel0, pax,                 &
      &                             pyj, pe2, pvvel, pvvel0, pay, pdt, pcfl_scale, &
      &                         ll_bouncedx, ll_bouncedy, ll_silted, &
      &                        zdrag_groundedxS,zdrag_groundedyS,unitWeightSilt,shearStrengthSilt,depthSilt,SiltScourWidth,muFric,&
      &                        zdrag_groundedx,zdrag_groundedy,FlagStopX,FlagStopY,&
      &       zdrag_groundedxG,zdrag_groundedyG,zFdisequil,zaxeBeforeFG, zayeBeforeFG,&
      &                        Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori,&
      &                        Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_accel  ***
      !!
      !! ** Purpose :   compute the iceberg acceleration.
      !!
      !! ** Method  : - sum the terms in the momentum budget
      !!----------------------------------------------------------------------
      TYPE(iceberg ), POINTER, INTENT(in   ) ::   berg             ! berg
      INTEGER                , INTENT(in   ) ::   kt               ! time step
      REAL(wp)               , INTENT(in   ) ::   pcfl_scale
      REAL(wp)               , INTENT(in   ) ::   pxi   , pyj      ! berg position in (i,j) referential
      REAL(wp)               , INTENT(in   ) ::   puvel , pvvel    ! berg velocity [m/s]
      REAL(wp)               , INTENT(in   ) ::   puvel0, pvvel0   ! initial berg velocity [m/s]
      REAL(wp)               , INTENT(  out) ::   pe1, pe2         ! horizontal scale factor at (xi,yj)
      REAL(wp)               , INTENT(inout) ::   pax, pay         ! berg acceleration
      REAL(wp)               , INTENT(in   ) ::   pdt              ! berg time step
      LOGICAL                , INTENT(in   ) ::   ll_bouncedx, ll_bouncedy, ll_silted ! whether we are grounded or not
      REAL(wp), INTENT(inout) ::   zdrag_groundedxS, zdrag_groundedyS   ! silt drag
      REAL(wp), INTENT(in   ) ::   muFric
      REAL(wp), INTENT(in   ) ::   unitWeightSilt, shearStrengthSilt, depthSilt, SiltScourWidth
      REAL(wp)               , INTENT(inout) ::   zdrag_groundedx, zdrag_groundedy
      REAL(wp)               , INTENT(  out) ::   zaxeBeforeFG, zayeBeforeFG
      INTEGER                , INTENT(out  ) ::   FlagStopX, FlagStopY
      REAL(wp)               , INTENT(in   ) ::   zdrag_groundedxG, zdrag_groundedyG
      REAL(wp)               , INTENT(in   ) ::   zFdisequil
      REAL(wp), INTENT(out  ) ::   Xdrag_ocn, Xdrag_atm, Xdrag_ice, Xgrad_press_sc, Xdrag_Cori
      REAL(wp), INTENT(out  ) ::   Ydrag_ocn, Ydrag_atm, Ydrag_ice, Ygrad_press_sc, Ydrag_Cori
      !
      REAL(wp), PARAMETER ::   pp_alpha     = 0._wp      !
      REAL(wp), PARAMETER ::   pp_beta      = 1._wp      !
      REAL(wp), PARAMETER ::   pp_vel_lim   =15._wp      ! max allowed berg speed
      REAL(wp), PARAMETER ::   pp_accel_lim = 1.e-2_wp   ! max allowed berg acceleration
      REAL(wp), PARAMETER ::   pp_Cr0       = 0.06_wp    !
      REAL(wp), PARAMETER ::   taperCoeffK  = 200000.0_wp
      REAL(wp), PARAMETER ::   Static2KineticRatio  = 100.0_wp

      !
      INTEGER  ::   itloop, ikb, jk
      REAL(wp) ::   zuo, zssu, zui, zua, zuwave, zssh_x, zcn, zhi
      REAL(wp) ::   zvo, zssv, zvi, zva, zvwave, zssh_y
      REAL(wp) ::   zff, zT, zD, zW, zL, zM, zF
      REAL(wp) ::   zaxe_ocn , zaye_ocn, zdrag_ocn, zdrag_atm, zdrag_ice, zwave_rad
      REAL(wp) ::   z_ocn, z_atm, z_ice, zdep
      REAL(wp) ::   zampl, zwmod, zCr, zLwavelength, zLcutoff, zLtop
      REAL(wp) ::   zlambda, zdetA, zA11, zA12, zaxe, zaye, zD_hi
      REAL(wp) ::   zuveln, zvveln, zus, zvs, zspeed, zloc_dx, zspeed_new
      REAL(wp) ::   press_grad_term_x, press_grad_term_y
      REAL(wp), DIMENSION(jpk) :: zaxe_ocn_profile, zaye_ocn_profile, zdrag_ocn_profile, zhpi_profile, zhpj_profile
      INTEGER  ::   DoNotDebugExplicitForcing
      REAL(wp), DIMENSION(jpk) :: zuoce, zvoce, ze3t, zdepw
      INTEGER  ::   berg_number
      REAL(wp) ::   pi0, pj0     ! previous iceberg position
      INTEGER  ::   ii, ii0, ii0px, ii0mx
      INTEGER  ::   ij, ij0, ij0py, ij0my
      REAL(wp) ::   pi, pj
      REAL(wp) :: pu, pv
      REAL(wp) :: distX1, distX2, distY1, distY2, dist
      REAL(wp) :: distXp, distYp, distXm, distYm
      REAL(wp) :: angleTopoMaxGradient, angleTopoTraject
      REAL(wp) :: ibathy_x0y0, ibathy_x1y1, zDequil
      REAL(wp) :: ibathy_x0y0pX, ibathy_x0y0mX
      REAL(wp) :: ibathy_x0y0pY, ibathy_x0y0mY
      REAL(wp) :: zdrag_groundedS, issh_x1y1

      !!----------------------------------------------------------------------
      

      ! Interpolate gridded fields to berg
      nknberg = berg%number(1)
      berg_number = berg%number(1)
      zdrag_groundedS = 0._wp
      FlagStopX=0
      FlagStopY=0
      zaxeBeforeFG=0._wp 
      zayeBeforeFG=0._wp
      CALL icb_utl_interp( pxi, pyj, pe1=pe1, pe2=pe2,     &   ! scale factor
         &                 pssu=zssu, pui=zui, pua=zua,    &   ! oce/ice/atm velocities
         &                 pssv=zssv, pvi=zvi, pva=zva,    &   ! oce/ice/atm velocities
         &                 pssh_i=zssh_x, pssh_j=zssh_y,   &   ! ssh gradient
         &                 phi=zhi, pff=zff)                   ! ice thickness and coriolis

      zM = berg%current_point%mass
      zT = berg%current_point%thickness               ! total thickness
      zD = rho_berg_1_oce * zT                        ! draught (keel depth) at equilibrium (free flotation)
      zDequil = zD
      zF = zT - zD                                    ! freeboard at equilibrium (free flotation)
      zW = berg%current_point%width
      zL = berg%current_point%length
      pi0 = berg%current_point%xi
      pj0 = berg%current_point%yj

      IF(ll_bouncedx .OR. ll_bouncedy) THEN
         zF = zF + zFdisequil
         zD = zT - zF
      ENDIF


      zhi   = MIN( zhi   , zD    )
      zD_hi = MAX( 0._wp, zD-zhi )
 
     ! Wave radiation
      zuwave = zua - zssu   ;   zvwave = zva - zssv   ! Use wind speed rel. to ocean for wave model
      zwmod  = zuwave*zuwave + zvwave*zvwave          ! The wave amplitude and length depend on the  current;
      !                                               ! wind speed relative to the ocean. Actually wmod is wmod**2 here.
      zampl        = 0.5_wp * 0.02025_wp * zwmod      ! This is "a", the wave amplitude
      zLwavelength =       0.32_wp    * zwmod         ! Surface wave length fitted to data in table at
      !                                               ! http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
      zLcutoff     = 0.125_wp * zLwavelength
      zLtop        = 0.25_wp  * zLwavelength
      zCr          = pp_Cr0 * MIN(  MAX( 0._wp, (zL-zLcutoff) / ((zLtop-zLcutoff)+1.e-30)) , 1._wp)  ! Wave radiation coefficient
      !                                               ! fitted to graph from Carrieres et al.,  POAC Drift Model.
      zwave_rad    = 0.5_wp * pp_rho_seawater / zM * zCr * grav * zampl * MIN( zampl,zF ) * (2._wp*zW*zL) / (zW+zL)
      zwmod        = SQRT( zua*zua + zva*zva )        ! Wind speed
      IF( zwmod /= 0._wp ) THEN
         zuwave = zua/zwmod   ! Wave radiation force acts in wind direction ...       !!gm  this should be the wind rel. to ocean ?
         zvwave = zva/zwmod
      ELSE
         zuwave = 0._wp   ;    zvwave=0._wp   ;    zwave_rad=0._wp ! ... and only when wind is present.     !!gm  wave_rad=0. is useless
      ENDIF

      ! Weighted drag coefficients
      z_ocn = pp_rho_seawater / zM * (0.5_wp*pp_Cd_wv*zW*(zD_hi)+pp_Cd_wh*zW*zL)
      z_atm = pp_rho_air      / zM * (0.5_wp*pp_Cd_av*zW*zF     +pp_Cd_ah*zW*zL)
      z_ice = pp_rho_ice      / zM * (0.5_wp*pp_Cd_iv*zW*zhi              )
      IF( abs(zui) + abs(zvi) == 0._wp )   z_ice = 0._wp

      ! lateral velocities
      ! default ssu and ssv
      ! ln_M2016: mean velocity along the profile
      IF ( ln_M2016 ) THEN
         ! interpol needed data
          CALL icb_utl_interp( pxi, pyj, puoce=zuoce, pvoce=zvoce, pe3t=ze3t )   ! 3d velocities
          CALL icb_utl_interp( pxi, pyj, izhpi=zhpi_profile, izhpj=zhpj_profile, pe3t=ze3t )        
         !compute bottom level
         CALL icb_utl_getkb( ikb, ze3t, zD )
         
         ! compute mean velocity 
         CALL icb_utl_zavg(zuo, zuoce, ze3t, zD, ikb)
         CALL icb_utl_zavg(zvo, zvoce, ze3t, zD, ikb)

         ! compute mean horizontal pressure gradients
         CALL icb_utl_zavg(press_grad_term_x, zhpi_profile, ze3t, zD, ikb)
         CALL icb_utl_zavg(press_grad_term_y, zhpj_profile, ze3t, zD, ikb)

         ! The pressure gradient term includes a component due to density gradients
         ! and a component due to SSH gradients
         press_grad_term_x= press_grad_term_x -grav * zssh_x
         press_grad_term_y= press_grad_term_y -grav * zssh_y

      ELSE
         zuo = zssu
         zvo = zssv
         ! If the Merino 2016 option is disabled, we approximate the pressure
         ! gradient using SSH gradients
         press_grad_term_x= -grav * zssh_x
         press_grad_term_y= -grav * zssh_y
      END IF

      zuveln = puvel   ;   zvveln = pvvel ! Copy starting uvel, vvel
      !
      DO itloop = 1, 2  ! Iterate on drag coefficients
         !
         zus = 0.5_wp * ( zuveln + puvel )
         zvs = 0.5_wp * ( zvveln + pvvel )
         zdrag_ocn = z_ocn * SQRT( (zus-zuo)*(zus-zuo) + (zvs-zvo)*(zvs-zvo) )
         zdrag_atm = z_atm * SQRT( (zus-zua)*(zus-zua) + (zvs-zva)*(zvs-zva) )
         zdrag_ice = z_ice * SQRT( (zus-zui)*(zus-zui) + (zvs-zvi)*(zvs-zvi) )
!        zdrag_atm = 0._wp

       IF ( ln_M2016 ) THEN

         DO jk=1,jpk
         zdrag_ocn_profile(jk) =  z_ocn * (( (zus-zuoce(jk))*(zus-zuoce(jk)) &
                 &                            + (zvs-zvoce(jk))*(zvs-zvoce(jk)))**0.5_wp)
         END DO
       END IF
         !
         ! Explicit accelerations

         ! The explicit acceleration includes the impact of horizontal pressure
         ! gradients rescaled by the density of water over the density of
         ! icebergs
         zaxe = press_grad_term_x *(pp_rho_seawater/rn_rho_bergs)*(zD/zT)
         zaye = press_grad_term_y *(pp_rho_seawater/rn_rho_bergs)*(zD/zT)
         ! Another explicit term is the accel. due to wave radiation
         zwave_rad=0._wp
         zaxe = zaxe  + zwave_rad * zuwave
         zaye = zaye  + zwave_rad * zvwave


         IF( pp_alpha > 0._wp ) THEN   ! If implicit, use time-level (n) rather than RK4 latest
            zaxe = zaxe + zff*pvvel0
            zaye = zaye - zff*puvel0
         ELSE
            zaxe = zaxe + zff*pvvel
            zaye = zaye - zff*puvel
         ENDIF
         IF( pp_beta > 0._wp ) THEN    ! If implicit, use time-level (n) rather than RK4 latest
            zaxe = zaxe - zdrag_atm*(puvel0-zua) -zdrag_ice*(puvel0-zui)
            zaye = zaye - zdrag_atm*(pvvel0-zva) -zdrag_ice*(pvvel0-zvi)
           IF ( ln_M2016 ) THEN
              ! We compute the drag due to the relative velocity at each level     
                      DO jk=1,jpk
                        zaxe_ocn_profile(jk) = -(zdrag_ocn_profile(jk))*(puvel0-zuoce(jk))
                        zaye_ocn_profile(jk) = -(zdrag_ocn_profile(jk))*(pvvel0-zvoce(jk))
                      END DO
              ! We average the ocean drag vertically over the keel depth        
             CALL icb_utl_zavg(zaxe_ocn, zaxe_ocn_profile, ze3t, zD, ikb)
             CALL icb_utl_zavg(zaye_ocn, zaye_ocn_profile, ze3t, zD, ikb)
             zaxe = zaxe + zaxe_ocn
             zaye = zaye + zaye_ocn
           ELSE
            zaxe = zaxe - zdrag_ocn*(puvel0-zuo)
            zaye = zaye - zdrag_ocn*(pvvel0-zvo)
           ENDIF
         ELSE     
            zaxe = zaxe - zdrag_atm*(puvel -zua) -zdrag_ice*(puvel -zui)
            zaye = zaye - zdrag_atm*(pvvel -zva) -zdrag_ice*(pvvel -zvi)
           IF ( ln_M2016 ) THEN
              ! We compute the drag due to the relative velocity at each level                   
                      DO jk=1,jpk
                        zaxe_ocn_profile(jk) = -(zdrag_ocn_profile(jk))*(puvel-zuoce(jk))
                        zaye_ocn_profile(jk) = -(zdrag_ocn_profile(jk))*(pvvel-zvoce(jk))
                      END DO
              ! We average the ocean drag vertically over the keel depth
             CALL icb_utl_zavg(zaxe_ocn, zaxe_ocn_profile, ze3t, zD, ikb)
             CALL icb_utl_zavg(zaye_ocn, zaye_ocn_profile, ze3t, zD, ikb)
             zaxe = zaxe + zaxe_ocn
             zaye = zaye + zaye_ocn
           ELSE
            zaxe = zaxe - zdrag_ocn*(puvel -zuo)
            zaye = zaye - zdrag_ocn*(pvvel -zvo)
           ENDIF
         ENDIF

      ! Output the acceleration components   
      zaxeBeforeFG=zaxe
      zayeBeforeFG=zaye
      Xdrag_ocn = zaxe_ocn !-zdrag_ocn*(puvel0-zuo)
      Xdrag_atm = -zdrag_atm*(puvel0-zua)
      Xdrag_ice = -zdrag_ice*(puvel0-zui)
      Xgrad_press_sc = press_grad_term_x*(pp_rho_seawater/rn_rho_bergs)*(zD/zT)
      Xdrag_Cori = +zff*pvvel
      Ydrag_ocn = zaye_ocn !-zdrag_ocn*(pvvel0-zvo)
      Ydrag_atm = -zdrag_atm*(pvvel0-zva)
      Ydrag_ice = -zdrag_ice*(pvvel0-zvi)
      Ygrad_press_sc = press_grad_term_y*(pp_rho_seawater/rn_rho_bergs)*(zD/zT)
      Ydrag_Cori = -zff*puvel

      ! Add gravitational acceleraton in the case of grounding along the sloping
      ! solid basement
      IF((ll_bouncedx).OR.(ll_bouncedx)) THEN
       zaxe = zaxe + zdrag_groundedxG
       zaye = zaye + zdrag_groundedyG
      ENDIF


      IF((ll_bouncedx).OR.(ll_bouncedx).OR.(ll_silted)) THEN
      ! What x,y location does the acceleration (excluding dissipative forces) 
      ! push the iceberg towards?      
        pi = pxi + 0.5*(zaxe/pe1)*((pdt)**2)
        pj = pyj + 0.5*(zaye/pe2)*((pdt)**2)
        ii  = INT(pi +0.5_wp) + (nn_hls-1)  ;  ij  = INT(pj +0.5_wp) + (nn_hls-1)
        ii = mi1( ii )
        ij = mj1( ij )


      ! The case of static icebergs when grounding and silt resistance (if present)
      ! have to be oriented against the acceleration, not agains the velocity vectors
       IF ((SQRT(puvel0**2 + pvvel0**2).LE.1.E-20_wp) .OR. &
           & (SQRT(puvel**2 + pvvel**2).LE.1.E-20_wp)) THEN

        ! What velocity would arise from the acceleration excluding dissipative forces   
        pu = puvel0 + zaxe*pdt
        pv = pvvel0 +  zaye*pdt


         CALL icb_utl_interp( pxi, pyj,   &
         &                 ibathy=ibathy_x0y0                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pi, pj,     &
         &                 ibathy=ibathy_x1y1                  &   ! bathy
         &                  )

         CALL icb_utl_interp( pxi+0.1_wp, pyj, &
         &                 ibathy=ibathy_x0y0pX                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pxi-0.1_wp, pyj, &
         &                 ibathy=ibathy_x0y0mX                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pxi, pyj+0.1_wp, &
         &                 ibathy=ibathy_x0y0pY                  &   ! bathy
         &                  )
         CALL icb_utl_interp( pxi, pyj-0.1_wp, &
         &                 ibathy=ibathy_x0y0mY                  &   ! bathy
         &                  )

         CALL icb_utl_interp( pxi, pyj,   &
         &                 issh=issh_x1y1                  &   ! ssh
         &                  )

      ! Compute the silt resistance in case of static grounding
      ! to orient it against acceleration vectors
      IF  (berg%current_point%mass .GT. 0._wp) THEN
            zdrag_groundedS = (0.5_wp*unitWeightSilt*SiltScourWidth &
      &                         *((MAX(0._wp,(MIN(depthSilt, &
      &                          zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
      &                          ))))**2._wp) &
      &                       +2._wp*shearStrengthSilt*SiltScourWidth &
      &                               *(MAX(0._wp,(MIN(depthSilt, &
      &                                zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
      &                                )))) &
      &                       +0.5_wp*SQRT(2._wp)*shearStrengthSilt &
      &                               *((MAX(0._wp,(MIN(depthSilt, &
      &                                zDequil-ABS(ibathy_x1y1)-issh_x1y1 &
      &                                ))))**2._wp) &
      &                      )/berg%current_point%mass
      ENDIF

         ii  = INT(pi +0.5_wp) + (nn_hls-1)  ;  ij  = INT(pj +0.5_wp) + (nn_hls-1) ! current
         ii0 = INT(pxi+0.5_wp) + (nn_hls-1)  ;  ij0 = INT(pyj+0.5_wp) +(nn_hls-1)
         ii0px= INT(pxi+0.5_wp+0.1_wp) + (nn_hls-1) ; ij0py=INT(pyj+0.5_wp+0.1_wp) +(nn_hls-1)
         ii0mx= INT(pxi+0.5_wp-0.1_wp) + (nn_hls-1) ; ij0my=INT(pyj+0.5_wp-0.1_wp) +(nn_hls-1)
         ii0 = mi1( ii0 )
         ij0 = mj1( ij0 )
         ii0px = mi1( ii0px )
         ij0py = mj1( ij0py )
         ii0mx = mi1( ii0mx )
         ij0my = mj1( ij0my )

         distX1=0.5_wp*(e1t(ii,ij0)+e1t(ii0,ij0))*ABS(pi-pxi)
         distX2=0.5_wp*(e1t(ii,ij)+e1t(ii0,ij))*ABS(pi-pxi)
         distY1=0.5_wp*(e2t(ii0,ij)+e2t(ii0,ij0))*ABS(pj-pyj)
         distY2=0.5_wp*(e2t(ii,ij)+e2t(ii,ij0))*ABS(pj-pyj)
         distXp=0.5_wp*(e1t(ii0,ij0)+e1t(ii0px,ij0))*0.1_wp
         distXm=0.5_wp*(e1t(ii0,ij0)+e1t(ii0mx,ij0))*0.1_wp
         distYp=0.5_wp*(e2t(ii0,ij0)+e2t(ii0,ij0py))*0.1_wp
         distYm=0.5_wp*(e2t(ii0,ij0)+e2t(ii0,ij0my))*0.1_wp

      ! The maximum topographic slope at the iceberg location
        angleTopoMaxGradient=atan(-(((0.5_wp*(ibathy_x0y0pX-ibathy_x0y0)/distXp+ &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mX)/distXm)**2._wp + &
      &                     (0.5_wp*(ibathy_x0y0pY-ibathy_x0y0)/distYp+ &
      &                      0.5_wp*(ibathy_x0y0-ibathy_x0y0mY)/distYm)**2._wp)**0.5_wp)) &
      &                      /(3.14_wp/180._wp)


         dist = ((0.5_wp*distX1+0.5_wp*distX2)**2 + &
      &  (0.5_wp*distY1+0.5_wp*distY2)**2)**0.5_wp

         IF (dist .LT. 0.01_wp) THEN
            angleTopoTraject=0._wp
         ELSE
           angleTopoTraject= atan(-(ibathy_x1y1-ibathy_x0y0)/dist) &
      &    /(3.14_wp/180._wp)
         ENDIF

       ! Static friction only if the drive for future motion is not negligible
       ! and if there is a freeboard disequilibrium
       ! pu and pv are the projected velocities       
         IF (SQRT(pu**2 + pv**2) .GT. 1.E-20_wp) THEN
            zdrag_groundedx = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*( Static2KineticRatio*muFric*(cos((3.14_wp/180._wp) &
            & * angleTopoTraject))*(cos((3.14_wp/180._wp) &
            & * angleTopoMaxGradient))*(abs(pu)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pu))
            zdrag_groundedy = &
            & -(pp_rho_seawater/rn_rho_bergs) &
            &  *(zFdisequil/zT) &
            &  *9.81_wp*( Static2KineticRatio*muFric*(cos((3.14_wp/180._wp) &
            & * angleTopoTraject))*(cos((3.14_wp/180._wp) &
            & * angleTopoMaxGradient))*(abs(pv)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pv))
        
        ! Static grounding: orient sediment resistance against the acceleration vectors
        ! pu and pv are in units of velocity but reflect the acceleration vectors
            zdrag_groundedxS = &
            &  -zdrag_groundedS*(abs(pu)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pu)
            zdrag_groundedyS = &
            &  -zdrag_groundedS*(abs(pv)/SQRT(pu**2 + pv**2)) &
            & * tanh(taperCoeffK*pv)
         ELSE
            FlagStopX=1
            FlagStopY=1
         ENDIF

         ! Note that in the comparisons below zaxe, zaye already include gravity but 
         ! none of the dissipative forces
              IF (ABS(zaxe) .LE. (ABS(zdrag_groundedx)+ABS(zdrag_groundedxS))) THEN
               zaxe =-puvel0/pdt
               FlagStopX=1 ! The iceberg stops or remains at rest
              ELSE
               zaxe = zaxe + zdrag_groundedx + zdrag_groundedxS
              ENDIF

              IF (ABS(zaye) .LE. (ABS(zdrag_groundedy)+ABS(zdrag_groundedyS))) THEN
               zaye =-pvvel0/pdt
               FlagStopY=1 ! The iceberg stops or remains at rest
              ELSE
              zaye = zaye + zdrag_groundedy + zdrag_groundedyS
              ENDIF


       ! Else, if the iceberg is not static:
       ELSE

              IF (ABS(puvel0+pdt*zaxe) .LE. (ABS(pdt*zdrag_groundedx)+ABS(pdt*zdrag_groundedxS))) THEN
               zaxe =-puvel0/pdt
               zuveln=0._wp
               FlagStopX=1 ! The iceberg stops or remains at rest
              ELSE
               zaxe = zaxe + zdrag_groundedx  + zdrag_groundedxS
              ENDIF

              IF (ABS(pvvel0+pdt*zaye) .LE. (ABS(pdt*zdrag_groundedy)+ABS(pdt*zdrag_groundedyS))) THEN
               zaye =-pvvel0/pdt
               zvveln=0._wp
               FlagStopY=1 ! The iceberg stops or remains at rest
              ELSE
              zaye = zaye + zdrag_groundedy  + zdrag_groundedyS
              ENDIF

       ENDIF

      ! If the iceberg encounters an ice-shelf, it stops.
      IF(tmask(ii,ij,1) .EQ. 0) THEN
               zaxe =-puvel0/pdt
               pax  =-puvel0/pdt
               zuveln=0._wp
               FlagStopX=1
               zaye =-pvvel0/pdt
               pay =-pvvel0/pdt
               zvveln=0._wp
               FlagStopY=1
      ENDIF

      ENDIF


!         ! Solve for implicit accelerations
!         IF( pp_alpha + pp_beta > 0._wp ) THEN
!            zlambda = zdrag_ocn + zdrag_atm + zdrag_ice !&
!!       &              + SQRT(zdrag_groundedX**2 + zdrag_groundedy**2)
!            zA11    = 1._wp + pp_beta *pdt*zlambda
!            zA12    =         pp_alpha*pdt*zff
!            zdetA   = 1._wp / ( zA11*zA11 + zA12*zA12 )
!            pax     = zdetA * ( zA11*zaxe + zA12*zaye )
!            pay     = zdetA * ( zA11*zaye - zA12*zaxe )
!         ELSE
!            pax = zaxe   ;   pay = zaye
!         ENDIF


         pax = zaxe   ;   pay = zaye
         zuveln = puvel0 + pdt*pax
         zvveln = pvvel0 + pdt*pay
         !
      END DO      ! itloop


      IF( rn_speed_limit > 0._wp ) THEN       ! Limit speed of bergs based on a CFL criteria (if asked)
         zspeed = SQRT( zuveln*zuveln + zvveln*zvveln )    ! Speed of berg
         IF( zspeed > 0._wp ) THEN
            zloc_dx = MIN( pe1, pe2 )                                ! minimum grid spacing
            ! cfl scale is function of the RK4 step
            zspeed_new = zloc_dx / pdt * rn_speed_limit * pcfl_scale ! Speed limit as a factor of dx / dt
            IF( zspeed_new < zspeed ) THEN
               zuveln = zuveln * ( zspeed_new / zspeed )             ! Scale velocity to reduce speed
               zvveln = zvveln * ( zspeed_new / zspeed )             ! without changing the direction
               pax = (zuveln - puvel0)/pdt
               pay = (zvveln - pvvel0)/pdt
               !
               ! print speeding ticket
               IF (nn_verbose_level > 0) THEN
                  !         ! WRITE(numicb, 9200) 'icb speeding : ',kt, nknberg, zspeed, &
                  !     &                pxi, pyj, zuo, zvo, zua, zva, zui, zvi
                !  9200 FORMAT(a,i9,i6,f9.2,1x,4(1x,2f9.2))
               END IF
               !
               CALL icb_dia_speed()
            ENDIF
         ENDIF
      ENDIF

         
      !                                      ! check the speed and acceleration limits
      IF (nn_verbose_level > 0) THEN
       !  IF( ABS( zuveln ) > pp_vel_lim   .OR. ABS( zvveln ) > pp_vel_lim   )   &
            !         ! WRITE(numicb,'("pe=",i3,x,a)') narea,'Dump triggered by excessive velocity'
       !  IF( ABS( pax    ) > pp_accel_lim .OR. ABS( pay    ) > pp_accel_lim )   &
            !         ! WRITE(numicb,'("pe=",i3,x,a)') narea,'Dump triggered by excessive acceleration'
      ENDIF
      !
      


   END SUBROUTINE icb_accel

   !!======================================================================
END MODULE icbdyn
