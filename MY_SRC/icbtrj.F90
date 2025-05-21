MODULE icbtrj
   !!======================================================================
   !!                       ***  MODULE  icbtrj  ***
   !! Ocean physics:  trajectory I/O routines
   !!======================================================================
   !! History :  3.3  !  2010-01  (Martin&Adcroft) Original code
   !!             -   !  2011-03  (Madec)          Part conversion to NEMO form
   !!             -   !                            Removal of mapping from another grid
   !!             -   !  2011-05  (Alderson)       New module to handle trajectory output
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   icb_trj_init  : initialise iceberg trajectory output files
   !!   icb_trj_write :
   !!   icb_trj_sync  :
   !!   icb_trj_end   :
   !!----------------------------------------------------------------------
   USE par_oce        ! NEMO parameters
   USE dom_oce        ! NEMO ocean domain
   USE phycst         ! NEMO physical constants
   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines
   !
   USE lib_mpp        ! NEMO MPI library, lk_mpp in particular
   USE in_out_manager ! NEMO IO, numout in particular
   USE ioipsl  , ONLY : ju2ymds    ! for calendar
   USE netcdf

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_trj_init    ! routine called in icbini.F90 module
   PUBLIC   icb_trj_write   ! routine called in icbstp.F90 module
   PUBLIC   icb_trj_sync    ! routine called in icbstp.F90 module
   PUBLIC   icb_trj_end     ! routine called in icbstp.F90 module

   INTEGER ::   num_traj = 0
   INTEGER ::   n_dim, m_dim
   INTEGER ::   ntrajid
   INTEGER ::   numberid, nstepid, nscaling_id
   INTEGER ::   nlonid, nlatid, nxid, nyid, nuvelid, nvvelid, nmassid
   INTEGER ::   nssuid, nssvid, nuaid, nvaid, nuiid, nviid
   INTEGER ::   nsshxid, nsshyid, nsstid, ncntid, nthkid
   INTEGER ::   nthicknessid, nwidthid, nlengthid
   INTEGER ::   nyearid, ndayid
   INTEGER ::   nmass_of_bits_id, nheat_density_id
   INTEGER ::   zdrag_groundedx_storedid, zdrag_groundedy_storedid
   INTEGER ::   zdrag_groundedxG_storedid, zdrag_groundedyG_storedid
   INTEGER ::   zdrag_groundedxS_storedid, zdrag_groundedyS_storedid
   INTEGER ::   Xdrag_ocn_storedid, Xdrag_atm_storedid, Xdrag_ice_storedid, Xgrad_press_sc_storedid, Xdrag_Cori_storedid
   INTEGER ::   Ydrag_ocn_storedid, Ydrag_atm_storedid, Ydrag_ice_storedid, Ygrad_press_sc_storedid, Ydrag_Cori_storedid
   INTEGER ::   zaxe_storedid, zaye_storedid
   INTEGER ::   zaxeBeforeFG_storedid, zayeBeforeFG_storedid
   INTEGER ::   zFequilid, zFdisequilid
   INTEGER ::   ibathy_x1y1id, issh_x1y1id


   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: icbtrj.F90 14030 2020-12-03 09:26:33Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_trj_init( ktend )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_trj_init  ***
      !!
      !! ** Purpose :   initialise iceberg trajectory output files
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   ktend   ! time step index
      !
      INTEGER                ::   iret, iyear, imonth, iday
      INTEGER                ::   idg  ! number of digits
      REAL(wp)               ::   zfjulday, zsec
      CHARACTER(len=80)      ::   cl_filename
      CHARACTER(LEN=12)      ::   clfmt            ! writing format
      CHARACTER(LEN=8 )      ::   cldate_ini, cldate_end
      TYPE(iceberg), POINTER ::   this
      TYPE(point)  , POINTER ::   pt
      !!----------------------------------------------------------------------

      ! compute initial time step date
      CALL ju2ymds( fjulday, iyear, imonth, iday, zsec )
      WRITE(cldate_ini, '(i4.4,2i2.2)') iyear, imonth, iday

      ! compute end time step date
      zfjulday = fjulday + rn_Dt / rday * REAL( nitend - nit000 + 1 , wp)
      IF( ABS(zfjulday - REAL(NINT(zfjulday),wp)) < 0.1 / rday )   zfjulday = REAL(NINT(zfjulday),wp)   ! avoid truncation error
      CALL ju2ymds( zfjulday, iyear, imonth, iday, zsec )
      WRITE(cldate_end, '(i4.4,2i2.2)') iyear, imonth, iday

      ! define trajectory output name
      cl_filename = 'trajectory_icebergs_'//cldate_ini//'-'//cldate_end
      IF ( lk_mpp ) THEN
         idg = MAX( INT(LOG10(REAL(MAX(1,jpnij-1),wp))) + 1, 4 )          ! how many digits to we need to write? min=4, max=9
         WRITE(clfmt, "('(a,a,i', i1, '.', i1, ',a)')") idg, idg          ! '(a,a,ix.x,a)'
         WRITE(cl_filename,  clfmt) TRIM(cl_filename), '_', narea-1, '.nc'
      ELSE
         WRITE(cl_filename,'(a,a)') TRIM(cl_filename),               '.nc'
      ENDIF
      IF( lwp .AND. nn_verbose_level >= 0 )   WRITE(numout,'(2a)') 'icebergs, icb_trj_init: creating ',TRIM(cl_filename)

      iret = NF90_CREATE( TRIM(cl_filename), NF90_CLOBBER, ntrajid )
      IF (iret .NE. NF90_NOERR)   CALL ctl_stop('icebergs, icb_trj_init: nf_create failed')

      ! Dimensions
      iret = NF90_DEF_DIM( ntrajid, 'n', NF90_UNLIMITED, n_dim )
      IF ( iret /= NF90_NOERR )   CALL ctl_stop('icebergs, icb_trj_init: nf_def_dim n failed')
      iret = NF90_DEF_DIM( ntrajid, 'k', nkounts, m_dim )
      IF ( iret /= NF90_NOERR )   CALL ctl_stop('icebergs, icb_trj_init: nf_def_dim k failed')

      ! Variables
      iret = NF90_DEF_VAR( ntrajid, 'iceberg_number', NF90_INT   , (/m_dim,n_dim/), numberid         )
      iret = NF90_DEF_VAR( ntrajid, 'timestep'      , NF90_INT   , n_dim          , nstepid          )
      iret = NF90_DEF_VAR( ntrajid, 'mass_scaling'  , NF90_DOUBLE, n_dim          , nscaling_id      )
      iret = NF90_DEF_VAR( ntrajid, 'lon'           , NF90_DOUBLE, n_dim          , nlonid           )
      iret = NF90_DEF_VAR( ntrajid, 'lat'           , NF90_DOUBLE, n_dim          , nlatid           )
      iret = NF90_DEF_VAR( ntrajid, 'xi'            , NF90_DOUBLE, n_dim          , nxid             )
      iret = NF90_DEF_VAR( ntrajid, 'yj'            , NF90_DOUBLE, n_dim          , nyid             )
      iret = NF90_DEF_VAR( ntrajid, 'uvel'          , NF90_DOUBLE, n_dim          , nuvelid          )
      iret = NF90_DEF_VAR( ntrajid, 'vvel'          , NF90_DOUBLE, n_dim          , nvvelid          )
      iret = NF90_DEF_VAR( ntrajid, 'ssu'           , NF90_DOUBLE, n_dim          , nssuid           )
      iret = NF90_DEF_VAR( ntrajid, 'ssv'           , NF90_DOUBLE, n_dim          , nssvid           )
      iret = NF90_DEF_VAR( ntrajid, 'uta'           , NF90_DOUBLE, n_dim          , nuaid            )
      iret = NF90_DEF_VAR( ntrajid, 'vta'           , NF90_DOUBLE, n_dim          , nvaid            )
      iret = NF90_DEF_VAR( ntrajid, 'uti'           , NF90_DOUBLE, n_dim          , nuiid            )
      iret = NF90_DEF_VAR( ntrajid, 'vti'           , NF90_DOUBLE, n_dim          , nviid            )
      iret = NF90_DEF_VAR( ntrajid, 'ssh_x'         , NF90_DOUBLE, n_dim          , nsshxid          )
      iret = NF90_DEF_VAR( ntrajid, 'ssh_y'         , NF90_DOUBLE, n_dim          , nsshyid          )
      iret = NF90_DEF_VAR( ntrajid, 'sst'           , NF90_DOUBLE, n_dim          , nsstid           )
      iret = NF90_DEF_VAR( ntrajid, 'icnt'          , NF90_DOUBLE, n_dim          , ncntid           )
      iret = NF90_DEF_VAR( ntrajid, 'ithk'          , NF90_DOUBLE, n_dim          , nthkid           )
      iret = NF90_DEF_VAR( ntrajid, 'mass'          , NF90_DOUBLE, n_dim          , nmassid          )
      iret = NF90_DEF_VAR( ntrajid, 'thickness'     , NF90_DOUBLE, n_dim          , nthicknessid     )
      iret = NF90_DEF_VAR( ntrajid, 'width'         , NF90_DOUBLE, n_dim          , nwidthid         )
      iret = NF90_DEF_VAR( ntrajid, 'length'        , NF90_DOUBLE, n_dim          , nlengthid        )
      iret = NF90_DEF_VAR( ntrajid, 'year'          , NF90_INT   , n_dim          , nyearid          )
      iret = NF90_DEF_VAR( ntrajid, 'day'           , NF90_DOUBLE, n_dim          , ndayid           )
      iret = NF90_DEF_VAR( ntrajid, 'mass_of_bits'  , NF90_DOUBLE, n_dim          , nmass_of_bits_id )
      iret = NF90_DEF_VAR( ntrajid, 'heat_density'  , NF90_DOUBLE, n_dim          , nheat_density_id )
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedx_stored'  , NF90_DOUBLE, n_dim          , zdrag_groundedx_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedy_stored'  ,NF90_DOUBLE, n_dim          , zdrag_groundedy_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedxG_stored'  ,NF90_DOUBLE, n_dim          , zdrag_groundedxG_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedyG_stored'  ,NF90_DOUBLE, n_dim          , zdrag_groundedyG_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedxS_stored'  ,NF90_DOUBLE, n_dim          , zdrag_groundedxS_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'zdrag_groundedyS_stored'  ,NF90_DOUBLE, n_dim          , zdrag_groundedyS_storedid)
      iret = NF90_DEF_VAR( ntrajid, 'Xdrag_ocn_stored'  , NF90_DOUBLE, n_dim, Xdrag_ocn_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Xdrag_atm_stored'  , NF90_DOUBLE, n_dim, Xdrag_atm_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Xdrag_ice_stored'  , NF90_DOUBLE, n_dim, Xdrag_ice_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Xgrad_press_sc_stored'  , NF90_DOUBLE, n_dim, Xgrad_press_sc_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Xdrag_Cori_stored'  , NF90_DOUBLE, n_dim, Xdrag_Cori_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Ydrag_ocn_stored'  , NF90_DOUBLE, n_dim,  Ydrag_ocn_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Ydrag_atm_stored'  , NF90_DOUBLE, n_dim,  Ydrag_atm_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Ydrag_ice_stored'  , NF90_DOUBLE, n_dim,  Ydrag_ice_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Ygrad_press_sc_stored'  , NF90_DOUBLE, n_dim,  Ygrad_press_sc_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'Ydrag_Cori_stored'  , NF90_DOUBLE, n_dim, Ydrag_Cori_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'zaxe_stored'  , NF90_DOUBLE, n_dim, zaxe_storedid ) 
      iret = NF90_DEF_VAR( ntrajid, 'zaye_stored'  , NF90_DOUBLE, n_dim, zaye_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'zaxeBeforeFG_stored'  , NF90_DOUBLE, n_dim, zaxeBeforeFG_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'zayeBeforeFG_stored'  , NF90_DOUBLE, n_dim, zayeBeforeFG_storedid )
      iret = NF90_DEF_VAR( ntrajid, 'zFequil'  , NF90_DOUBLE, n_dim, zFequilid )
      iret = NF90_DEF_VAR( ntrajid, 'zFdisequil'  , NF90_DOUBLE, n_dim, zFdisequilid )
      iret = NF90_DEF_VAR( ntrajid, 'ibathy_x1y1'  , NF90_DOUBLE, n_dim, ibathy_x1y1id )
      iret = NF90_DEF_VAR( ntrajid, 'issh_x1y1'  , NF90_DOUBLE, n_dim, issh_x1y1id )


      ! Attributes
      iret = NF90_PUT_ATT( ntrajid, numberid        , 'long_name', 'iceberg number on this processor' )
      iret = NF90_PUT_ATT( ntrajid, numberid        , 'units'    , 'count' )
      iret = NF90_PUT_ATT( ntrajid, nstepid         , 'long_name', 'timestep number kt' )
      iret = NF90_PUT_ATT( ntrajid, nstepid         , 'units'    , 'count' )
      iret = NF90_PUT_ATT( ntrajid, nlonid          , 'long_name', 'longitude' )
      iret = NF90_PUT_ATT( ntrajid, nlonid          , 'units'    , 'degrees_E')
      iret = NF90_PUT_ATT( ntrajid, nlatid          , 'long_name', 'latitude' )
      iret = NF90_PUT_ATT( ntrajid, nlatid          , 'units'    , 'degrees_N' )
      iret = NF90_PUT_ATT( ntrajid, nxid            , 'long_name', 'x grid box position' )
      iret = NF90_PUT_ATT( ntrajid, nxid            , 'units'    , 'fractional' )
      iret = NF90_PUT_ATT( ntrajid, nyid            , 'long_name', 'y grid box position' )
      iret = NF90_PUT_ATT( ntrajid, nyid            , 'units'    , 'fractional' )
      iret = NF90_PUT_ATT( ntrajid, nuvelid         , 'long_name', 'zonal velocity' )
      iret = NF90_PUT_ATT( ntrajid, nuvelid         , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nvvelid         , 'long_name', 'meridional velocity' )
      iret = NF90_PUT_ATT( ntrajid, nvvelid         , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nssuid          , 'long_name', 'ocean u component' )
      iret = NF90_PUT_ATT( ntrajid, nssuid          , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nssvid          , 'long_name', 'ocean v component' )
      iret = NF90_PUT_ATT( ntrajid, nssvid          , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nuaid           , 'long_name', 'atmosphere u component' )
      iret = NF90_PUT_ATT( ntrajid, nuaid           , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nvaid           , 'long_name', 'atmosphere v component' )
      iret = NF90_PUT_ATT( ntrajid, nvaid           , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nuiid           , 'long_name', 'sea ice u component' )
      iret = NF90_PUT_ATT( ntrajid, nuiid           , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nviid           , 'long_name', 'sea ice v component' )
      iret = NF90_PUT_ATT( ntrajid, nviid           , 'units'    , 'm/s' )
      iret = NF90_PUT_ATT( ntrajid, nsshxid         , 'long_name', 'sea surface height gradient from x points' )
      iret = NF90_PUT_ATT( ntrajid, nsshxid         , 'units'    , 'm/m' )
      iret = NF90_PUT_ATT( ntrajid, nsshyid         , 'long_name', 'sea surface height gradient from y points' )
      iret = NF90_PUT_ATT( ntrajid, nsshyid         , 'units'    , 'm/m' )
      iret = NF90_PUT_ATT( ntrajid, nsstid          , 'long_name', 'sea surface temperature' )
      iret = NF90_PUT_ATT( ntrajid, nsstid          , 'units'    , 'degC')
      iret = NF90_PUT_ATT( ntrajid, ncntid          , 'long_name', 'sea ice concentration' )
      iret = NF90_PUT_ATT( ntrajid, ncntid          , 'units'    , 'degC')
      iret = NF90_PUT_ATT( ntrajid, nthkid          , 'long_name', 'sea ice thickness' )
      iret = NF90_PUT_ATT( ntrajid, nthkid          , 'units'    , 'm'   )
      iret = NF90_PUT_ATT( ntrajid, nmassid         , 'long_name', 'mass')
      iret = NF90_PUT_ATT( ntrajid, nmassid         , 'units'    , 'kg'  )
      iret = NF90_PUT_ATT( ntrajid, nthicknessid    , 'long_name', 'thickness' )
      iret = NF90_PUT_ATT( ntrajid, nthicknessid    , 'units'    , 'm'   )
      iret = NF90_PUT_ATT( ntrajid, nwidthid        , 'long_name', 'width' )
      iret = NF90_PUT_ATT( ntrajid, nwidthid        , 'units'    , 'm'    )
      iret = NF90_PUT_ATT( ntrajid, nlengthid       , 'long_name', 'length' )
      iret = NF90_PUT_ATT( ntrajid, nlengthid       , 'units'    , 'm'    )
      iret = NF90_PUT_ATT( ntrajid, nyearid         , 'long_name', 'calendar year' )
      iret = NF90_PUT_ATT( ntrajid, nyearid         , 'units'    , 'years' )
      iret = NF90_PUT_ATT( ntrajid, ndayid          , 'long_name', 'day of year' )
      iret = NF90_PUT_ATT( ntrajid, ndayid          , 'units'    , 'days' )
      iret = NF90_PUT_ATT( ntrajid, nscaling_id     , 'long_name', 'scaling factor for mass of berg' )
      iret = NF90_PUT_ATT( ntrajid, nscaling_id     , 'units'    , 'none' )
      iret = NF90_PUT_ATT( ntrajid, nmass_of_bits_id, 'long_name', 'mass of bergy bits' )
      iret = NF90_PUT_ATT( ntrajid, nmass_of_bits_id, 'units'    , 'kg'   )
      iret = NF90_PUT_ATT( ntrajid, nheat_density_id, 'long_name', 'heat density' )
      iret = NF90_PUT_ATT( ntrajid, nheat_density_id, 'units'    , 'J/kg' )
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedx_storedid         ,'long_name', 'Xfric')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedx_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedy_storedid         ,'long_name', 'Yfric')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedy_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedxG_storedid         ,'long_name', 'Xgrav')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedxG_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedyG_storedid         ,'long_name', 'Ygrav')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedyG_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedxS_storedid         ,'long_name', 'Xsedi')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedxS_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedyS_storedid         ,'long_name', 'Ysedi')
      iret = NF90_PUT_ATT( ntrajid, zdrag_groundedyS_storedid         ,'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_ocn_storedid         , 'long_name','XOcnAccel')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_ocn_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_atm_storedid         , 'long_name','XAtmAccel')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_atm_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_ice_storedid         , 'long_name','XIceAccel')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_ice_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Xgrad_press_sc_storedid         , 'long_name','XSSHAccel')
      iret = NF90_PUT_ATT( ntrajid, Xgrad_press_sc_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_Cori_storedid        , 'long_name','XCoriAccel')
      iret = NF90_PUT_ATT( ntrajid, Xdrag_Cori_storedid        , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_ocn_storedid         , 'long_name','YOcnAccel')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_ocn_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_atm_storedid         , 'long_name','YAtmAccel')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_atm_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_ice_storedid         , 'long_name','YIceAccel')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_ice_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Ygrad_press_sc_storedid         , 'long_name','YSSHAccel')
      iret = NF90_PUT_ATT( ntrajid, Ygrad_press_sc_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_Cori_storedid        , 'long_name','YCoriAccel')
      iret = NF90_PUT_ATT( ntrajid, Ydrag_Cori_storedid        , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zaxe_storedid          ,'long_name','XNetAccel')
      iret = NF90_PUT_ATT( ntrajid, zaxe_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zaye_storedid         , 'long_name','YNetAccel')
      iret = NF90_PUT_ATT( ntrajid, zaye_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zaxeBeforeFG_storedid         , 'long_name','XAccelBeforeFG')
      iret = NF90_PUT_ATT( ntrajid, zaxeBeforeFG_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zayeBeforeFG_storedid         , 'long_name','YAccelBeforeFG')
      iret = NF90_PUT_ATT( ntrajid, zayeBeforeFG_storedid         , 'units', 'm/s2')
      iret = NF90_PUT_ATT( ntrajid, zFequilid         , 'long_name','FreeBequilib')
      iret = NF90_PUT_ATT( ntrajid, zFequilid         , 'units', 'm')
      iret = NF90_PUT_ATT( ntrajid, zFdisequilid         , 'long_name','FreeBdisequil')
      iret = NF90_PUT_ATT( ntrajid, zFdisequilid         , 'units', 'm')
      iret = NF90_PUT_ATT( ntrajid, ibathy_x1y1id         , 'long_name','Bathy')
      iret = NF90_PUT_ATT( ntrajid, ibathy_x1y1id         , 'units', 'm')
      iret = NF90_PUT_ATT( ntrajid, issh_x1y1id         , 'long_name','SSH')
      iret = NF90_PUT_ATT( ntrajid, issh_x1y1id         , 'units', 'm')


      !
      ! End define mode
      iret = NF90_ENDDEF( ntrajid )
      !
   END SUBROUTINE icb_trj_init


   SUBROUTINE icb_trj_write( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_trj_write  ***
      !!
      !! ** Purpose :   write out iceberg trajectories
      !!
      !! ** Method  : - for the moment write out each snapshot of positions later
      !!                can rewrite so that it is buffered and written out more efficiently
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! time-step index
      !
      INTEGER                ::   iret, jn
      CHARACTER(len=80)      ::   cl_filename
      TYPE(iceberg), POINTER ::   this
      TYPE(point  ), POINTER ::   pt
      !!----------------------------------------------------------------------

      ! Write variables
      ! sga - just write out the current point of the trajectory

      this => first_berg
      jn = num_traj
      DO WHILE( ASSOCIATED(this) )
         pt => this%current_point
         jn = jn + 1
         !
         iret = NF90_PUT_VAR( ntrajid, numberid        , this%number      , (/1,jn/) , (/nkounts,1/) )
         iret = NF90_PUT_VAR( ntrajid, nstepid         , kt               , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nscaling_id     , this%mass_scaling, (/ jn /) )
         !
         iret = NF90_PUT_VAR( ntrajid, nlonid          , pt%lon           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nlatid          , pt%lat           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nxid            , pt%xi            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nyid            , pt%yj            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nuvelid         , pt%uvel          , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nvvelid         , pt%vvel          , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nssuid          , pt%ssu           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nssvid          , pt%ssv           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nuaid           , pt%ua            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nvaid           , pt%va            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nuiid           , pt%ui            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nviid           , pt%vi            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nsshxid         , pt%ssh_x         , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nsshyid         , pt%ssh_y         , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nsstid          , pt%sst           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, ncntid          , pt%cn            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nthkid          , pt%hi            , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nmassid         , pt%mass          , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nthicknessid    , pt%thickness     , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nwidthid        , pt%width         , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nlengthid       , pt%length        , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nyearid         , pt%year          , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, ndayid          , pt%day           , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nmass_of_bits_id, pt%mass_of_bits  , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, nheat_density_id, pt%heat_density  , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedx_storedid, pt%zdrag_groundedx_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedy_storedid, pt%zdrag_groundedy_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedxG_storedid, pt%zdrag_groundedxG_stored  , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedyG_storedid, pt%zdrag_groundedyG_stored  ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedxS_storedid, pt%zdrag_groundedxS_stored  , (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zdrag_groundedyS_storedid, pt%zdrag_groundedyS_stored  ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zaxe_storedid, pt%zaxe_stored   ,  (/ jn /) ) 
         iret = NF90_PUT_VAR( ntrajid, zaye_storedid, pt%zaye_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Xdrag_ocn_storedid, pt%Xdrag_ocn_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Xdrag_atm_storedid, pt%Xdrag_atm_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Xdrag_ice_storedid, pt%Xdrag_ice_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Xgrad_press_sc_storedid, pt%Xgrad_press_sc_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Xdrag_Cori_storedid, pt%Xdrag_Cori_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Ydrag_ocn_storedid, pt%Ydrag_ocn_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Ydrag_atm_storedid, pt%Ydrag_atm_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Ydrag_ice_storedid, pt%Ydrag_ice_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Ygrad_press_sc_storedid, pt%Ygrad_press_sc_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, Ydrag_Cori_storedid, pt%Ydrag_Cori_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zaxeBeforeFG_storedid, pt%zaxeBeforeFG_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zayeBeforeFG_storedid, pt%zayeBeforeFG_stored   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zFequilid, pt%zFequil   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, zFdisequilid, pt%zFdisequil    ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, ibathy_x1y1id, pt%ibathy_x1y1   ,  (/ jn /) )
         iret = NF90_PUT_VAR( ntrajid, issh_x1y1id, pt%issh_x1y1   ,  (/ jn /) )
         !
         this => this%next
      END DO
      IF( lwp .AND. nn_verbose_level > 0 )   WRITE(numout,*) 'trajectory write to frame ', jn
      num_traj = jn
      !
   END SUBROUTINE icb_trj_write


   SUBROUTINE icb_trj_sync()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_trj_sync  ***
      !!
      !! ** Purpose :   
      !!----------------------------------------------------------------------
      INTEGER ::   iret
      !!----------------------------------------------------------------------
      ! flush to file
      iret = NF90_SYNC( ntrajid )
      IF ( iret /= NF90_NOERR )   CALL ctl_stop( 'icebergs, icb_trj_sync: nf_sync failed' )
      !
   END SUBROUTINE icb_trj_sync


   SUBROUTINE icb_trj_end()
      !!----------------------------------------------------------------------
      INTEGER ::   iret
      !!----------------------------------------------------------------------
      ! Finish up
      iret = NF90_CLOSE( ntrajid )
      IF ( iret /= NF90_NOERR )   CALL ctl_stop( 'icebergs, icb_trj_end: nf_close failed' )
      !
   END SUBROUTINE icb_trj_end

   !!======================================================================
END MODULE icbtrj
