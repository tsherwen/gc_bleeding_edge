!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bromocarb_mod.F
!
! !DESCRIPTION: Module BROMOCARB\_MOD contains variables and routines for the
!  GEOS-Chem bromocarbon simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE BROMOCARB_MOD
!
! !USES:
!
      USE PRECISION_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: SET_TROP_CONC
      PUBLIC :: SET_CH3Br
      PUBLIC :: SET_BrO
!
! !REMARKS:
!  HEMCO has made most of this routine obsolete.  HEMCO reads in the
!  data files in netCDF format and tracks them in its data structure.
!
! !REVISION HISTORY:
!  13 Aug 2007 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  27 Aug 2012 - M. Payer    - Now parallelize DO loops
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Jun 2014 - R. Yantosca - Removed code made obsolete by HEMCO
!  24 Jun 2014 - R. Yantosca - Removed INIT_BROMOCARB, CLEANUP_BROMOCARB
!  06 Nov 2015 - M. Yannetti - Added PRECISION_MOD
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of CHBr3  [kg/mol]
      REAL(fp), PARAMETER :: MWT_CHBr3  = 2.53e-1_fp

      ! Molecular weight of CH2Br2 [kg/mol]
      REAL(fp), PARAMETER :: MWT_CH2Br2 = 1.74e-1_fp

      ! Number of elements for TOT_OC_AREA
      INTEGER,  PARAMETER :: NUM_REGIONS = 2
!
! PRIVATE TYPES:
!
      ! Surface area of the box; only depends on latitude
      REAL(fp)            :: A_M2

      ! Species ID flags
      INTEGER             :: id_CH3Br
      INTEGER             :: id_BrO

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ch3br
!
! !DESCRIPTION: Subroutine SET\_CH3Br set CH3Br Concentrations in the planetary
!  boundary layer. Based on latitude bands (1) 90-55N, (2) 55N-0, (3) 0-55S,
!  (4) 55-90S. Values for setting pbl flux were determined by surface
!  measurements from NOAA 2006 data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CH3Br( am_I_Root, Input_Opt, State_Met,
     &                      State_Chm, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD
      USE GC_GRID_MOD,        ONLY : GET_YMID
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PhysConstants,      ONLY : AIRMW
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!   1) Spc is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_CH3Br in main.f (it's
!      with the emissions) means that it should be in [kg].
!
! !REVISION HISTORY:
!  12 Feb 2008 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  27 Aug 2012 - M. Payer    - Added parallel DO loop
!  28 Aug 2012 - M. Payer    - Add error check for CH3Br to avoid OOB error
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  24 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Met,
!                              State_Chm.  Get fields from these objects.
!  10 Sep 2015 - E. Lundgren - Modification for State_Chm%TRACERS units
!                              now in kg/kg dry air (previously kg/box)
!  10 Sep 2015 - E. Lundgren - Remove passed argument 'unitflag'
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  25 Jul 2016 - M. Yannetti - Replaced TCVV with MW from spec db and phys const
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3Br values ( from pptv )
      REAL(fp), PARAMETER :: gt55N     = 8.35e-12_fp
      REAL(fp), PARAMETER :: gt0_lt55N = 8.27e-12_fp
      REAL(fp), PARAMETER :: lt0_gt55S = 6.94e-12_fp
      REAL(fp), PARAMETER :: lt55S     = 6.522e-12_fp
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL             :: FIRST = .TRUE.

      ! Scalars
      INTEGER             :: I,  J,  L
      REAL(fp)            :: CH3Br_conc
      REAL(fp)            :: LAT_MID

      ! Pointers
      REAL(fp),   POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! Initialization
      !=================================================================

      ! Locate the species ID of CH3Br on the first call
      IF ( FIRST ) THEN
         id_CH3Br = Ind_('CH3Br')
         FIRST    = .FALSE.
      ENDIF

      ! Return if CH3Br is not found
      IF ( id_CH3Br <= 0 .and. am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_CH3Br: CH3Br not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in Spc'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN
      ENDIF

      ! Check that species units are in kg/kg dry air (ewl, 9/10/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
         CALL GC_Error( 'Incorrect species units: ' //
     &                  State_Chm%Spc_Units, RC,
     &                  'Routine SET_CH3Br in bromocarb_mod.F' )
         RETURN
      ENDIF

      ! Point to the chemical species array [kg/kg dry air]
      Spc        => State_Chm%Species

      ! Zero
      LAT_MID    = 0e+0_fp
      CH3Br_conc = 0e+0_fp

      !=================================================================
      ! SET_CH3Br begins here!
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, LAT_MID, CH3Br_conc )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0e+0_fp ) THEN

            ! base lat band selection on midpoint
            ! latitude of the box
            LAT_MID = GET_YMID( I, J, L )

            ! Selecting the latitude bands:
            IF ( LAT_MID > 55e+0_fp ) THEN
               CH3Br_conc = gt55N
            ELSEIF ( (LAT_MID >= 0e+0_fp) .and. (LAT_MID <= 55e+0_fp) )
     &       THEN
               CH3Br_conc = gt0_lt55N
            ELSEIF ( (LAT_MID < 0e+0_fp) .and. (LAT_MID >= -55e+0_fp) )
     &       THEN
               CH3Br_conc = lt0_gt55S
            ELSEIF ( LAT_MID < -55e+0_fp ) THEN
               CH3Br_conc = lt55S
            ENDIF

            ! Convert the [v/v] CH3Br conc units (dry air) to [kg/kg]
            ! when setting species concentration (ewl, 9/10/15)
!            Spc(I,J,L,id_CH3Br) = CH3Br_conc / ( AIRMW /
!     &                         State_Chm%SpcData(id_CH3Br)%Info%emMW_g )
            ! Set CH3Br to a tiny number (1E-20)
            Spc(I,J,L,id_CH3Br) = 1e-20_fp

         ENDIF  ! end selection of PBL boxes

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE SET_CH3Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_trop_conc
!
! !DESCRIPTION: Subroutine SET\_TROP\_CONC sets concentrations in troposhere.
!            This subroutine is used to make sure halogen concetrations in
!            troposphere are set to zero for no halogen troposphere only runs.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TROP_CONC( am_I_Root, Input_Opt, State_Met,
     &                      State_Chm, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!   1) Spc is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_CH3Br in main.f (it's
!      with the emissions) means that it should be in [kg].
!
! !REVISION HISTORY:
!  29 Oct 2015 - T. Sherwen - Initial version
!  31 Jul 2018 - T. Sherwen - Updated to use v12.0.0 nomenclature
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Small number instead of zero to prevent integrator issues
      REAL(fp), PARAMETER :: ZERO_EFF   = 1.00e-20_fp
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL             :: FIRST = .TRUE.

      ! Scalars
      INTEGER             :: I,  J,  L
      INTEGER             :: AS, COUNT, N, NN, nSpecies, id_HOI

      ! Setup a list of halogen species and allocatable array for their indices
      INTEGER, ALLOCATABLE :: INDICES(:)
      CHARACTER  (len=10), DIMENSION(47) :: HAL_Species

      ! Pointers
      REAL(fp),   POINTER :: Spc(:,:,:,:)
      REAL(fp)            :: Val2Set

      ! Objects
      TYPE(Species), POINTER  :: SpcInfo

      ! Number of total species
      nSpecies  = State_Chm%nSpecies

      ! Set the value to set concentrations to
      Val2Set   = ZERO_EFF

      !=================================================================
      ! Initialization
      !=================================================================

      ! tms bebug
      WRITE(6,*) 'SET_TROP_CONC called - tms - trop [halogens] =  0'

      ! list halogen species
      ! for development reasons ICl, IBr, CH2IBr, CH2ICl are considered
      ! as iodine species. BrCl is considered as a bromine species
      ! chlorine species, bromine species, and iodine species
      ! The first 47 species below were set to zero in v10-01
      HAL_Species = (/
     & 'CH2Cl2    ', 'HOCl      ', 'ClNO3     ', 'ClNO2     ',
     & 'HCl       ', 'ClOO      ', 'Cl2O2     ', 'ClO       ',
     & 'CHCl3     ', 'Cl        ', 'OClO      ', 'CH3Cl     ',
     & 'Cl2       ' ,
     & 'Br2       ', 'BrNO3     ', 'Br        ', 'HBr       ',
     & 'BrO       ', 'CH3Br     ', 'CH2Br2    ', 'BrCl      ',
     & 'BrNO2     ', 'BrSALC    ', 'BrSALA    ', 'HOBr      ',
     & 'CHBr3     ' ,
     & 'IONO      ', 'OIO       ', 'HI        ', 'CH2I2     ',
     & 'IONO2     ', 'IO        ', 'CH2IBr    ', 'CH3IT     ',
     & 'CH2ICl    ','ISALC      ', 'I2O2      ', 'I2O3      ',
     & 'I2O4      ', 'IBr       ', 'I2        ', 'INO       ',
     & 'I         ', 'AERI      ', 'HOI       ', 'ICl       ',
     &  'ISALA     ' /)


      ! Index array for halogen tracers indices
      ALLOCATE( INDICES( SIZE( HAL_Species ) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'INDICES' )
      INDICES = 0

      ! get index of species in nSpecies
      COUNT = 1
      DO N = 1, nSpecies
         ! NOTE: requires 1:1 tracer index to species index mapping
         SpcInfo => State_Chm%SpcData(N)%Info

         ! if N_TRACER(N) in species add indices
         DO NN = 1, SIZE( HAL_Species )
            IF( TRIM( SpcInfo%Name ) == TRIM( HAL_Species(NN) ) ) THEN
               INDICES( COUNT ) = Ind_( TRIM(HAL_Species(NN)) )
               COUNT = COUNT + 1
            ENDIF
         ENDDO
      ENDDO

      ! debug
      WRITE(6,*) 'tms debug 3 - SET_TROP_CONC - ', COUNT, INDICES,
     & nSpecies, HAL_Species

      ! Locate the species ID of HOI on the first call
      IF ( FIRST ) THEN
         id_HOI = Ind_('HOI')
         FIRST    = .FALSE.
      ENDIF

      ! Return if HOI is not found
      IF ( id_HOI <= 0 .and. am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_CH3Br: CH3Br not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in Spc'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN
      ENDIF

      ! Check that species units are in kg/kg dry air (ewl, 9/10/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
         CALL GC_Error( 'Incorrect species units: ' //
     &                  State_Chm%Spc_Units, RC,
     &                  'Routine SET_TROP_CONC in bromocarb_mod.F' )
         RETURN
      ENDIF

      ! Point to the chemical species array [kg/kg dry air]
      Spc        => State_Chm%Species

      !=================================================================
      ! SET_TROP_CONC begins here!
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, Val2Set )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( State_Met%InTropoSphere( I, J, L ) ) THEN
            ! if N_TRACER(N) in species add indices
            DO N = 1, SIZE( HAL_Species )
               ! get the Spec array index for the HAL species
               NN = INDICES( N )
               ! Set the halogen conc value to ZERO_EFF
               Spc(I,J,L, NN) = Val2Set
            ENDDO
         ENDIF  ! end selection of Troposphere boxes

      ENDDO ! End I
      ENDDO ! End J
      ENDDO ! End L
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE SET_TROP_CONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_bro
!
! !DESCRIPTION: Subroutine SET\_BRO sets BrO concentrations in the planetary
!  boundary layer. Based on latitude bands (1) 90-55N, (2) 55N-0, (3) 0-55S,
!  (4) 55-90S. Values for setting pbl flux were determined by surface
!  measurements from NOAA 2006 data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_BRO( am_I_Root, Input_Opt, State_Met,
     &                    State_Chm, RC         )

!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : IS_WATER
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PhysConstants,      ONLY : AIRMW
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!   1) Spc is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_Bro in main.f (it's
!      with the emissions) means that it should be in [kg].
!
! !REVISION HISTORY:
!  12 Feb 2008 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  27 Aug 2012 - M. Payer    - Added parallel DO loop
!  28 Aug 2012 - M. Payer    - Add error check for BrO to avoid OOB error
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  24 Jun 2014 - R. Yantosca - Remove IJLOOP, it's not used anymroe
!  06 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp)
!  10 Sep 2015 - E. Lundgren - Modification for State_Chm%TRACERS units
!                              now in kg/kg dry air (previously kg/box)
!  10 Sep 2015 - E. Lundgren - Remove passed argument 'unitflag'
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  25 Jul 2016 - M. Yannetti - Replaced TCVV with MW from spec db and phys const
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL           :: FIRST = .TRUE.

      ! Scalars
      INTEGER           :: I,  J,  L,  N
      REAL(fp)          :: BrO_conc

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !====================================================================
      ! SET_BrO begins here!
      !====================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Get the species ID for BrO only on the first call
      IF ( FIRST ) THEN
         id_BrO = Ind_('BrO')
         FIRST  = .FALSE.
      ENDIF

      ! Return if BrO is not found
      IF ( id_BrO <= 0 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_BrO: BrO not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in Spc'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN
      ENDIF

      ! -----------------------------------------
      ! If we aren't using this 1pptv experiment
      ! then return without updating Spc array.
      ! -----------------------------------------
      IF ( .not. Input_Opt%LFIX_PBL_BRO ) THEN
         RETURN
      ENDIF

      ! Point to the chemical species array [kg/kg dry air]
      ! Reference fields in derived-type objects
      Spc => State_Chm%Species

      ! ------------------------------------------------------
      ! Set the BrO concentration to 1 pptv inside the PBL.
      ! ------------------------------------------------------
      BrO_conc = 1.0e-12_fp

      ! Check that species units are in kg/kg dry air (ewl, 9/10/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
         CALL GC_Error( 'Incorrect species units: ' //
     &                  State_Chm%Spc_Units, RC,
     &                  'Routine SET_BrO in bromocarb_mod.F' )
         RETURN
      ENDIF

      !=================================================================
      ! SET_BrO begins here!
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, BrO_conc )
      DO J = 1, JJPAR
         lon_loop: DO I = 1, IIPAR

            ! -----------------------------------------
            ! 1. Determine if we're in the marine
            !    boundary layer. If so, procede,
            !    otherwise, skip. (note, we should NOT
            !    0 the concentration... it can be
            !    contributed from other timesteps and
            !    sources.
            ! -----------------------------------------
            IF ( .not. IS_WATER( I, J, State_Met) ) THEN
               CYCLE lon_loop
            ENDIF

            ! -----------------------------------------
            ! 2. Get the cosine of the SZA to determine
            !    if there's available sunlight for
            !    activation of bromine-chemistry.
            !    If so, set the PBL BrO to 1ppt.
            ! -----------------------------------------

            IF ( State_Met%SUNCOS(I,J) > 0.d0 ) THEN
               BrO_conc = 1.0e-12_fp ! 1pptv if daytime
            ELSE
               BrO_conc = 0.e+0_fp    ! 0 otherwise
            ENDIF

            DO L = 1, LLPAR
               IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

                  ! Convert the [v/v] BrO conc units (dry air) to [kg/kg]
                  ! when setting species concentration (ewl, 9/10/15)
                  Spc(I,J,L,id_BrO) = BrO_conc / ( AIRMW
     &                       / State_Chm%SpcData(id_BrO)%Info%emMW_g )

               ENDIF            ! end selection of PBL boxes

            END DO              !END l-loop over altitudes
         END DO lon_loop        !END i-loop over longitudes
      END DO                    !END j-loop over latitudes
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE SET_BRO
!EOC
      END MODULE BROMOCARB_MOD
