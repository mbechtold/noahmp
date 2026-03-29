module SoilWaterMainMod

!!! Main soil water module including all soil water processes & update soil moisture
!!! surface runoff, infiltration, soil water diffusion, subsurface runoff, tile drainage
!!! Peatland-specific options (Chakraborty et al., 2025; revised Bechtold, 2026)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use RunoffSurfaceTopModelGrdMod,       only : RunoffSurfaceTopModelGrd
  use RunoffSurfaceTopModelEquiMod,      only : RunoffSurfaceTopModelEqui
  use RunoffSurfaceFreeDrainMod,         only : RunoffSurfaceFreeDrain
  use RunoffSurfaceBatsMod,              only : RunoffSurfaceBATS
  use RunoffSurfaceTopModelMmfMod,       only : RunoffSurfaceTopModelMMF
  use RunoffSurfaceVicMod,               only : RunoffSurfaceVIC
  use RunoffSurfaceXinAnJiangMod,        only : RunoffSurfaceXinAnJiang
  use RunoffSurfaceDynamicVicMod,        only : RunoffSurfaceDynamicVic
  use RunoffSubSurfaceEquiWaterTableMod, only : RunoffSubSurfaceEquiWaterTable
  use RunoffSubSurfaceGroundWaterMod,    only : RunoffSubSurfaceGroundWater
  use RunoffSubSurfaceDrainageMod,       only : RunoffSubSurfaceDrainage
  use RunoffSubSurfaceShallowMmfMod,     only : RunoffSubSurfaceShallowWaterMMF
  use RunoffSubSurfacePeatlandMod,       only : RunoffSubSurfacePeatland
  use MicroTopoCorrectionMod,            only : MicroTopoCorrection
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics
  use PeatMicroTopoMod,                  only : FloodedFrac, FsoilMicroTopo, &
                                                 EquilibriumSMFlat,           &
                                                 EquilibriumSMMicroTopo,      &
                                                 HeadShiftFromThetaFlat,      &
                                                 HeadShiftFromThetaMicro,     &
                                                 FindWaterTable,              &
                                                 FindWaterTableTotal,         &
                                                 SurfaceWaterStorage,         &
                                                 SoilWaterStorageMicroTopoLite, &
                                                 ThetaFromHeadShiftFlat,      &
                                                 ThetaFromHeadShiftMicro,     &
                                                 z_trunc
  use SoilWaterDiffusionRichardsMod,     only : SoilWaterDiffusionRichards
  use SoilMoistureSolverMod,             only : SoilMoistureSolver
  use TileDrainageSimpleMod,             only : TileDrainageSimple
  use TileDrainageHooghoudtMod,          only : TileDrainageHooghoudt
  use WaterTableEquilibriumMod,          only : WaterTableEquilibrium

  implicit none

contains

  subroutine SoilWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SOILWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Option=9; Peatland runoff scheme (Chakraborty & Bechtold, 2025; revised Bechtold, 2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout)  :: noahmp

! local variables
    integer                           :: LoopInd1, LoopInd2           ! loop index
    integer                           :: IndIter                      ! iteration index
    integer                           :: NumIterSoilWat               ! iteration times soil moisture
    real(kind=kind_noahmp)            :: TimeStepFine                 ! fine time step [s]
    real(kind=kind_noahmp)            :: SoilSatExcAcc                ! accumulation of soil saturation excess [m]
    real(kind=kind_noahmp)            :: SoilWatConductAcc            ! sum of SoilWatConductivity*ThicknessSnowSoilLayer
    real(kind=kind_noahmp)            :: WaterRemove                  ! water mass removal [mm]
    real(kind=kind_noahmp)            :: SoilWatRem                   ! temporary remaining soil water [mm]
    real(kind=kind_noahmp)            :: SoilWaterMin                 ! minimum soil water [mm]
    real(kind=kind_noahmp)            :: DrainSoilBotAcc              ! accumulated drainage water [mm] at fine time step
    real(kind=kind_noahmp)            :: RunoffSurfaceAcc             ! accumulated surface runoff [mm] at fine time step
    real(kind=kind_noahmp)            :: InfilSfcAcc                  ! accumulated infiltration rate [m/s]
    real(kind=kind_noahmp)            :: InfilRateSfc_FSW_change      ! surface water infiltration flux
    real(kind=kind_noahmp)            :: W_soil_peat                  ! total soil water content [m]
    real(kind=kind_noahmp)            :: z_wt_peat                    ! water table elevation (z positive up) [m]
    real(kind=kind_noahmp)            :: z_wt_begin                   ! z_wt at timestep start [m]
    real(kind=kind_noahmp)            :: z_wt_end                     ! z_wt at timestep end [m]
    real(kind=kind_noahmp)            :: WTD_begin                    ! WTD at timestep start [m], positive downward
    real(kind=kind_noahmp)            :: FSW_change_flux              ! FSW_change from flux accumulation [mm]
    real(kind=kind_noahmp), parameter :: WTD_equil_threshold = 0.3_kind_noahmp  ! WTD threshold for equilibrium bypass [m]
    real(kind=kind_noahmp), parameter :: SoilImpPara = 4.0            ! soil impervious fraction parameter
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft1     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft2     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft3     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqTmp   ! temporary soil liquid water [mm]
    real(kind=kind_noahmp)            :: d_top_peat, d_bot_peat       ! layer depth bounds for peat [m]
    real(kind=kind_noahmp)            :: z_col_bot_peat               ! column bottom depth [m]
    real(kind=kind_noahmp)            :: thetas_peat, ae_peat, bb_peat! Campbell peat parameters
    real(kind=kind_noahmp)            :: SM_eq_flat_tmp                ! temporary flat equilibrium SM
    real(kind=kind_noahmp)            :: SM_eq_micro_tmp               ! temporary microtopo equilibrium SM
    real(kind=kind_noahmp)            :: SM_excess_tmp                  ! excess (departure from equil) [m3/m3]
    real(kind=kind_noahmp)            :: SM_old_peat                   ! SoilLiqWater before transfer [m3/m3]
    real(kind=kind_noahmp)            :: excess_vol                    ! overflow excess volume [m]
    real(kind=kind_noahmp)            :: space_avail                   ! available pore space [m]
    real(kind=kind_noahmp)            :: transfer_vol                  ! overflow transfer volume [m]
    real(kind=kind_noahmp)            :: InfilSoil_peat                ! soil-portion infiltration [m/s]
    real(kind=kind_noahmp)            :: EvapSoil_peat                 ! soil-portion evaporation [mm/s]
    real(kind=kind_noahmp)            :: TranspSoil_peat               ! soil-portion transpiration [mm/s]
    real(kind=kind_noahmp)            :: RunoffSoil_peat               ! soil-portion runoff [mm]
    real(kind=kind_noahmp)            :: W_soil_check                   ! normalization check [m]
    real(kind=kind_noahmp)            :: W_total_peat                   ! total water (soil + surface) [m]
    real(kind=kind_noahmp)            :: Q_net_total                    ! total net flux [m]
    real(kind=kind_noahmp)            :: W_fsw_begin                    ! surface water at timestep start [m]
    real(kind=kind_noahmp)            :: W_soil_eq_end                  ! equilibrium soil water at z_wt_end [m]
    real(kind=kind_noahmp)            :: mean_delta_peat                ! mean Richards delta for conservation [m3/m3]
    real(kind=kind_noahmp)            :: delta_richards                 ! per-layer Richards SM change [m3/m3]
   real(kind=kind_noahmp)            :: head_shift_tmp                 ! equivalent pressure-head anomaly [m]
   real(kind=kind_noahmp)            :: head_shift_cap                 ! cap for head anomaly magnitude [m]
   real(kind=kind_noahmp)            :: lambda_peat                    ! linear scaling of head anomalies [-]
   real(kind=kind_noahmp)            :: W_target_peat                  ! conserved soil water target after Richards [m]
   real(kind=kind_noahmp)            :: WTD_lo_peat                    ! lower bracket for WTD solve [m]
   real(kind=kind_noahmp)            :: WTD_hi_peat                    ! upper bracket for WTD solve [m]
   real(kind=kind_noahmp)            :: WTD_mid_peat                   ! midpoint for WTD solve [m]
   real(kind=kind_noahmp)            :: W_lo_peat                      ! wet-side storage bracket [m]
   real(kind=kind_noahmp)            :: W_hi_peat                      ! dry-side storage bracket [m]
   real(kind=kind_noahmp)            :: W_mid_peat                     ! midpoint storage in WTD solve [m]
    integer                           :: LoopJ                         ! overflow cascade index
    integer                           :: SatTopInd                      ! topmost fully-saturated layer index
   integer                           :: IterWTD                        ! iteration index for WTD closure
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqWaterOrig   ! original SoilLiqWater before forward transfer [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqWater1D_bef ! 1D profile before Richards [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqGap    ! per-layer change in forward transfer [m3/m3]
   real(kind=kind_noahmp), allocatable, dimension(:) :: HeadShiftMicro ! microtopography head anomalies [m]
   real(kind=kind_noahmp), allocatable, dimension(:) :: HeadShiftFlat  ! flat-column head anomalies [m]

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    noahmp soil time step [s]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              FlagUrban              => noahmp%config%domain%FlagUrban              ,& ! in,    logical flag for urban grid
              OptRunoffSurface       => noahmp%config%nmlist%OptRunoffSurface       ,& ! in,    options for surface runoff
              OptRunoffSubsurface    => noahmp%config%nmlist%OptRunoffSubsurface    ,& ! in,    options for subsurface runoff
              OptPeatlandPhysics     => noahmp%config%nmlist%OptPeatlandPhysics     ,& ! in,    options for peatland physics
              OptTileDrainage        => noahmp%config%nmlist%OptTileDrainage        ,& ! in,    options for tile drainage
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! in,    soil ice content [m3/m3]
              TileDrainFrac          => noahmp%water%state%TileDrainFrac            ,& ! in,    tile drainage map (fraction)
              SoilSfcInflowMean      => noahmp%water%flux%SoilSfcInflowMean         ,& ! in,    mean water input on soil surface [m/s]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,    saturated soil matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,    Campbell b exponent [-]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil water content [m3/m3]
              RechargeGwDeepWT       => noahmp%water%state%RechargeGwDeepWT         ,& ! inout, recharge to or from the water table when deep [m]
              DrainSoilBot           => noahmp%water%flux%DrainSoilBot              ,& ! out,   soil bottom drainage [m/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface             ,& ! out,   surface runoff [mm per soil timestep]
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface          ,& ! out,   subsurface runoff [mm per soil timestep] 
              InfilRateSfc           => noahmp%water%flux%InfilRateSfc              ,& ! out,   infiltration rate at surface [m/s]
              TileDrain              => noahmp%water%flux%TileDrain                 ,& ! out,   tile drainage [mm per soil timestep]
              Transpiration          => noahmp%water%flux%Transpiration             ,& ! in,    transpiration rate [mm/s]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet             ,& ! in,    net ground (soil/snow) evaporation [mm/s]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! inout,   water table depth [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax        ,& ! out,   maximum soil imperviousness fraction
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity      ,& ! out,   soil hydraulic conductivity [m/s]
              SoilEffPorosity        => noahmp%water%state%SoilEffPorosity          ,& ! out,   soil effective porosity [m3/m3]
              SoilImpervFrac         => noahmp%water%state%SoilImpervFrac           ,& ! out,   impervious fraction due to frozen soil
              SoilIceFrac            => noahmp%water%state%SoilIceFrac              ,& ! out,   ice fraction in frozen soil
              SoilSaturationExcess   => noahmp%water%state%SoilSaturationExcess     ,& ! out,   saturation excess of the total soil [m]
              SoilIceMax             => noahmp%water%state%SoilIceMax               ,& ! out,   maximum soil ice content [m3/m3]
              FSW_change             => noahmp%water%state%FSW_change               ,& ! inout,   surface storage change [mm]
              FSW_peat_error         => noahmp%water%state%FSW_peat_error           ,& ! inout,   peatland numerical water balance error [mm]
              FloodedFraction        => noahmp%water%state%FloodedFraction          ,& ! inout,   flooded fraction [-]
              f_soil                 => noahmp%water%state%f_soil                   ,& ! inout, fraction of flux in and out of soil [-]
              SoilLiqWaterMin        => noahmp%water%state%SoilLiqWaterMin         ,& ! out,   minimum soil liquid water content [m3/m3]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth [m] of layer-bottom from soil surface
              TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean      & ! inout, mean transpiration water loss from soil layers [m/s]
             )
! ----------------------------------------------------------------------

    ! initialization
    if (.not. allocated(MatRight)  ) allocate(MatRight  (1:NumSoilLayer))
    if (.not. allocated(MatLeft1)  ) allocate(MatLeft1  (1:NumSoilLayer))
    if (.not. allocated(MatLeft2)  ) allocate(MatLeft2  (1:NumSoilLayer))
    if (.not. allocated(MatLeft3)  ) allocate(MatLeft3  (1:NumSoilLayer))
    if (.not. allocated(SoilLiqTmp)) allocate(SoilLiqTmp(1:NumSoilLayer))
    MatRight         = 0.0
    MatLeft1         = 0.0
    MatLeft2         = 0.0
    MatLeft3         = 0.0
    SoilLiqTmp       = 0.0
    RunoffSurface    = 0.0
    RunoffSubsurface = 0.0
    InfilRateSfc     = 0.0
    SoilSatExcAcc    = 0.0
    InfilSfcAcc      = 1.0e-06
    
    if ( OptPeatlandPhysics == 1 ) then
            call ApplyPeatlandPhysics(noahmp) 
    endif

    ! for the case when snowmelt water is too large
    do LoopInd1 = 1, NumSoilLayer
       SoilEffPorosity(LoopInd1) = max(1.0e-4, (SoilMoistureSat(LoopInd1) - SoilIce(LoopInd1)))
       SoilSatExcAcc             = SoilSatExcAcc + max(0.0, SoilLiqWater(LoopInd1) - SoilEffPorosity(LoopInd1)) * &
                                                   ThicknessSnowSoilLayer(LoopInd1)
       SoilLiqWater(LoopInd1)    = min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1))
    enddo
    write(*,*) 'DEBUG: SoilSatExcAcc = ', SoilSatExcAcc

    ! impermeable fraction due to frozen soil
    do LoopInd1 = 1, NumSoilLayer
       SoilIceFrac(LoopInd1)    = min(1.0, SoilIce(LoopInd1) / SoilMoistureSat(LoopInd1))
       SoilImpervFrac(LoopInd1) = max(0.0, exp(-SoilImpPara*(1.0-SoilIceFrac(LoopInd1))) - exp(-SoilImpPara)) / &
                                  (1.0 - exp(-SoilImpPara))
    enddo

    ! maximum soil ice content and minimum liquid water of all layers
    SoilIceMax        = 0.0
    SoilImpervFracMax = 0.0
    SoilLiqWaterMin   = SoilMoistureSat(1)
    do LoopInd1 = 1, NumSoilLayer
       if ( SoilIce(LoopInd1) > SoilIceMax )               SoilIceMax        = SoilIce(LoopInd1)
       if ( SoilImpervFrac(LoopInd1) > SoilImpervFracMax ) SoilImpervFracMax = SoilImpervFrac(LoopInd1)
       if ( SoilLiqWater(LoopInd1) < SoilLiqWaterMin )     SoilLiqWaterMin   = SoilLiqWater(LoopInd1)
    enddo

    ! subsurface runoff for runoff scheme option 2
    if ( OptRunoffSubsurface == 2 ) call RunoffSubSurfaceEquiWaterTable(noahmp)
    
    ! ================================================================
    ! Peatland: Set up parameters and z_wt_begin.
    ! Use stored WaterTableDepth directly (from restart or previous
    ! timestep) rather than re-diagnosing from soil water alone.
    ! This ensures cross-timestep consistency of combined
    ! soil + surface water tracking.  (Bechtold, 2026)
    ! ================================================================
    if ( OptPeatlandPhysics == 1 ) then
       thetas_peat    = SoilMoistureSat(1)
       ae_peat        = abs(SoilMatPotentialSat(1))
       bb_peat        = SoilExpCoeffB(1)
       z_col_bot_peat = abs(DepthSoilLayer(NumSoilLayer))

       ! Compute total soil water [m]
       W_soil_peat = 0.0_kind_noahmp
       do LoopInd1 = 1, NumSoilLayer
          W_soil_peat = W_soil_peat + &
              SoilLiqWater(LoopInd1) * abs(ThicknessSnowSoilLayer(LoopInd1))
       enddo

       ! Use stored WTD (continuous across timesteps)

       z_wt_begin = -WaterTableDepth
       WTD_begin  = WaterTableDepth
       !z_wt_begin = FindWaterTable(W_soil_peat, thetas_peat, ae_peat, &
       !bb_peat, z_col_bot_peat, -WaterTableDepth)
       !WTD_begin  = -z_wt_begin
    endif

    ! Peatland: Ivanov runoff (uses diagnosed WTD)
    if ( OptRunoffSubsurface == 9 ) then
            call RunoffSubSurfacePeatland(noahmp)
    endif


    ! jref impermable surface at urban
    if ( FlagUrban .eqv. .true. ) SoilImpervFrac(1) = 0.95

    ! ================================================================
    ! Peatland: Two-regime algorithm (Bechtold, 2026)
    ! Equilibrium bypass for WTD < 0.3 m, Richards path for WTD >= 0.3 m
    ! ================================================================
    if ( OptPeatlandPhysics == 1 ) then

       ! --- Common initial steps ---
       call MicroTopoCorrection(noahmp)   ! Compute f_soil and FloodedFraction

       ! surface runoff and infiltration rate
       if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
       if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
       if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
       if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
       if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
       if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
       if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
       if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)

       ! Add RunoffSurface back to infiltration for peatland
       InfilRateSfc = InfilRateSfc + RunoffSurface
       RunoffSurface = 0.0

       ! Partition infiltration flux
       InfilRateSfc_FSW_change = (1.0_kind_noahmp - f_soil) * InfilRateSfc
       InfilRateSfc = f_soil * InfilRateSfc

       ! ================================================================
       ! EQUILIBRIUM PATH: WTD < 0.3 m (Shallow WT, bypass Richards)
       ! Uses TOTAL water (soil + surface) to find z_wt consistently.
       ! ================================================================
       if ( WTD_begin < WTD_equil_threshold ) then

          ! Surface water storage at timestep begin [m]
          W_fsw_begin = SurfaceWaterStorage(z_wt_begin)

          ! Total water (soil + surface) at begin [m]
          W_total_peat = W_soil_peat + W_fsw_begin

          ! Total net flux [m]:  InfilRateSfc + InfilRateSfc_FSW_change = full (pre-split) infiltration
          Q_net_total = (InfilRateSfc + InfilRateSfc_FSW_change) * SoilTimeStep           &
                      - EvapGroundNet * SoilTimeStep / 1000.0_kind_noahmp                 &
                      - Transpiration * SoilTimeStep / 1000.0_kind_noahmp                 &
                      - RunoffSubsurface * SoilTimeStep / 1000.0_kind_noahmp

          W_total_peat = W_total_peat + Q_net_total

          ! Find z_wt from combined soil+surface storage (exact conservation)
          z_wt_end = FindWaterTableTotal(W_total_peat, thetas_peat, ae_peat, &
              bb_peat, z_col_bot_peat, z_wt_begin)
          WaterTableDepth = -z_wt_end

          ! Soil water from the new z_wt [m]
          W_soil_peat = SoilWaterStorageMicroTopoLite(z_wt_end, thetas_peat, ae_peat, &
              bb_peat, z_col_bot_peat)

          ! Set soil moisture to column-averaged hydrostatic equilibrium
          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_peat = 0.0_kind_noahmp
             else
                d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
             endif
             d_bot_peat = abs(DepthSoilLayer(LoopInd1))
             SoilLiqWater(LoopInd1) = EquilibriumSMMicroTopo(d_top_peat, d_bot_peat, &
                 WaterTableDepth, thetas_peat, ae_peat, bb_peat)
          enddo

          ! Normalize equilibrium profile to match W_soil_peat exactly
          W_soil_check = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             W_soil_check = W_soil_check + &
                 SoilLiqWater(LoopInd1) * abs(ThicknessSnowSoilLayer(LoopInd1))
          enddo
          if (abs(W_soil_check) > 1.0e-12_kind_noahmp) then
             do LoopInd1 = 1, NumSoilLayer
                SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) * (W_soil_peat / W_soil_check)
             enddo
          endif

          ! FSW_change: exact from SurfaceWaterStorage [mm]
          FSW_change = (SurfaceWaterStorage(z_wt_end) - W_fsw_begin) * 1000.0_kind_noahmp

          ! Flux-based diagnostic reference
          FSW_change_flux = InfilRateSfc_FSW_change * SoilTimeStep * 1000.0_kind_noahmp &
                          - (1.0_kind_noahmp - f_soil) * EvapGroundNet * SoilTimeStep   &
                          - (1.0_kind_noahmp - f_soil) * Transpiration * SoilTimeStep   &
                          - (1.0_kind_noahmp - f_soil) * RunoffSubsurface * SoilTimeStep

          ! Peatland numerical water balance error diagnostic
          FSW_peat_error = FSW_change - FSW_change_flux

          ! Update FloodedFraction from final WTD
          FloodedFraction = FloodedFrac(z_wt_end)

          ! No Richards iterations needed
          DrainSoilBot = 0.0
          ! Route any pre-clipped saturation excess as surface runoff [mm/s]
          RunoffSurface = SoilSatExcAcc * 1000.0 / SoilTimeStep
          write(*,*) 'DEBUG: RunoffSurface1 = ', RunoffSurface

       ! ================================================================
       ! RICHARDS PATH: WTD >= 0.3 m (Deep WT, full Richards)
       ! ================================================================
       else

          ! --- Forward transfer: map disequilibrium in head space ---
          if (.not. allocated(SoilLiqWaterOrig))   allocate(SoilLiqWaterOrig(1:NumSoilLayer))
          if (.not. allocated(SoilLiqWater1D_bef)) allocate(SoilLiqWater1D_bef(1:NumSoilLayer))
          if (.not. allocated(HeadShiftMicro))     allocate(HeadShiftMicro(1:NumSoilLayer))
          if (.not. allocated(HeadShiftFlat))      allocate(HeadShiftFlat(1:NumSoilLayer))

          HeadShiftMicro(:) = 0.0_kind_noahmp
          HeadShiftFlat(:)  = 0.0_kind_noahmp
          head_shift_cap    = 4.0_kind_noahmp * ae_peat

          SatTopInd = NumSoilLayer + 1
          do LoopInd1 = NumSoilLayer, 2, -1
             if ( abs(DepthSoilLayer(LoopInd1-1)) + ae_peat >= WaterTableDepth ) then
                SatTopInd = LoopInd1
             else
                exit
             endif
          enddo
          if ( SatTopInd == 2 .and. WaterTableDepth <= 0.0_kind_noahmp ) then
             SatTopInd = 1
          endif

          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_peat = 0.0_kind_noahmp
             else
                d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
             endif
             d_bot_peat = abs(DepthSoilLayer(LoopInd1))

             ! Store original column-averaged profile
             SoilLiqWaterOrig(LoopInd1) = SoilLiqWater(LoopInd1)

             if ( LoopInd1 >= SatTopInd .or. SoilLiqWaterOrig(LoopInd1) >= SoilEffPorosity(LoopInd1) - 1.0e-8_kind_noahmp ) then
                HeadShiftMicro(LoopInd1) = 0.0_kind_noahmp
             else
                head_shift_tmp = HeadShiftFromThetaMicro(SoilLiqWaterOrig(LoopInd1), d_top_peat, d_bot_peat, &
                    WaterTableDepth, thetas_peat, ae_peat, bb_peat)
                HeadShiftMicro(LoopInd1) = max(-head_shift_cap, min(head_shift_cap, head_shift_tmp))
             endif

             SoilLiqWater(LoopInd1) = ThetaFromHeadShiftFlat(d_top_peat, d_bot_peat, WaterTableDepth, &
                 HeadShiftMicro(LoopInd1), thetas_peat, ae_peat, bb_peat)
             SoilLiqWater(LoopInd1) = max(0.001_kind_noahmp, min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
          enddo

         ! Debug: Write out SoilLiqWaterOrig and z_wt_begin
         write(*,*) 'DEBUG: z_wt_begin = ', z_wt_begin
         write(*,*) 'DEBUG: SoilLiqWaterOrig:'
         do LoopInd1 = 1, NumSoilLayer
            write(*,*) '  Layer', LoopInd1, SoilLiqWaterOrig(LoopInd1)
         enddo

         ! Redistribute water that exceeds porosity to neighbors
         do LoopInd1 = 1, NumSoilLayer
            if (SoilLiqWater(LoopInd1) > SoilEffPorosity(LoopInd1)) then

               excess_vol = (SoilLiqWater(LoopInd1) - SoilEffPorosity(LoopInd1)) * &
                            abs(ThicknessSnowSoilLayer(LoopInd1))

               SoilLiqWater(LoopInd1) = SoilEffPorosity(LoopInd1)

               ! Try all layers above, starting from the nearest one
               if (LoopInd1 > 1) then
                  LoopJ = LoopInd1 - 1
                  do while (LoopJ >= 1 .and. excess_vol > 0.0_kind_noahmp)
                     space_avail = (SoilEffPorosity(LoopJ) - SoilLiqWater(LoopJ)) * &
                                   abs(ThicknessSnowSoilLayer(LoopJ))
                     if (space_avail > 0.0_kind_noahmp) then
                        transfer_vol = min(excess_vol, space_avail)
                        SoilLiqWater(LoopJ) = SoilLiqWater(LoopJ) + &
                            transfer_vol / abs(ThicknessSnowSoilLayer(LoopJ))
                        excess_vol = excess_vol - transfer_vol
                     endif
                     LoopJ = LoopJ - 1
                  enddo
               endif

               ! Try layers below
               if (excess_vol > 0.0_kind_noahmp .and. LoopInd1 < NumSoilLayer) then
                  LoopJ = LoopInd1 + 1
                  do while (LoopJ <= NumSoilLayer .and. excess_vol > 0.0_kind_noahmp)
                     space_avail = (SoilEffPorosity(LoopJ) - SoilLiqWater(LoopJ)) * &
                                   abs(ThicknessSnowSoilLayer(LoopJ))
                     if (space_avail > 0.0_kind_noahmp) then
                        transfer_vol = min(excess_vol, space_avail)
                        SoilLiqWater(LoopJ) = SoilLiqWater(LoopJ) + &
                            transfer_vol / abs(ThicknessSnowSoilLayer(LoopJ))
                        excess_vol = excess_vol - transfer_vol
                     endif
                     LoopJ = LoopJ + 1
                  enddo
               endif

            endif
         enddo

          ! Store 1D profile before Richards for delta computation
          do LoopInd1 = 1, NumSoilLayer
             SoilLiqWater1D_bef(LoopInd1) = SoilLiqWater(LoopInd1)
          enddo

          ! --- Determine iteration times ---
          NumIterSoilWat = 3
          if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
             NumIterSoilWat = NumIterSoilWat*2
          endif
          TimeStepFine = SoilTimeStep / NumIterSoilWat

          ! --- Solve soil moisture via Richards ---
          InfilSfcAcc      = 1.0e-06
          DrainSoilBotAcc  = 0.0
          RunoffSurfaceAcc = 0.0

          do IndIter = 1, NumIterSoilWat
             if ( SoilSfcInflowMean > 0.0 ) then
                if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,TimeStepFine)
                if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,TimeStepFine)
                if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,TimeStepFine)
                if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,TimeStepFine,InfilSfcAcc)
             endif
             write(*,*) 'SoilLiqWater, IndIter (before): ',IndIter, SoilLiqWater
             call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
             call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
             write(*,*) 'SoilLiqWater, IndIter (after): ',IndIter, SoilLiqWater
                 write(*,*) 'DEBUG: SoilSatExcAccBef = ', SoilSatExcAcc
             SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess
                 write(*,*) 'DEBUG: SoilSatExcAccAfter = ', SoilSatExcAcc
             DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot
             RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface
          enddo

          DrainSoilBot  = DrainSoilBotAcc / NumIterSoilWat
          RunoffSurface = RunoffSurfaceAcc / NumIterSoilWat
          write(*,*) 'DEBUG: RunoffSurface3 = ', RunoffSurface
          RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep
          write(*,*) 'DEBUG: RunoffSurface4 = ', RunoffSurface

          ! --- Remove f_soil fraction of subsurface runoff from soil ---
          ! Exclude fully-saturated layers (below WT + capillary fringe) from K-weighted removal
          SatTopInd = NumSoilLayer + 1
          do LoopInd1 = NumSoilLayer, 2, -1
             if ( abs(DepthSoilLayer(LoopInd1-1)) + ae_peat >= WaterTableDepth ) then
                SatTopInd = LoopInd1
             else
                exit
             endif
          enddo
          if ( SatTopInd == 2 .and. WaterTableDepth <= 0.0_kind_noahmp ) then
             SatTopInd = 1
          endif

          SoilWatConductAcc = 0.0
          do LoopInd1 = 1, min(SatTopInd - 1, NumSoilLayer)
             SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
          enddo
          if (SoilWatConductAcc > 0.0) then
             do LoopInd1 = 1, min(SatTopInd - 1, NumSoilLayer)
                WaterRemove = f_soil * RunoffSubsurface * SoilTimeStep * &
                             (SoilWatConductivity(LoopInd1)*ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
                SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1)*1000.0)
                write(*,*) 'DEBUG PEAT WaterRemove [mm] Layer', LoopInd1, WaterRemove
             enddo
          endif

          ! --- Backward transfer target: conserved soil water from Richards path ---

          ! Flux-based FSW_change diagnostic reference
          FSW_change_flux = InfilRateSfc_FSW_change * SoilTimeStep * 1000.0_kind_noahmp &
                          - (1.0_kind_noahmp - f_soil) * EvapGroundNet * SoilTimeStep   &
                          - (1.0_kind_noahmp - f_soil) * Transpiration * SoilTimeStep   &
                          - (1.0_kind_noahmp - f_soil) * RunoffSubsurface * SoilTimeStep

         ! Debug: Write out all fluxes used in FSW_change_flux
         write(*,*) 'DEBUG: FSW_change_flux terms:'
         write(*,*) '  InfilRateSfc_FSW_change = ', InfilRateSfc_FSW_change
         write(*,*) '  SoilTimeStep = ', SoilTimeStep
         write(*,*) '  f_soil = ', f_soil
         write(*,*) '  EvapGroundNet = ', EvapGroundNet
         write(*,*) '  Transpiration = ', Transpiration
         write(*,*) '  RunoffSubsurface = ', RunoffSubsurface

          ! Conserved soil water target after Richards and runoff removal [m]
          W_target_peat = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             W_target_peat = W_target_peat + &
                 (SoilLiqWaterOrig(LoopInd1) + &
                  (SoilLiqWater(LoopInd1) - SoilLiqWater1D_bef(LoopInd1))) * &
                 abs(ThicknessSnowSoilLayer(LoopInd1))
          enddo

          ! Convert the post-Richards flat profile into equivalent head anomalies.
          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_peat = 0.0_kind_noahmp
             else
                d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
             endif
             d_bot_peat = abs(DepthSoilLayer(LoopInd1))

             if ( SoilLiqWater(LoopInd1) >= SoilEffPorosity(LoopInd1) - 1.0e-8_kind_noahmp ) then
                HeadShiftFlat(LoopInd1) = 0.0_kind_noahmp
             else
                head_shift_tmp = HeadShiftFromThetaFlat(SoilLiqWater(LoopInd1), d_top_peat, d_bot_peat, &
                    WTD_begin, thetas_peat, ae_peat, bb_peat)
                HeadShiftFlat(LoopInd1) = max(-head_shift_cap, min(head_shift_cap, head_shift_tmp))
             endif
          enddo

          ! Solve final WTD so the head-shifted microtopography profile
          ! matches the conserved soil water after the Richards step.
          lambda_peat = 1.0_kind_noahmp
          WTD_lo_peat = -z_trunc
          WTD_hi_peat = z_col_bot_peat + z_trunc + head_shift_cap

          W_lo_peat = 0.0_kind_noahmp
          W_hi_peat = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_peat = 0.0_kind_noahmp
             else
                d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
             endif
             d_bot_peat = abs(DepthSoilLayer(LoopInd1))
             W_lo_peat = W_lo_peat + ThetaFromHeadShiftMicro(d_top_peat, d_bot_peat, WTD_lo_peat, &
                 lambda_peat * HeadShiftFlat(LoopInd1), thetas_peat, ae_peat, bb_peat) * &
                 abs(ThicknessSnowSoilLayer(LoopInd1))
             W_hi_peat = W_hi_peat + ThetaFromHeadShiftMicro(d_top_peat, d_bot_peat, WTD_hi_peat, &
                 lambda_peat * HeadShiftFlat(LoopInd1), thetas_peat, ae_peat, bb_peat) * &
                 abs(ThicknessSnowSoilLayer(LoopInd1))
          enddo

          if ( W_target_peat >= W_lo_peat ) then
             WaterTableDepth = WTD_lo_peat
          else if ( W_target_peat <= W_hi_peat ) then
             WaterTableDepth = WTD_hi_peat
          else
             do IterWTD = 1, 50
                WTD_mid_peat = 0.5_kind_noahmp * (WTD_lo_peat + WTD_hi_peat)
                W_mid_peat = 0.0_kind_noahmp
                do LoopInd1 = 1, NumSoilLayer
                   if (LoopInd1 == 1) then
                      d_top_peat = 0.0_kind_noahmp
                   else
                      d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
                   endif
                   d_bot_peat = abs(DepthSoilLayer(LoopInd1))
                   W_mid_peat = W_mid_peat + ThetaFromHeadShiftMicro(d_top_peat, d_bot_peat, WTD_mid_peat, &
                       lambda_peat * HeadShiftFlat(LoopInd1), thetas_peat, ae_peat, bb_peat) * &
                       abs(ThicknessSnowSoilLayer(LoopInd1))
                enddo
                if ( abs(W_mid_peat - W_target_peat) < 1.0e-8_kind_noahmp ) exit
                if ( W_mid_peat > W_target_peat ) then
                   WTD_lo_peat = WTD_mid_peat
                else
                   WTD_hi_peat = WTD_mid_peat
                endif
             enddo
             WaterTableDepth = 0.5_kind_noahmp * (WTD_lo_peat + WTD_hi_peat)
          endif

          z_wt_end = -WaterTableDepth

          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_peat = 0.0_kind_noahmp
             else
                d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
             endif
             d_bot_peat = abs(DepthSoilLayer(LoopInd1))
             SoilLiqWater(LoopInd1) = ThetaFromHeadShiftMicro(d_top_peat, d_bot_peat, WaterTableDepth, &
                 lambda_peat * HeadShiftFlat(LoopInd1), thetas_peat, ae_peat, bb_peat)
          enddo

          W_soil_check = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             W_soil_check = W_soil_check + SoilLiqWater(LoopInd1) * abs(ThicknessSnowSoilLayer(LoopInd1))
          enddo

          ! --- Saturation overflow cascade ---
          do LoopInd1 = 1, NumSoilLayer
             if (SoilLiqWater(LoopInd1) > SoilEffPorosity(LoopInd1)) then
                excess_vol = (SoilLiqWater(LoopInd1) - SoilEffPorosity(LoopInd1)) * &
                             abs(ThicknessSnowSoilLayer(LoopInd1))
                SoilLiqWater(LoopInd1) = SoilEffPorosity(LoopInd1)
                ! Try layer above
                if (LoopInd1 > 1 .and. SoilLiqWater(LoopInd1-1) < SoilEffPorosity(LoopInd1-1)) then
                   space_avail = (SoilEffPorosity(LoopInd1-1) - SoilLiqWater(LoopInd1-1)) * &
                                 abs(ThicknessSnowSoilLayer(LoopInd1-1))
                   transfer_vol = min(excess_vol, space_avail)
                   SoilLiqWater(LoopInd1-1) = SoilLiqWater(LoopInd1-1) + &
                       transfer_vol / abs(ThicknessSnowSoilLayer(LoopInd1-1))
                   excess_vol = excess_vol - transfer_vol
                endif
                ! Try layers below
                if (excess_vol > 0.0_kind_noahmp .and. LoopInd1 < NumSoilLayer) then
                   LoopJ = LoopInd1 + 1
                   do while (LoopJ <= NumSoilLayer .and. excess_vol > 0.0_kind_noahmp)
                      space_avail = (SoilEffPorosity(LoopJ) - SoilLiqWater(LoopJ)) * &
                                    abs(ThicknessSnowSoilLayer(LoopJ))
                      if (space_avail > 0.0_kind_noahmp) then
                         transfer_vol = min(excess_vol, space_avail)
                         SoilLiqWater(LoopJ) = SoilLiqWater(LoopJ) + &
                             transfer_vol / abs(ThicknessSnowSoilLayer(LoopJ))
                         excess_vol = excess_vol - transfer_vol
                      endif
                      LoopJ = LoopJ + 1
                   enddo
                endif
                ! Remaining excess → surface runoff
                if (excess_vol > 0.0_kind_noahmp) then
                   RunoffSurface = RunoffSurface + excess_vol * 1000.0 / SoilTimeStep
                   write(*,*) 'DEBUG: RunoffSurface5 = ', RunoffSurface
                endif
             endif
          enddo

          ! Safety floor
          do LoopInd1 = 1, NumSoilLayer
             SoilLiqWater(LoopInd1) = max(0.001_kind_noahmp, SoilLiqWater(LoopInd1))
          enddo

          ! FSW_change stays flux-based because the final microtopography
          ! profile is solved to the conserved Richards soil water target.
          FSW_change = FSW_change_flux

          FSW_peat_error = (W_soil_check - W_target_peat) * 1000.0_kind_noahmp

          ! Update FloodedFraction from final WTD
          FloodedFraction = FloodedFrac(z_wt_end)

          DrainSoilBot = DrainSoilBot * 1000.0  ! m/s -> mm/s

       endif   ! end Richards path

       ! ================================================================
       ! Common final steps: unit conversion, SoilMoisture, deallocation
       ! (FSW_change, FloodedFraction, WTD already set by each path)
       ! ================================================================

       ! Accumulated RunoffSurface and RunoffSubsurface [mm per soil timestep]
       RunoffSurface    = RunoffSurface    * SoilTimeStep
       write(*,*) 'DEBUG: RunoffSurface6 = ', RunoffSurface
       RunoffSubsurface = RunoffSubsurface * SoilTimeStep
       TileDrain        = 0.0

       ! Update soil moisture
       do LoopInd1 = 1, NumSoilLayer
           SoilMoisture(LoopInd1) = SoilLiqWater(LoopInd1) + SoilIce(LoopInd1)
       enddo

       ! Deallocate peatland local arrays
       if (allocated(SoilLiqWaterOrig))   deallocate(SoilLiqWaterOrig)
       if (allocated(SoilLiqWater1D_bef)) deallocate(SoilLiqWater1D_bef)
      if (allocated(HeadShiftMicro))     deallocate(HeadShiftMicro)
      if (allocated(HeadShiftFlat))      deallocate(HeadShiftFlat)
       if (allocated(SoilLiqGap))         deallocate(SoilLiqGap)

    else   ! non-peatland path

    ! surface runoff and infiltration rate
    if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
    if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
    if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
    if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
    if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)

    ! determine iteration times
    NumIterSoilWat = 3
    if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
       NumIterSoilWat = NumIterSoilWat*2
    endif
    TimeStepFine = SoilTimeStep / NumIterSoilWat

    ! solve soil moisture via Richards equation
    InfilSfcAcc      = 1.0e-06
    DrainSoilBotAcc  = 0.0
    RunoffSurfaceAcc = 0.0

    do IndIter = 1, NumIterSoilWat
       if ( SoilSfcInflowMean > 0.0 ) then
          if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,TimeStepFine,InfilSfcAcc)
       endif
       call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
       call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
       SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess
       DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot
       RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface
    enddo

    DrainSoilBot  = DrainSoilBotAcc / NumIterSoilWat
    RunoffSurface = RunoffSurfaceAcc / NumIterSoilWat
    RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep  ! m/s -> mm/s
    write(*,*) 'DEBUG: RunoffSurface8 = ', RunoffSurface
    DrainSoilBot  = DrainSoilBot * 1000.0  ! m/s -> mm/s

    ! compute tile drainage
    if ( (OptTileDrainage == 1) .and. (TileDrainFrac > 0.3) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageSimple(noahmp)
    endif
    if ( (OptTileDrainage == 2) .and. (TileDrainFrac > 0.1) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageHooghoudt(noahmp)
    END IF

    ! removal of soil water due to subsurface runoff (option 2)
    if ( OptRunoffSubsurface == 2 ) then
       SoilWatConductAcc = 0.0
       do LoopInd1 = 1, NumSoilLayer
          SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
       enddo
       do LoopInd1 = 1, NumSoilLayer
          WaterRemove            = RunoffSubsurface * SoilTimeStep * &
                                  (SoilWatConductivity(LoopInd1)*ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
          SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1)*1000.0)
       enddo
    endif

    ! Limit SoilLiqTmp safety
    if ( OptRunoffSubsurface /= 1 ) then
       do LoopInd2 = 1, NumSoilLayer
          SoilLiqTmp(LoopInd2) = SoilLiqWater(LoopInd2) * ThicknessSnowSoilLayer(LoopInd2) * 1000.0
       enddo

       SoilWaterMin = 0.01   ! mm
       do LoopInd2 = 1, NumSoilLayer-1
          if ( SoilLiqTmp(LoopInd2) < 0.0 ) then
             SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2)
          else
             SoilWatRem = 0.0
          endif
          SoilLiqTmp(LoopInd2  ) = SoilLiqTmp(LoopInd2  ) + SoilWatRem
          SoilLiqTmp(LoopInd2+1) = SoilLiqTmp(LoopInd2+1) - SoilWatRem
       enddo
       LoopInd2 = NumSoilLayer
       if ( SoilLiqTmp(LoopInd2) < SoilWaterMin ) then
           SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2)
       else
           SoilWatRem = 0.0
       endif
       SoilLiqTmp(LoopInd2) = SoilLiqTmp(LoopInd2) + SoilWatRem
       RunoffSubsurface     = RunoffSubsurface - SoilWatRem/SoilTimeStep

       if ( OptRunoffSubsurface == 5 ) RechargeGwDeepWT = RechargeGwDeepWT - SoilWatRem * 1.0e-3

       do LoopInd2 = 1, NumSoilLayer
          SoilLiqWater(LoopInd2) = SoilLiqTmp(LoopInd2) / (ThicknessSnowSoilLayer(LoopInd2)*1000.0)
       enddo
    endif

    ! compute groundwater and subsurface runoff
    if ( OptRunoffSubsurface == 1 ) call RunoffSubSurfaceGroundWater(noahmp)

    ! compute subsurface runoff based on drainage rate
    if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 4) .or. (OptRunoffSubsurface == 6) .or. &
         (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
         call RunoffSubSurfaceDrainage(noahmp)
    endif

    ! update soil moisture
    do LoopInd2 = 1, NumSoilLayer
        SoilMoisture(LoopInd2) = SoilLiqWater(LoopInd2) + SoilIce(LoopInd2)
    enddo

    ! compute subsurface runoff and shallow water table for MMF scheme
    if ( OptRunoffSubsurface == 5 ) call RunoffSubSurfaceShallowWaterMMF(noahmp)

    ! accumulated water flux over soil timestep [mm]
    RunoffSurface    = RunoffSurface    * SoilTimeStep
    RunoffSubsurface = RunoffSubsurface * SoilTimeStep
    TileDrain        = TileDrain        * SoilTimeStep

    ! deallocate local arrays to avoid memory leaks
    deallocate(MatRight  )
    deallocate(MatLeft1  )
    deallocate(MatLeft2  )
    deallocate(MatLeft3  )
    deallocate(SoilLiqTmp)

    endif   ! end if (OptPeatlandPhysics == 1) / else

    end associate

  end subroutine SoilWaterMain

end module SoilWaterMainMod
