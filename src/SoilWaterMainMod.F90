module SoilWaterMainMod

!!! Main soil water module including all soil water processes & update soil moisture
!!! surface runoff, infiltration, soil water diffusion, subsurface runoff, tile drainage
!!! Peatland Gaussian microtopography (Dettmann & Bechtold 2015) with multi-column
!!! ensemble, Sy-weighted flux partitioning, and post-Richards equilibration.

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
  use MicroTopoCorrectionMod,            only : MicroTopoCorrection, SurfaceWaterStorage_mm, &
                                                 MicroTopoStorageDeficitExtended, &
                                                 MicroTopoSpecificYield, sigma_z
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics
  use SoilWaterDiffusionRichardsMod,     only : SoilWaterDiffusionRichards
  use SoilMoistureSolverMod,             only : SoilMoistureSolver
  use TileDrainageSimpleMod,             only : TileDrainageSimple
  use TileDrainageHooghoudtMod,          only : TileDrainageHooghoudt
  use WaterTableEquilibriumMod,          only : WaterTableEquilibrium
  use WaterTableEquilibriumPeatMod,      only : WaterTableEquilibriumPeat

  implicit none

contains

  subroutine SoilWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SOILWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Peatland Gaussian microtopo: Dettmann & Bechtold (2015), implemented 2025
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
    ! Peatland-specific local variables
    real(kind=kind_noahmp)            :: W_surface_begin              ! surface water storage at timestep start [mm]
    real(kind=kind_noahmp)            :: W_surface_end                ! surface water storage at timestep end [mm]
    real(kind=kind_noahmp)            :: D_soil_actual                ! actual soil moisture deficit [m]
    real(kind=kind_noahmp)            :: D_soil_equil                 ! equilibrium deficit at target WTD [m]
    real(kind=kind_noahmp)            :: WTD_soil                     ! WTD diagnosed from soil layers [m]
    real(kind=kind_noahmp)            :: WTD_surface                  ! WTD diagnosed from surface water [m]
    real(kind=kind_noahmp)            :: WTD_target                   ! equilibrated target WTD [m]
    real(kind=kind_noahmp)            :: W_total                      ! total storage (soil+surface) for conservation [mm]
    real(kind=kind_noahmp)            :: W_soil_actual                ! actual soil water [mm]
    real(kind=kind_noahmp)            :: W_soil_equil                 ! equilibrium soil water at WTD_target [mm]
    real(kind=kind_noahmp)            :: DeltaW_transfer              ! water transfer soil↔surface [mm]
    real(kind=kind_noahmp)            :: Sy_total                     ! total specific yield for distribution
    real(kind=kind_noahmp)            :: thetas, ae, bb               ! peat soil params
    real(kind=kind_noahmp)            :: zwt_lo, zwt_hi, zwt_mid      ! bisection variables
    real(kind=kind_noahmp)            :: func_lo, func_hi, func_mid   ! bisection function values
    real(kind=kind_noahmp)            :: InfilRateSfcTotal             ! total infiltration before partitioning [m/s]
    integer                           :: IterEquil                     ! equilibration iteration counter

    real(kind=kind_noahmp), parameter :: WaterTableDepthMinPeat = -0.2449_kind_noahmp
    real(kind=kind_noahmp), parameter :: SoilImpPara = 4.0            ! soil impervious fraction parameter
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft1     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft2     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft3     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqTmp   ! temporary soil liquid water [mm]

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    noahmp soil time step [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth of layer-bottom from surface [m]
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
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,    saturated matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,    Campbell b exponent
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
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! inout, water table depth [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax        ,& ! out,   maximum soil imperviousness fraction
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity      ,& ! out,   soil hydraulic conductivity [m/s]
              SoilEffPorosity        => noahmp%water%state%SoilEffPorosity          ,& ! out,   soil effective porosity [m3/m3]
              SoilImpervFrac         => noahmp%water%state%SoilImpervFrac           ,& ! out,   impervious fraction due to frozen soil
              SoilIceFrac            => noahmp%water%state%SoilIceFrac              ,& ! out,   ice fraction in frozen soil
              SoilSaturationExcess   => noahmp%water%state%SoilSaturationExcess     ,& ! out,   saturation excess of the total soil [m]
              SoilIceMax             => noahmp%water%state%SoilIceMax               ,& ! out,   maximum soil ice content [m3/m3]
              FSW_change             => noahmp%water%state%FSW_change               ,& ! out,   surface water storage change [mm]
              FloodedFrac            => noahmp%water%state%FloodedFrac              ,& ! inout, flooded fraction [-]
              f_soil                 => noahmp%water%state%f_soil                   ,& ! inout, Sy-weighted flux partition fraction [-]
              SoilLiqWaterMin        => noahmp%water%state%SoilLiqWaterMin           & ! out,   minimum soil liquid water content [m3/m3]
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
    
    ! Peatland-specific subsurface runoff (Ivanov formulation, unchanged)
    if ( OptRunoffSubsurface == 9 ) then
            call RunoffSubSurfacePeatland(noahmp) 
    endif

    ! jref impermable surface at urban
    if ( FlagUrban .eqv. .true. ) SoilImpervFrac(1) = 0.95
    
    if ( OptPeatlandPhysics == 1 ) then
       WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)
    endif

    ! surface runoff and infiltration rate using different schemes
    if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
    if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
    if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
    if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
    if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)
    
    ! For peatlands: add RunoffSurface to infiltration (no overland flow)
    if ( OptPeatlandPhysics == 1 ) then
       InfilRateSfc = InfilRateSfc + RunoffSurface
       RunoffSurface = 0.0
    endif

    ! determine iteration times to solve soil water diffusion and moisture
    NumIterSoilWat = 3
    if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
       NumIterSoilWat = NumIterSoilWat*2
    endif
    TimeStepFine = SoilTimeStep / NumIterSoilWat

    ! solve soil moisture
    InfilSfcAcc      = 1.0e-06
    DrainSoilBotAcc  = 0.0
    RunoffSurfaceAcc = 0.0

    !=========================================================================
    ! PEATLAND BRANCH: Gaussian microtopo with Sy-weighted partitioning
    !=========================================================================
    if ( OptPeatlandPhysics == 1 ) then

       ! --- Step 1: Diagnose WTD from current soil moisture ---
       ! Must come BEFORE W_surface_begin so that WTD is consistent
       ! with the current soil moisture state (not stale from previous timestep).
       call WaterTableEquilibriumPeat(noahmp)
       WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

       ! --- Step 2: Save beginning-of-timestep surface water storage ---
       W_surface_begin = SurfaceWaterStorage_mm(WaterTableDepth, sigma_z)

       ! --- Step 3: Compute FloodedFrac and f_soil (= f_part) ---
       call MicroTopoCorrection(noahmp)

       ! --- Step 4-6: Partition fluxes between soil and surface water ---
       InfilRateSfcTotal = InfilRateSfc
       InfilRateSfc = f_soil * InfilRateSfcTotal     ! soil portion → enters Richards

       ! --- Step 7: Standard Richards (4 layers, unmodified) ---
       ! InfilRateSfc already contains only the soil portion
       ! Transpiration and EvapGroundNet are handled inside Richards as usual
       ! NOTE: RunoffSurface routines are NOT called inside the loop for peatlands
       ! because we handle infiltration partitioning ourselves (f_soil scaling).
       ! Calling them would overwrite InfilRateSfc with SoilSfcInflowMean.
       do IndIter = 1, NumIterSoilWat
          call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
          call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
          SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess
          DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot
          RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface
       enddo

       ! --- Step 8: Post-Richards equilibration ---
       ! Peat soil parameters
       thetas = SoilMoistureSat(1)
       ae     = abs(SoilMatPotentialSat(1))
       bb     = SoilExpCoeffB(1)

       ! Compute current soil storage [mm]
       W_soil_actual = 0.0_kind_noahmp
       do LoopInd1 = 1, NumSoilLayer
          W_soil_actual = W_soil_actual + SoilLiqWater(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1) * 1000.0_kind_noahmp
       enddo

       ! Surface water budget: add infiltration to surface, subtract ET and runoff shares
       ! W_surface after fluxes applied
       W_surface_end = W_surface_begin + &
                        (1.0_kind_noahmp - f_soil) * InfilRateSfcTotal * SoilTimeStep * 1000.0_kind_noahmp - &
                        (1.0_kind_noahmp - f_soil) * EvapGroundNet * SoilTimeStep - &
                        (1.0_kind_noahmp - f_soil) * Transpiration * SoilTimeStep
       W_surface_end = max(0.0_kind_noahmp, W_surface_end)

       ! Total storage for conservation
       W_total = W_soil_actual + W_surface_end

       ! Equilibrate: find WTD_target that is consistent with both domains
       ! Bisect for WTD_target where soil_storage(WTD) + surface_storage(WTD) = W_total
       do IterEquil = 1, 3
          ! Diagnose WTD_soil from soil layers
          call WaterTableEquilibriumPeat(noahmp)
          WTD_soil = max(WaterTableDepth, WaterTableDepthMinPeat)

          ! Diagnose WTD_surface from surface water via inverse of SurfaceWaterStorage_mm
          ! Bisect: find WTD where SurfaceWaterStorage_mm(WTD) = W_surface_end
          if (W_surface_end <= 0.0_kind_noahmp) then
             WTD_surface = 1.0_kind_noahmp   ! no surface water → deep WTD
          else
             zwt_lo = WaterTableDepthMinPeat
             zwt_hi = 1.0_kind_noahmp
             do LoopInd1 = 1, 60
                zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
                func_mid = SurfaceWaterStorage_mm(zwt_mid, sigma_z) - W_surface_end
                if (abs(func_mid) < 0.001_kind_noahmp .or. (zwt_hi - zwt_lo) < 1.0e-4_kind_noahmp) exit
                if (func_mid > 0.0_kind_noahmp) then
                   zwt_lo = zwt_mid   ! storage too high → deeper WTD
                else
                   zwt_hi = zwt_mid   ! storage too low → shallower WTD
                endif
             enddo
             WTD_surface = zwt_mid
          endif

          ! If soil and surface WTDs already agree, done
          if (abs(WTD_soil - WTD_surface) < 1.0e-4_kind_noahmp) exit

          ! Find WTD_target conserving total storage
          ! f(WTD) = soil_storage(WTD) + SurfaceWaterStorage(WTD) - W_total = 0
          ! soil_storage(WTD) = thetas*total_depth*1000 - deficit(WTD)*1000
          zwt_lo = WaterTableDepthMinPeat
          zwt_hi = 3.0_kind_noahmp * (-DepthSoilLayer(NumSoilLayer))
          do LoopInd1 = 1, 60
             zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
             ! Soil storage at this WTD
             D_soil_equil = MicroTopoStorageDeficitExtended(zwt_mid, NumSoilLayer, DepthSoilLayer, &
                                                             thetas, ae, bb, sigma_z)
             W_soil_equil = thetas * (-DepthSoilLayer(NumSoilLayer)) * 1000.0_kind_noahmp - &
                            D_soil_equil * 1000.0_kind_noahmp
             func_mid = W_soil_equil + SurfaceWaterStorage_mm(zwt_mid, sigma_z) - W_total
             if (abs(func_mid) < 0.01_kind_noahmp .or. (zwt_hi - zwt_lo) < 1.0e-4_kind_noahmp) exit
             if (func_mid > 0.0_kind_noahmp) then
                zwt_lo = zwt_mid   ! too much total storage → deeper WTD
             else
                zwt_hi = zwt_mid   ! too little → shallower WTD
             endif
          enddo
          WTD_target = max(zwt_mid, WaterTableDepthMinPeat)

          ! Transfer water between domains
          D_soil_equil = MicroTopoStorageDeficitExtended(WTD_target, NumSoilLayer, DepthSoilLayer, &
                                                          thetas, ae, bb, sigma_z)
          W_soil_equil = thetas * (-DepthSoilLayer(NumSoilLayer)) * 1000.0_kind_noahmp - &
                         D_soil_equil * 1000.0_kind_noahmp
          DeltaW_transfer = W_soil_equil - W_soil_actual   ! positive = surface→soil

          ! Distribute transfer across layers proportional to thickness
          ! (simplified: uniform distribution weighted by layer thickness)
          do LoopInd1 = 1, NumSoilLayer
             SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) + &
                (DeltaW_transfer / 1000.0_kind_noahmp) * &
                (ThicknessSnowSoilLayer(LoopInd1) / (-DepthSoilLayer(NumSoilLayer)))  / &
                ThicknessSnowSoilLayer(LoopInd1)
             ! Clamp to valid range
             SoilLiqWater(LoopInd1) = max(0.0_kind_noahmp, &
                                      min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
          enddo

          ! Update surface water
          W_surface_end = W_total - 0.0_kind_noahmp
          W_soil_actual = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             W_soil_actual = W_soil_actual + SoilLiqWater(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1) * 1000.0_kind_noahmp
          enddo
          W_surface_end = W_total - W_soil_actual
          W_surface_end = max(0.0_kind_noahmp, W_surface_end)

          WaterTableDepth = WTD_target
       enddo ! IterEquil

       ! --- Step 9: RunoffSubsurface already computed (Ivanov, capped at 0.0002) ---
       ! --- Step 10: Remove RunoffSubsurface from budget ---
       ! Split between soil and surface by f_soil
       SoilWatConductAcc = 0.0_kind_noahmp
       do LoopInd1 = 1, NumSoilLayer
          SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
       enddo
       if (SoilWatConductAcc > 0.0_kind_noahmp) then
          do LoopInd1 = 1, NumSoilLayer
             WaterRemove = f_soil * RunoffSubsurface * SoilTimeStep * &
                          (SoilWatConductivity(LoopInd1)*ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
             SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1)*1000.0_kind_noahmp)
             SoilLiqWater(LoopInd1) = max(0.0_kind_noahmp, SoilLiqWater(LoopInd1))
          enddo
       endif

       ! --- Step 11: Diagnose final WTD ---
       call WaterTableEquilibriumPeat(noahmp)
       WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

       ! --- Step 12: FSW_change from explicit surface water budget ---
       ! Remove (1-f_soil) share of RunoffSubsurface from surface water
       W_surface_end = W_surface_end - (1.0_kind_noahmp - f_soil) * RunoffSubsurface * SoilTimeStep
       W_surface_end = max(0.0_kind_noahmp, W_surface_end)
       FSW_change = W_surface_end - W_surface_begin

       ! Finalize Richards outputs
       DrainSoilBot  = DrainSoilBotAcc / NumIterSoilWat
       RunoffSurface = RunoffSurfaceAcc / NumIterSoilWat
       RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep
       DrainSoilBot  = DrainSoilBot * 1000.0

    !=========================================================================
    ! NON-PEATLAND BRANCH: standard Noah-MP soil water
    !=========================================================================
    else

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
       RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep
       DrainSoilBot  = DrainSoilBot * 1000.0

       ! compute tile drainage
       if ( (OptTileDrainage == 1) .and. (TileDrainFrac > 0.3) .and. (OptRunoffSurface == 3) ) then
          call TileDrainageSimple(noahmp)
       endif
       if ( (OptTileDrainage == 2) .and. (TileDrainFrac > 0.1) .and. (OptRunoffSurface == 3) ) then
          call TileDrainageHooghoudt(noahmp)
       endif

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

       ! Limit SoilLiqTmp to be greater than or equal to watmin
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

    endif ! OptPeatlandPhysics

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

    end associate

  end subroutine SoilWaterMain

end module SoilWaterMainMod
