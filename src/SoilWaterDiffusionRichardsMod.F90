module SoilWaterDiffusionRichardsMod

!!! Solve Richards equation for soil water movement/diffusion
!!! Compute the right hand side of the time tendency term of the soil
!!! water diffusion equation.  also to compute (prepare) the matrix
!!! coefficients for the tri-diagonal matrix of the implicit time scheme.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilHydraulicPropertyMod


  implicit none

contains

  subroutine SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SRT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft1     ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft2     ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft3     ! left-hand side term of the matrix

! local variable
    integer                                           :: LoopInd                     ! loop index
    integer                                           :: SatTopInd                   ! topmost fully-saturated layer index
    integer                                           :: TransInd                    ! transitional layer index (contains WT)
    real(kind=kind_noahmp)                            :: DepthSnowSoilTmp            ! temporary snow/soil layer depth [m]
    real(kind=kind_noahmp)                            :: SoilMoistTmpToWT            ! temporary soil moisture between bottom of the soil and water table
    real(kind=kind_noahmp)                            :: SoilMoistBotTmp             ! temporary soil moisture below bottom to calculate flux
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepthSnowSoilInv            ! inverse of snow/soil layer depth [1/m]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilThickTmp                ! temporary soil thickness
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilWaterGrad               ! temporary soil moisture vertical gradient
    real(kind=kind_noahmp), allocatable, dimension(:) :: WaterExcess                 ! temporary excess water flux
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilMoistureTmp             ! temporary soil moisture

! --------------------------------------------------------------------
    associate(                                                                             &
              NumSoilLayer              => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of soil layers
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,& ! in,  depth [m] of layer-bottom from soil surface
              OptSoilPermeabilityFrozen => noahmp%config%nmlist%OptSoilPermeabilityFrozen ,& ! in,  options for frozen soil permeability
              OptRunoffSubsurface       => noahmp%config%nmlist%OptRunoffSubsurface       ,& ! in,  options for drainage and subsurface runoff
              OptPeatlandPhysics        => noahmp%config%nmlist%OptPeatlandPhysics        ,& ! in,  options for peatland physics
              SoilDrainSlope            => noahmp%water%param%SoilDrainSlope              ,& ! in,  slope index for soil drainage
              SoilMatPotentialSat       => noahmp%water%param%SoilMatPotentialSat          ,& ! in,  saturated matric potential [m]
              InfilRateSfc              => noahmp%water%flux%InfilRateSfc                 ,& ! in,  infiltration rate at surface [m/s]
              EvapSoilSfcLiqMean        => noahmp%water%flux%EvapSoilSfcLiqMean           ,& ! in,  mean evaporation from soil surface [m/s]
              TranspWatLossSoilMean     => noahmp%water%flux%TranspWatLossSoilMean        ,& ! in,  mean transpiration water loss from soil layers [m/s]
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,& ! in,  soil water content [m3/m3]
              SoilMoisture              => noahmp%water%state%SoilMoisture                ,& ! in,  total soil moisture [m3/m3]
              WaterTableDepth           => noahmp%water%state%WaterTableDepth             ,& ! in,  water table depth [m]
              SoilImpervFrac            => noahmp%water%state%SoilImpervFrac              ,& ! in,  fraction of imperviousness due to frozen soil
              SoilImpervFracMax         => noahmp%water%state%SoilImpervFracMax           ,& ! in,  maximum soil imperviousness fraction
              SoilIceMax                => noahmp%water%state%SoilIceMax                  ,& ! in,  maximum soil ice content [m3/m3]
              SoilMoistureToWT          => noahmp%water%state%SoilMoistureToWT            ,& ! in,  soil moisture between bottom of the soil and the water table
              SoilWatConductivity       => noahmp%water%state%SoilWatConductivity         ,& ! out, soil hydraulic conductivity [m/s]
              SoilWatDiffusivity        => noahmp%water%state%SoilWatDiffusivity          ,& ! out, soil water diffusivity [m2/s]
              FSW_change                => noahmp%water%state%FSW_change                  ,& ! inout,   surface storage change [mm]
              f_soil                    => noahmp%water%state%f_soil                      ,& ! in, fraction of flux to soil [-] (Sy_soil/Sy_total)
              DrainSoilBot              => noahmp%water%flux%DrainSoilBot                  & ! out, soil bottom drainage [m/s]
             )
! ----------------------------------------------------------------------

    ! initialization
    if (.not. allocated(DepthSnowSoilInv)) allocate(DepthSnowSoilInv(1:NumSoilLayer))
    if (.not. allocated(SoilThickTmp)    ) allocate(SoilThickTmp    (1:NumSoilLayer))
    if (.not. allocated(SoilWaterGrad)   ) allocate(SoilWaterGrad   (1:NumSoilLayer))
    if (.not. allocated(WaterExcess)     ) allocate(WaterExcess     (1:NumSoilLayer))
    if (.not. allocated(SoilMoistureTmp) ) allocate(SoilMoistureTmp (1:NumSoilLayer))
    MatRight(:)         = 0.0
    MatLeft1(:)         = 0.0
    MatLeft2(:)         = 0.0
    MatLeft3(:)         = 0.0
    DepthSnowSoilInv(:) = 0.0
    SoilThickTmp(:)     = 0.0
    SoilWaterGrad(:)    = 0.0
    WaterExcess(:)      = 0.0
    SoilMoistureTmp(:)  = 0.0

    ! compute soil hydraulic conductivity and diffusivity
    if ( OptSoilPermeabilityFrozen == 1 ) then
       do LoopInd = 1, NumSoilLayer
          call SoilDiffusivityConductivityOpt1(noahmp,SoilWatDiffusivity(LoopInd),SoilWatConductivity(LoopInd),&
                                               SoilMoisture(LoopInd),SoilImpervFrac(LoopInd),LoopInd) 
          SoilMoistureTmp(LoopInd) = SoilMoisture(LoopInd)
       enddo
       if ( OptRunoffSubsurface == 5 ) SoilMoistTmpToWT = SoilMoistureToWT
    endif

    if ( OptSoilPermeabilityFrozen == 2 ) then
       do LoopInd = 1, NumSoilLayer
          call SoilDiffusivityConductivityOpt2(noahmp,SoilWatDiffusivity(LoopInd),SoilWatConductivity(LoopInd),&
                                               SoilLiqWater(LoopInd),SoilIceMax,LoopInd)
          SoilMoistureTmp(LoopInd) = SoilLiqWater(LoopInd)
       enddo
       if ( OptRunoffSubsurface == 5 ) &
          SoilMoistTmpToWT = SoilMoistureToWT * SoilLiqWater(NumSoilLayer) / SoilMoisture(NumSoilLayer)  !same liquid fraction as in the bottom layer
    endif

    ! compute gradient and flux of soil water diffusion terms
    do LoopInd = 1, NumSoilLayer
       if ( LoopInd == 1 ) then
          SoilThickTmp(LoopInd)     = - DepthSoilLayer(LoopInd)
          DepthSnowSoilTmp          = - DepthSoilLayer(LoopInd+1)
          DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
          SoilWaterGrad(LoopInd)    = 2.0 * (SoilMoistureTmp(LoopInd)-SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp
          WaterExcess(LoopInd)      = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                      InfilRateSfc + TranspWatLossSoilMean(LoopInd) + EvapSoilSfcLiqMean
          !if (OptRunoffSubsurface == 9) then
          if ( OptPeatlandPhysics == 1 ) then
             if (f_soil < 0.000001) then
                WaterExcess(LoopInd) = 0.0
             else
                WaterExcess(LoopInd) = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                       InfilRateSfc + f_soil*TranspWatLossSoilMean(LoopInd) + f_soil*EvapSoilSfcLiqMean
             endif
          endif
       else if ( LoopInd < NumSoilLayer ) then
          SoilThickTmp(LoopInd)     = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd))
          DepthSnowSoilTmp          = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd+1))
          DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
          SoilWaterGrad(LoopInd)    = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp
          WaterExcess(LoopInd)      = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                      SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1) - SoilWatConductivity(LoopInd-1) + &
                                      TranspWatLossSoilMean(LoopInd)
          !if (OptRunoffSubsurface == 9) then
          if ( OptPeatlandPhysics == 1 ) then
             if (f_soil < 0.000001) then
                WaterExcess(LoopInd) = 0.0
             else
                WaterExcess(LoopInd) = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                       SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1) - SoilWatConductivity(LoopInd-1) + &
                                       f_soil*TranspWatLossSoilMean(LoopInd)
             endif
          endif
       else
          SoilThickTmp(LoopInd) = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd))
          ! MB: For peatlands we don't want to lose water through the bottom ... instead it should raise the water level
          ! using the equilibrium approach that is also used in RunoffSubsurfaceOption 2
          !if ( (OptRunoffSubsurface == 1) .or. (OptRunoffSubsurface == 2) .or. (OptRunoffSubsurface == 9)) then
          if ( (OptRunoffSubsurface == 1) .or. (OptRunoffSubsurface == 2) .or. (OptPeatlandPhysics == 1)) then
             DrainSoilBot = 0.0
          endif
          if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 6) .or. &
               (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
             DrainSoilBot = SoilDrainSlope * SoilWatConductivity(LoopInd)
          endif
          if ( OptRunoffSubsurface == 4 ) then
             DrainSoilBot = (1.0 - SoilImpervFracMax) * SoilWatConductivity(LoopInd)
          endif
          if ( OptRunoffSubsurface == 5 ) then   ! gmm new m-m&f water table dynamics formulation
             DepthSnowSoilTmp  = 2.0 * SoilThickTmp(LoopInd)
             if ( WaterTableDepth < (DepthSoilLayer(NumSoilLayer)-SoilThickTmp(NumSoilLayer)) ) then
                ! gmm interpolate from below, midway to the water table, 
                ! to the middle of the auxiliary layer below the soil bottom
                SoilMoistBotTmp = SoilMoistureTmp(LoopInd) - (SoilMoistureTmp(LoopInd)-SoilMoistTmpToWT) * &
                                  SoilThickTmp(LoopInd)*2.0 / (SoilThickTmp(LoopInd)+DepthSoilLayer(LoopInd)-WaterTableDepth)
             else
                SoilMoistBotTmp = SoilMoistTmpToWT
             endif
             SoilWaterGrad(LoopInd) = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistBotTmp) / DepthSnowSoilTmp
             DrainSoilBot           = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd)
          endif
          WaterExcess(LoopInd) = -(SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1)) - SoilWatConductivity(LoopInd-1) + &
                                 TranspWatLossSoilMean(LoopInd) + DrainSoilBot
          !if (OptRunoffSubsurface == 9) then
          if ( OptPeatlandPhysics == 1 ) then
             if (f_soil < 0.000001) then
                WaterExcess(LoopInd) = 0.0
             else
                WaterExcess(LoopInd) = -(SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1)) &
                                 - SoilWatConductivity(LoopInd-1) + &
                                 f_soil*TranspWatLossSoilMean(LoopInd) + DrainSoilBot
             endif
          endif
       endif
    enddo

    ! prepare the matrix coefficients for the tri-diagonal matrix
    do LoopInd = 1, NumSoilLayer
       if ( LoopInd == 1 ) then
          MatLeft1(LoopInd) =   0.0
          MatLeft2(LoopInd) =   SoilWatDiffusivity(LoopInd  ) * DepthSnowSoilInv(LoopInd  ) / SoilThickTmp(LoopInd)
          MatLeft3(LoopInd) = - MatLeft2(LoopInd)
       else if ( LoopInd < NumSoilLayer ) then
          MatLeft1(LoopInd) = - SoilWatDiffusivity(LoopInd-1) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd)
          MatLeft3(LoopInd) = - SoilWatDiffusivity(LoopInd  ) * DepthSnowSoilInv(LoopInd  ) / SoilThickTmp(LoopInd)
          MatLeft2(LoopInd) = - (MatLeft1(LoopInd) + MatLeft3(LoopInd))
       else
          MatLeft1(LoopInd) = - SoilWatDiffusivity(LoopInd-1) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd)
          MatLeft3(LoopInd) =   0.0
          MatLeft2(LoopInd) = - (MatLeft1(LoopInd) + MatLeft3(LoopInd))
       endif
       MatRight(LoopInd) = WaterExcess(LoopInd) / (-SoilThickTmp(LoopInd))
    enddo

    ! Peatland: decouple inactive layers from Richards domain.
    ! A layer is inactive (decoupled) unless the WT + capillary fringe
    ! is entirely below the layer bottom. This prevents the theta-based
    ! discretization from generating spurious gravitational flux across
    ! the poorly-resolved WT boundary.
    if ( OptPeatlandPhysics == 1 ) then
       ! Find topmost inactive layer using layer-bottom criterion.
       SatTopInd = NumSoilLayer + 1
       do LoopInd = NumSoilLayer, 1, -1
          if ( abs(DepthSoilLayer(LoopInd)) + abs(SoilMatPotentialSat(1)) >= WaterTableDepth ) then
             SatTopInd = LoopInd
          else
             exit
          endif
       enddo

       if ( SatTopInd <= NumSoilLayer ) then
          TransInd = SatTopInd - 1

          ! --- Transitional layer: remove bottom flux, keep top flux + sinks ---
          if ( TransInd >= 1 ) then
             if ( f_soil < 0.000001_kind_noahmp ) then
                WaterExcess(TransInd) = 0.0
             else if ( TransInd == 1 ) then
                WaterExcess(TransInd) = -InfilRateSfc &
                    + f_soil*TranspWatLossSoilMean(TransInd) + f_soil*EvapSoilSfcLiqMean
             else
                WaterExcess(TransInd) = -(SoilWatDiffusivity(TransInd-1)*SoilWaterGrad(TransInd-1)) &
                    - SoilWatConductivity(TransInd-1) &
                    + f_soil*TranspWatLossSoilMean(TransInd)
             endif
             MatRight(TransInd) = WaterExcess(TransInd) / (-SoilThickTmp(TransInd))
             MatLeft3(TransInd) = 0.0
             MatLeft2(TransInd) = -(MatLeft1(TransInd) + MatLeft3(TransInd))
          endif

          ! --- Saturated layers: completely inert (no flux, no transpiration) ---
          do LoopInd = SatTopInd, NumSoilLayer
             WaterExcess(LoopInd) = 0.0
             MatRight(LoopInd)    = 0.0
             MatLeft1(LoopInd)    = 0.0
             MatLeft2(LoopInd)    = 0.0
             MatLeft3(LoopInd)    = 0.0
          enddo
       endif
    endif

    ! deallocate local arrays to avoid memory leaks
    deallocate(DepthSnowSoilInv)
    deallocate(SoilThickTmp    )
    deallocate(SoilWaterGrad   )
    deallocate(WaterExcess     )
    deallocate(SoilMoistureTmp )

    end associate

  end subroutine SoilWaterDiffusionRichards

end module SoilWaterDiffusionRichardsMod
