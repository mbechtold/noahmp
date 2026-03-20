module SoilWaterDiffusionRichardsMod

!!! Solve Richards equation for soil water movement/diffusion
!!! Compute the right hand side of the time tendency term of the soil
!!! water diffusion equation.  also to compute (prepare) the matrix
!!! coefficients for the tri-diagonal matrix of the implicit time scheme.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilHydraulicPropertyMod
  use MicroTopoCorrectionMod, only: SoilFractionAtDepth

  implicit none

contains

  subroutine SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SRT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Peatland microtopo: Bechtold (2025)
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
    real(kind=kind_noahmp)                            :: DepthSnowSoilTmp            ! temporary snow/soil layer depth [m]
    real(kind=kind_noahmp)                            :: SoilMoistTmpToWT            ! temporary soil moisture between bottom of the soil and water table
    real(kind=kind_noahmp)                            :: SoilMoistBotTmp             ! temporary soil moisture below bottom to calculate flux
    real(kind=kind_noahmp)                            :: theta_hyd_tmp               ! per-soil-area moisture for peatland hydraulic properties
    real(kind=kind_noahmp)                            :: f_soil_iface                ! soil fraction at interface for peatland
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
              FSW_change                => noahmp%water%state%FSW_change                  ,& ! inout, surface storage change [mm]
              f_part                    => noahmp%water%state%f_part                      ,& ! inout, soil flux partition fraction
              FloodedFrac               => noahmp%water%state%FloodedFrac                 ,& ! inout, flooded fraction
              f_soil_k                  => noahmp%water%state%f_soil_k                    ,& ! in,  per-layer soil fraction
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
          SoilMoistTmpToWT = SoilMoistureToWT * SoilLiqWater(NumSoilLayer) / SoilMoisture(NumSoilLayer)
    endif

    ! Peatland: recompute hydraulic properties at per-soil-area moisture
    ! and scale D, K by f_soil at layer interfaces
    if ( OptPeatlandPhysics == 1 ) then
       do LoopInd = 1, NumSoilLayer
          theta_hyd_tmp = SoilMoistureTmp(LoopInd) / max(f_soil_k(LoopInd), 0.01_kind_noahmp)
          theta_hyd_tmp = min(theta_hyd_tmp, noahmp%water%param%SoilMoistureSat(LoopInd))
          if ( OptSoilPermeabilityFrozen == 1 ) then
             call SoilDiffusivityConductivityOpt1(noahmp,SoilWatDiffusivity(LoopInd),SoilWatConductivity(LoopInd),&
                                                  theta_hyd_tmp,SoilImpervFrac(LoopInd),LoopInd)
          else
             call SoilDiffusivityConductivityOpt2(noahmp,SoilWatDiffusivity(LoopInd),SoilWatConductivity(LoopInd),&
                                                  theta_hyd_tmp,SoilIceMax,LoopInd)
          endif
          ! Scale by f_soil at the interface (bottom of layer k)
          f_soil_iface = SoilFractionAtDepth(abs(DepthSoilLayer(LoopInd)))
          SoilWatDiffusivity(LoopInd)  = SoilWatDiffusivity(LoopInd) * f_soil_iface
          SoilWatConductivity(LoopInd) = SoilWatConductivity(LoopInd) * f_soil_iface
       enddo
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
          if ( OptPeatlandPhysics == 1 ) then
             ! D, K already scaled by f_soil_interface; source terms scaled by f_soil_k
             WaterExcess(LoopInd)   = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                      InfilRateSfc + f_soil_k(LoopInd)*TranspWatLossSoilMean(LoopInd) + &
                                      f_soil_k(LoopInd)*EvapSoilSfcLiqMean
          endif
       else if ( LoopInd < NumSoilLayer ) then
          SoilThickTmp(LoopInd)     = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd))
          DepthSnowSoilTmp          = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd+1))
          DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
          SoilWaterGrad(LoopInd)    = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp
          WaterExcess(LoopInd)      = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                      SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1) - SoilWatConductivity(LoopInd-1) + &
                                      TranspWatLossSoilMean(LoopInd)
          if ( OptPeatlandPhysics == 1 ) then
             WaterExcess(LoopInd)   = SoilWatDiffusivity(LoopInd)*SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - &
                                      SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1) - SoilWatConductivity(LoopInd-1) + &
                                      f_soil_k(LoopInd)*TranspWatLossSoilMean(LoopInd)
          endif
       else
          SoilThickTmp(LoopInd) = (DepthSoilLayer(LoopInd-1) - DepthSoilLayer(LoopInd))
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
          if ( OptRunoffSubsurface == 5 ) then
             DepthSnowSoilTmp  = 2.0 * SoilThickTmp(LoopInd)
             if ( WaterTableDepth < (DepthSoilLayer(NumSoilLayer)-SoilThickTmp(NumSoilLayer)) ) then
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
          if ( OptPeatlandPhysics == 1 ) then
             WaterExcess(LoopInd) = -(SoilWatDiffusivity(LoopInd-1)*SoilWaterGrad(LoopInd-1)) - SoilWatConductivity(LoopInd-1) + &
                                    f_soil_k(LoopInd)*TranspWatLossSoilMean(LoopInd) + DrainSoilBot
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

    ! deallocate local arrays to avoid memory leaks
    deallocate(DepthSnowSoilInv)
    deallocate(SoilThickTmp    )
    deallocate(SoilWaterGrad   )
    deallocate(WaterExcess     )
    deallocate(SoilMoistureTmp )

    end associate

  end subroutine SoilWaterDiffusionRichards

end module SoilWaterDiffusionRichardsMod
