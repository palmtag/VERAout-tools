   module mod_geom
      implicit none

!  hard-wired geometry information for Watts Bar Reactor
!  this information written to CORE block of HDF file

      integer, parameter :: kdfuel=49  ! number of axial levels
      integer, parameter :: iafull=15  ! number of assemblies across

      integer :: maxasm=56          ! number of assemblies quarter core
      integer :: npin=17            ! 17x17 pins across each assembly

      integer :: icoremap(iafull,iafull)

      real(8) :: apitch=21.50d0     ! assembly pitch (cm)
      real(8) :: prated=3411.0d0    ! rated power (MW) full-core
      real(8) :: crated=16591.3524d0 ! rated flow (kg/s) full-core  (131.68 Mlb/hr)

      real(8) :: zmesh(0:kdfuel)    ! axial mesh bounds (cm)

      data zmesh / &
           11.951d0,  15.817d0,  24.028d0,  32.239d0,  40.450d0, &
           48.662d0,  56.873d0,  65.084d0,  73.295d0,  77.105d0, &
           85.170d0,  93.235d0, 101.300d0, 109.365d0, 117.430d0, &
          125.495d0, 129.305d0, 137.370d0, 145.435d0, 153.500d0, &
          161.565d0, 169.630d0, 177.695d0, 181.505d0, 189.570d0, &
          197.635d0, 205.700d0, 213.765d0, 221.830d0, 229.895d0, &
          233.705d0, 241.770d0, 249.835d0, 257.900d0, 265.965d0, &
          274.030d0, 282.095d0, 285.905d0, 293.970d0, 302.035d0, &
          310.100d0, 318.165d0, 326.230d0, 334.295d0, 338.105d0, &
          346.0262d0,353.9474d0,361.8686d0,369.7898d0,377.711d0 /

      data icoremap / &
          0, 0, 0, 0,56,55,54,53,54,55,56, 0, 0, 0, 0, &
          0, 0,52,51,50,49,48,47,48,49,50,51,52, 0, 0, &
          0,46,45,44,43,42,41,40,41,42,43,44,45,46, 0, &
          0,39,38,37,36,35,34,33,34,35,36,37,38,39, 0, &
         32,31,30,29,28,27,26,25,26,27,28,29,30,31,32, &
         24,23,22,21,20,19,18,17,18,19,20,21,22,23,24, &
         16,15,14,13,12,11,10, 9,10,11,12,13,14,15,16, &
          8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, &
         16,15,14,13,12,11,10, 9,10,11,12,13,14,15,16, &
         24,23,22,21,20,19,18,17,18,19,20,21,22,23,24, &
         32,31,30,29,28,27,26,25,26,27,28,29,30,31,32, &
          0,39,38,37,36,35,34,33,34,35,36,37,38,39, 0, &
          0,46,45,44,43,42,41,40,41,42,43,44,45,46, 0, &
          0, 0,52,51,50,49,48,47,48,49,50,51,52, 0, 0, &
          0, 0, 0, 0,56,55,54,53,54,55,56, 0, 0, 0, 0  /   ! core map for qtr-symmetry

   end module mod_geom

