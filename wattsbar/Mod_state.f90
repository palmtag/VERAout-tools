   module mod_state

!  hard-wired state information
!  this information written to STATE block of HDF file

!  this information should be updated at every statepoint

      real(8) :: power=100.0d0      ! power (%)
      real(8) :: flow=100.0d0       ! flow (%)
      real(8) :: xkeff=1.0d0
      real(8) :: boron=1200.0d0
      real(8) :: pressure=2250.0d0  ! pressure (psia)

      real(8) :: ave_exp=0.0d0      ! cycle exposure (GWD/MT)
      real(8) :: efpd=0.0d0         ! cycle exposure (EFPD)

      real(8), allocatable :: pinpow(:,:,:,:)   ! 4-D pin-by-pin powers
      real(8), allocatable :: pinexp(:,:,:,:)   ! 4-D pin-by-pin exposures

   contains

!=======================================================================
!
!  Subroutine to allocate and initially fill distributions
!
!=======================================================================

      subroutine mem_state(npin, maxasm, kdfuel)
      implicit none
      integer, intent(in) :: npin
      integer, intent(in) :: maxasm
      integer, intent(in) :: kdfuel

      allocate(pinpow(npin,npin,maxasm,kdfuel))
      allocate(pinexp(npin,npin,maxasm,kdfuel))

      pinpow=0.0    ! reset each statepoint
      pinexp=0.0

      return
      end subroutine mem_state

!=======================================================================

   end module mod_state

