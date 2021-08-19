   program create_wb_hdf
!=======================================================================
!
!  Program to create an VERA HDF file using user supplied data
!
!  Geometry information is hard-wired for the WB reactor.
!  This information is set in the module "Mod_geom"
!
!  State information is hard-wired in module "Mod_state"
!  This information should be updated every statepoint from user data
!
!  Dependencies:  This program uses the HDF data library that can be
!  downloaded from: https://hdfgroup.org/
!  
!  Scott Palmtag
!  August 18, 2021
!
!=======================================================================
      use hdf5
      use mod_geom,  only : npin, maxasm, kdfuel
      use mod_state, only : mem_state, pinpow
      implicit none

      character(len=100) :: hdfname   ! HDF output file name
      character(len=256) :: powfile   ! input file containing pin powers

      integer :: ierror    ! error code
      integer :: ista      ! statepoint number
      integer :: iout      ! output unit number

      integer :: k, na

      integer :: maxsta=1  ! ***** only write one state for now

      logical :: ifxst     ! logical to check if input file exists

!--- hdf data

      integer(hid_t) :: file_id  ! HDF File identifier

!--- initialize

      powfile='power.inp'
      hdfname='output.h5'
      iout=8

!--- read command line

      k=iargc()
      if (k.eq.0) then
        stop 'usage: reader [powfile] {h5 file}'
      endif
      if (k.ge.1) then
        call getarg(1,powfile)
      endif
      if (k.ge.2) then
        call getarg(2,hdfname)
      else
        hdfname=trim(powfile)//'.h5'
      endif

      if (powfile.eq.hdfname) stop 'input and output file names cannot be the same'

!--- allocate memory for pin distributions

      call mem_state(npin, maxasm, kdfuel)

!--- Open text input file

      inquire (file=powfile,exist=ifxst)
      if (.not.ifxst) then
        write (*,'(2a)') 'reading: ', trim(powfile)
        stop 'powfile is not found'
      endif

!-----------------------------------------
!  Create HDF file and write CORE block
!-----------------------------------------

      call write_hdf_core(hdfname, file_id, iout)

!-----------------------------------------
!  Read user data
!-----------------------------------------

!!    open(81,file=powfile,status='old',action='read')

      pinpow=0.1d0

!!    pinpow(na,k,i,j)

      na=10   ! assembly
      do k=1, npin   ! j
        pinpow(na,:,1,k)=1.0d0+dble(k)*0.1d0
        pinpow(na,:,2,k)=1.0d0+dble(k)*0.1d0
        pinpow(na,:,3,k)=1.0d0+dble(k)*0.1d0
        pinpow(na,:,4,k)=1.0d0+dble(k)*0.1d0
        pinpow(na,:,5,k)=1.0d0+dble(k)*0.1d0
        pinpow(na,:,10,k)=1.0d0+dble(k)*0.1d0
      enddo

      na=44   ! assembly
      do k=1, npin   ! j
        pinpow(na,:,k,4)=2.0d0+dble(k)*0.1d0
        pinpow(na,:,k,5)=2.0d0+dble(k)*0.1d0
        pinpow(na,:,k,6)=2.0d0+dble(k)*0.1d0
        pinpow(na,:,k,7)=2.0d0+dble(k)*0.1d0
      enddo

      na=43   ! assembly
      pinpow(na,:, 1, 1)=1.1d0
      pinpow(na,:, 1,17)=1.2d0
      pinpow(na,:,17, 1)=1.3d0
      pinpow(na,:,17,17)=1.7d0

 ! test k with assembly 52

!!    do k=1, kdfuel
!!      pinpow(:,:,52,k)=dble(k)
!!    enddo


!-------------------------
!  Write statepoints
!-------------------------

!--- write statepoint data to HDF file

      do ista=1, maxsta

        call write_hdf_state(file_id, ista, iout)

      enddo

!--- close powfile input

!!    close(81)   ! close pin file

!--- close HDF file

      call h5fclose_f(file_id, ierror)

      write (0,*)    'finished'
      write (iout,*) 'finished'

      end program create_wb_hdf
