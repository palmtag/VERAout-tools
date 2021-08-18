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

      integer :: k

      integer :: maxsta=2  ! ******** temp

      logical :: ifxst     ! logical to check if input file exists

!--- hdf data

      integer(hid_t)   :: file_id       ! File identifier

!--- initialize

      powfile='power.inp'
      hdfname='output.h5'
      iout=8

!--- check to see if the name of the powfile was specified on the command line

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

!!    open(81,file=powfile,status='old',action='read')

!-----------------------------------------
!  Create HDF file and write CORE block
!-----------------------------------------

      call create_hdf(hdfname, file_id, iout)

!-------------------------
!  Loop over statepoints
!-------------------------

      do ista=1, maxsta

        pinpow=0.6d0*ista

!--- write statepoint data to HDF file

        call write_hdf_state(file_id, ista, iout)

      enddo

!--- close powfile input

!!    close(81)   ! close pin file

!--- close HDF file

      call h5fclose_f(file_id, ierror)

      write (0,*)    'finished'
      write (iout,*) 'finished'

      end program create_wb_hdf
