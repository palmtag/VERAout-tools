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
      use mod_state, only : mem_state
      implicit none

      character(len=100) :: hdfname   ! HDF output file name
      character(len=256) :: powfile   ! input file containing pin powers

      integer :: ierror    ! error code
      integer :: ista      ! statepoint number
      integer :: iout      ! output unit number

      integer :: k

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

      call readpow(powfile)

!-------------------------
!  Write statepoints
!-------------------------

!--- write statepoint data to HDF file

      do ista=1, maxsta

        call write_hdf_state(file_id, ista, iout)

      enddo

!--- close HDF file

      call h5fclose_f(file_id, ierror)

      write (0,*)    'finished'
      write (iout,*) 'finished'

      end program create_wb_hdf

!=======================================================================
!
!  Read user data from file
!
!=======================================================================
      subroutine readpow(powfile)
      use mod_state
      implicit none

      character(len=*), intent(in) :: powfile   ! input file name

      character(len=200) :: line
      character(len=10)  :: key
      integer :: na, i, j, k

      pinpow=0.0d0

      open(81,file=powfile,status='old',action='read')

      do
        read (81,'(a)',end=100) line
        if (line.eq.' ') cycle   ! skip blank line
        read (line,*) key        ! read first keyword, might be a number
        if     (key.eq.'boron') then
           read (line,*) key, boron
        elseif (key.eq.'keff') then
           read (line,*) key, xkeff
        elseif (key.eq.'power') then
           read (line,*) key, power
        elseif (key.eq.'flow') then
           read (line,*) key, flow
        elseif (key.eq.'efpd') then
           read (line,*) key, efpd
        elseif (key.eq.'exposure') then
           read (line,*) key, ave_exp
        else          ! read power
          read (line,*) na, k, i, j, pinpow(na,k,i,j)
        endif
      enddo
  100 continue

      close (81)

! ???? check normalization

      return
      end subroutine readpow

