   program omctally
   use  hdf5
   use  mod_hdftools
   implicit none
!=======================================================================
!
!  Program to read OpenMC statepoint file and print tallies
!
!  Copyright (c) 2024 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2024/06/12 - original version based on MPACTread
!
!-----------------------------------------------------------------------

      character(len=80)  :: inputfile
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j, k
      integer            :: ierror

      integer            :: ids
      integer            :: ncells
      integer            :: n_particles
      integer            :: n_batches
      integer            :: n_inactive
      integer            :: n_tallies
      integer            :: n_filters
      integer, allocatable :: cell_list(:)   ! (ncells)
      real(8), allocatable :: results(:,:,:) ! (2,1,ncells)

      logical            :: ifxst

      real(8)            :: kcombined(2)

      character(len=80)  :: dataset          ! HDF dataset name
      character(len=80)  :: group_name

      integer(hid_t)     :: file_id

      integer, parameter :: npin=9       ! *** hard-wired mask for 17x17
      integer :: mask(npin,npin)
      real(8) :: power(npin,npin)
      real(8) :: xsum
      data mask /  &
         0, 1, 1, 0, 1, 1, 0, 1, 1, &
         1, 1, 1, 1, 1, 1, 1, 1, 1, &
         1, 1, 1, 1, 1, 1, 1, 1, 1, &
         0, 1, 1, 0, 1, 1, 0, 1, 1, &
         1, 1, 1, 1, 1, 1, 1, 1, 1, &
         1, 1, 1, 1, 1, 0, 1, 1, 1, &
         0, 1, 1, 0, 1, 1, 1, 1, 1, &
         1, 1, 1, 1, 1, 1, 1, 1, 1, &
         1, 1, 1, 1, 1, 1, 1, 1, 1 /

!--- initialize

      inputfile=' '

      write (*,*) 'Reading OpenMC HDF Cross Section file'

!----------------------------------------------------------------------
!  Read arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  omctally.exe [hdf5_file] {options}'
        write (*,*)
        stop
      endif

! parse command line arguments

      do i=1, iargs
        call get_command_argument(i,carg)
!       if     (carg.eq.'3D' .or. carg.eq.'3d') then
!         if3d=.true.
!       elseif (carg.eq.'2D' .or. carg.eq.'2d') then
!         if2d=.true.
!       elseif (carg.eq.'2DA' .or. carg.eq.'2da') then
!         if2da=.true.
!       elseif (carg.eq.'-help' .or. carg.eq.'--help') then  ! add for Ben
!         write (*,*) 'run omctally with no command line arguments for help'
!       else
          inputfile=carg
!       endif
      enddo

      if (inputfile.eq.' ') then
        write (*,*) 'usage: omctally.exe [file.h5]'
        stop 'no input file specified on command line'
      endif

!--- initialize fortran interface

      call h5open_f(ierror)

!--- open HDF file

      write (*,'(2a)') 'reading h5 file: ', trim(inputfile)
      inquire(file=inputfile, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(inputfile),' does not exist'
        stop 'input file does not exist'
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop 'H5 input file could not be opened'
      endif

!--- Read HDF file data


!  keff

      dataset='k_combined'
      write (*,*) 'debug: reading dataset ', trim(dataset)
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (.not.ifxst) then
        write (*,'(/,1x,a)') 'this does not appear to be an OpenMC statepoint file - exiting'
        goto 800
      endif
      call hdf5_read_double(file_id, dataset, 2, kcombined)
      write (0,*) 'k_combined ', kcombined(:)

!  run info

      dataset='n_batches'
      call hdf5_read_integer(file_id, dataset, n_batches)
      dataset='n_inactive'
      call hdf5_read_integer(file_id, dataset, n_inactive)
      dataset='n_particles'
      call hdf5_read_integer(file_id, dataset, n_particles)

      write (0,*) 'n_particles', n_particles
      write (0,*) 'n_batches  ', n_batches
      write (0,*) 'n_inactive ', n_inactive

!  tally info

      group_name='/tallies'

      dataset='n_tallies'   ! actually an attribute
      call read_grpatt_int(file_id, group_name, dataset, n_tallies)
      write (0,*) 'n_tallies  ', n_tallies
      if (n_tallies.eq.0) then
        write (*,'(/,a)') ' tallies do not exist on this file'
        goto 800
      endif

      dataset='ids'   ! actually an attribute
      call read_grpatt_int(file_id, group_name, dataset, ids)
      write (0,*) 'ids       ', ids

      group_name='/tallies/filters'

      dataset='n_filters'  ! actually an attribute
      call read_grpatt_int(file_id, group_name, dataset, n_filters)
      write (0,*) 'n_filters ', n_filters

!  read tallies for this filter

        group_name='/tallies/filters/filter 1'

        dataset=trim(group_name)//'/n_bins'
        call hdf5_read_integer(file_id, dataset, ncells)
        write (0,*) 'nbins     ', ncells

        allocate (cell_list(ncells))
        allocate (results(2,1,ncells))

        dataset=trim(group_name)//'/bins'
        call hdf5_read_integer(file_id, dataset, ncells, cell_list)
!d      do i=1, ncells
!d        write (*,'(i4,i6)') i, cell_list(i)
!d      enddo

        group_name='/tallies/tally 1'

        dataset=trim(group_name)//'/results'
        call hdf5_read_double(file_id, dataset, 2, 1, ncells, results)
        do i=1, ncells
          write (*,'(i4,i4,1p,2e14.5)') i, cell_list(i), results(1,1,i), results(2,1,i)
        enddo

!--- fill 2D array

        if (sum(mask).eq.ncells) then    ! 17x17  264 power pins
          k=0
          xsum=0.0d0
          power=0.0d0
          do j=1, npin
            do i=1, npin
              if (mask(i,j).gt.0) then
                k=k+1
                power(i,j)=results(1,1,k)
                xsum=xsum+power(i,j)
              endif
            enddo
          enddo

          write (*,*) 'tally sum ', xsum
          if (xsum.gt.0.0d0) power=power*264.0d0*0.25d0/xsum  ! normalize

          write (*,*) 'pin powers (normalized)'
          do j=1, npin
            write (*,'(i3,20f8.5)') j, (power(i,j),i=1,npin)
          enddo
        endif

        deallocate (results)
        deallocate (cell_list)

!  done with filter

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------
  800 continue

      call h5fclose_f(file_id, ierror)

      write (*,'(/,a)') 'done'

      end program omctally
