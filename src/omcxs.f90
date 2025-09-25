   program omcread
   use  hdf5
   use  mod_hdftools
   implicit none
!=======================================================================
!
!  Program to read OpenMC cross section HDF file and create 
!  text input for Lupine (and other codes)
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
!  User may have to modify program manually to set which
!  cross section types are present (if necessary)
!
!-----------------------------------------------------------------------

      character(len=80)  :: inputfile
      character(len=80)  :: carg            ! command line argument
      integer            :: iargs           ! number of command line arguments
      integer            :: i, j
      integer            :: ierror
      integer            :: icell=0         ! cell number

      integer            :: ngroups
      integer            :: maxcell

      integer, parameter :: max_cell_list=100    ! maximum cell number
      integer            :: cell_list(max_cell_list)

      logical            :: ifxst
      logical            :: ifdebug=.false. ! debug flag
      logical            :: ifrmsmg=.false. ! special monte carlo edits

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=12)  :: group_name

      integer(hid_t)     :: file_id

      real(8)            :: xtemp
      real(8)            :: chisum
      real(8)            :: d1, siga1

! cross section data

      type type_cell
        sequence
        real(8), allocatable :: xsabs(:)    ! (ngroup)
        real(8), allocatable :: xscap(:)    ! (ngroup)
        real(8), allocatable :: xschi(:)    ! (ngroup)
        real(8), allocatable :: xsfis(:)    ! (ngroup)
        real(8), allocatable :: xskapfis(:) ! (ngroup)
        real(8), allocatable :: xsnufis(:)  ! (ngroup)
        real(8), allocatable :: xstran(:)   ! (ngroup)
        real(8), allocatable :: xstot(:)    ! (ngroup)
        real(8), allocatable :: xsscat(:,:) ! (ngroup,ngroup)
        integer  :: id
        logical  :: iffis
      end type type_cell
      type(type_cell), allocatable :: cell(:)

!--- initialize

      inputfile=' '

      cell_list(:)=0

      write (*,*) 'Reading OpenMC HDF Cross Section file'

!----------------------------------------------------------------------
!  Read arguments from command line
!----------------------------------------------------------------------

      iargs = command_argument_count()
      if (iargs.lt.1) then
        write (*,*) 'usage:  omcread.exe [hdf5_file] {options}'
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
!         write (*,*) 'run omcread with no command line arguments for help'
!       else
          inputfile=carg
!       endif
      enddo

      if (inputfile.eq.' ') then
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

      ngroups=0

!   number of groups - read group attribute (I had to write a special routine for this)

      group_name='/'
      dataset='# groups'   ! actually an attribute
      call read_grpatt_int(file_id, group_name, dataset, ngroups)

!--- count number of cells

      maxcell=0
      do icell=1, max_cell_list
        write (group_name,'(a,i0,a)') '/cell/', icell, '/'
        if (ifdebug) write (*,*) 'debug: checking state= ', group_name
        call h5lexists_f(file_id, group_name, ifxst, ierror)
        if (ifxst) then
          maxcell=maxcell+1
          cell_list(maxcell)=icell
        endif
      enddo

      write (0,*) 'number of groups ', ngroups
      write (0,*) 'number of cells  ', maxcell
      write (*,*) 'number of groups ', ngroups
      write (*,*) 'number of cells  ', maxcell

!--- allocate cross section arrays

! Possible groups on HDF5:
!       GROUP "absorption"
!       GROUP "capture"
!       GROUP "chi"
!       GROUP "fission"
!       GROUP "kappa-fission"
!       GROUP "nu-fission"      ! required for LUPINE
!       GROUP "scatter matrix"  ! required for LUPINE
!       GROUP "total"
!       GROUP "transport"    ! required for LUPINE

      allocate (cell(maxcell))
      do icell=1, maxcell
        allocate (cell(icell)%xsabs(ngroups))
        allocate (cell(icell)%xscap(ngroups))
        allocate (cell(icell)%xschi(ngroups))
        allocate (cell(icell)%xsfis(ngroups))
        allocate (cell(icell)%xskapfis(ngroups))
        allocate (cell(icell)%xsnufis(ngroups))
        allocate (cell(icell)%xstran(ngroups))
        allocate (cell(icell)%xstot(ngroups))
        allocate (cell(icell)%xsscat(ngroups,ngroups))
        cell(icell)%id=cell_list(icell)
        cell(icell)%iffis=.false.
        cell(icell)%xsabs    =0.0d0
        cell(icell)%xscap    =0.0d0
        cell(icell)%xschi    =0.0d0
        cell(icell)%xsfis    =0.0d0
        cell(icell)%xskapfis =0.0d0
        cell(icell)%xsnufis  =0.0d0
        cell(icell)%xstran   =0.0d0
        cell(icell)%xstot    =0.0d0
        cell(icell)%xsscat   =0.0d0
      enddo

!--- read cell data

      do icell=1, maxcell
        write (*,'(/,a,i0)') 'reading cell ', cell(icell)%id

        write (group_name,'(a,i0,a)') '/cell/', cell(icell)%id, '/'

!  absorption

        dataset=trim(group_name)//'absorption/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cell(icell)%xsabs=0.0d0
        else
          call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xsabs)
          do i=1, ngroups
            write (*,160) i, cell(icell)%xsabs(i)
          enddo
        endif

!  capture

        dataset=trim(group_name)//'capture/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cell(icell)%xscap=0.0d0
        else
          call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xscap)
          do i=1, ngroups
            write (*,160) i, cell(icell)%xscap(i)
          enddo
        endif

!  fission (put this before other fission-like cross sections)

        dataset=trim(group_name)//'fission/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cycle
        endif

        call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xsfis)

        xtemp=0.0d0   ! check if fissionable
        do i=1, ngroups
          xtemp=xtemp+cell(icell)%xsfis(i)
        enddo
        if (xtemp.eq.0.0d0) then
           cell(icell)%iffis=.false.
        else
           cell(icell)%iffis=.true.
        endif

        if (cell(icell)%iffis) then
          do i=1, ngroups
            write (*,160) i, cell(icell)%xsfis(i)
          enddo
        endif

!  chi

        dataset=trim(group_name)//'chi/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cycle
        endif

        call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xschi)

        chisum=0.0d0
        do i=1, ngroups
          chisum=chisum+cell(icell)%xschi(i)
        enddo
        if (chisum.gt.0.0d0 .and. .not.cell(icell)%iffis) then
           write (*,*) 'ERROR: chi found in non-fissionable cell'
           cell(icell)%iffis=.true.
        endif

        if (cell(icell)%iffis) then
          do i=1, ngroups
            write (*,160) i, cell(icell)%xschi(i)
          enddo
          write (*,*) 'chi sum ', chisum
        endif

!  kappa fission

        dataset=trim(group_name)//'kappa-fission/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cycle
        endif

        call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xskapfis)

        xtemp=0.0d0
        do i=1, ngroups
          xtemp=xtemp+cell(icell)%xskapfis(i)
        enddo
        if (xtemp.gt.0.0d0 .and. .not.cell(icell)%iffis) then
           write (*,*) 'ERROR: kappa-fission found in non-fissionable cell'
           cell(icell)%iffis=.true.
        endif

        if (cell(icell)%iffis) then
          do i=1, ngroups
            if (cell(icell)%xsfis(i).gt.0.0d0) then
               xtemp=cell(icell)%xskapfis(i)/cell(icell)%xsfis(i)
            else
               xtemp=0.0d0
            endif
            write (*,160) i, cell(icell)%xskapfis(i), xtemp
          enddo
        endif

!  nu fission

        dataset=trim(group_name)//'nu-fission/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
          cycle
        endif

        call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xsnufis)

        xtemp=0.0d0
        do i=1, ngroups
          xtemp=xtemp+cell(icell)%xsnufis(i)
        enddo
        if (xtemp.gt.0.0d0 .and. .not.cell(icell)%iffis) then
           write (*,*) 'ERROR: nu-fission found in non-fissionable cell'
           cell(icell)%iffis=.true.
        endif

        if (cell(icell)%iffis) then
          do i=1, ngroups
            if (cell(icell)%xsfis(i).gt.0.0d0) then
               xtemp=cell(icell)%xsnufis(i)/cell(icell)%xsfis(i)
            else
               xtemp=0.0d0
            endif
            write (*,160) i, cell(icell)%xsnufis(i), xtemp
          enddo
        endif

!  total

        dataset=trim(group_name)//'total/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
        else
          call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xstot)
          do i=1, ngroups
            write (*,160) i, cell(icell)%xstot(i)
          enddo
        endif

!  transport

        dataset=trim(group_name)//'transport/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
        else
          write (*,*) 'transport and D'
          call hdf5_read_double(file_id, dataset, ngroups, cell(icell)%xstran)
          do i=1, ngroups
            write (*,160) i, cell(icell)%xstran(i), 1.0d0/(3.0d0*cell(icell)%xstran(i))
          enddo
        endif

!  scatter matrix

        dataset=trim(group_name)//'scatter matrix/average'
        call h5lexists_f(file_id, dataset, ifxst, ierror)
        if (.not.ifxst) then
          write (*,'(2a)') '>>', trim(dataset)
          write (*,'(2a)') 'WARNING: dataset does not exist for this file'
        else
          call hdf5_read_double(file_id, dataset, ngroups, ngroups, cell(icell)%xsscat)
          do j=1, ngroups
            write (*,180) i, cell(icell)%xsscat(:,j)
          enddo
        endif

      enddo
  160 format (i4,1p,5e15.8)
  180 format (i4,1p,40e14.6)

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------

      call h5fclose_f(file_id, ierror)

!--- one group edits

      if (ngroups.eq.1) then
         write (*,*)
         write (*,*) 'one group edits'
         icell=2
         i=1    ! group
         d1=1.0d0/(3.0d0*cell(icell)%xstran(i))
         siga1=cell(icell)%xsfis(i)+cell(icell)%xscap(i)
         write (*,'(1x,a,f10.4)')    'd1    ', d1
         write (*,'(1x,a,1p,e12.4)') 'siga1 ', siga1
         write (*,'(1x,a,f10.4)')    'L^2   ', d1/siga1
         write (*,'(1x,a,f10.4)')    'L     ', sqrt(d1/siga1)
         write (*,'(f10.4, 1p, e12.4, 0p, f10.4)') d1, siga1, sqrt(d1/siga1)
      endif

!--- balance
!      abs = fis + cap
!      tot = abs + scat

      write (*,*)
      write (*,*) 'Cross Section Balance - check abs'
      write (*,*) '  abs,  fis+cap, diff'
      do icell=1, maxcell
        write (*,'(1x,a,i0,1x,i0)') 'cell', icell, cell(icell)%id
        do i=1, ngroups
          xtemp=cell(icell)%xsfis(i)+cell(icell)%xscap(i)   ! abs
          write (*,180) i, cell(icell)%xsabs(i), xtemp, cell(icell)%xsabs(i)-xtemp
        enddo
      enddo

      write (*,*)
      write (*,*) 'Cross Section Balance - check transport'
      write (*,*) 'Is the in-group cross section corrected for????'
      write (*,*) ' trans, tot, fis+cap+scat, trans-sum'
      do icell=1, maxcell
        write (*,'(1x,a,i0,1x,i0)') 'cell', icell, cell(icell)%id
        do i=1, ngroups
          xtemp=cell(icell)%xsfis(i)+cell(icell)%xscap(i)   ! abs
          do j=1, ngroups
             xtemp=xtemp+cell(icell)%xsscat(j,i)  !   j <-- i
          enddo
          write (*,180) i, cell(icell)%xstran(i), cell(icell)%xstot(i), xtemp, &
                  cell(icell)%xstran(i)-xtemp
        enddo
      enddo

!--- convert total to removal (transport - in-group scattering)

!r    write (*,*) 'setting total to transport - in-group scatter  (i.e. removal)'
!r    write (*,*) 'setting in-group scatter to zero'
!r    do icell=1, maxcell
!r      do i=1, ngroups
!r        cell(icell)%xstot(i)=cell(icell)%xstran(i)-cell(icell)%xsscat(i,i)
!r        cell(icell)%xsscat(i,i)=0.0d0
!r      enddo
!r    enddo

!--- edits for Monte Carlo - hardcode to last cell number

      if (ifrmsmg) then
        icell=maxcell
        write (*,*)
        write (*,*) 'Monte Carlo edits'
        write (*,"(' ngrp ', i0,' /')") ngroups
        write (*,220,advance='no') (cell(icell)%xsabs(i),i=1,ngroups)
        write (*,*) '/ abs'
  220   format (' siga ', 1p, 100e12.5)
        write (*,*) 'sigs / start input on next line'
        do j=1, ngroups
          write (*,'(1p,100e12.5)',advance='no') (cell(icell)%xsscat(i,j),i=1,ngroups)
          write (*,'(a,i0)') ' / from ', j
        enddo
      endif

!--- Final edits to text file

      write (*,*) 'creating file: xs.txt for Lupine style edits'

      open (12,file='xs.txt')
      write (12,'(a,i0)') 'ngroup ', ngroups
      write (12,'(a,i0)') 'niso   ', maxcell
      write (12,*)
      do icell=1, maxcell
         write (12,'(a,i4,a,i0)') '#',icell,' MAT',cell(icell)%id
      enddo
      do icell=1, maxcell
        write (12,*)
        write (12,'(a,i0)') 'name MAT', cell(icell)%id

        write (12,'(a)') 'transport'
        do i=1, ngroups
          write (12,'(1p,e15.8)') cell(icell)%xstran(i)
        enddo

        write (12,'(a)') 'total'     ! now the removal cross section
        do i=1, ngroups
          write (12,'(1p,e15.8)') cell(icell)%xstot(i)
        enddo

        write (12,'(a)') 'scatter'
        do i=1, ngroups
          write (12,'(1p,100e15.8)') (cell(icell)%xsscat(i,j),j=1,ngroups)
        enddo

        if (cell(icell)%iffis) then

          write (12,'(a)') 'nusf'
          do i=1, ngroups
            write (12,'(1p,e15.8)') cell(icell)%xsnufis(i)
          enddo

          write (12,'(a)') 'chi'
          do i=1, ngroups
            write (12,'(1p,e15.8)') cell(icell)%xschi(i)
          enddo
        endif
      enddo
      close (12)

!--- edit finemesh edits to xs2

      write (*,*) 'creating file: xs2.txt'

      open (12,file='xs2.txt')
      write (12,'(a,i0)') 'group ', ngroups
      write (12,'(a,i0)') 'niso   ', maxcell
      write (12,'(a)') '# ============================================'
      write (12,'(a)') '#   Cross sections'
      write (12,*)
      do icell=1, maxcell
         write (12,'(a,i4,a,i0,a)') 'mat ',icell,' MAT',cell(icell)%id,' /'
      enddo

      do icell=1, maxcell
        write (12,*)
        write (12,'(a,i0)') '# name MAT', cell(icell)%id

        write (12,'(a,i0)') 'diff ', icell
        do i=1, ngroups
          write (12,'(1p,e15.8)') 1.0d0/(3.0d0*cell(icell)%xstran(i))
        enddo

        write (12,'(a,i0)') 'abs ', icell     ! now the removal cross section
        do i=1, ngroups
          write (12,'(1p,e15.8)') cell(icell)%xsabs(i)
        enddo

        write (12,'(a,i0)') 'scat ', icell
        do i=1, ngroups
          write (12,'(1p,100e16.8)') (cell(icell)%xsscat(j,i),j=1,ngroups)
        enddo

        if (cell(icell)%iffis) then

          write (12,'(a,i0)') 'nufis ', icell
          do i=1, ngroups
            write (12,'(1p,e15.8)') cell(icell)%xsnufis(i)
          enddo

          write (12,'(a,i0)') 'chi ', icell
          do i=1, ngroups
            write (12,'(1p,e15.8)') cell(icell)%xschi(i)
          enddo
        endif
      enddo
      write (12,'(a)') 'sta'
      close (12)

      write (*,'(/,a)') 'done'

      end program omcread
