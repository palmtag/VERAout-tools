!=======================================================================
!
!  Program to read two MPACT Restart file and perform general edits
!
!  Copyright (c) 2016 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  ********* NOTICE *************************
!  This utility doesn't really do much, it's purpose is to serve as an example
!  on how to read the MPACT restart file
!  ******************************************
!
!  2016/03/02 - original file
!
!  Use utility "h5dump -H restart.h5" to see format of restart file
!
!-----------------------------------------------------------------------
      program mpactres
      use  hdf5
      use  mod_hdftools
      implicit none

      integer            :: iargs           ! number of command line arguments
      integer            :: k, n
      integer            :: ierror

      integer            :: itype        ! temp variable
      integer            :: ndim         ! temp variable
      integer            :: idim(10)     ! temp variable

      logical            :: ifxst

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=80)  :: group_name
      character(len=20)  :: assm_name       ! assembly label  Assembly_01_01
!!    character(len=10)  :: symmetry
      character(len=80)  :: inputfile
      character(len=80)  :: statename

      integer(hid_t)     :: file_id        ! HDF file ID

! restart file data

      integer  :: irot                      ! rotation
      integer  :: ia                        ! assembly coordinates
      integer  :: ja                        ! assembly coordinates
      integer  :: nxpin                     ! number of pins across assembly
      integer  :: nypin                     ! number of pins across assembly
      integer  :: nxs
      integer  :: nz                        ! number of axial planes
      integer  :: ndat
      integer  :: nxasy                     ! number of assemblies across core
      integer  :: nyasy                     ! number of assemblies across core
      integer  :: nzaid                     ! number of zaids in statepoint
      integer, allocatable :: zaids(:)      ! list of zaids

      real(8)  :: amass
      real(8)  :: exposure

      real(8), allocatable :: axial(:)      ! axial heights
      real(8), allocatable :: axialmesh(:)  ! axial mesh
      real(8), allocatable :: axmass(:)     ! axial mass

!--- read restart file name from command line

      inputfile=' '
      statename=' '

      iargs = command_argument_count()
      if (iargs.lt.2) then
        write (*,*) 'usage:  mpactres.exe [restart file] [statepoint name]'
        write (*,*)
        write (*,*) 'For example:'
        write (*,*) '> mpactres.exe restart.h5 EXP310'
        write (*,*)
        stop
      endif

      call get_command_argument(1,inputfile)
      call get_command_argument(2,statename)

      write (*,'(2a)') 'reading restart file: ', trim(inputfile)
      write (*,'(2a)') 'statepoint: ', trim(statename)

!--- initialize fortran interface (required)

      call h5open_f(ierror)

!--- open HDF file

      inquire(file=inputfile, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(inputfile),' does not exist'
        stop
      endif
      call h5fopen_f (inputfile, H5F_ACC_RDONLY_F, file_id, ierror)   ! read only
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(inputfile), &
                   ' could not be opened'
        stop
      endif

!-------------------------
! Read statepoint data
!-------------------------

      group_name='/CORE_Core/RESTART_'//trim(statename)
      call h5lexists_f(file_id, group_name, ifxst, ierror)
      if (.not.ifxst) then
        write (*,*) trim(group_name)
        write (*,*) 'Statepoint not found on file'
        goto 900
      endif

      dataset=trim(group_name)//'/'//'exposure'
      call hdf5_read_double(file_id, dataset, exposure)
      write (*,*) 'Statepoint exposure ', exposure

      dataset=trim(group_name)//'/'//'nxasy'
      call hdf5_read_integer(file_id, dataset, nxasy)

      dataset=trim(group_name)//'/'//'nyasy'
      call hdf5_read_integer(file_id, dataset, nyasy)

      dataset=trim(group_name)//'/'//'ZAIDs'
      call h5info(file_id, dataset, itype, ndim, idim)
      if (ndim.ne.1) then
        write (*,*) '*** ndim = ', ndim
        write (*,*) '*** idim = ', idim
        stop 'invalid dimensions in zaids'
      endif
      nzaid=idim(1)
      allocate (zaids(nzaid))
      call hdf5_read_integer(file_id, dataset, nzaid, k, zaids)

      write (*,*)
      write (*,*) trim(statename)
      write (*,*) '  exposure ', exposure
      write (*,*) '  nxasy    ', nxasy
      write (*,*) '  nyasy    ', nyasy
      write (*,*) '  nzaid    ', nzaid

!!    write (*,*) zaids(:)

!-------------------------
!  Loop over assemblies
!-------------------------

      assm_name='Assembly_00_00'       ! assembly name template

!!   do ja=1, nxasy
!!     do ia=1, nyasy

     ia=1    ! assembly coordinates
     ja=1    ! assembly coordinates

!--- build assembly label

         write (assm_name(10:11),'(i2.2)') ia
         write (assm_name(13:14),'(i2.2)') ja

!--- read dimensions

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'nxpin'
         call hdf5_read_integer(file_id, dataset, nxpin)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'nypin'
         call hdf5_read_integer(file_id, dataset, nypin)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'nz'
         call hdf5_read_integer(file_id, dataset, nz)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'nxs'
         call hdf5_read_integer(file_id, dataset, nxs)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'ndat'
         call hdf5_read_integer(file_id, dataset, ndat)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'Rotation'
         call hdf5_read_integer(file_id, dataset, irot)

!!       dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'symmetry'
!!       call read_string(file_id, dataset, symmetry)

         allocate (axial(nz))
         allocate (axmass(nz))
         allocate (axialmesh(nz))

!--- distributions

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'Initial Heavy Metal Mass'
         call hdf5_read_double(file_id, dataset, amass)
         write (*,*) 'Initial Heavy Metal Mass', amass

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'Axial Heights'
         call hdf5_read_double(file_id, dataset, nz, n, axial)      ! n will be set to naxial

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'Axial Heavy Metal Mass'
         call hdf5_read_double(file_id, dataset, nz, k, axmass)     ! n will be set to naxial

         dataset=trim(group_name)//'/'//trim(assm_name)//'/'//'Axial Mesh'
         call hdf5_read_double(file_id, dataset, nz, k, axialmesh)

!--- edits

         write (*,*)
         write (*,*) trim(assm_name)
         write (*,*) '   nxpin= ', nxpin
         write (*,*) '   nypin= ', nypin
         write (*,*) '   nz   = ', nz
         write (*,*) '   nxs  = ', nz
         write (*,*) '   ndat = ', ndat
         write (*,*) '   rotation = ', irot
         write (*,*) '    K, axial heights, axial mesh, axial mass'
         do k=nz, 1, -1
           write (*,'(i5,3f10.4)') k, axial(k), axialmesh(k), axmass(k)
         enddo
       


         deallocate (axialmesh)
         deallocate (axmass)
         deallocate (axial)

!!      enddo   ! ia
!!    enddo     ! ja

      deallocate (zaids)

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------
  900 continue

      call h5fclose_f(file_id, ierror)

      write (*,'(/,a)') 'done'

      end program


