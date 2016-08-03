   program mpactres
   use  hdf5
   use  mod_hdftools
   implicit none
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
!  *** Note that the MPACT restart file uses copressed HDF
!
!  2016/03/02 - original file
!  2016/07/13 - add average isotopic edits (intendend for a pincell)
!  2016/07/23 - add isotopic edits of first 3 rings
!
!  Use utility "h5dump -H restart.h5" to see format of restart file
!
!-----------------------------------------------------------------------

      logical  :: ifdata=.false.         ! print long data

      integer  :: iargs        ! number of command line arguments
      integer  :: k, n
      integer  :: i, ii
      integer  :: ierror
      integer  :: izsum

      integer  :: itype        ! temp variable
      integer  :: ndim         ! temp variable
      integer  :: idim(10)     ! temp variable

      logical  :: ifxst

      character(len=80)  :: dataset         ! HDF dataset name
      character(len=80)  :: group_name
      character(len=20)  :: assm_name       ! assembly label  Assembly_01_01
!!    character(len=10)  :: symmetry
      character(len=80)  :: inputfile
      character(len=80)  :: statename

      integer(hid_t)     :: file_id         ! HDF file ID

      real(8)  :: vsum

! restart file data

      integer  :: irot                      ! rotation
      integer  :: iver                      ! version
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
      integer  :: nzone
      integer, allocatable :: zaids(:)      ! list of zaids
      integer, allocatable :: nucid(:)      ! list of nuclear id's by assembly
      real(8), allocatable :: fcomp(:)      ! fuel composition     by assembly

      real(8)  :: amass
      real(8)  :: exposure
      real(8)  :: dephours                  ! depletion hours

      real(8), allocatable :: axial(:)      ! axial heights
      real(8), allocatable :: axialmesh(:)  ! axial mesh
      real(8), allocatable :: axmass(:)     ! axial mass

      integer, allocatable :: zoneindx(:)   ! nzone
      real(8), allocatable :: zoneburn(:)   ! nzone
      real(8), allocatable :: zonemass(:)   ! nzone
      real(8), allocatable :: zonevolm(:)   ! nzone

!--- read restart file name from command line

      inputfile=' '
      statename=' '

      iargs = command_argument_count()
      if (iargs.lt.2) then
        write (*,*) 'usage:  mpactres.exe [restart file] [statepoint name] {data}'
        write (*,*)
        write (*,*) 'For example:'
        write (*,*) '> mpactres.exe restart.h5 EXP310'
        write (*,*)
        stop
      endif

      call get_command_argument(1,inputfile)
      call get_command_argument(2,statename)

      if (iargs.ge.3) ifdata=.true.

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
!
! Strings:
!   cycle label
!   op date
!   serials
!   unit number
!   xlabel
!   ylabel
!
! Non-strings:
!   exposure           double
!   depletion hours    double
!   nxasy              integer
!   nyasy              integer
!   version            integer
!   ZAIDs              integer array
!
      dataset=trim(group_name)//'/exposure'
      call hdf5_read_double(file_id, dataset, exposure)

      dataset=trim(group_name)//'/depletion hours'
      call hdf5_read_double(file_id, dataset, dephours)

      dataset=trim(group_name)//'/nxasy'
      call hdf5_read_integer(file_id, dataset, nxasy)

      dataset=trim(group_name)//'/nyasy'
      call hdf5_read_integer(file_id, dataset, nyasy)

      dataset=trim(group_name)//'/version'
      call hdf5_read_integer(file_id, dataset, iver)

!  if program crashes on next read, it may be that the
!  HDF library does not have compression installed

      dataset=trim(group_name)//'/ZAIDs'
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
      write (*,*) '  statename ', trim(statename)
      write (*,*) '  exposure ', exposure
      write (*,*) '  dep hours', dephours
      write (*,130) 'nxasy    ', nxasy
      write (*,130) 'nyasy    ', nyasy
      write (*,130) 'version  ', iver
      write (*,130) 'nzaid    ', nzaid

  130 format (3x,a,i6)

      if (ifdata) then
        write (*,*)
        write (*,*) 'Zaids (by statepoint)'
        write (*,*) zaids(:)
      endif

!-------------------------
!  Loop over assemblies
!-------------------------

      assm_name='Assembly_00_00'       ! assembly name template

!d   ia=1    ! assembly coordinates
!d   ja=1    ! assembly coordinates

      do ja=1, nxasy
        do ia=1, nyasy


!--- build assembly label

         write (assm_name(10:11),'(i2.2)') ia
         write (assm_name(13:14),'(i2.2)') ja

!--- read dimensions

!
!  Scalars:
! * Rotation                 H5T_STD_I32LE
! * ndat                     H5T_STD_I32LE
! * nxpin                    H5T_STD_I32LE
! * nxs                      H5T_STD_I32LE
! * nypin                    H5T_STD_I32LE
! * nz                       H5T_STD_I32LE
! * Initial Heavy Metal Mass H5T_IEEE_F64LE
!
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

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Initial Heavy Metal Mass'
         call hdf5_read_double(file_id, dataset, amass)

!!       dataset=trim(group_name)//'/'//trim(assm_name)//'/symmetry'
!!       call read_string(file_id, dataset, symmetry)

!  the following 3D array can be used to determine which zones belong to each pin
!  **** need to check the order of dimensions (nxpin, nypin, nz)
!!       dataset=trim(group_name)//'/'//trim(assm_name)//'/nxsreg pin'
!!       call hdf5_read_double(file_id, dataset, nxpin, nypin, nz, nxsreg)

         write (*,'(/,2a)') ' Assembly ', trim(assm_name)
         write (*,130) 'nxpin   ', nxpin
         write (*,130) 'nypin   ', nypin
         write (*,130) 'nz      ', nz
         write (*,130) 'nxs     ', nz
         write (*,130) 'ndat    ', ndat
         write (*,130) 'rotation', irot
         write (*,*) '  Initial Heavy Metal Mass', amass

         allocate (axial(nz))
         allocate (axmass(nz))
         allocate (axialmesh(nz))

         allocate (nucid(ndat))
         allocate (fcomp(ndat))

! Strings:
!   symmetry                 H5T_STRING
!   Serial Number            H5T_STRING

!--- axial distributions

         write (*,*)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Axial Heights'
         call hdf5_read_double(file_id, dataset, nz, n, axial)      ! n will be set to naxial

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Axial Heavy Metal Mass'
         call hdf5_read_double(file_id, dataset, nz, k, axmass)     ! n will be set to naxial

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Axial Mesh'
         call hdf5_read_double(file_id, dataset, nz, k, axialmesh)

         write (*,'(2a)') ' Assembly ', trim(assm_name)
         write (*,'(3x,a)') 'K, axial heights, axial mesh, axial mass'
         do k=nz, 1, -1
           write (*,'(i4,2f10.4,f12.6)') k, axial(k), axialmesh(k), axmass(k)
         enddo


!--- depletion region data
!
!  Example: pincell calculation has 12 zones per rod, nzone=12
!
!  'Isotope Index Map'                  H5T_STD_I32LE  (nzone)
!  'XS Region Burnup'                   H5T_IEEE_F64LE (nzone)
!  'XS Region Initial Heavy Metal Mass' H5T_IEEE_F64LE (nzone)
!  'XS Region Volume'                   H5T_IEEE_F64LE (nzone)

         write (*,*)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Isotope Index Map'
         call h5info(file_id, dataset, itype, ndim, idim)
         nzone=idim(1)

         allocate (zoneindx(nzone))
         allocate (zoneburn(nzone))
         allocate (zonemass(nzone))
         allocate (zonevolm(nzone))

         call hdf5_read_integer(file_id, dataset, nzone, zoneindx)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/XS Region Burnup'
         call hdf5_read_double(file_id, dataset, nzone, zoneburn)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/XS Region Initial Heavy Metal Mass'
         call hdf5_read_double(file_id, dataset, nzone, zonemass)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/XS Region Volume'
         call hdf5_read_double(file_id, dataset, nzone, zonevolm)

         izsum=0
         vsum=0.0d0

         write (*,'(/,2a)') ' Assembly ', trim(assm_name)
         write (*,*) '   Zone Index,   Burnup,  Initial Mass, Volume'
         do n=1, nzone
           write (*,540) n, zoneindx(n), zoneburn(n), zonemass(n), zonevolm(n)
           izsum=izsum+zoneindx(n)
           vsum =vsum +zonevolm(n)
         enddo
     540 format (i6, i6, f12.3, 1p, 3e14.5)

         write (*,*) 'sum of zone indexes = ', izsum
         if (izsum.ne.ndat) stop 'sum of zone indexes does not match ndat'

         write (*,*) 'total volume = ', vsum

         if (nxpin.eq.17 .and. nypin.eq.17) then
           write (*,*) 'debug: temp vol  ', 3.141592653589d0*0.4096d0*0.4096d0*264.0d0
         endif
         if (nxpin.eq.1 .and. nypin.eq.1) then
           write (*,*) 'debug: temp vol  ', 3.141592653589d0*0.4096d0*0.4096d0
         endif

!--- Nuclide distributions

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Nuclide ID'
         call hdf5_read_integer(file_id, dataset, ndat, nucid)

         dataset=trim(group_name)//'/'//trim(assm_name)//'/Fuel Composition'
         call hdf5_read_double (file_id, dataset, ndat, fcomp)

         if (ifdata) then
           ii=0
           do n=1, nzone
             write (*,*) 'Zone ', n
             do i=1, zoneindx(n)
                ii=ii+1
                write (*,550) ii, zaids(nucid(ii)), fcomp(ii), zonevolm(n)
             enddo
           enddo
           if (ii.ne.ndat) stop 'ndat sum error'
         endif
     550 format (i10, i8, 1p, 2e14.6)

!   calculate average isotopics over assembly

         if (nz.eq.1) call ave_pin(nzone, ndat, nzaid, zoneindx, zonevolm, nucid, fcomp, zaids)

!--- end of assembly

         deallocate (zonevolm)
         deallocate (zonemass)
         deallocate (zoneburn)
         deallocate (zoneindx)

         deallocate (axialmesh)
         deallocate (axmass)
         deallocate (axial)
         deallocate (nucid)
         deallocate (fcomp)

        enddo   ! ia
      enddo     ! ja

      deallocate (zaids)

!--------------------------------------------------------------------------------
! Finish
!--------------------------------------------------------------------------------
  900 continue

      call h5fclose_f(file_id, ierror)

      write (*,'(/,a)') 'done'

      end program



!=======================================================================
!
!   Subroutine to calculate average isotopes in a fuel pin
!   and print a list of mcnp-compatible isotopes
!
!=======================================================================
      subroutine ave_pin(nzone, ndat, nzaid, zoneindx, zonevolm, nucid, fcomp, zaids)
      implicit none

      integer, intent(in) :: nzone
      integer, intent(in) :: ndat
      integer, intent(in) :: nzaid

      integer, intent(in) :: zoneindx(nzone)
      real(8), intent(in) :: zonevolm(nzone)
      integer, intent(in) :: nucid(ndat)
      real(8), intent(in) :: fcomp(ndat)
      integer, intent(in) :: zaids(nzaid)

!--- local

      integer :: i, j, n, ii
      integer :: j5, j4
      integer :: iz, isom
      integer :: niso
      real(8) :: vsum
      real(8) :: xtmp

      real(8), allocatable :: ring(:,:)  ! long list
      real(8), allocatable :: avepin(:)  ! long list
      real(8), allocatable :: ave2(:)    ! short list
      real(8), allocatable :: ave3(:,:)  ! short list in first three rings
      integer, allocatable :: izaid(:)   ! short list
      integer, allocatable :: isave(:)   ! short list


!--- start

      allocate (ring(3,nzaid))
      allocate (avepin(nzaid))
      allocate (ave2(nzaid))       ! short list
      allocate (ave3(3,nzaid))     ! short list
      allocate (izaid(nzaid))      ! short list
      allocate (isave(nzaid))      ! short list
      avepin=0.0d0
      ring=0.0d0

!--- average over all zones

      vsum=0.0d0
      ii=0
      do n=1, nzone
        do i=1, zoneindx(n)
           ii=ii+1
           avepin(nucid(ii))=avepin(nucid(ii)) + fcomp(ii)*zonevolm(n)
        enddo
        vsum=vsum+zonevolm(n)
      enddo
      if (ii.ne.ndat) stop 'ndat sum error'

  550 format (i6, i8, 1p, 8e14.6)
  560 format (i6, 2i8, 1p, 8e14.6)

      ii=0
      do n=1, 3    ! **** only first 3 rings
        do i=1, zoneindx(n)
           ii=ii+1
           ring(n,nucid(ii))=ring(n,nucid(ii)) + fcomp(ii)
        enddo
      enddo

!--- find non-zero values and create short list

      write (*,*)
      write (*,*) 'Assembly average number densities in fuel'
      ii=0
      do n=1, nzaid
        if (avepin(n).ne.0.0d0) then
           ii=ii+1    ! number of non-zero values
           izaid(ii)=zaids(n)
           ave2 (ii)=avepin(n)/vsum
           ave3 (1,ii)=ring(1,n)
           ave3 (2,ii)=ring(2,n)
           ave3 (3,ii)=ring(3,n)
        endif
      enddo
      niso=ii      ! short list

!--- convert to MCNP isotopes

      isave=izaid     ! save copy of original MPACT zaids

      do i=1, niso

        iz=izaid(i)/10
        isom=izaid(i)-10*iz   ! isomer flag

        j5=iz/100      ! remove 500 from mpact fission product isotopes
        j4=iz/1000
        j5=j5-j4*10    ! third digit
        if (j5.ge.5) iz=iz-500

!  MPACT isomers:
!  224   47610 1    Ag-110M  add 400 to isomer
!  233   52627 1    Te-127M  add 400 to isomer
!  234   52629 1    Te-129M  add 400 to isomer
!  137   61648 1    Pm-148M  add 400 to isomer
!  178   95242 1    Pm-242M  add 400 to stable nuclide!

        if (isom.gt.0) then 
!d        write (0,*) 'debug: isomer ', iz
          if (iz.eq.47110 .or. iz.eq.52127 .or. iz.eq.52129 .or. iz.eq.61148) then
             iz=iz+400    ! put in mcnp notation
          elseif (iz.eq.95242) then   ! do nothing
             continue
          else
             write (*,*) 'zaid = izaid(i)'
             stop 'unknown isomer in mpact isotope list'
          endif
        endif

        if (isom.eq.0 .and. iz.eq.95242) then
          iz=iz+400   ! add 400 to stable isomer
        endif

        if (iz.eq.8001) iz=8016

        izaid(i)=iz

      enddo

!--- sort

      do i=1, niso-1
        do j=i+1, niso
          if (izaid(i).gt.izaid(j)) then
            n=izaid(i)
            izaid(i)=izaid(j)
            izaid(j)=n

            n=isave(i)
            isave(i)=isave(j)
            isave(j)=n

            xtmp=ave2(i)
            ave2(i)=ave2(j)
            ave2(j)=xtmp

            xtmp=ave3(1,i)
            ave3(1,i)=ave3(1,j)
            ave3(1,j)=xtmp

            xtmp=ave3(2,i)
            ave3(2,i)=ave3(2,j)
            ave3(2,j)=xtmp

            xtmp=ave3(3,i)
            ave3(3,i)=ave3(3,j)
            ave3(3,j)=xtmp

          endif
        enddo
      enddo

!--- look for duplicates (list has already been sorted)

      do i=1, niso-1
        if (izaid(i).eq.izaid(i+1)) then
           write (*,*) 'zaid = ', izaid(i)
           stop 'duplicate zaid found'
        endif
      enddo

!--- print short list

!  ring 1 is outside ring -> ring 3 is inside ring

      write (*,'(/,9x,a)')  'mpact    mcnp    ave           reg 1 (out)   reg 2         reg 3 (in)'
      do i=1, niso
         write (*,560) i, isave(i), izaid(i), ave2(i), ave3(:,i)
      enddo

      write (*,*) 'Warning: may have to manually remove 51127 - not in MCNP library'
      write (*,*) 'Warning: may have to manually remove 65161 - not in MCNP library'

      open (24,file='isolist3')
      do i=1, niso
        if (izaid(i).eq.51127) cycle
        if (izaid(i).eq.65161) cycle
        write (24,550) i, izaid(i), ave2(i), ave3(:,i)
      enddo
      close(24)

      deallocate(ring)
      deallocate(avepin)
      deallocate(ave2)
      deallocate(izaid)
      deallocate(isave)

      return
      end subroutine ave_pin

