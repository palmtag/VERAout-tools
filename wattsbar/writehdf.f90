   subroutine create_hdf(hdfname, file_id, iout)
!=======================================================================
!
!  Subroutine to create HDF file and write CORE block
!
!=======================================================================
      use hdf5
      use mod_hdftools
      use mod_geom
      implicit none

      character(len=*), intent(in) :: hdfname
      integer(hid_t),   intent(out):: file_id    ! HDF file identifier
      integer,          intent(in) :: iout       ! output unit

!--- local

      integer :: itmp
      integer :: ierror
      integer :: idims(10)           ! dataset dimensions

      integer(hid_t) :: group_id     ! HDF group identifier

      real(8), allocatable :: pinload(:,:,:,:)  ! (maxasm,kdfuel,npin,npin)

      character(len=60) :: dsetname     ! HDF dataset name
      character(len=20) :: group_name   ! HDF group name

!--- Initialize fortran interface.

      call h5open_f(ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

!--- Create a new HDF file using default properties.

      write (*,'(/,2a,i10)') ' creating HDF file: ', trim(hdfname)

      call h5fcreate_f(hdfname, H5F_ACC_TRUNC_F, file_id, ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

!--- edits

        write (iout,*)
        write (iout,55) 'Number of Assemblies', maxasm
        write (iout,55) 'Number of axial     ', kdfuel
        write (iout,59) 'apitch            ', apitch
        write (iout,65) 'Rated Power       ', prated, 'MW'
        write (iout,65) 'Rated Flow        ', crated, 'Mlb/hr'
        write (iout,*)

   55   format(3X,A,I8)
   59   format(3X,A,F8.3)
   65   format(3X,A,F8.2,1X,A,' (full core)')
!  70   format(3X,A,F8.3)

!--- create CORE group

      group_name='/CORE/'
      call h5gcreate_f (file_id, group_name, group_id, ierror)
      if (ierror.ne.0) stop 'error creating hdf group'

      dsetname = trim(group_name)//'apitch'
      call hwrite_double_scalar(file_id, dsetname, apitch)

      dsetname = trim(group_name)//'rated_power'
      call hwrite_double_scalar(file_id, dsetname, prated)

      dsetname = trim(group_name)//'rated_flow'
      call hwrite_double_scalar(file_id, dsetname, crated)

      if (maxasm.gt.100) then
        itmp=1   ! full core
        write (*,*) 'WARNING: assuming full-core problem'
      else
        itmp=4
        write (*,*) 'WARNING: assuming qtr-core problem'
      endif
      dsetname = trim(group_name)//'core_sym'
      call hwrite_integer_scalar(file_id, dsetname, itmp)

! axial mesh

      idims(:)=0     ! clear
      idims(1)=kdfuel+1

      dsetname = trim(group_name)//'axial_mesh'
      call hwrite_double(file_id, dsetname, idims, zmesh)

!  core map

      idims(:)=0     ! clear
      idims(1)=iafull
      idims(2)=iafull

      dsetname = trim(group_name)//'core_map'
      call hwrite_integer(file_id, dsetname, idims, icoremap)

!   pin loadings

      allocate (pinload(maxasm,kdfuel,npin,npin))
      pinload=1.0d0

      idims(:)=0     ! clear
      idims(1)=maxasm
      idims(2)=kdfuel
      idims(3)=npin
      idims(4)=npin

      dsetname = trim(group_name)//'initial_mass'
      call hwrite_double(file_id, dsetname, idims, pinload)

      deallocate (pinload)

! close CORE group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

      return
      end subroutine create_hdf
!=======================================================================
!
!  Write STATE block to HDF file
!
!=======================================================================
      subroutine write_hdf_state(file_id, nstate, iout)
      use hdf5
      use mod_hdftools
      use mod_geom, only : npin, maxasm, kdfuel
      use mod_state
      implicit none

      integer(hid_t), intent(in) :: file_id    ! HDF file identifier
      integer,        intent(in) :: nstate     ! statepoint number
      integer,        intent(in) :: iout       ! output unit

!--- local

      integer :: ia, ja, na, k

      integer :: ierror
      integer :: idims(10)           ! dataset dimensions

      integer(hid_t) :: group_id     ! HDF group identifier

      real(8), allocatable :: dpow(:,:,:,:)

      character(len=60) :: dsetname    ! HDF dataset name
      character(len=20) :: group_name  ! HDF group name

!--- edit statepoint values

      write (iout,*)
      write (iout,50) nstate
      write (iout,70) 'Percent Power     ', power
      write (iout,70) 'Percent Flow      ', flow
      write (iout,70) 'keff              ', xkeff
      write (iout,70) 'boron             ', boron
      write (iout,70) 'pressure          ', pressure,'psia'
      write (iout,70) 'exposure          ', ave_exp, 'GWD/MT'
      write (iout,70) 'exposure          ', efpd,    'EFPD'
      write (iout,*)

   50 format(' Writing Statepoint',i4)
   70 format(3x,a,f8.3,1x,a)

!--- create Statepoint

      group_name = '/STATE_0000/'
      write (group_name(8:11),'(i4.4)') nstate

      call h5gcreate_f (file_id, group_name, group_id, ierror)
      if (ierror.ne.0) stop 'error creating hdf group'

      dsetname = trim(group_name)//'keff'
      call hwrite_double_scalar(file_id, dsetname, xkeff)

      dsetname = trim(group_name)//'boron'
      call hwrite_double_scalar(file_id, dsetname, boron)

      dsetname = trim(group_name)//'pressure'
      call hwrite_double_scalar(file_id, dsetname, pressure)

      dsetname = trim(group_name)//'power'
      call hwrite_double_scalar(file_id, dsetname, power)

      dsetname = trim(group_name)//'flow'
      call hwrite_double_scalar(file_id, dsetname, flow)

      dsetname = trim(group_name)//'exposure'
      call hwrite_double_scalar(file_id, dsetname, ave_exp)

      dsetname = trim(group_name)//'exposure_efpd'
      call hwrite_double_scalar(file_id, dsetname, efpd)

   ! pin distributions

      allocate (dpow(maxasm,kdfuel,npin,npin))

   ! pin power

      idims(:)=0     ! clear
      idims(1)=maxasm
      idims(2)=kdfuel
      idims(3)=npin
      idims(4)=npin

!**** WARNING: need to confirm pin power order

!d    k=kdfuel/2
!d    na=maxasm/2
!d    write (*,*) 'sim debug pin power k=', k,'  na=', na
!d    do ja=1, npin
!d      write (*,245) ja, (pinpow(ia,ja,na,k),ia=1, npin)
!d    enddo
!d245 format (i3,20f8.3)

      do k=1, kdfuel
        do na=1, maxasm
          do ia=1, npin
            do ja=1, npin
               dpow(na,k,ja,ia)=pinpow(ia,ja,na,k)  ! transpose
            enddo
          enddo
        enddo
      enddo

!d    k=kdfuel/2
!d    na=maxasm/2
!d    write (*,*) 'hdf debug pin power k=', k,'  na=', na
!d    do ja=1, npin
!d      write (*,245) ja, (dpow(na,k,ja,ia),ia=1, npin)
!d    enddo

      dsetname = trim(group_name)//'pin_powers'
      call hwrite_double(file_id, dsetname, idims, dpow)

! debug
!d          write (*,*)
!d          write (*,*)  'HDF Pin Power Edits - Assembly ', na
!d          k=13
!d          na=maxasm/2
!d          do ja=1, npin
!d            write (*,320) ja, (dpow(na,k,ja,ia),ia=1,npin)
!d          enddo
!d    320   format(i3, (20f8.4))

   ! pin exposure

      idims(:)=0     ! clear
      idims(1)=maxasm
      idims(2)=kdfuel
      idims(3)=npin
      idims(4)=npin

!**** WARNING: need to confirm pin power order

      do k=1, kdfuel
        do na=1, maxasm
          do ia=1, npin
            do ja=1, npin
               dpow(na,k,ja,ia)=pinexp(ia,ja,na,k)  ! transpose
            enddo
          enddo
        enddo
      enddo

      dsetname = trim(group_name)//'pin_exposures'
      call hwrite_double(file_id, dsetname, idims, dpow)

      deallocate (dpow)

     ! close group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

      return
      end subroutine write_hdf_state
!=======================================================================
