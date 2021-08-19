   subroutine write_hdf_core(hdfname, file_id, iout)
!=======================================================================
!
!  Subroutine to create HDF file and write VERA CORE block
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

   55   format(3x,a,i8)
   59   format(3x,a,f8.3)
   65   format(3x,a,f8.2,1x,a,' (full core)')
!  70   format(3x,a,f8.3)

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

!  pin loadings

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

!--- close CORE group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

      return
      end subroutine write_hdf_core
