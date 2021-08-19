   subroutine write_hdf_state(file_id, nstate, iout)
!=======================================================================
!
!  Write VERA STATE block to HDF file
!
!=======================================================================
      use hdf5
      use mod_hdftools
      use mod_geom, only : npin, maxasm, kdfuel
      use mod_state
      implicit none

      integer(hid_t), intent(in) :: file_id    ! HDF file identifier
      integer,        intent(in) :: nstate     ! statepoint number
      integer,        intent(in) :: iout       ! output unit

!--- local

      integer :: ierror
      integer :: idims(10)           ! dataset dimensions

      integer(hid_t) :: group_id     ! HDF group identifier

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

   ! pin power

      idims(:)=0     ! clear
      idims(1)=maxasm
      idims(2)=kdfuel
      idims(3)=npin
      idims(4)=npin

      dsetname = trim(group_name)//'pin_powers'
      call hwrite_double(file_id, dsetname, idims, pinpow)

   ! pin exposure

      idims(:)=0     ! clear
      idims(1)=maxasm
      idims(2)=kdfuel
      idims(3)=npin
      idims(4)=npin

      dsetname = trim(group_name)//'pin_exposures'
      call hwrite_double(file_id, dsetname, idims, pinexp)

!  Other possible distributions:
!    'pin_fuel_temp'
!    'pin_max_clad_surface_temp'
!    'pin_mod_temps'
!    'pin_mod_dens'
!    'pin_steamrate'
!    'pin_cool_dens'
!    'pin_cool_temp'


!--- close group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

      return
      end subroutine write_hdf_state
!=======================================================================
