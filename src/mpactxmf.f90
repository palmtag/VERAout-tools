!=======================================================================
!
!  Program to read an existing MPACT HDF output file and
!   generate a mesh file and XMF file for plotting with VISIT or PARAVIEW
!
!  Copyright (c) 2014-2015 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2014/07/18 - original file based on ctfxmf
!  2014/07/24 - moved grid to separate HDF file
!  2014/07/25 - added support for multiple statepoints
!  2014/07/26 - create multiple xmf files when multiple statepoints
!                are found so we can generate movies
!
!  Usage:
!    mpactxmf.exe [file.h5]
!
!  This will create two new files:
!     file.grid.h5
!     file.h5.xmf
!
!  If there are more than one statepoint, multiple xmf files will be generated
!
!  Open the xmf file in either visit or paraview
!
!=======================================================================
      program mpactxmf

      use hdf5
      use mod_hdftools

      implicit none

      character(len=120):: filename      ! HDF File name
      character(len=120):: filectf       ! HDF CTF File name
      character(len=120):: gridfile      ! HDF grid File name
      character(len=40) :: dataset       ! dataset name

      integer(hid_t)    :: file_id       ! File identifier
      integer           :: idims(10)     ! Dataset dimensions (use with h5info)
      integer           :: ndim

      integer :: ll
      integer :: ia
      integer :: itype
      integer :: ierror

      logical :: ifxst

      integer :: naxial
      integer :: nassm
      integer :: npin
      integer :: nnrod    ! number of nodes in rod mesh
      integer :: nstate   ! number of statepoints

      real(8) :: xexp

      character(len=12) :: state_name

!--- Initialize HDF fortran interface.

      call h5open_f(ierror)
      if (ierror.ne.0) stop 'ierror: h5open_f'

!--- read existing HDF data filename from command line

      ia=iargc()
      if (ia.ne.1) stop 'usage: mpactxmf [file.h5]'
      call getarg(1,filename)

      ll=len_trim(filename)
      if (filename(ll-2:ll).ne.'.h5') then
        write (0,*) 'filename must end with .h5'
        stop 'invalid HDF filename'
      endif

      gridfile=filename(1:ll-3)//'.grid.h5'

!--- look for corresponding CTF filename

      filectf=' '
      ll=len_trim(filename)
      if (filename(ll-2:ll).eq.'.h5') filectf=filename(1:ll-3)//'.ctf.h5'
      inquire (file=filectf,exist=ifxst)
      if (ifxst) then
         write (*,*) '  *** found ctf file ', trim(filectf)
      else
         write (*,*) '  *** NO ctf file found ', trim(filectf)
         filectf=' '      ! reset if not found
      endif

!--- open HDF data file

      write (*,'(2a)') 'reading h5 file: ', trim(filename)
      inquire(file=filename, exist=ifxst)
      if (.not.ifxst) then
        write (*,'(3a)') 'error: input file ',trim(filename),' does not exist'
        stop
      endif
      call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(filename), ' could not be opened'
        stop
      endif

!--- check if first statepoint exists (it should always exist)

      nstate=1
      call make_state_name(state_name, nstate)

      call h5lexists_f(file_id, state_name, ifxst, ierror)
      if (.not.ifxst) then
         stop 'statepoint data does not exist on file'
      endif

!--- read pin dimensions from pin power array

      dataset=trim(state_name)//'pin_powers'
      call h5info(file_id, dataset, itype, ndim, idims)

!d    write (*,*) 'idim: ', idims(1:ndim)

      nassm =idims(1)         ! order (nassm, naxial, npin, npin)
      naxial=idims(2)
      npin  =idims(4)

      write (*,20) 'ndim  ', ndim
      write (*,20) 'nassm ', nassm
      write (*,20) 'naxial', naxial
      write (*,20) 'npin  ', idims(3)
      write (*,20) 'npin  ', npin
   20 format (' debug pin: ', a, '=', i6)

      if (ndim.ne.4) stop 'invalid dimensions in pin data'
      if (npin.eq.0) stop 'invalid number of pins in pin data'
      if (idims(3).ne.idims(4)) stop 'invalid npin in pin data'

!--- write XMFMESH data to a separate grid file

      call coordinates(file_id, gridfile, nassm, naxial, npin, nnrod)

!--- find number of statepoints on file

      do

! read exposure value at this statepoint

        dataset=trim(state_name)//'exposure_efpd'
        call hdf5_read_double  (file_id, dataset, xexp)
        write (*,'(2a,f12.5)') trim(dataset), ' exposure ', xexp

! write xmf file

        call write_xmf(npin, npin, naxial, nassm, nnrod, filename, gridfile, filectf, nstate, xexp)

! check if the next state exists

        call make_state_name(state_name, nstate+1)
        call h5lexists_f(file_id, state_name, ifxst, ierror)
        if (.not.ifxst) exit
        nstate=nstate+1

      enddo
      write (*,*) 'nstate=', nstate

!--- close HDF file

      call h5fclose_f(file_id, ierror)
      write (*,'(2a)') ' finished reading HDF file: ', trim(filename)

!--- close fortran interface.

      call h5close_f(ierror)

!-----------------------------

      end
!=======================================================================
!
!  Subroutine to create a statepoint string in the form '/STATE_0000/'
!
!=======================================================================
      subroutine make_state_name(state_name, nstate)
      implicit none
      integer,          intent(in)  :: nstate
      character(len=*), intent(out) :: state_name

      state_name='/STATE_0000/'
      write (state_name(8:11),'(i4.4)') nstate

      return
      end subroutine make_state_name
!=======================================================================
!
!  Subroutine to write grid coordinates to XMFMESH group
!  Needed for xdmf plotting
!
!=======================================================================

      subroutine coordinates(file_id, gridfile, nassm, naxial, npin, nnrod)
      use hdf5
      use mod_hdftools
      implicit none
      integer(hid_t),   intent(in)  :: file_id       ! File identifier
      character(len=*), intent(in)  :: gridfile      ! grid file name
      integer,          intent(in)  :: nassm, naxial, npin
      integer,          intent(out) :: nnrod         ! output number of rod grids

!--- local

      integer :: i, j, k
      integer :: na              ! assembly number
      integer :: nn
      integer :: iasm, jasm      ! number of assemblies in x and y direction
      integer :: nkoff           ! axial node offset
      integer :: npl             ! number of nodes in a plane
      integer :: nnodes
      integer :: ia, ip, ja, jp
      integer :: ierror

      real(8) :: delta

      real(8), allocatable :: zbounds(:)      ! axial bounds (naxial+1)

      real,    allocatable :: xyz(:,:)        ! single precision
      real(8), allocatable :: xi(:)
      real(8), allocatable :: yi(:)
      integer, allocatable :: gg(:,:,:,:,:)   ! map nodes to cells
      integer, allocatable :: ijk(:,:,:,:)
      integer, allocatable :: mapassm(:,:)    ! assembly map
      integer, allocatable :: itemp(:,:)
      integer, allocatable :: nodemap(:,:)

      character(len=12) :: group_name
      character(len=40) :: dataset       ! Dataset name
      integer           :: idims(10)     ! Dataset dimensions
      integer(hid_t)    :: grid_id       ! Grid File identifier
      integer(hid_t)    :: group_id      ! Group identifier

      integer :: ndim
      integer :: itype
      integer :: itmp   ! temporary storage
      integer :: isym   ! core symmetry
      integer :: ioff   ! core symmetry

      logical :: ifxst

      nnodes=0
      nnrod=0

      delta=1.26d0  ! ************ temp fix: assumed pin pitch

!--- check if CORE exists, if not, do special fixup

      dataset='/CORE'
      call h5lexists_f(file_id, dataset, ifxst, ierror)
      if (.not.ifxst) then
        write (0,*) '*** CORE group does not exist - performing fixup ***'

        isym=1
        iasm=1
        jasm=1

        if (nassm.ne.1) then
          write (*,*) 'CORE fixup will only work with one assembly'
          stop 'CORE fixup will only work with one assembly'
        endif
        if (naxial.ne.1) then
          write (*,*) 'CORE fixup will only work with 2D problems'
          stop 'CORE fixup will only work with 2D problems'
        endif

        allocate (mapassm(iasm,jasm))    ! define assembly map
        mapassm(1,1)=1

        allocate (zbounds(0:naxial))
        zbounds(0)=0.0d0
        zbounds(1)=1.0d0

        goto 30
      endif

!--- read core symmetry

      dataset='/CORE/core_sym'
      call hdf5_read_integer (file_id, dataset, isym)

      write (*,*) 'core symmetry = ', isym

!--- read core map

      dataset='/CORE/core_map'
      call h5info(file_id, dataset, itype, ndim, idims)

      if (ndim.ne.2) then
        write (*,*) 'debug: ndim=', ndim
        write (*,*) 'debug: core_map dims ', idims(1:ndim)
        stop 'invalid core map dimensions'
      endif

      jasm=idims(1)   ! C-order
      iasm=idims(2)

      if (iasm.ne.jasm) stop 'only square cores are supported'

      allocate (mapassm(iasm,jasm))    ! define assembly map
      mapassm=0
      call hdf5_read_integer (file_id, dataset, jasm, iasm, mapassm)

!  map is in "C-order" so we need to transpose into "Fortran-order"

      na=0
      do j=1, jasm
        do i=1, iasm
           if (i.ne.j) then
             itmp=mapassm(i,j)      ! transpose
             mapassm(i,j)=mapassm(j,i)
             mapassm(j,i)=itmp
           endif
           na=max(na,mapassm(i,j))
        enddo
      enddo

      if (na.ne.nassm) stop 'error in assembly map numbers'

      write (*,*) 'full-core assembly map'
      do j=1, jasm
        write (*,'(3x,50i4)') (mapassm(i,j),i=1, iasm)
      enddo

!--- assembly map on file is always full-core
!--- decrease map size if qtr-symmetry

      if (isym.eq.4) then

        jasm=(iasm+1)/2    ! core must be square - reduce j but keep i the same
        ioff=iasm-jasm

        allocate (itemp(jasm,jasm))
        do j=1, jasm
          do i=1, jasm
            itemp(i,j)=mapassm(i+ioff,j+ioff)
          enddo
        enddo
        deallocate (mapassm)
        iasm=jasm   ! decrease i also
        allocate (mapassm(iasm,jasm))
        mapassm=itemp
        deallocate (itemp)

        write (*,*) 'qtr-core assembly map'
        do j=1, jasm
          write (*,'(3x,50i4)') (mapassm(i,j),i=1, iasm)
        enddo

      endif

!--- allocate and read axial cell heights

      allocate (zbounds(0:naxial))

      dataset='/CORE/axial_mesh'
      call hdf5_read_double(file_id, dataset, naxial+1, k, zbounds)
      if (k.ne.naxial+1) stop 'invalid axial_mesh size'

      write (*,*) 'debug: axial elevations [cm]'
      do k=naxial, 0, -1
         write (*,'(i4,f12.6)') k, zbounds(k)
      enddo

  30 continue   ! end of fix-up if CORE does not exist

!-----------------------------------------------------------------------
!    Create Rod Mesh
!    (CTF files will have both rod mesh and channel mesh)
!-----------------------------------------------------------------------

!--- create 2D map of nodes - needed for ragged cores with missing assemblies on edges

      allocate (nodemap(0:npin*jasm,0:npin*iasm))   ! needs to be in C-order
      nodemap=0

      j=0
      do ja=1, jasm
        do jp=1, npin
          j=j+1

          i=0
          do ia=1, iasm
            do ip=1, npin
              i=i+1

              na=mapassm(ia,ja)
              if (na.gt.0) then
                nodemap(j-1,i-1)=1     ! flag node as used, we will number later
                nodemap(j-1,i  )=1
                nodemap(j  ,i-1)=1
                nodemap(j  ,i  )=1
              endif

            enddo  ! ip
          enddo    ! ia
        enddo  ! jp
      enddo    ! ja

!d    write (*,*) 'rod node map: (before numbering)'
!d    do j=0, npin*jasm
!d      write (*,'(2x,500i1)') (nodemap(j,i),i=0,npin*iasm)
!d    enddo

!  Number the nodes
!  nodes are across, then down

      nn=0
      do i=0, npin*iasm
        do j=0, npin*jasm
          if (nodemap(j,i).gt.0) then
            nn=nn+1
            nodemap(j,i)=nn
          endif
        enddo
      enddo

      write (*,*) 'number of rod nodes in one plane = ', nn

      npl=nn    ! save number of nodes in a plane (needed for offset later)

      if (nn.lt.350) then
        write (*,*) 'rod node map: (after numbering)'
        do j=0, npin*jasm
          write (*,'(2x,100i3)') (nodemap(j,i),i=0,npin*iasm)
        enddo
      endif

!--- write pin dimensions

!
!  (i,j) origin is in top left corner, i goes across, j goes down
!  (x,y) origin is in bottom left corner,  x goes across left to right, y goes up
!
!  mapping:  j goes down from top to bottom   (negative y direction)
!            i goes across from left to right (positive x direction)
!

! allocate pointer arrays

      nnodes=npl*(naxial+1)
      write (*,*) 'number of 3D nodes = ', nnodes

      nnrod=nnodes    ! save for output

      allocate (xyz(3,nnodes))
      allocate (gg(8,nassm,naxial,npin,npin))   ! 8 node numbers per cell
      xyz=-1000.0d0
      gg=-100

!--- fill xyz values for each node

      allocate (xi(0:npin*iasm))
      allocate (yi(0:npin*jasm))
      do i=0, npin*iasm
        xi(i)=i*delta              ! x aligns with i
      enddo
      do j=0, npin*iasm
        yi(j)=(npin*iasm-j)*delta   ! y aligns with j, but is the reverse direction
      enddo

      nkoff=0
      do k=0, naxial
        do j=0, npin*jasm
          do i=0, npin*iasm
            nn=nodemap(j,i)
            if (nn.gt.0) then
              xyz(1,nn+nkoff)=real(xi(i))
              xyz(2,nn+nkoff)=real(yi(j))
              xyz(3,nn+nkoff)=real(zbounds(k))
            endif
          enddo
        enddo
        nkoff=nkoff+npl
      enddo

      deallocate (xi,yi)

!--- fill element to node map (gg)

      nkoff=0

      do k=1, naxial

        j=0
        do ja=1, jasm
          do jp=1, npin
            j=j+1

            i=0
            do ia=1, iasm
              do ip=1, npin
                i=i+1

                na=mapassm(ia,ja)
                if (na.gt.0) then

          ! I'm using (ip,jp) instead of (jp,ip) in order to get the (i,j)
          ! coordinates to show up right in the plots

                  gg(1,na,k,ip,jp)=nodemap(j,  i-1)+nkoff      ! 1,0,0
                  gg(2,na,k,ip,jp)=nodemap(j,  i  )+nkoff      ! 1,1,0
                  gg(3,na,k,ip,jp)=nodemap(j-1,i  )+nkoff      ! 0,1,0
                  gg(4,na,k,ip,jp)=nodemap(j-1,i-1)+nkoff      ! 0,0,0

                  gg(5,na,k,ip,jp)=gg(1,na,k,ip,jp)+npl        ! 1,0,1
                  gg(6,na,k,ip,jp)=gg(2,na,k,ip,jp)+npl        ! 1,1,1
                  gg(7,na,k,ip,jp)=gg(3,na,k,ip,jp)+npl        ! 0,1,1
                  gg(8,na,k,ip,jp)=gg(4,na,k,ip,jp)+npl        ! 0,0,1

                endif  ! na

              enddo  ! ip
            enddo    ! ia

          enddo   ! jp
        enddo     ! ja

        nkoff=nkoff+npl    ! increment axial offset

      enddo  ! k

!--- sanity check coordinates and convert to zero based numbering

      do ip=1, npin
        do jp=1, npin
          do k=1, naxial
            do na=1, nassm
              do i=1, 8
                if (gg(i,na,k,ip,jp).lt.1)      stop 'invalid node number .lt.1'
                if (gg(i,na,k,ip,jp).gt.nnodes) stop 'invalid node number .gt.max'
                gg(i,na,k,ip,jp)=gg(i,na,k,ip,jp)-1
              enddo
            enddo
          enddo
        enddo
      enddo

!--- open HDF grid file

      call h5fcreate_f (gridfile, H5F_ACC_TRUNC_F, grid_id, ierror)
      if (ierror.lt.0) then
        write (*,'(3a)') 'error: H5 input file ',trim(gridfile), ' could not be created'
        stop
      endif

!--- create HDF group for XMFMESH data

      group_name='/XMFMESH/'
      write (*,*) 'create group: ', group_name

      call h5gcreate_f (grid_id, group_name, group_id, ierror)
      if (ierror.ne.0) stop 'error creating group'

!--- write mapping to HDF

      idims(:)=0     ! clear
      idims(1)=8
      idims(2)=nassm
      idims(3)=naxial
      idims(4)=npin
      idims(5)=npin

      dataset = trim(group_name)//'ggrod'
      call hwrite_integer(grid_id, dataset, idims, gg)


      idims(:)=0     ! clear
      idims(1)=3
      idims(2)=nnodes

      dataset = trim(group_name)//'xyzrod'
      call hwrite_real (grid_id, dataset, idims, xyz)

      deallocate (gg)
      deallocate (nodemap)
      deallocate (xyz)

!--- write (i,j,k) data for debug

      allocate (ijk(nassm,naxial,npin,npin))
      ijk=0

      do i=1, npin
        do j=1, npin
          do k=1, naxial
            do na=1, nassm
              ijk(na,k,j,i)=i
            enddo
          enddo
        enddo
      enddo

      idims(:)=0     ! clear
      idims(1)=nassm
      idims(2)=naxial
      idims(3)=npin
      idims(4)=npin

      dataset = trim(group_name)//'ipin'
      call hwrite_integer (grid_id, dataset, idims, ijk)

      do i=1, npin
        do j=1, npin
          do k=1, naxial
            do na=1, nassm
              ijk(na,k,j,i)=j
            enddo
          enddo
        enddo
      enddo

      dataset = trim(group_name)//'jpin'
      call hwrite_integer (grid_id, dataset, idims, ijk)

      do i=1, npin
        do j=1, npin
          do k=1, naxial
            do na=1, nassm
              ijk(na,k,j,i)=k
            enddo
          enddo
        enddo
      enddo

      dataset = trim(group_name)//'kcoord'
      call hwrite_integer (grid_id, dataset, idims, ijk)

      do i=1, npin
        do j=1, npin
          do k=1, naxial
            do na=1, nassm
              ijk(na,k,j,i)=na
            enddo
          enddo
        enddo
      enddo

      dataset = trim(group_name)//'nassm'
      call hwrite_integer (grid_id, dataset, idims, ijk)

      deallocate (ijk)

!-----------------------------------------------------------------------

      deallocate (mapassm)

!--- close XMFMESH group

      call h5gclose_f (group_id, ierror)
      if (ierror.ne.0) stop 'error closing group'

!--- close grid HDF file

      call h5fclose_f(grid_id, ierror)

      write (*,'(2a)') ' finished writing HDF grid file: ', trim(gridfile)

!--- return

      return
      end subroutine coordinates

!=======================================================================
!
!  Subroutine to write XMF file
!
!=======================================================================
      subroutine write_xmf(nx, ny, nz, nassm, nnrod, filename, gridname, filectf, nstate, xexp)
      implicit none
      integer, intent(in) :: nx, ny, nz, nassm, nnrod
      integer, intent(in) :: nstate    ! statepoint number
      real(8), intent(in) :: xexp      ! exposure
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: filectf
      character(len=*), intent(in) :: gridname

      character(len=12) :: state_name
      character(len=10) :: label_xmf

      integer :: io=10
      integer :: idis

      integer, parameter :: maxdist=6    ! maximum number of distributions
      character(len=20) :: dist_label(maxdist)

      integer, parameter :: maxctf=10    ! maximum number of CTF distributions
      character(len=40) :: ctf_label(maxctf)

      dist_label(1)='pin_powers'
      dist_label(2)='pin_fueltemps'
      dist_label(3)='pin_cladtemps'
      dist_label(4)='pin_modtemps'
      dist_label(5)='pin_moddens'
      dist_label(6)='pin_exposures'

      ctf_label( 1)="Rod_Surface_Temp_NE_Quad [C]"
      ctf_label( 2)="Rod_Surface_Temp_NW_Quad [C]"
      ctf_label( 3)="Rod_Surface_Temp_SE_Quad [C]"
      ctf_label( 4)="Rod_Surface_Temp_SW_Quad [C]"
      ctf_label( 5)="Steaming_Rate_NE_Quad [kg_per_s]"
      ctf_label( 6)="Steaming_Rate_NW_Quad [kg_per_s]"
      ctf_label( 7)="Steaming_Rate_SE_Quad [kg_per_s]"
      ctf_label( 8)="Steaming_Rate_SW_Quad [kg_per_s]"
      ctf_label( 9)="pin_fueltemps [C]"
      ctf_label(10)="pin_powers"     ! no units

      write (label_xmf,'(".",i0,".xmf")') nstate

      write (0,'(3a)') ' writing XMF file: ',trim(filename)//trim(label_xmf)

      open (io,file=trim(filename)//trim(label_xmf))

!--- write xml header

      write (io,20) '<?xml version="1.0" ?>'
      write (io,20) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write (io,20) '<Xdmf Version="2.0">'
      write (io,20) ' <Domain>'

!-----------------------------------------------------------------------
!--- write rod grid - hexahedron

      write (io,120) 'rod_mesh', nx, ny, nz, nassm, nx, ny, nz, nassm, trim(gridname), '/XMFMESH/ggrod'

  120 format ( '  <Grid Name="',a,'" GridType="Uniform">', &
            /,5x,'<Topology TopologyType="Hexahedron" Dimensions="',i0,1x,i0,1x,i0,1x,i0,'">', &
            /,6x,'  <DataItem Dimensions="',i0,1x,i0,1x,i0,1x,i0,' 8" NumberType="Int" Format="HDF">', &
            /,8x,a,':',a, &
            /,6x,'  </DataItem>', &
            /,5x,'</Topology>')

      write (io,125) xexp
  125 format (5x,'<Time Value="',f0.3,'" />')

!--- rod node coordinates

      write (io,130) nnrod, trim(gridname), '/XMFMESH/xyzrod'

  130 format (5x,'<Geometry Type="XYZ">', &
            /,8x,'<DataStructure Name="xyz" DataType="Float" Precision="4" ', &
                      'Dimensions="',i0,' 3" Format="HDF">', &
            /,8x,a,':',a, &
            /,8x,'</DataStructure>', &
            /,5x,'</Geometry>')

!--- rod cell data

!  the following four datasets are for debug purposes and can be removed

      if (nstate.eq.1) then

        write (io,80) 'ipin'
        write (io,32) nx, ny, nz, nassm        ! 32=integer array
        write (io,10) trim(gridname),'/XMFMESH/','ipin'
        write (io,20) '        </DataItem>'
        write (io,20) '     </Attribute>'

        write (io,80) 'jpin'
        write (io,32) nx, ny, nz, nassm        ! 32=integer array
        write (io,10) trim(gridname),'/XMFMESH/','jpin'
        write (io,20) '        </DataItem>'
        write (io,20) '     </Attribute>'

        write (io,80) 'kcoord'
        write (io,32) nx, ny, nz, nassm        ! 32=integer array
        write (io,10) trim(gridname),'/XMFMESH/','kcoord'
        write (io,20) '        </DataItem>'
        write (io,20) '     </Attribute>'

        write (io,80) 'nassm'
        write (io,32) nx, ny, nz, nassm        ! 32=integer array
        write (io,10) trim(gridname),'/XMFMESH/','nassm'
        write (io,20) '        </DataItem>'
        write (io,20) '     </Attribute>'

      endif

!  pin data

      call make_state_name(state_name, nstate)

      do idis=1, maxdist

        write (io,80) trim(dist_label(idis))
        write (io,30) nx, ny, nz, nassm
        write (io,10) trim(filename),trim(state_name),trim(dist_label(idis))
        write (io,20) '       </DataItem>'
        write (io,20) '     </Attribute>'

      enddo

!--- CTF data

      if (filectf.ne.' ') then
         write (*,*) ' adding CTF data'
         do idis=1, maxctf
           write (io,80) trim(ctf_label(idis))
           write (io,30) nx, ny, nz, nassm
           write (io,10) trim(filectf),trim(state_name),trim(ctf_label(idis))
           write (io,20) '       </DataItem>'
           write (io,20) '     </Attribute>'
         enddo
      endif



!--- close grid

      write (io,20) '  </Grid>'

  10  format (8x,a,':',2a)
  20  format (a)

  30  format (8x,'<DataItem Dimensions="',i0,1x,i0,1x,i0,1x,i0,'" NumberType="Float" ', &
                       'Precision="8" Format="HDF">')

  32  format (8x,'<DataItem Dimensions="',i0,1x,i0,1x,i0,1x,i0,'" NumberType="Int" Format="HDF">')

  80  format (5x,'<Attribute Name="',a,'" Center="Cell">')

!-----------------------------------------------------------------------
!--- finished

      write (io,20) ' </Domain>'
      write (io,20) '</Xdmf>'
      close (io)

      return
      end
!=======================================================================
!=======================================================================
