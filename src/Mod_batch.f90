   module mod_batch
   use  hdf5
   use  mod_hdftools
   implicit none
!=======================================================================
!
!  Module to store batch edit information
!
!  Copyright (c) 2017 Core Physics, Inc.
!
!  Distributed under the MIT license.
!  See the LICENSE file in the main directory for details.
!
!=======================================================================
!
!  2017/04/25 - initial version to print pin exposure batch edits
!
!  Batch nbatch+1 is the core average
!
!  Remaining issues:
!
!  ** string read assumes that the string length is size 2
!  ** need to fix HDF reader
!
!  ** the average skips zero exposure pins, really need to mask with power, not exposure
!
!-----------------------------------------------------------------------

      integer :: nbatch

      integer, parameter :: maxbatch=20

      character(len=6) :: batchlabel(maxbatch)
      integer          :: batchcnt(maxbatch+1)
      real(8)          :: batch1(maxbatch+1)    ! batch ave
      real(8)          :: batch2(maxbatch+1)    ! batch min
      real(8)          :: batch3(maxbatch+1)    ! batch max

      integer, allocatable :: mapbatch(:,:)  ! (icore,jcore)

   contains

!=======================================================================
!
!  Subroutine to read assm_map and initialize list of unique assemblies (batches)
!
!=======================================================================

      subroutine batch_init(file_id, mapcore, icore, jcore)
      implicit none
      integer(hid_t), intent(in) :: file_id
      integer,        intent(in) :: icore, jcore
      integer,        intent(in) :: mapcore(icore,jcore)

      integer            :: i, j, k, kk, n
      integer            :: nassm
      integer            :: itype
      integer            :: ndim
      integer            :: idim(10)

      character(len=30) :: dataset    ! HDF dataset name
      character(len=6)  :: btmp
      character(len=2), allocatable :: assm_map(:) ! *** guess length 2

      nbatch=0
      batchlabel(:)=' '
      batchcnt  (:)=0

      batch1(:)=0.0d0
      batch2(:)=1.0d20
      batch3(:)=0.0d0

!--- read assembly map (i.e. batch labels)

      write (0,*) 'debug: reading assm_map'

      dataset='/INPUT/CASEID/CORE/assm_map'
      call h5info(file_id, dataset, itype, ndim, idim)
      write (*,*) 'ndim = ', ndim
      write (*,*) 'idim = ', idim(1:ndim)
      write (*,*) 'itype= ', itype
      if (ndim.ne.1) stop 'invalid dimensions in assm_map'
      nassm=idim(1)

      allocate (assm_map(nassm))
      assm_map=' '

      call read_string1d(file_id, dataset, assm_map, nassm)

!--- check for blank assembly names
!--- **** indicates that assm_map character length is not 2 ***

      do i=1, nassm
        if (assm_map(i).eq.' ') then
           write (*,*) 'missing assm_name - skipping batch edits'
           write (*,*) 'read_string1d assumes 2 character labels***'
           goto 999
        endif
      enddo

!--- calculate number of batches in core

      do i=1, nassm
        kk=0
        do k=1, nbatch
          if (assm_map(i).eq.batchlabel(k)) kk=k
        enddo
        if (kk.eq.0) then
          nbatch=nbatch+1
          kk=nbatch
          if (kk.gt.maxbatch) stop 'maxbatch exceeded'
          batchlabel(kk)=assm_map(i)
        endif
        batchcnt(kk)=batchcnt(kk)+1
      enddo

      do i=1, nbatch-1       ! sort
        do j=i+1, nbatch
          if (batchlabel(i).gt.batchlabel(j)) then
            btmp=batchlabel(i)
            batchlabel(i)=batchlabel(j)
            batchlabel(j)=btmp
            kk=batchcnt(i)
            batchcnt(i)=batchcnt(j)
            batchcnt(j)=kk
          endif
        enddo
      enddo

      write (*,*)
      write (*,*) 'Batch labels:'
      kk=0   ! sum
      do k=1, nbatch
        write (*,40) k, batchlabel(k), batchcnt(k)
        kk=kk+batchcnt(k)
      enddo
      batchcnt(nbatch+1)=kk   ! core average
      write (*,42) kk
   40 format (i6,1x,a,i4)
   42 format (7x,'sum   ',i4)

!--- create core map of batches for later use

      allocate (mapbatch(icore,jcore))
      mapbatch=0
      n=0
      do j=1, jcore
        do i=1, icore
          if (mapcore(i,j).gt.0) then
             n=n+1
             kk=0
             do k=1, nbatch
               if (assm_map(n).eq.batchlabel(k)) kk=k
             enddo
             if (kk.eq.0) stop 'error creating batch type map'
             mapbatch(i,j)=kk
          endif
        enddo
      enddo
      if (n.ne.nassm) stop 'error creating batch type map'

      write (*,*)
      write (*,*) 'map of batch numbers'
      do j=1, jcore
        write (*,'(30i3)') (mapbatch(i,j),i=1,icore)
      enddo

!--- finished

  999 continue
      deallocate (assm_map)

      return
      end subroutine batch_init
!=======================================================================
!
!  Subroutines to calculate batch statistics for a given distribution
!
!  Loop over full-core in order to get full-core statistics correct
!
!=======================================================================
      subroutine batchstat (npin, kd, nassm, icore, jcore, mapcore, axial, power)
      implicit none

      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: power(npin,npin,kd,nassm)

!--- local

      integer  :: i, j, k
      integer  :: ia, ja, na, nb
      real(8)  :: pp
      real(8)  :: zlen, zave      ! axial values
      real(8)  :: c3min, c3max    ! 3D values

!--- calculate 3D statistics

      batch1(:)=0.0d0
      batch2(:)=1.0d20
      batch3(:)=0.0d0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        zave=0.0d0       ! assembly values
        zlen=0.0d0
        c3min=1.0d20
        c3max=0.0d0

        do k=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                zave=zave+axial(k)*pp
                zlen=zlen+axial(k)
                if (pp.lt.c3min) c3min=pp
                if (pp.gt.c3max) c3max=pp
              endif
            enddo
          enddo
        enddo

        if (zlen.gt.0.0d0) then
          zave=zave/zlen
        else
          c3min=0.0d0    ! avoid overflow
        endif

        nb=mapbatch(ia,ja)
        batch1(nb)=batch1(nb)+zave         ! batch ave
        batch2(nb)=min(batch2(nb),c3min)   ! batch min
        batch3(nb)=max(batch3(nb),c3max)   ! batch max

      enddo    ! ia
      enddo    ! ja

      nb=nbatch+1   ! core average
      do i=1, nbatch
        batch1(nb)=batch1(nb)    +batch1(i)    ! core ave *** average of non-zero points
        batch2(nb)=min(batch2(nb),batch2(i))   ! core min
        batch3(nb)=max(batch3(nb),batch3(i))   ! core max
        batch1(i)=batch1(i)/batchcnt(i)
      enddo
      batch1(nb)=batch1(nb)/batchcnt(nb)

      return
      end subroutine batchstat
!=======================================================================
!
!  Subroutines to print final batch edits
!
!=======================================================================
      subroutine batchedit
      implicit none

!--- local

      integer  :: k, kk

!--- print batch statistics

      write (*,30)
      kk=0   ! sum
      do k=1, nbatch
        write (*,40) k, batchlabel(k), batchcnt(k), batch1(k), batch2(k), batch3(k)
        kk=kk+batchcnt(k)
      enddo

      k=nbatch+1
      batch1(k)=batch1(k)/batchcnt(k)
      write (*,42)  batchcnt(k), batch2(k), batch3(k)
  ! *** don't print core average batch 1 - it is average of non-zero points

   30 format (/,' Batch edits:           pin_exposure', &
              /,'   Batch    nassm     ave    minpin   maxpin')

   40 format (i6,1x,a,i4, 3f9.4)
   42 format ('  core ave',3x,i4,9x, 3f9.4)

      return
      end subroutine batchedit
!=======================================================================


   end module mod_batch


