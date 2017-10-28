   module mod_batch
!=======================================================================
!
!  Module to store batch edit information
!  *** still in development
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
!-----------------------------------------------------------------------
   use  hdf5
   use  mod_hdftools
   implicit none

      integer :: nbatch

      integer, parameter :: maxbatch=20

      character(len=6) :: batchlabel(maxbatch)
      integer          :: batchcnt(maxbatch+1)

      real(8)          :: batch1(maxbatch+1)    ! pinexp batch ave
      real(8)          :: batch2(maxbatch+1)    ! pinexp batch min
      real(8)          :: batch3(maxbatch+1)    ! pinexp batch max

      real(8)          :: batch4(maxbatch+1)    ! pinpow batch ave
      real(8)          :: batch5(maxbatch+1)    ! pinpow batch min
      real(8)          :: batch6(maxbatch+1)    ! pinpow batch max

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
      batch4(:)=0.0d0
      batch5(:)=1.0d20
      batch6(:)=0.0d0

!--- read assembly map (i.e. batch labels)

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
      subroutine batchstat (npin, kd, nassm, icore, jcore, mapcore, axial, pinexp, power)
      implicit none

      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: icore, jcore
      integer, intent(in) :: mapcore(icore,jcore)
      real(8), intent(in) :: axial(kd)
      real(8), intent(in) :: pinexp(npin,npin,kd,nassm)
      real(8), intent(in) :: power(npin,npin,kd,nassm)

!--- local

      integer  :: i, j, k, kpin, kpinsave
      integer  :: ia, ja, na, nb
      real(8)  :: ee, pp
      real(8)  :: zlen                  ! axial length
      real(8)  :: eave, e3min, e3max    ! 3D exp values
      real(8)  :: pave, p3min, p3max    ! 3D pow values

!--- calculate 3D statistics

      if (nbatch.eq.0) return

      kpinsave=0

      batch1(:)=0.0d0
      batch2(:)=1.0d20
      batch3(:)=0.0d0
      batch4(:)=0.0d0
      batch5(:)=1.0d20
      batch6(:)=0.0d0

      do ja=1, jcore     ! loop over assemblies in full-core
        do ia=1, icore   ! loop over assemblies in full-core
        na=mapcore(ia,ja)
        if (na.eq.0) cycle

        zlen=0.0d0
        eave=0.0d0       ! assembly values
        e3min=1.0d20
        e3max=0.0d0
        pave=0.0d0       ! assembly values
        p3min=1.0d20
        p3max=0.0d0

        kpin=0
        do k=1, kd
          do j=1, npin
            do i=1, npin
              pp=power(i,j,k,na)
              if (pp.gt.0.0d0) then
                zlen=zlen+axial(k)
                ee=pinexp(i,j,k,na)
                eave=eave+axial(k)*ee
                if (ee.lt.e3min) e3min=ee
                if (ee.gt.e3max) e3max=ee
                pave=pave+axial(k)*pp
                if (pp.lt.p3min) p3min=pp
                if (pp.gt.p3max) p3max=pp
                kpin=kpin+1
              endif
            enddo
          enddo
        enddo

        if (kpinsave.eq.0) kpinsave=kpin
        if (kpin.ne.kpinsave) then
           write (*,*) kpin
           write (*,*) kpinsave
           write (*,*) 'WARNING: number of 3D pins in a bundle is changing'
        endif
!d      write (*,*) 'debug: ia, ja, kpin', ia, ja, kpin, real(kpin)/kd

        if (zlen.gt.0.0d0) then
          eave=eave/zlen
          pave=pave/zlen
        else
          e3min=0.0d0    ! avoid overflow
          p3min=0.0d0    ! avoid overflow
        endif

        nb=mapbatch(ia,ja)
        if (nb.eq.0) stop 'zero batch number found'
        batch1(nb)=    batch1(nb)+eave     ! batch ave
        batch2(nb)=min(batch2(nb),e3min)   ! batch min
        batch3(nb)=max(batch3(nb),e3max)   ! batch max
        batch4(nb)=    batch4(nb)+pave     ! batch ave
        batch5(nb)=min(batch5(nb),p3min)   ! batch min
        batch6(nb)=max(batch6(nb),p3max)   ! batch max

      enddo    ! ia
      enddo    ! ja

      nb=nbatch+1   ! core average
      do i=1, nbatch
        batch1(nb)=    batch1(nb)+batch1(i)    ! core ave
        batch2(nb)=min(batch2(nb),batch2(i))   ! core min
        batch3(nb)=max(batch3(nb),batch3(i))   ! core max
        batch4(nb)=    batch4(nb)+batch4(i)    ! core ave
        batch5(nb)=min(batch5(nb),batch5(i))   ! core min
        batch6(nb)=max(batch6(nb),batch6(i))   ! core max
        batch1(i)=batch1(i)/batchcnt(i)
        batch4(i)=batch4(i)/batchcnt(i)
      enddo
      batch1(nb)=batch1(nb)/batchcnt(nb)
      batch4(nb)=batch4(nb)/batchcnt(nb)

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

      if (nbatch.eq.0) return

      write (*,30)
      kk=0   ! sum
      do k=1, nbatch
        write (*,40) k, batchlabel(k), batchcnt(k), batch1(k), batch2(k), batch3(k), batch4(k), batch5(k), batch6(k)
        kk=kk+batchcnt(k)
      enddo

      k=nbatch+1    ! core average values
      write (*,42)  batchcnt(k), batch1(k), batch2(k), batch3(k), batch4(k), batch5(k), batch6(k)

   30 format (/,' Batch edits:      ------ pin_exposure -----      ------- pin_power -------', &
              /,'   Batch    nassm     ave    minpin   maxpin         ave    minpin   maxpin')

   40 format (i6,1x,a,i4, 3f9.4,4x,3f9.4)
   42 format ('  core ave',3x,i4,3f9.4,4x,3f9.4)

      write (*,*)
      write (*,*) '*** pin exposure averages are not weighted by mass'

      return
      end subroutine batchedit
!=======================================================================

   end module mod_batch


