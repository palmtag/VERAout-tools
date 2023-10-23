   module mod_nea
   implicit none
!=======================================================================
!
!  Special edits to write pin powers to CSV files which can be
!  imported into NEA Benchmark spreadsheet
!
!  Open one CSV file per assembly
!  Write pin powers for multiple exposures
!
!=======================================================================
  
   contains

!=======================================================================
!
!  Subroutine to open CSV files - one per detector location
!
!=======================================================================
      subroutine nea_open()
      use mod_coregeom, only : icore, jcore, xlabel, detector_map
      implicit none

      integer :: i, j
      integer :: iunit
      integer :: ndet
      character(len=10) :: fname

      write (*,*)
      write (*,*) 'Special NEA benchmark edits'

      write (*,'(/,a)') ' Detector Map:'
      do j=1, jcore
         write (*,'(2x,20i3)') (detector_map(i,j),i=1,icore)
      enddo

      ndet=0
      do j=1, jcore
        do i=1, icore
!          if (mapcore(i,j).gt.0) then
           if (detector_map(i,j).gt.0) then
             ndet=ndet+1
             write (fname,32) xlabel(i), j
             write (*,20) ndet, i, j, trim(fname)
             iunit=100+ndet
             open (iunit,file=trim(fname)//'.csv')
           endif
        enddo
      enddo
  32  format ('Assy_',a1,i0)
  20  format (i5,' (',i2,',',i2,')  ',a)

      write (*,*) 'number of detectors = ', ndet

      if (ndet.ne.58) then
         write (*,*) '****** WARNING: number of detectors does not appear correct *******'
         write (*,*) '****** WARNING: if this is qtr-core, the map is rotated     *******'
      endif

      return
      end subroutine nea_open
!=======================================================================
!
!  Subroutine to close CSV files - one per detector location
!
!=======================================================================
      subroutine nea_close()
      use mod_coregeom, only : icore, jcore, detector_map
      implicit none

      integer :: i, j
      integer :: iunit
      integer :: ndet

      ndet=0
      do j=1, jcore
        do i=1, icore
!          if (mapcore(i,j).gt.0) then
           if (detector_map(i,j).gt.0) then
             ndet=ndet+1
             iunit=100+ndet
             close (iunit)
           endif
        enddo
      enddo

      write (*,*)
      write (*,*) 'closing NEA Benchmark CSV files'
      write (*,*) 'number of detectors = ', ndet

      if (ndet.ne.58) then
         write (*,*) '****** WARNING: number of detectors does not appear correct *******'
         write (*,*) '****** WARNING: if this is qtr-core, the map is rotated     *******'
      endif

      return
      end subroutine nea_close
!=======================================================================
!
!  Subroutine to print 3D pin powers by assembly using pretty output
!
!=======================================================================
      subroutine print_nea(efpd, title, nstate, npin, kd, nassm, power)
      use mod_coregeom, only : icore, jcore, mapcore, detector_map, axial, axial0
      implicit none
      integer, intent(in) :: npin, kd, nassm
      integer, intent(in) :: nstate
      real(8), intent(in) :: efpd
      real(8), intent(in) :: power(npin,npin,kd,nassm)
      character(len=*), intent(in) :: title

!--- local

      integer :: ia, ja, k, na
      integer :: i, j
      integer :: ndet
      integer :: iunit

      real(8) :: z1, z2, zmid

      character(len=1) :: clet

!--- start

      ndet=0
      do ja=1, jcore
        do ia=1, icore
           na=mapcore(ia,ja)
           if (na.eq.0) cycle
           if (detector_map(ia,ja).eq.0) cycle
           ndet=ndet+1
           iunit=100+ndet
           write (iunit,'(/,a,i0,a,f8.2,a)') 'Pin Powers, Case ',nstate,': ', efpd, ' EFPD'
           write (iunit,'(a)') ',,,,Pin'

!  write row of labels

           write (iunit,'(a)',advance='no') 'Axial Elevation,,,,'
           do j=1, npin
             clet=char(ichar('A')+j-1)
             do i=1, npin
                write (iunit,'(",",a1,i0)',advance='no') clet,i
             enddo
           enddo
           write (iunit,*)
           write (iunit,'(a)') 'Bottom,Midpoint,Top,Height'


!***** check are letters the row or column?

!  write row of pin powers 

           z1=axial0
           do k=1, kd
             z2=z1+axial(k) 
             zmid=0.5d0*(z1+z2)
             write (iunit,'(f8.3,3(",",f8.3))',advance='no') z1, zmid, z2, axial(k) ! bottom, mid, top, height
             do j=1, npin
               do i=1, npin
                 write (iunit,'(",",f8.4)',advance='no') power(i,j,k,na)
               enddo
             enddo
             write (iunit,*)
             z1=z2
           enddo    ! k
           write (iunit,*)    ! trailing blank line

        enddo
      enddo


      return
      end subroutine print_nea
!=======================================================================
   end module mod_nea
