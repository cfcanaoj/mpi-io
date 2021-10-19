      module binmod
      implicit none
      private
      integer::SAG1D,SAG2D,SAG3D,SAD3D
      public bindmp
      contains

      subroutine bindmp
      use root, nfldhyd=>nfld
      use mpipara
      implicit NONE

      integer :: i, j, k, l, m, n
      integer :: nx1, nx2, nx3

      integer strtoi,incr
      character(10) :: datadir
      data datadir /"bindata/"/ 
      character(15) :: unffile
      character(40)::usrfile
      character(30) :: fpathbin,fpathunf
      integer, parameter :: unitunf=560
      integer,save:: unitd3d,unitg1d, unitg2d, unitg3d, unitg0d
      data unitd3d / 512 /
      data unitg1d / 513 /
      data unitg2d / 514 /
      data unitg3d / 515 /
      data unitg0d / 516 /

      logical :: fileflag

      logical,save :: is_inited
      data is_inited / .false. /

      integer::iss,jss,kss
      integer::iee,jee,kee
      integer::itot,jtot,ktot
      integer:: nfld
      integer::lastnum
      real(8),dimension(:,:),allocatable,save :: grid1D,grid2D,grid3D &
     &                                          ,grid0D
      real(8),dimension(:,:,:,:),allocatable,save :: data3D
      integer,dimension(4)::Asize,Ssize,Start
      integer(kind=MPI_OFFSET_KIND) idisp
      data idisp / 0 /

!==============

      iss = 1; jss = 1; kss = 1
      iee = npart(1); jee = npart(2); kee = npart(3)
      nfld=nfldhyd
      initdata: if(.not. is_inited )then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Asize(1) = nfld
         Ssize(1) = nfld
         Start(1) = 0 
         Asize(2) = ntotal(1) ! total zones for 1D 
         Asize(3) = ntotal(2) ! total zones for 2D
         Asize(4) = ntotal(3) ! total zones for 3D
         Ssize(2) = (iee-iss+1) ! partial zones in 1 process 
         Ssize(3) = (jee-jss+1) ! partial zones in 1 process 
         Ssize(4) = (kee-kss+1) ! partial zones in 1 process 
         Start(2) = (iee-iss+1) * coords(1)
         Start(3) = (jee-jss+1) * coords(2)
         Start(4) = (kee-kss+1) * coords(3)

         call MPI_TYPE_CREATE_SUBARRAY(&
     & 4, &! dimension of array
     & Asize,Ssize,Start,&
     & MPI_ORDER_FORTRAN,&
     & MPI_DOUBLE_PRECISION,&
     & SAD3D,& ! Data type of Subarray for data 3D
     & ierr)
         call MPI_TYPE_COMMIT(SAD3D,ierr)

         allocate(data3D(nfld,iss:iee,jss:jee,kss:kee))
      endif initdata
      
      data3D(1:nfld,iss:iee,jss:jee,kss:kee) = vpart(1:nfld,iss:iee,jss:jee,kss:kee)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      incr=nhy
      write(usrfile,"(a3,a2,a1,i5.5)")'d3d',id,'.',incr
      fpathbin = trim(datadir)//usrfile

      call MPI_FILE_OPEN(MPI_COMM_WORLD, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitd3d,ierr)
      call MPI_FILE_SET_VIEW(&
     &  unitd3d,  &! file path
     &     idisp, & ! 
     & MPI_DOUBLE_PRECISION,& 
     &     SAD3D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL(&
     &   unitd3d,  &! file path
     &    data3D,  &! the data
     & (iee-iss+1)*(jee-jss+1)*(kee-kss+1)*nfld,& ! total data number
     & MPI_DOUBLE_PRECISION,&  
     &      stat,&
     &      ierr)
      call MPI_FILE_CLOSE(unitd3d,ierr)

      if (myid_w .eq. 0) then
         write(6,"('Binary data dump written at time=',1pe12.5,&
     &   ' cycle=',i0)") time,nhy
      endif

      is_inited = .true.

      return
      end subroutine bindmp
      end module binmod
