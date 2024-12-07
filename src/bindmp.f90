      module binmod
      implicit none
      private
      integer::SAG1D,SAG2D,SAG3D,SAD3D
      integer :: commG1D,commG2D,commG3D
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
      integer,save:: unitd3d,unitg1d, unitg2d, unitg3d
      data unitd3d / 512 /
      data unitg1d / 513 /
      data unitg2d / 514 /
      data unitg3d / 515 /

      logical :: fileflag

      logical,save :: is_inited
      data is_inited / .false. /

      integer::iss,jss,kss
      integer::iee,jee,kee
      integer:: nfld
      integer,parameter::ngrid=2
      integer,parameter::xdir=1,ydir=2,zdir=3
      integer::lastnum
      real(8),dimension(:,:),allocatable,save :: grid1D,grid2D,grid3D
      real(8),dimension(:,:,:,:),allocatable,save :: data3D
      integer,dimension(4)::Asize,Ssize,Start
      integer(kind=MPI_OFFSET_KIND) idisp
      data idisp / 0 /
      
      integer::color,key

!====================================================
      init1D: if(.not. is_inited)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         iss = 1
         if(coords(xdir) < ntiles(xdir)-1)then
            iee = npart(xdir)
         else
            iee = npart(xdir)+1
         endif
         
         allocate(grid1D(ngrid,iss:iee))
         grid1D(1,iss:iee) =   xe(iss:iee)
         grid1D(2,iss:iee) =   xc(iss:iee)

! Grid
         Asize(1) = ngrid
         Ssize(1) = ngrid
         Start(1) = 0
         Asize(2) = ntotal(xdir)+1             ! total izones
         Ssize(2) = (iee-iss+1)                ! izones in 1 process
         Start(2) = npart(xdir) * coords(xdir) ! Start point
         call MPI_TYPE_CREATE_SUBARRAY( &
     & 2,  &! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION, &
     & SAG1D,  &! Data type of Subarray for Grid 1D
     & ierr)
         call MPI_TYPE_COMMIT(SAG1D,ierr)

      color =  coords(2)*ntiles(3)+coords(3)
      key   =  coords(1)
      call MPI_COMM_SPLIT(comm3d,color,key,commG1D,ierr)
      if(color == 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D GRID WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(unffile,"(a3,a2)")'g1d',id
      fpathbin = trim(datadir)//unffile
      
      call MPI_FILE_OPEN(commG1D, &
     &                         fpathbin, &! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg1d,ierr)

      call MPI_FILE_SET_VIEW( &
     &   unitg1d,   &! file path
     &     idisp,   &! 
     & MPI_DOUBLE_PRECISION,   &
     &     SAG1D,   &! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg1d,         &! file path
     &    grid1D,         &! the data
     & (iee-iss+1)*ngrid, &! total data number
     & MPI_DOUBLE_PRECISION, &
     & stat, ierr)
      call MPI_FILE_CLOSE(unitg1d,ierr)
      endif
      endif init1D

      init2D: if(.not. is_inited)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         jss = 1
         if(coords(ydir) < ntiles(ydir)-1)then
            jee = npart(ydir)
         else
            jee = npart(ydir)+1
         endif
         
         allocate(grid2D(ngrid,jss:jee))
         grid2D(1,jss:jee) =   ye(jss:jee)
         grid2D(2,jss:jee) =   yc(jss:jee)

! Grid
         Asize(1) = ngrid
         Ssize(1) = ngrid
         Start(1) = 0
         Asize(2) = ntotal(ydir)+1     ! total jzones
         Ssize(2) = (jee-jss+1)        ! jzones in 1 process
         Start(2) = npart(ydir) * coords(ydir) ! Start point
         call MPI_TYPE_CREATE_SUBARRAY( &
     & 2,  &! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION,  &! This is changed to double by cpp
     & SAG2D,  &! Data type of Subarray for Grid 2D
     & ierr)
         call MPI_TYPE_COMMIT(SAG2D,ierr)

      color =  coords(3)*ntiles(1)+coords(1)
      key   =  coords(2)
      call MPI_COMM_SPLIT(comm3d,color,key,commG2D,ierr)
      if(color == 0) then
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D GRID WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(unffile,"(a3,a2)")'g2d',id
      fpathbin = trim(datadir)//unffile
      
      call MPI_FILE_OPEN(commG2D, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg2d,ierr)
      call MPI_FILE_SET_VIEW( &
     &   unitg2d,   &! file ID
     &     idisp,   &! 
     & MPI_DOUBLE_PRECISION,   &! that is replaced to double by cpp
     &     SAG2D,   &! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg2d,         &! file ID
     &    grid2D,         &! the data
     & (jee-jss+1)*ngrid, &! total data number
     & MPI_DOUBLE_PRECISION,         &! that is replaced to double by cpp
     & stat, ierr)
      
      call MPI_FILE_CLOSE(unitg2d,ierr)
      endif
      endif init2D

      init3D: if(.not. is_inited)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         kss = 1
         if(coords(zdir) < ntiles(zdir)-1)then
            kee = npart(zdir)
         else
            kee = npart(zdir)+1
         endif
         
         allocate(grid3D(ngrid,kss:kee))
         grid3D(1,kss:kee) =   ze(kss:kee)
         grid3D(2,kss:kee) =   zc(kss:kee)

! Grid
         Asize(1) = ngrid
         Ssize(1) = ngrid
         Start(1) = 0
         Asize(2) = ntotal(zdir)+1     ! total kzones
         Ssize(2) = (kee-kss+1)        ! kzones in 1 process
         Start(2) = npart(zdir) * coords(3) ! Start point
         call MPI_TYPE_CREATE_SUBARRAY( &
     & 2,  &! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION, &
     & SAG3D,  &! Data type of Subarray for Grid 3D
     & ierr)
         call MPI_TYPE_COMMIT(SAG3D,ierr)

      color =  coords(1)*ntiles(2)+coords(2)
      key   =  coords(3)
      call MPI_COMM_SPLIT(comm3d,color,key,commG3D,ierr)
      if(color == 0) then
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D GRID WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(unffile,"(a3,a2)")'g3d',id
      fpathbin = trim(datadir)//unffile

      call MPI_FILE_OPEN(commG3D, &
     &                         fpathbin, &! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg3d,ierr)
      call MPI_FILE_SET_VIEW( &
     &  unitg3d,   &! file path
     &     idisp,  &! 
     & MPI_DOUBLE_PRECISION,  &
     &     SAG3D,  &! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &  unitg3d,          &! file path
     &    grid3D,         &! the data
     & (kee-kss+1)*ngrid, &! total data number
     & MPI_DOUBLE_PRECISION, &
     & stat, ierr)
      call MPI_FILE_CLOSE(unitg3d,ierr)
      endif
      endif init3D

!====================================================
      iss = 1; jss = 1; kss = 1
      iee = npart(xdir); jee = npart(ydir); kee = npart(zdir)
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
         write(6,"('Binary data dump written at time=',1pe12.5,' cycle=',i0)") time,nhy
      endif

      is_inited = .true.

      return
      end subroutine bindmp
      end module binmod
