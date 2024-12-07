      module root
      implicit none
      integer,dimension(3):: ntotal
      integer,dimension(3):: npart
      integer:: nfld

      character(2),parameter::id="Sc"
      integer::nhy
      real(8)::time
      real(8),dimension(:,:,:,:),allocatable:: vpart
      real(8),dimension(:),allocatable:: xc,xe
      real(8),dimension(:),allocatable:: yc,ye
      real(8),dimension(:),allocatable:: zc,ze

      end module  root

      module mpipara
      implicit none
      include "mpif.h"
      integer, parameter :: mreq  = 300
      integer :: stat(MPI_STATUS_SIZE,mreq)                     
      integer :: req(mreq)

      integer :: ierr,myid_w, nprocs_w
      integer :: mpi_comm_hyd,myid_hyd, nprocs_hyd
      integer :: comm3d,myid, nprocs
      logical :: periodic(3)
      integer :: ntiles(3), coords(3)
      logical :: reorder
      integer :: n1m, n1p, n2m, n2p, n3m, n3p
      integer :: nreq, nsub

      end module  mpipara

      subroutine InitializeMPI
      use mpipara
      implicit none
! Initialize MPI
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )

      return
      end subroutine InitializeMPI

      subroutine ReadParameters
      use root
      use mpipara
      implicit none
      character(20)::filename
      data filename /"inputparam"/
      logical::flag

      namelist /totalgrid/ ntotal,nfld
      namelist /paragrid/ ntiles,periodic

      nfld=200

      ntotal(1)=512
      ntotal(2)=64
      ntotal(3)=128

      ntiles(1)=16
      ntiles(2)=8
      ntiles(3)=16
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.false.

       INQUIRE(FILE = filename,EXIST = flag)
       if(flag) then
          open(unit=1,file=filename,status='old')
       else   ! flag       
          write(6,*) filename,' not found'
          call MPI_FINALIZE ( ierr )
          stop 'Cannot Open input parameter file'
       end if ! flag 

      if(myid_w .eq. 0) then
       read(1,totalgrid)
       read(1,paragrid)
      endif ! myid_w

      call MPI_BCAST(ntotal,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ntiles,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(periodic,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine ReadParameters

      subroutine SetMPIParameters
      use mpipara
      implicit none
      integer::key,color
      integer::np_hyd

! Making 3D strucure
      np_hyd = ntiles(1)*ntiles(2)*ntiles(3)
      color = int(myid_w/np_hyd)
      key   = myid_w   
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,mpi_comm_hyd,ierr)
      call MPI_COMM_SIZE( mpi_comm_hyd, nprocs_hyd, ierr )
      call MPI_COMM_RANK( mpi_comm_hyd, myid_hyd , ierr )     

! Create a virtual Cartesian topology for the domain decomposition.
!
      call MPI_CART_CREATE( mpi_comm_hyd, 3, ntiles, periodic &
     &                    , reorder, comm3d, ierr )
      call MPI_COMM_RANK( comm3d, myid,     ierr )
      call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
!
! Find the ranks of my neighbors; find my virtual Cartesian coords.
!
      call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
      call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
      call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
!
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )
      
      return
      end subroutine SetMPIParameters

      subroutine AllocateVariable
      use root
      use mpipara
      implicit none
      integer:: i,j,k

      npart(:) = ntotal(:)/ntiles(:)

      if(myid_w .eq. 1 )write(6,*) "toral grid       ",ntotal(:)
      if(myid_w .eq. 1 )write(6,*) "parallel tiles   ",ntiles(:)
      if(myid_w .eq. 1 )write(6,*) "grid for 1 procss",npart(:)

!========================
      allocate(vpart(nfld,npart(1),npart(2),npart(3)))

      vpart(:,:,:,:) = 0.0d0
!========================
! here xc, yc, zc do not need the last mesh but allocate it for the later convenience
      allocate(xc(npart(1)+1),xe(npart(1)+1))
      allocate(yc(npart(2)+1),ye(npart(2)+1))
      allocate(zc(npart(3)+1),ze(npart(3)+1))

!========================
      do i=1,npart(1)+1
         xe(i) = i-1 + coords(1)*npart(1)
      enddo
      do i=1,npart(1)
         xc(i) = 0.5d0*(xe(i)+xe(i+1))
      enddo
!========================      
      do j=1,npart(2)+1
         ye(j) = j-1 + coords(2)*npart(2)
      enddo
      do j=1,npart(2)
         yc(j) = 0.5d0*(ye(j)+ye(j+1))
      enddo
!========================
      do k=1,npart(3)+1
         ze(k) = k-1 + coords(3)*npart(3)
      enddo
      do k=1,npart(3)
         zc(k) = 0.5d0*(ze(k)+ze(k+1))
      enddo
!========================
 
      return
      end subroutine AllocateVariable
