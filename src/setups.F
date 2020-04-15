      module root
      implicit none
      integer,dimension(3):: ntotal
      integer,dimension(3):: npart
      integer:: nfld

      character(2),parameter::id="Sc"
      integer::nhy
      real(8)::time
      real(8),dimension(:,:,:,:),allocatable:: vpart

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

      subroutine SetParameters
      use root
      implicit none

      ntotal(1)=512
      ntotal(2)=64
      ntotal(3)=128

      return
      end subroutine SetParameters

      subroutine SetMPIParameters
      use mpipara
      implicit none
      integer::key,color
      integer::np_hyd

      ntiles(1)=16
      ntiles(2)=8
      ntiles(3)=16
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.false.

c Initialize MPI
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )

c Making 3D strucure
      np_hyd = ntiles(1)*ntiles(2)*ntiles(3)
      color = int(myid_w/np_hyd)
      key   = myid_w   
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,mpi_comm_hyd,ierr)
      call MPI_COMM_SIZE( mpi_comm_hyd, nprocs_hyd, ierr )
      call MPI_COMM_RANK( mpi_comm_hyd, myid_hyd , ierr )     

c Create a virtual Cartesian topology for the domain decomposition.
c
      call MPI_CART_CREATE( mpi_comm_hyd, 3, ntiles, periodic
     &                    , reorder, comm3d, ierr )
      call MPI_COMM_RANK( comm3d, myid,     ierr )
      call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
c
c Find the ranks of my neighbors; find my virtual Cartesian coords.
c
      call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
      call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
      call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
c
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )
      
      return
      end subroutine SetMPIParameters

      subroutine AllocateVariable
      use root
      use mpipara
      implicit none

      nfld=200
      npart(:) = ntotal(:)/ntiles(:)

      if(myid_w .eq. 0 )write(6,*) "grid for 1 procss",npart(:)

      allocate(vpart(nfld,npart(1),npart(2),npart(3)))

      vpart(:,:,:,:) = 0.0d0

      return
      end subroutine AllocateVariable
