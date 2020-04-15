      program iotest
      use root
      use binmod
      implicit none
      
      call SetParameters
      call SetMPIParameters
      call AllocateVariable

      do nhy=1,10
         call bindmp
      enddo

      stop
      end program iotest
