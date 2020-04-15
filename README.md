# A sample code for Parallel IO

## Proceadure
    
    mkdir bindata 
    cd src
    make
    cd ..
    qsub pbs_xc.csh
   
## Change Parameters
    
    vim inputparam

### Meaning
Grid number for 3D data
    
    ntotal(1)=512,ntotal(2)=64,ntotal(3)=128
    
Number of hydrodynamic variables
    
    nfld=200
    
Number of stenciles for the parallelization
    
    ntiles(1)=2,ntiles(2)=2,ntiles(3)=2
    
The number of the MPI process should match ntiles(1)*ntiles(2)*ntiles(3)


