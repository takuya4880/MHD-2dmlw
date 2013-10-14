module defstruct
    implicit none
    integer,parameter :: nx=300
    integer,parameter :: nz=200
    integer,parameter :: marg=4
    integer,parameter :: ix=nx+2*marg
    integer,parameter :: iz=nz+2*marg
    

    type constants 
        integer nx, nz, ix, iz, marg
        double precision dx, dz, dt, wid, hig
        double precision gam, q, a
        double precision gx, gy, gz 
    end type
    
    type cell
        type(constants) con
        double precision x(ix), z(iz)
        double precision ro(ix,iz), rovx(ix,iz), rovy(ix,iz), rovz(ix,iz)
        double precision bx(ix,iz), by(ix,iz), bz(ix,iz), e(ix,iz) 
        double precision pr(ix,iz)
        double precision bpot(ix,iz)       
    end type    


    
end module

