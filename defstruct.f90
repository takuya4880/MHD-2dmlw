module defstruct
    implicit none
    integer,parameter :: nx=256
    integer,parameter :: ny=256
    integer,parameter :: marg=2
    integer,parameter :: ix=nx+2*marg
    integer,parameter :: iy=ny+2*marg

    type constants 
        integer nx, ny, ix, iy, marg
        double precision dx, dy, dt, wid, hig
        double precision gam, q, a
        double precision gx, gy, gz 
    end type
    
    type cell
        type(constants) con
        double precision x(ix), y(iy)
        double precision ro(ix,iy), rovx(ix,iy), rovy(ix,iy), rovz(ix,iy)
        double precision bx(ix,iy), by(ix,iy), bz(ix,iy), e(ix,iy) 
        double precision pr(ix,iy)       
    end type    

end module

