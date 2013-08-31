subroutine flux(box, fx, fy)
    use defstruct
    implicit none
    type(cell) :: box, fx, fy
    
    double precision, allocatable :: b2(:,:),vx(:,:),vy(:,:),vz(:,:),roh(:,:)
    allocate(b2(ix,iy))
    allocate(vx(ix,iy),vy(ix,iy),vz(ix,iy))
    allocate(roh(ix,iy))

    b2 = box%bx**2 + box%by**2 + box%bz**2
    vx = box%rovx/box%ro
    vy = box%rovy/box%ro
    vz = box%rovz/box%ro
    roh = 0.5*(box%rovx**2+box%rovy**2+box%rovz**2)*roi &
            + box%pr*box%con%gam/(box%con%gam-1.) + b2

    fx%ro = box%rovx
    fx%rovx = vx*box%rovx - box%bx*box%bx + box%pr + 0.5*b2
    fx%rovx = vx*box%rovy - box%bx*box%by
    fx%rovz = vx*box%rovz - box%bx*box%bz
    fx%bx = 0.
    fx%by = vx*box%by - vy*box%bx
    fx%bz = vx*box%bz - vz*box%bx
    fx%e = roh*vx - box%bx*(vx*box%bx + vy*box%by + vz*box%bz)

    fy%ro = box%rovy
    fy%rovx = vy*box%rovx - box%by*box%bx
    fy%rovx = vy*box%rovy - box%by*box%by + box%pr + 0.5*b2
    fy%rovz = vy*box%rovz - box%by*box%bz
    fy%bx = vy*box%bx - vx*box%by
    fy%by = 0.
    fy%bz = vy*box%bz - vz*box%by
    fy%e = roh*vy - box%by*(vx*box%bx + vy*box%by + vz*box%bz) 
    
    deallocate(b2,vx,vy,vz,roh)

end subroutine
