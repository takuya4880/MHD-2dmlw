subroutine flux(box, fx, fy)
    use defstruct
    implicit none
    type(cell) :: box, fx, fy
    
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision, allocatable :: b2(:,:),roi(:,:), roh(:,:)
    double precision, allocatable :: eta(:,:), ex(:,:), ey(:,:), ez(:,:)
    double precision, allocatable :: jx(:,:), jy(:,:), jz(:,:)
    allocate(b2(ix,iy))
    allocate(roi(ix,iy))
    allocate(roh(ix,iy))
    allocate(eta(ix,iy))
    allocate(ex(ix,iy), ey(ix,iy), ez(ix,iy))
    allocate(jx(ix,iy), jy(ix,iy), jz(ix,iy))

    b2 = box%bx**2 + box%by**2 + box%bz**2
    roi = 1./box%ro
    roh = 0.5*(box%rovx**2+box%rovy**2+box%rovz**2)*roi &
            + box%pr*box%con%gam/(box%con%gam-1.) + b2

    fx%ro = box%rovx
    fx%rovx = box%rovx*box%rovx*roi - box%bx*box%bx + box%pr + 0.5*b2
    fx%rovy = box%rovx*box%rovy*roi - box%bx*box%by
    fx%rovz = box%rovx*box%rovz*roi - box%bx*box%bz
    fx%bx = 0.
    fx%by = (box%rovx*box%by - box%rovy*box%bx) * roi
    fx%bz = (box%rovx*box%bz - box%rovz*box%bx) * roi
    fx%e = (roh*box%rovx - box%bx*(box%rovx*box%bx + box%rovy*box%by &
            + box%rovz*box%bz) )*roi 
    fx%bpot = 0

    fy%ro = box%rovy
    fy%rovx = box%rovy*box%rovx*roi - box%by*box%bx
    fy%rovy = box%rovy*box%rovy*roi - box%by*box%by + box%pr + 0.5*b2
    fy%rovz = box%rovy*box%rovz*roi - box%by*box%bz
    fy%bx = (box%rovy*box%bx - box%rovx*box%by) * roi
    fy%by = 0.
    fy%bz = (box%rovy*box%bz - box%rovz*box%by) * roi 
    fy%e = (roh*box%rovy - box%by*(box%rovx*box%bx + box%rovy*box%by &
            + box%rovz*box%bz) )*roi 
    fy%bpot = 0
    
    deallocate(b2,roi,roh,eta,ex,ey,ez,jx,jy,jz)

end subroutine
