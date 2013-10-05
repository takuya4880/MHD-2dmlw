subroutine flux(box, fx, fz)
    use defstruct
    implicit none
    type(cell) :: box, fx, fz
    
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision, allocatable :: b2(:,:),roi(:,:), roh(:,:)
    double precision, allocatable :: eta(:,:), ex(:,:), ey(:,:), ez(:,:)
    double precision, allocatable :: jx(:,:), jy(:,:), jz(:,:)
    allocate(b2(ix,iz))
    allocate(roi(ix,iz))
    allocate(roh(ix,iz))
    allocate(eta(ix,iz))
    allocate(ex(ix,iz), ey(ix,iz), ez(ix,iz))
    allocate(jx(ix,iz), jy(ix,iz), jz(ix,iz))

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

    fz%ro = box%rovz
    fz%rovx = box%rovz*box%rovx*roi - box%bz*box%bx
    fz%rovy = box%rovz*box%rovy*roi - box%bz*box%by
    fz%rovz = box%rovz*box%rovz*roi - box%bz*box%bz + box%pr + 0.5*b2
    fz%bx = (box%rovz*box%bx - box%rovx*box%bz) * roi
    fz%by = (box%rovz*box%by - box%rovy*box%bz) * roi 
    fz%bz = 0.
    fz%e = (roh*box%rovz - box%bz*(box%rovx*box%bx + box%rovy*box%by &
            + box%rovz*box%bz) )*roi 
    fz%bpot = 0
    
    deallocate(b2,roi,roh,eta,ex,ey,ez,jx,jy,jz)

end subroutine
