subroutine flux(box, fx, fz)
    use defstruct
    implicit none
    type(cell) :: box, fx, fz
    
    integer :: i,j
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

    jx(2:ix-1,2:iz-1) = -(box%by(2:ix-1,3:iz)-box%by(2:ix-1,1:iz-2))/(2.*box%con%dz)
    jy(2:ix-1,2:iz-1) = (box%bx(2:ix-1,3:iz)-box%bx(2:ix-1,1:iz-2))/(2.*box%con%dz) &
                        -(box%bz(3:ix,2:iz-1)-box%bz(1:ix-2,2:iz-1))/(2.*box%con%dx)
    jx(2:ix-1,2:iz-1) = (box%by(3:ix,2:iz-1)-box%by(1:ix-2,2:iz-1))/(2.*box%con%dx)

    eta = sqrt((jx**2+jy**2+jz**2)*16.*atan(1.0))*roi  !calculate vd
    do j=1,iz
        do i=1,ix
            if (eta(i,j)<vc) then
                eta(i,j) = 0.
            else 
                eta(i,j) = alp*(eta(i,j)-1.)**2
                if (eta(i,j)>etamax) then
                    eta(i,j) = etamax
                end if
            end if
        end do
    end do

    ex = eta*jx + (-box%rovy*box%bz+box%rovz*box%by)*roi
    ey = eta*jy + (-box%rovz*box%bx+box%rovx*box%bz)*roi
    ez = eta*jz + (-box%rovx*box%by+box%rovy*box%bx)*roi

    fx%ro = box%rovx
    fx%rovx = box%rovx*box%rovx*roi - box%bx*box%bx + box%pr + 0.5*b2
    fx%rovy = box%rovx*box%rovy*roi - box%bx*box%by
    fx%rovz = box%rovx*box%rovz*roi - box%bx*box%bz
    fx%bx = 0.
    fx%by = -ez
    fx%bz = ey
    fx%e = (roh*box%rovx - box%bx*(box%rovx*box%bx + box%rovy*box%by &
            + box%rovz*box%bz) )*roi + (ey*box%bz - ez*box%by)
    fx%bpot = 0

    fz%ro = box%rovz
    fz%rovx = box%rovz*box%rovx*roi - box%bz*box%bx
    fz%rovy = box%rovz*box%rovy*roi - box%bz*box%by
    fz%rovz = box%rovz*box%rovz*roi - box%bz*box%bz + box%pr + 0.5*b2
    fz%bx = -ey
    fz%by = ex
    fz%bz = 0.
    fz%e = (roh*box%rovz - box%bz*(box%rovx*box%bx + box%rovy*box%by &
            + box%rovz*box%bz) )*roi + (ex*box%by - ey*box%bx)
    fz%bpot = 0
    
    deallocate(b2,roi,roh,eta,ex,ey,ez,jx,jy,jz)

end subroutine
