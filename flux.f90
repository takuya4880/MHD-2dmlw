module flx
implicit none
contains
subroutine flux(box, fx, fz)
    use defstruct
    implicit none
    type(cell) :: box, fx, fz
    
    integer :: i,j
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision :: b2,roi,roh,vx,vy,vz,bx,by,bz,pr
    double precision :: eta, ex, ey, ez
    double precision :: jx, jy, jz

    !$omp parallel do private(i,roi,vx,vy,vz,bx,by,bz,pr,b2,roh,jx,jy,jz,eta,ex,ey,ez)
    do j=1,iz
        do i=1,ix
            roi = 1./box%ro(i,j)
            vx = box%rovx(i,j) * roi
            vy = box%rovy(i,j) * roi
            vz = box%rovz(i,j) * roi
            bx = box%bx(i,j)
            by = box%by(i,j)
            bz = box%bz(i,j)
            pr = box%pr(i,j)
            b2 = bx**2 + by**2 + bz**2
            roh = 0.5*(vx**2+vy**2+vz**2) &
                        + pr*box%con%gam/(box%con%gam-1.) + b2

            jx = -(box%by(i,j+1)-box%by(i,j-1))/(2.*box%con%dz)
            jy = (box%bx(i,j+1)-box%bx(i,j-1))/(2.*box%con%dz) &
                        -(box%bz(i+1,j)-box%bz(i-1,j))/(2.*box%con%dx)
            jz = (box%by(i+1,j)-box%by(i-1,j))/(2.*box%con%dx)

            eta = sqrt((jx**2+jy**2+jz**2)*16.*atan(1.0))*roi  !calculate vd
            
            if (eta<vc) then
                eta = 0.
            else 
                eta = alp*(eta/vc-1.)**2
                if (eta>etamax) then
                    eta = etamax
                end if
            end if
            ex = eta*jx + (-vy*bz+vz*by)
            ey = eta*jy + (-vz*bx+vx*bz)
            ez = eta*jz + (-vx*by+vy*bx)

            fx%ro(i,j) = box%rovx(i,j)
            fx%rovx(i,j) = vx*box%rovx(i,j) - bx*bx + pr + 0.5*b2
            fx%rovy(i,j) = vx*box%rovy(i,j) - bx*by
            fx%rovz(i,j) = vx*box%rovz(i,j) - bx*bz
            fx%bx(i,j) = 0.
            fx%by(i,j) = -ez
            fx%bz(i,j) = ey
            fx%e(i,j) = (roh*vx - bx*(vx*bx + vy*by + vz*bz) ) + (ey*bz - ez*by) 
            fx%bpot(i,j) = 0

            fz%ro(i,j) = box%rovz(i,j)
            fz%rovx(i,j) = vz*box%rovx(i,j) - bz*bx
            fz%rovy(i,j) = vz*box%rovy(i,j) - bz*by
            fz%rovz(i,j) = vz*box%rovz(i,j) - bz*bz + pr + 0.5*b2
            fz%bx(i,j) = -ey
            fz%by(i,j) = ex
            fz%bz(i,j) = 0.
            fz%e(i,j) = (roh*vz - bz*(vx*bx + vy*by + vz*bz) ) + (ex*by - ey*bx)
            fz%bpot(i,j) = 0
        end do
    end do
    !$omp end parallel do 
    

end subroutine
end module 
