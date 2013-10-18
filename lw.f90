module lw
implicit none 
contains
subroutine lw1(box, h, fx, fz, s)
    use defstruct
    implicit none
    type(cell) :: box, h, fx, fz, s
    
    double precision :: ddx, ddz, dt
    ddx = box%con%dt / box%con%dx
    ddz = box%con%dt / box%con%dz
    dt  = box%con%dt

    call each1(box%ro,h%ro,fx%ro,fz%ro,s%ro,ddx,ddz,dt)
    call each1(box%rovx,h%rovx,fx%rovx,fz%rovx,s%rovx,ddx,ddz,dt)
    call each1(box%rovy,h%rovy,fx%rovy,fz%rovy,s%rovy,ddx,ddz,dt)
    call each1(box%rovz,h%rovz,fx%rovz,fz%rovz,s%rovz,ddx,ddz,dt)
    call each1(box%bx,h%bx,fx%bx,fz%bx,s%bx,ddx,ddz,dt)
    call each1(box%by,h%by,fx%by,fz%by,s%by,ddx,ddz,dt)
    call each1(box%bz,h%bz,fx%bz,fz%bz,s%bz,ddx,ddz,dt)
    call each1(box%e,h%e,fx%e,fz%e,s%e,ddx,ddz,dt)
    call each1(box%bpot,h%bpot,fx%bpot,fz%bpot,s%bpot,ddx,ddz,dt)
    
    
end subroutine

subroutine each1(u,h,fx,fz,s,ddx,ddz,dt)
    use defstruct
    !$use omp_lib
    implicit none
    double precision :: u(ix,iz),h(ix,iz),fx(ix,iz),fz(ix,iz),s(ix,iz)
    double precision :: ddx, ddz, dt
    
    integer i,j
    double precision fffx, fffz, ss

    !$omp parallel do private(j,fffx,fffz,ss) 
    do i=1,iz-1
        do j=1,ix-1
            fffx = 0.5 * ( fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i) )
            fffz = 0.5 * ( fz(j+1,i+1)-fz(j+1,i)+fz(j,i+1)-fz(j,i) )
            ss   = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            h(j,i) = 0.25 * ( u(j+1,i+1)+u(j+1,i)+u(j,i+1)+u(j,i) ) &
                        - 0.5*(ddx*fffx + ddz*fffz - dt*ss)
        end do
    end do
    !$omp end parallel do
end subroutine

subroutine lw2(box, h, f, fx, fz, s)
    use defstruct
    implicit none
    type(cell) :: box, h, f, fx, fz, s
    
    double precision :: ddx, ddz, dt
    ddx = box%con%dt / box%con%dx
    ddz = box%con%dt / box%con%dz
    dt  = box%con%dt

    call each2(box%ro,h%ro,f%ro,fx%ro,fz%ro,s%ro,ddx,ddz,dt)
    call each2(box%rovx,h%rovx,f%rovx,fx%rovx,fz%rovx,s%rovx,ddx,ddz,dt)
    call each2(box%rovy,h%rovy,f%rovy,fx%rovy,fz%rovy,s%rovy,ddx,ddz,dt)
    call each2(box%rovz,h%rovz,f%rovz,fx%rovz,fz%rovz,s%rovz,ddx,ddz,dt)
    call each2(box%bx,h%bx,f%bx,fx%bx,fz%bx,s%bx,ddx,ddz,dt)
    call each2(box%by,h%by,f%by,fx%by,fz%by,s%by,ddx,ddz,dt)
    call each2(box%bz,h%bz,f%bz,fx%bz,fz%bz,s%bz,ddx,ddz,dt)
    call each2(box%e,h%e,f%e,fx%e,fz%e,s%e,ddx,ddz,dt)
    call each2(box%bpot,h%bpot,f%bpot,fx%bpot,fz%bpot,s%bpot,ddx,ddz,dt)
   

end subroutine

subroutine each2(box,h,f,fx,fz,s,ddx,ddz,dt)
    use defstruct
    !$use omp_lib
    implicit none
    double precision :: box(ix,iz),f(ix,iz),h(ix,iz)
    double precision :: fx(ix,iz),fz(ix,iz),s(ix,iz)
    double precision :: ddx, ddz, dt
    
    integer i,j
    double precision fffx,fffz,ss

    !$omp parallel do private(j,fffx,fffz,ss) 
    do i=1,iz-2
        do j=1,ix-2
            fffx = 0.5 * (fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i))
            fffz = 0.5 * (fz(j+1,i+1)-fz(j+1,i)+fz(j,i+1)-fz(j,i))
            ss  = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            f(j+1,i+1) = box(j+1,i+1) - ddx*fffx - ddz*fffz + dt*ss
        end do
    end do
    !$omp end parallel do

end subroutine 

subroutine artvis(box, f)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box, f
    double precision, allocatable :: kapx(:,:),kapz(:,:)
    allocate(kapx(ix,iz),kapz(ix,iz))
    
    !$omp parallel workshare
    kapx(2:ix,:) = box%con%q * abs(box%rovx(2:ix,:)/box%ro(2:ix,:) &
                                    - box%rovx(1:ix-1,:)/box%ro(1:ix-1,:))
    kapz(:,2:iz) = box%con%q * abs(box%rovz(:,2:iz)/box%ro(:,2:iz) &
                                    - box%rovz(:,1:iz-1)/box%ro(:,1:iz-1)) 
    !$omp end parallel workshare

    call eachav(box%ro, f%ro, kapx, kapz, box%con)
    call eachav(box%rovx, f%rovx, kapx, kapz, box%con)
    call eachav(box%rovy, f%rovy, kapx, kapz, box%con)
    call eachav(box%rovz, f%rovz, kapx, kapz, box%con)
    call eachav(box%bx, f%bx, kapx, kapz, box%con)
    call eachav(box%by, f%by, kapx, kapz, box%con)
    call eachav(box%bz, f%bz, kapx, kapz, box%con)
    call eachav(box%e, f%e, kapx, kapz, box%con)
    call eachav(box%bpot, f%bpot, kapx, kapz, box%con)

    deallocate(kapx,kapz)

end subroutine

subroutine eachav(box,f,kapx,kapz,con)
    use defstruct
    !$use omp_lib
    double precision :: box(ix,iz),f(ix,iz)
    double precision :: kapx(ix,iz),kapz(ix,iz)
    type(constants) con
    
    double precision :: ddx, ddz
    double precision, allocatable :: difx(:,:),difz(:,:)
    allocate(difx(ix,iz),difz(ix,iz))
    ddx = con%dt / con%dx
    ddz = con%dt / con%dz
    
    !$omp parallel
    !$omp workshare
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    difz(:,2:iz) = f(:,2:iz) - f(:,1:iz-1)
    !$omp end workshare
    !$omp workshare
    box(3:ix-2,3:iz-2) = f(3:ix-2,3:iz-2) &
             + ddx * ( kapx(4:ix-1,3:iz-2)*difx(4:ix-1,3:iz-2)     &
                       - kapx(3:ix-2,3:iz-2)*difx(3:ix-2,3:iz-2) ) &
             + ddz * ( kapz(3:ix-2,4:iz-1)*difz(3:ix-2,4:iz-1)     &
                       - kapz(3:ix-2,3:iz-2)*difz(3:ix-2,3:iz-2) ) 
    !$omp end workshare
    !$omp end parallel 
    deallocate(difx,difz)
end subroutine 

subroutine flux(box, fx, fz)
    use defstruct
    implicit none
    type(cell) :: box, fx, fz
    
    integer :: i,j
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision :: b2,roi,h,vx,vy,vz,bx,by,bz,pr
    double precision :: eta, ex, ey, ez
    double precision :: jx, jy, jz

    !$omp parallel do private(i,roi,vx,vy,vz,bx,by,bz,pr,b2,h,jx,jy,jz,eta,ex,ey,ez)
    do j=2,iz-1
        do i=2,ix-1
            roi = 1./box%ro(i,j)
            vx = box%rovx(i,j) * roi
            vy = box%rovy(i,j) * roi
            vz = box%rovz(i,j) * roi
            bx = box%bx(i,j)
            by = box%by(i,j)
            bz = box%bz(i,j)
            pr = box%pr(i,j)
            b2 = bx**2 + by**2 + bz**2
            h = 0.5*(vx**2+vy**2+vz**2)*box%ro(i,j) &
                        + pr*box%con%gam/(box%con%gam-1.)

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
            fx%e(i,j) = h*vx + (ey*bz - ez*by) 
            fx%bpot(i,j) = 0

            fz%ro(i,j) = box%rovz(i,j)
            fz%rovx(i,j) = vz*box%rovx(i,j) - bz*bx
            fz%rovy(i,j) = vz*box%rovy(i,j) - bz*by
            fz%rovz(i,j) = vz*box%rovz(i,j) - bz*bz + pr + 0.5*b2
            fz%bx(i,j) = -ey
            fz%by(i,j) = ex
            fz%bz(i,j) = 0.
            fz%e(i,j) = h*vz + (ex*by - ey*bx)
            fz%bpot(i,j) = 0
        end do
    end do
    !$omp end parallel do 
    

end subroutine

subroutine source(box, s)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box, s
    integer :: i
    double precision :: fugou(iz)
    
    do i=1,iz
        if (box%z(i)>0) then 
            fugou = 1.
        else if (box%z(i)<0) then
            fugou = -1.
        else 
            fugou = 0.
        end if
    end do

    !$omp parallel workshare
    s%ro = 0.
    s%bx = 0.
    s%by = 0.
    s%bz = 0.
    s%rovx = box%ro*box%con%gx
    !s%rovy = box%ro*box%con%gy
    s%rovy = box%ro*box%con%gy
    s%rovz = box%ro*box%con%gz*spread(fugou,1,ix)
    s%e = box%rovx*box%con%gx + box%rovy*box%con%gy + box%rovz*box%con%gz

    s%bpot = (box%rovx*box%bz - box%rovz*box%bx)/box%ro
    !$omp end parallel workshare
end subroutine

end module 
