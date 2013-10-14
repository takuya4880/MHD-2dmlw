module l1
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
!contains
subroutine each1(u,h,fx,fz,s,ddx,ddz,dt)
    use defstruct
    implicit none
    double precision :: u(ix,iz),h(ix,iz),fx(ix,iz),fz(ix,iz),s(ix,iz)
    double precision :: ddx, ddz, dt
    
    integer i,j
    double precision fffx, fffz, ss

    do i=1,iz-1
        do j=1,ix-1
            fffx = 0.5 * ( fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i) )
            fffz = 0.5 * ( fz(j+1,i+1)-fz(j+1,i)+fz(j,i+1)-fz(j,i) )
            ss   = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            h(j,i) = 0.25 * ( u(j+1,i+1)+u(j+1,i)+u(j,i+1)+u(j,i) ) &
                        - 0.5*(ddx*fffx + ddz*fffz - dt*ss)
        end do
    end do

end subroutine

end module 
