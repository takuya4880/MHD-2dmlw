subroutine lw1(box, h, fx, fy, s)
    use defstruct
    implicit none
    type(cell) :: box, h, fx, fy, s
    
    double precision :: ddx, ddy, dt
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy
    dt  = box%con%dt

    call each1(box%ro,h%ro,fx%ro,fy%ro,s%ro,ddx,ddy,dt)
    call each1(box%rovx,h%rovx,fx%rovx,fy%rovx,s%rovx,ddx,ddy,dt)
    call each1(box%rovy,h%rovy,fx%rovy,fy%rovy,s%rovy,ddx,ddy,dt)
    call each1(box%rovz,h%rovz,fx%rovz,fy%rovz,s%rovz,ddx,ddy,dt)
    call each1(box%bx,h%bx,fx%bx,fy%bx,s%bx,ddx,ddy,dt)
    call each1(box%by,h%by,fx%by,fy%by,s%by,ddx,ddy,dt)
    call each1(box%bz,h%bz,fx%bz,fy%bz,s%bz,ddx,ddy,dt)
    call each1(box%e,h%e,fx%e,fy%e,s%e,ddx,ddy,dt)
    call each1(box%bpot,h%bpot,fx%bpot,fy%bpot,s%bpot,ddx,ddy,dt)
    
    
end subroutine
!contains
subroutine each1(u,h,fx,fy,s,ddx,ddy,dt)
    use defstruct
    implicit none
    double precision :: u(ix,iy),h(ix,iy),fx(ix,iy),fy(ix,iy),s(ix,iy)
    double precision :: ddx, ddy, dt
    
    integer i,j
    double precision fffx, fffy, ss

    do i=1,iy-1
        do j=1,ix-1
            fffx = 0.5 * ( fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i) )
            fffy = 0.5 * ( fy(j+1,i+1)-fy(j+1,i)+fy(j,i+1)-fy(j,i) )
            ss   = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            h(j,i) = 0.25 * ( u(j+1,i+1)+u(j+1,i)+u(j,i+1)+u(j,i) ) &
                        - 0.5*(ddx*fffx + ddy*fffy - dt*ss)
        end do
    end do

end subroutine

!end subroutine
