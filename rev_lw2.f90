subroutine lw2(box, h, f, fx, fy, s)
    use defstruct
    implicit none
    type(cell) :: box, h, f, fx, fy, s
    
    double precision :: ddx, ddy, dt
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy
    dt  = box%con%dt

    call each2(box%ro,h%ro,f%ro,fx%ro,fy%ro,s%ro,ddx,ddy,dt)
    call each2(box%rovx,h%rovx,f%rovx,fx%rovx,fy%rovx,s%rovx,ddx,ddy,dt)
    call each2(box%rovy,h%rovy,f%rovy,fx%rovy,fy%rovy,s%rovy,ddx,ddy,dt)
    call each2(box%rovz,h%rovz,f%rovz,fx%rovz,fy%rovz,s%rovz,ddx,ddy,dt)
    call each2(box%bx,h%bx,f%bx,fx%bx,fy%bx,s%bx,ddx,ddy,dt)
    call each2(box%by,h%by,f%by,fx%by,fy%by,s%by,ddx,ddy,dt)
    call each2(box%bz,h%bz,f%bz,fx%bz,fy%bz,s%bz,ddx,ddy,dt)
    call each2(box%e,h%e,f%e,fx%e,fy%e,s%e,ddx,ddy,dt)
   

end subroutine
!contains    
subroutine each2(box,h,f,fx,fy,s,ddx,ddy,dt)
    use defstruct
    implicit none
    double precision :: box(ix,iy),f(ix,iy),h(ix,iy)
    double precision :: fx(ix,iy),fy(ix,iy),s(ix,iy)
    double precision :: ddx, ddy, dt
    
    integer i,j
    double precision fffx,fffy,ss

    do i=1,iy-2
        do j=1,ix-2
            fffx = 0.5 * (fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i))
            fffy = 0.5 * (fy(j+1,i+1)-fy(j+1,i)+fy(j,i+1)-fy(j,i))
            ss  = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            f(j+1,i+1) = box(j+1,i+1) - ddx*fffx - ddy*fffy + dt*ss
        end do
    end do

end subroutine 
!end subroutine 

