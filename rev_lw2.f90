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
!contains    
subroutine each2(box,h,f,fx,fz,s,ddx,ddz,dt)
    use defstruct
    implicit none
    double precision :: box(ix,iz),f(ix,iz),h(ix,iz)
    double precision :: fx(ix,iz),fz(ix,iz),s(ix,iz)
    double precision :: ddx, ddz, dt
    
    integer i,j
    double precision fffx,fffz,ss

    do i=1,iz-2
        do j=1,ix-2
            fffx = 0.5 * (fx(j+1,i+1)+fx(j+1,i)-fx(j,i+1)-fx(j,i))
            fffz = 0.5 * (fz(j+1,i+1)-fz(j+1,i)+fz(j,i+1)-fz(j,i))
            ss  = 0.25 * ( s(j+1,i+1) +s(j+1,i) +s(j,i+1) +s(j,i) )
            f(j+1,i+1) = box(j+1,i+1) - ddx*fffx - ddz*fffz + dt*ss
        end do
    end do

end subroutine 
!end subroutine 

