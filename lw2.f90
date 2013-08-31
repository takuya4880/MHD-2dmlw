subroutine lw2(box, h, f, fx, fy)
    use defstruct
    implicit none
    type(cell) :: box, h, f, fx, fy
    
    double precision, allocatable :: ffx(:,:), ffy(:,:)
    double precision :: ddx, ddy
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy
    allocate(ffx(ix, iy), ffy(ix, iy))

    call each2(box%ro,h%ro,f%ro,fx%ro,fy%ro,ffx,ffy,ddx,ddy)
    call each2(box%rovx,h%rovx,f%rovx,fx%rovx,fy%rovx,ffx,ffy,ddx,ddy)
    call each2(box%rovy,h%rovy,f%rovy,fx%rovy,fy%rovy,ffx,ffy,ddx,ddy)
    call each2(box%rovz,h%rovz,f%rovz,fx%rovz,fy%rovz,ffx,ffy,ddx,ddy)
    call each2(box%bx,h%bx,f%bx,fx%bx,fy%bx,ffx,ffy,ddx,ddy)
    call each2(box%by,h%by,f%by,fx%by,fy%by,ffx,ffy,ddx,ddy)
    call each2(box%bz,h%bz,f%bz,fx%bz,fy%bz,ffx,ffy,ddx,ddy)
    call each2(box%e,h%e,f%e,fx%e,fy%e,ffx,ffy,ddx,ddy)
   
    deallocate(ffx,ffy) 

end subroutine
!contains    
subroutine each2(box,h,f,fx,fy,ffx,ffy,ddx,ddy)
    use defstruct
    implicit none
    double precision :: box(ix,iy),f(ix,iy),h(ix,iy),fx(ix,iy),fy(ix,iy)
    double precision :: ffx(ix,iy),ffy(ix,iy)
    double precision :: ddx, ddy

    ffx(1:ix-2,1:iy-2) = 0.5 * ( fx(2:ix-1,2:iy-1) + fx(2:ix-1,1:iy-2) &
                               - fx(1:ix-2,2:iy-1) - fx(1:ix-2,1:iy-2) ) 
    ffy(1:ix-2,1:iy-2) = 0.5 * ( fy(2:ix-1,2:iy-1) - fy(2:ix-1,1:iy-2) &
                               + fy(1:ix-2,2:iy-1) - fy(1:ix-2,1:iy-2) ) 
    f(2:ix-1,2:iy-1) = box(2:ix-1,2:iy-1) - ffx(1:ix-2,1:iy-2) * ddx &
                                  - ffy(1:ix-2,1:iy-2) * ddy

end subroutine 
!end subroutine 

