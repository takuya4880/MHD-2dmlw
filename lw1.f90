subroutine lw1(box, h, fx, fy)
    use defstruct
    implicit none
    type(cell) :: box, h, fx, fy
    
    double precision, allocatable :: ffx(:,:), ffy(:,:)
    double precision :: ddx, ddy
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy
    allocate(ffx(ix, iy), ffy(ix, iy))

    call each1(box%ro,h%ro,fx%ro,fy%ro,ffx,ffy,ddx,ddy)
    call each1(box%rovx,h%rovx,fx%rovx,fy%rovx,ffx,ffy,ddx,ddy)
    call each1(box%rovy,h%rovy,fx%rovy,fy%rovy,ffx,ffy,ddx,ddy)
    call each1(box%rovz,h%rovz,fx%rovz,fy%rovz,ffx,ffy,ddx,ddy)
    call each1(box%bx,h%bx,fx%bx,fy%bx,ffx,ffy,ddx,ddy)
    call each1(box%by,h%by,fx%by,fy%by,ffx,ffy,ddx,ddy)
    call each1(box%bz,h%bz,fx%bz,fy%bz,ffx,ffy,ddx,ddy)
    call each1(box%e,h%e,fx%e,fy%e,ffx,ffy,ddx,ddy)
    
    deallocate(ffx,ffy)
end subroutine

!contains
subroutine each1(u,h,fx,fy,ffx,ffy,ddx,ddy)
    use defstruct
    implicit none
    double precision :: u(ix,iy),h(ix,iy),fx(ix,iy),fy(ix,iy)
    double precision :: ffx(ix,iy),ffy(ix,iy)
    double precision :: ddx, ddy

    ffx(1:ix-1,1:iy-1) =  ( fx(2:ix,2:iy)   + fx(2:ix,1:iy-1) &
                                - fx(1:ix-1,2:iy) - fx(1:ix-1,1:iy-1) ) 
    ffy(1:ix-1,1:iy-1) =  ( fy(2:ix,2:iy)   - fy(2:ix,1:iy-1) &
                                + fy(1:ix-1,2:iy) - fy(1:ix-1,1:iy-1) ) 
    h(1:ix-1,1:iy-1) = 0.25 * ( u(2:ix,2:iy)   + u(2:ix,1:iy-1) &
                              + u(1:ix-1,2:iy) + u(1:ix-1,1:iy-1) &
                      - ddx*ffx(1:ix-1,1:iy-1) - ddy*ffy(1:ix-1,1:iy-1))

end subroutine

!end subroutine
