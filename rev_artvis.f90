subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f
    
    double precision :: ddx, ddy
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy

    call eachav(box%ro,f%ro,f%ro,box%con,ddx,ddy)
    call eachav(box%rovx,f%rovx,f%ro,box%con,ddx,ddy)
    call eachav(box%rovy,f%rovy,f%ro,box%con,ddx,ddy)
    call eachav(box%rovz,f%rovz,f%ro,box%con,ddx,ddy)
    call eachav(box%bx,f%bx,f%ro,box%con,ddx,ddy)
    call eachav(box%by,f%by,f%ro,box%con,ddx,ddy)
    call eachav(box%bz,f%bz,f%ro,box%con,ddx,ddy)
    call eachav(box%e,f%e,f%ro,box%con,ddx,ddy)
    

end subroutine
!contains    
subroutine eachav(box,f,ro,con,ddx,ddy)
    use defstruct
    implicit none
    double precision :: box(ix,iy),f(ix,iy),ro(ix,iy)
    type(constants) :: con
    double precision :: ddx,ddy

    double precision :: difx1,difx2,dify1,dify2,kapx1,kapx2,kapy1,kapy2
    integer :: i,j
    
    do j=3,iy-2
        do i=3,ix-2

            difx1 = f(i+1,j) - f(i,j)
            kapx1 = con%q * abs(difx1) / ro(i,j)
            difx2 = f(i,j) - f(i-1,j)
            kapx2 = con%q * abs(difx2) / ro(i,j)
            dify1 = f(i,j+1) - f(i,j)
            kapy1 = con%q * abs(dify1) / ro(i,j)
            dify2 = f(i,j) - f(i,j-1)
            kapy2 = con%q * abs(dify2) / ro(i,j)

            box(i,j) = f(i,j) + ddx*(kapx1*difx1-kapx2*difx2) &
                              + ddy*(kapy1*dify1-kapy2*dify2)
        end do
    end do

end subroutine
!end subroutine
