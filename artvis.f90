subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f

    call eachav(box%ro,f%ro,box%con)
    call eachav(box%rovx,f%rovx,box%con)
    call eachav(box%rovy,f%rovy,box%con)
    call eachav(box%rovz,f%rovz,box%con)
    call eachav(box%bx,f%bx,box%con)
    call eachav(box%by,f%by,box%con)
    call eachav(box%bz,f%bz,box%con)
    call eachav(box%e,f%e,box%con)

end subroutine
!contains    
subroutine eachav(box,f,con)
    use defstruct
    double precision :: box(ix,iy),f(ix,iy)
    type(constants) con
    
    double precision :: ddx, ddy
    double precision, allocatable :: difx(:,:),dify(:,:),kapx(:,:),kapy(:,:)
    allocate(difx(ix,iy),dify(ix,iy),kapx(ix,iy),kapy(ix,iy))
    
    ddx = con%dt / con%dx
    ddy = con%dt / con%dy
    
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    dify(:,2:iy) = f(:,2:iy) - f(:,1:iy-1)
    kapx = con%q * abs(difx)
    kapy = con%q * abs(dify)

    box(3:ix-2,3:iy-2) = f(3:ix-2,3:iy-2) &
             + ddx * ( kapx(4:ix-1,3:iy-2)*difx(4:ix-1,3:iy-2)     &
                       - kapx(3:ix-2,3:iy-2)*difx(3:ix-2,3:iy-2) ) &
             + ddy * ( kapy(3:ix-2,4:iy-1)*dify(3:ix-2,4:iy-1)     &
                       - kapy(3:ix-2,3:iy-2)*dify(3:ix-2,3:iy-2) ) 
    deallocate(difx,dify,kapx,kapy)
end subroutine 
!end subroutine 