subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f
    
    double precision, allocatable :: difx(:,:),dify(:,:),kapx(:,:),kapy(:,:)
    double precision :: ddx, ddy
    allocate(difx(ix,iy),dify(ix,iy),kapx(ix,iy),kapy(ix,iy))
    ddx = box%con%dt / box%con%dx
    ddy = box%con%dt / box%con%dy

    call eachav(box%ro,f%ro,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%rovx,f%rovx,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%rovy,f%rovy,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%rovz,f%rovz,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%bx,f%bx,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%by,f%by,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%bz,f%bz,box%con,ddx,ddy,difx,dify,kapx,kapy)
    call eachav(box%e,f%e,box%con,ddx,ddy,difx,dify,kapx,kapy)
    
    deallocate(difx,dify,kapx,kapy)

end subroutine
!contains    
subroutine eachav(box,f,con,ddx,ddy,difx,dify,kapx,kapy)
    use defstruct
    double precision :: box(ix,iy),f(ix,iy)
    type(constants) con
    double precision :: ddx,ddy
    double precision :: difx(ix,iy),dify(ix,iy),kapx(ix,iy),kapy(ix,iy) 
    
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    dify(:,2:iy) = f(:,2:iy) - f(:,1:iy-1)
    kapx = con%q * abs(difx)
    kapy = con%q * abs(dify)

    box(3:ix-2,3:iy-2) = f(3:ix-2,3:iy-2) &
             + ddx * ( kapx(4:ix-1,3:iy-2)*difx(4:ix-1,3:iy-2)     &
                       - kapx(3:ix-2,3:iy-2)*difx(3:ix-2,3:iy-2) ) &
             + ddy * ( kapy(3:ix-2,4:iy-1)*dify(3:ix-2,4:iy-1)     &
                       - kapy(3:ix-2,3:iy-2)*dify(3:ix-2,3:iy-2) ) 
end subroutine 
!end subroutine 
