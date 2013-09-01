subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f

    call eachav(box%ro, f%ro, f%ro, box%con)
    call eachav(box%rovx, f%rovx, f%ro, box%con)
    call eachav(box%rovy, f%rovy, f%ro, box%con)
    call eachav(box%rovz, f%rovz, f%ro, box%con)
    call eachav(box%bx, f%bx, f%ro, box%con)
    call eachav(box%by, f%by, f%ro, box%con)
    call eachav(box%bz, f%bz, f%ro, box%con)
    call eachav(box%e, f%e, f%ro, box%con)

end subroutine
!contains    
subroutine eachav(box,f,ro,con)
    use defstruct
    double precision :: box(ix,iy),f(ix,iy),ro(ix,iy)
    type(constants) con
    
    double precision :: ddx, ddy
    double precision, allocatable :: difx(:,:),dify(:,:),kapx(:,:),kapy(:,:)
    allocate(difx(ix,iy),dify(ix,iy),kapx(ix,iy),kapy(ix,iy))
    ddx = con%dt / con%dx
    ddy = con%dt / con%dy
    
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    dify(:,2:iy) = f(:,2:iy) - f(:,1:iy-1)
    kapx = con%q * abs(difx) / ro
    kapy = con%q * abs(dify) / ro

    box(3:ix-2,3:iy-2) = f(3:ix-2,3:iy-2) &
             + ddx * ( kapx(4:ix-1,3:iy-2)*difx(4:ix-1,3:iy-2)     &
                       - kapx(3:ix-2,3:iy-2)*difx(3:ix-2,3:iy-2) ) &
             + ddy * ( kapy(3:ix-2,4:iy-1)*dify(3:ix-2,4:iy-1)     &
                       - kapy(3:ix-2,3:iy-2)*dify(3:ix-2,3:iy-2) ) 
    deallocate(difx,dify,kapx,kapy)
end subroutine 
!end subroutine 
