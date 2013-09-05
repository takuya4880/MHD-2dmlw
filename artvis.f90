subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f
    double precision, allocatable :: kapx(:,:),kapy(:,:)
    allocate(kapx(ix,iy),kapy(ix,iy))
    
    kapx(2:ix,:) = box%con%q * abs(box%rovx(2:ix,:)/box%ro(2:ix,:) &
                                    - box%rovx(1:ix-1,:)/box%ro(1:ix-1,:))
    kapy(:,2:iy) = box%con%q * abs(box%rovy(:,2:iy)/box%ro(:,2:iy) &
                                    - box%rovy(:,1:iy-1)/box%ro(:,1:iy-1)) 

    call eachav(box%ro, f%ro, kapx, kapy, box%con)
    call eachav(box%rovx, f%rovx, kapx, kapy, box%con)
    call eachav(box%rovy, f%rovy, kapx, kapy, box%con)
    call eachav(box%rovz, f%rovz, kapx, kapy, box%con)
    call eachav(box%bx, f%bx, kapx, kapy, box%con)
    call eachav(box%by, f%by, kapx, kapy, box%con)
    call eachav(box%bz, f%bz, kapx, kapy, box%con)
    call eachav(box%e, f%e, kapx, kapy, box%con)

    deallocate(kapx,kapy)

end subroutine
!contains    
subroutine eachav(box,f,kapx,kapy,con)
    use defstruct
    double precision :: box(ix,iy),f(ix,iy)
    double precision :: kapx(ix,iy),kapy(ix,iy)
    type(constants) con
    
    double precision :: ddx, ddy
    double precision, allocatable :: difx(:,:),dify(:,:)
    allocate(difx(ix,iy),dify(ix,iy))
    ddx = con%dt / con%dx
    ddy = con%dt / con%dy
    
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    dify(:,2:iy) = f(:,2:iy) - f(:,1:iy-1)

    box(3:ix-2,3:iy-2) = f(3:ix-2,3:iy-2) &
             + ddx * ( kapx(4:ix-1,3:iy-2)*difx(4:ix-1,3:iy-2)     &
                       - kapx(3:ix-2,3:iy-2)*difx(3:ix-2,3:iy-2) ) &
             + ddy * ( kapy(3:ix-2,4:iy-1)*dify(3:ix-2,4:iy-1)     &
                       - kapy(3:ix-2,3:iy-2)*dify(3:ix-2,3:iy-2) ) 
    deallocate(difx,dify)
end subroutine 
!end subroutine 
