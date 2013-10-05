subroutine artvis(box, f)
    use defstruct
    implicit none
    type(cell) :: box, f
    double precision, allocatable :: kapx(:,:),kapz(:,:)
    allocate(kapx(ix,iz),kapz(ix,iz))
    
    kapx(2:ix,:) = box%con%q * abs(box%rovx(2:ix,:)/box%ro(2:ix,:) &
                                    - box%rovx(1:ix-1,:)/box%ro(1:ix-1,:))
    kapz(:,2:iz) = box%con%q * abs(box%rovz(:,2:iz)/box%ro(:,2:iz) &
                                    - box%rovz(:,1:iz-1)/box%ro(:,1:iz-1)) 

    call eachav(box%ro, f%ro, kapx, kapz, box%con)
    call eachav(box%rovx, f%rovx, kapx, kapz, box%con)
    call eachav(box%rovy, f%rovy, kapx, kapz, box%con)
    call eachav(box%rovz, f%rovz, kapx, kapz, box%con)
    call eachav(box%bx, f%bx, kapx, kapz, box%con)
    call eachav(box%by, f%by, kapx, kapz, box%con)
    call eachav(box%bz, f%bz, kapx, kapz, box%con)
    call eachav(box%e, f%e, kapx, kapz, box%con)
    call eachav(box%bpot, f%bpot, kapx, kapz, box%con)

    deallocate(kapx,kapz)

end subroutine
!contains    
subroutine eachav(box,f,kapx,kapz,con)
    use defstruct
    double precision :: box(ix,iz),f(ix,iz)
    double precision :: kapx(ix,iz),kapz(ix,iz)
    type(constants) con
    
    double precision :: ddx, ddz
    double precision, allocatable :: difx(:,:),difz(:,:)
    allocate(difx(ix,iz),difz(ix,iz))
    ddx = con%dt / con%dx
    ddz = con%dt / con%dz
    
    difx(2:ix,:) = f(2:ix,:) - f(1:ix-1,:)
    difz(:,2:iz) = f(:,2:iz) - f(:,1:iz-1)

    box(3:ix-2,3:iz-2) = f(3:ix-2,3:iz-2) &
             + ddx * ( kapx(4:ix-1,3:iz-2)*difx(4:ix-1,3:iz-2)     &
                       - kapx(3:ix-2,3:iz-2)*difx(3:ix-2,3:iz-2) ) &
             + ddz * ( kapz(3:ix-2,4:iz-1)*difz(3:ix-2,4:iz-1)     &
                       - kapz(3:ix-2,3:iz-2)*difz(3:ix-2,3:iz-2) ) 
    deallocate(difx,difz)
end subroutine 
!end subroutine 
