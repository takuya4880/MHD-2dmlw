subroutine source(box, s)
    use defstruct
    implicit none
    type(cell) :: box, s
    integer :: i
    double precision :: fugou(iy)
    
    do i=1,iy
        if (box%y(i)>0) then 
            fugou = 1.
        else if (box%y(i)<0) then
            fugou = -1.
        else 
            fugou = 0.
        end if
    end do

    s%ro = 0.
    s%bx = 0.
    s%by = 0.
    s%bz = 0.
    s%rovx = box%ro*box%con%gx
    !s%rovy = box%ro*box%con%gy
    s%rovy = box%ro*box%con%gy*spread(fugou,1,ix)
    s%rovz = box%ro*box%con%gz
    s%e = box%rovx*box%con%gx + box%rovy*box%con%gy + box%rovz*box%con%gz
end subroutine
