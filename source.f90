module sc
implicit none
contains
subroutine source(box, s)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box, s
    integer :: i
    double precision :: fugou(iz)
    
    do i=1,iz
        if (box%z(i)>0) then 
            fugou = 1.
        else if (box%z(i)<0) then
            fugou = -1.
        else 
            fugou = 0.
        end if
    end do

    !$omp parallel workshare
    s%ro = 0.
    s%bx = 0.
    s%by = 0.
    s%bz = 0.
    s%rovx = box%ro*box%con%gx
    !s%rovy = box%ro*box%con%gy
    s%rovy = box%ro*box%con%gy
    s%rovz = box%ro*box%con%gz*spread(fugou,1,ix)
    s%e = box%rovx*box%con%gx + box%rovy*box%con%gy + box%rovz*box%con%gz

    s%bpot = (box%rovx*box%bz - box%rovz*box%bx)/box%ro
    !$omp end parallel workshare
end subroutine
end module 
