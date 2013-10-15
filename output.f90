subroutine output(box)
    use defstruct
    implicit none
    type(cell) :: box

    integer :: i,j,m
    double precision :: v2(ix,iz), dt(ix,iz)
    m = box%con%marg

    do j=1,iz
        do i=1,ix
            write(23, *)  box%x(i), box%z(j), box%ro(i,j), box%pr(i,j), box%rovx(i,j), box%rovz(i,j), box%bpot(i,j)
        end do
        write(23,*) " "
    end do
    write(23,*) " "

    !i=100
    !do j=1,iz
    !    write(23, *)  box%z(j), box%ro(i,j), box%pr(i,j), box%rovx(i,j), box%rovz(i,j), box%bpot(i,j)
    !end do
    !write(23,*) " "
    !write(23,*) " "

end subroutine

