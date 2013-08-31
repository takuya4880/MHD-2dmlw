subroutine output(box)
    use defstruct
    implicit none
    type(cell) :: box

    integer :: i,j,m
    double precision :: v2(ix,iy), dt(ix,iy)
    m = box%con%marg

    do j=1,iy
        do i=1,ix
            write(23, *)  box%x(i), box%y(j), box%ro(i,j), box%pr(i,j), box%rovx(i,j), box%rovy(i,j)
        end do
        write(23,*) " "
    end do
    write(23,*) " "

end subroutine
