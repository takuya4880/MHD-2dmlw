subroutine step(box)
    use defstruct
    implicit none
    type(cell) :: box
    type(cell), pointer :: fx, fy, s, h, f
    integer :: i
    allocate(fx, fy, s, h, f)
    fx = box
    fy = box
    s = box
    h = box
    f = box
        
    forall(i=1:ix) h%x(i)=h%con%dx*(i-h%con%marg-1.)
    forall(i=1:iy) h%y(i)=h%con%dy*(i-h%con%marg-1.)
    
    call flux(box, fx, fy)
    call source(box, s)
    call lw1(box, h, fx, fy, s)
    call pressure(h)
    
    call flux(h, fx, fy)
    call source(h, s)
    call lw2(box, h, f, fx, fy, s)

    call artvis(box, f)
    call pressure(box)

    deallocate(fx, fy, s, h, f)

end subroutine
