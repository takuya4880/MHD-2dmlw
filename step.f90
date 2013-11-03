module st
implicit none 
contains
subroutine step(box)
    use defstruct
    use lw
    use pr

    implicit none
    type(cell) :: box
    type(cell), pointer :: fx, fz, s, h, f
    integer :: i
    allocate(fx, fz, s, h, f)
    fx = box
    fz = box
    s = box
    h = box
    h%con%a = -1.   !use a to knowtify this is half step value for source term
    f = box
        
    !forall(i=1:ix) h%x(i)=h%con%dx*(i-h%con%marg-1.)
    !forall(i=1:iz) h%z(i)=h%con%dz*(i-h%con%marg-1.)
    h%x = box%x + 0.5*box%con%dx
    h%z = box%z + 0.5*box%con%dz

    call flux(box, fx, fz)
    call source(box, s)
    call lw1(box, h, fx, fz, s)
    call pressure(h)
    
    call flux(h, fx, fz)
    call source(h, s)
    call lw2(box, h, f, fx, fz, s)

    call artvis(box, f)
    call pressure(box)

    deallocate(fx, fz, s, h, f)

end subroutine
end module
