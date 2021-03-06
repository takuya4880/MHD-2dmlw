
module dt
implicit none 
contains
subroutine detdt(box)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box

    double precision, allocatable :: v2(:,:)
    double precision d
    allocate(v2(ix,iz))
    
    if (box%con%dx<box%con%dz) then
        d = box%con%dx
    else
        d = box%con%dz
    end if

    !$omp parallel workshare
    v2 = (box%rovx**2 + box%rovy**2 + box%rovz**2)/(box%ro**2)
    v2 = v2 + (box%bx**2 + box%bz**2 + box%bz**2)/box%ro
    v2 = v2 + box%con%gam * box%pr / box%ro
    !$omp end parallel workshare

    box%con%dt = box%con%a * d / sqrt(maxval(v2))
    deallocate(v2)
end subroutine

end module
