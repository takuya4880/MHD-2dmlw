program main
    use defstruct    
    implicit none
    
    type(cell),pointer :: box
    double precision :: t, tint, tend, tnxt
    double precision :: uboundary(9,2)

    allocate(box)

    open(23,file="result.dat",status="replace")

    box%con%nx = nx
    box%con%ny = ny
    box%con%ix = ix
    box%con%iy = iy
    box%con%marg = marg
    box%con%wid = 150.
    box%con%hig = 60.
    box%con%dx = box%con%wid/dble(box%con%nx-1)
    box%con%dy = box%con%hig/dble(box%con%ny-1)
    box%con%a = 0.4
    box%con%q = 3.
    box%con%gam = 5./3.

    t = 0.
    tint = 10.
    tnxt = tint
    tend = 200.

    call initial(box, uboundary)
    call boundary(box, uboundary)
    call output(box)
    call pressure(box)

    do
        call detdt(box)    
        call step(box)
        call boundary(box, uboundary)
        t = t + box%con%dt
        print *, t, box%con%dt
        if (t>=tnxt) then
            call output(box)
            tnxt = tnxt + tint
        endif
        if (t>tend) exit
    end do

    close(23)

end program main
        
