program main
    use defstruct   
    use ic 
    use bc 
    use op 
    use pr 
    use dt 
    use st 
    use l1 
    use l2 
    use sc 
    use flx 
    use av
    implicit none
    
    type(cell),pointer :: box
    double precision :: uboundary(9,marg)
    double precision :: t, tint, tend, tnxt

    allocate(box)

    box%con%nx = nx
    box%con%nz = nz
    box%con%ix = ix
    box%con%iz = iz
    box%con%marg = marg
    box%con%wid = 150.
    box%con%hig = 60.
    box%con%dx = box%con%wid/dble(box%con%nx-1)
    box%con%dz = box%con%hig/dble(box%con%nz-1)
    box%con%a = 0.4
    box%con%q = 3.
    box%con%gam = 5./3.

    t = 0.
    tint = 1.
    tnxt = tint
    tend = 1.


    open(23,file="result.dat",status="replace")

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
        
