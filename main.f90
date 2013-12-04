program main
    use defstruct   
    use ic
    use bc
    use op 
    use pr
    use dt
    use st
    !$use omp_lib
    implicit none
    
    type(cell),pointer :: box
    double precision :: uboundary(9,marg)
    double precision :: t, tint, tend, tnxt
    integer :: start(8), time(8), minits, timelimit
    character*10 :: tmp
    integer :: mcont

    call omp_set_num_threads(2)
    allocate(box)
    !open(23,file="result.dat",status="replace")

    call date_and_time(tmp,tmp,tmp,start)

    mcont = 0
    timelimit = 210  ! 210:3.5 hours
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
    box%op%mf_params=9

    t = 0.
    tint = 1.
    tnxt = tint
    tend = 110.

    call initial(box, uboundary)
    call boundary(box, uboundary)
    call outpinit(box)
    if (mcont==1) then
        call outputread(box,t)
        tnxt = dint(t) + tint
    end if 
    call outp(box,t)
    call pressure(box)

    do
        call detdt(box)    
        call step(box)
        call boundary(box, uboundary)
        t = t + box%con%dt
        print *, t, box%con%dt
        if (t>=tnxt) then
            call outp(box,t)
            tnxt = tnxt + tint
        endif
        if (t>tend) exit
        if (box%con%dt<1.e-10) exit
    
        call date_and_time(tmp,tmp,tmp,time)
        time = time - start
        minits = time(3)*24*60+time(5)*60+time(6)
        if (minits>timelimit) exit
        
    end do

end program main
        
