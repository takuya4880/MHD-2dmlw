subroutine initial(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box
    double precision :: uboundary(9,2)

    integer :: i,j,m
    double precision :: gami               !inberse of gamma
    double precision :: wid
    double precision :: amp, tpt, tch, tcor, ycor, ym, x, y
    double precision :: wtr, w0, w1, y0, y1, betast
    double precision :: den(iy), pre(iy),temp(iy),beta(iy),f(iy),b(iy),ee(iy) 
    m = box%con%marg
    wid = box%con%wid

    amp = 0.05 
    tpt = 25.
    tch = 1.
    tcor = tch * tpt
    ycor = 18.
    wtr = 0.6
    w0 = 0.5
    w1 = 0.5
    y0 = 4.
    y1 = 8.
    ym = y1-y0
    betast=1.
    
    forall(i=1:ix) box%x(i)=box%con%dx*(i-m-0.5)
    forall(i=1:iy) box%y(i)=box%con%dy*(i-m-0.5)
    
    gami = 1./box%con%gam
    box%con%gx = 0.
    box%con%gy = -gami
    box%con%gz = 0.

    den = 1.
    den(m+1) = 1.
    pre(m+1) = gami
    
    temp = tch + 0.5 * (tcor-tch)*(tanh((box%y-ycor)/wtr) + 1.)
    f = 0.25 * (tanh((box%y-y0)/w0) + 1.) * (-tanh((box%y-y1)/w1) + 1.)
    beta = betast/f
    do i=m+2,iy
        den(i) = den(i-1) * ((1.+1./beta(i-1))*temp(i-1) + 0.5*box%con%gam*box%con%dy*box%con%gy)&
                          / ((1.+1./beta(i))*temp(i) - 0.5*box%con%gam*box%con%dy*box%con%gy)
    end do
    pre = pre(m+1) * (den/den(m+1)) * (temp/temp(m+1))
    b = sqrt(2.*pre/beta)
    
    open(24,file="initial.dat",status="replace")
    do i=1,iy
       write (24,*) box%y(i), den(i), pre(i), temp(i), f(i), beta(i), b(i)
    end do
    close(24)

    box%ro = spread(den,1,ix)
    box%rovx = 0.
    box%rovy = 0.
    box%rovz = 0.
    do i=1,ix
        do j=1,iy
            x = box%x(i)
            y = box%y(j)
            if (y>y0 .and. y<y1 .and. x>0.25*wid .and. x<0.75*wid) then
                box%rovx(i,j) = box%ro(i,j)*amp*f(j)*sin(16.*atan(1.d0)*(x-0.5*wid)/wid)
            end if
        end do
    end do
    box%bx = spread(b,1,ix)
    box%by = 0.
    box%bz = 0.
    box%pr = spread(pre,1,ix)  
    box%e = 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
            + box%pr/(box%con%gam-1.) &
            + 0.5*(box%bx**2 + box%by**2 + box%bz**2)

    uboundary(1,1:2) = den(iy-1:iy)
    uboundary(2:7,1:2) = 0
    uboundary(8,1:2) = box%e(1,iy-1:iy)
    uboundary(9,1:2) = pre(iy-1:iy)

end subroutine
