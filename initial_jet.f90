subroutine initial(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box
    double precision :: uboundary(9,marg)

    integer :: i,j,m,origin
    double precision :: gami               !inberse of gamma
    double precision :: wid
    double precision :: amp, tpt, tpho, tcor, x, y, a, ad, lp, phicor
    double precision :: w, ytr, yfsl, yfsu, ymgc, yinv, betafs, betacor
    double precision :: den(iy), pre(iy),temp(iy),phi(iy)
    double precision :: beta(iy),beta1i(iy),beta2i(iy),b(iy),ee(iy) 
    m = box%con%marg
    wid = box%con%wid
    gami = box%con%gam

    amp = 0.05
    lp = 20. 
    phicor = 4.*atan(1.d0) 
    tpt = 25.
    tpho = 1.
    tcor = tpho * tpt
    w = 0.5
    ytr = 8.
    yfsl = -4.
    yfsu = -2.
    ymgc = 9.5
    yinv = 5.
    betafs=4.
    betacor=0.2
    a = 2.
    ad = a*(box%con%gam-1.)/box%con%gam

    origin = int(5./box%con%hig*ny)+1+m
    
    forall(i=1:ix) box%x(i)=box%con%dx*(i-m)
    forall(i=1:iy) box%y(i)=box%con%dy*(i-origin)
    
    box%con%gx = 0.
    box%con%gy = -gami
    box%con%gz = 0.

    do i=origin+1,iy
        temp(i) = tpho + 0.5 * (tcor-tpho)*(tanh((box%y(i)-ytr)/w) + 1.)
    end do
    do i=origin,1,-1
        temp(i) = tpho - ad*box%y(i)
    end do

    den(origin) = 1.
    pre(origin) = gami*temp(origin)
    
    beta1i = 0.25/betafs*(tanh((box%y-yfsl)/w) + 1.) * (-tanh((box%y-yfsu)/w) + 1.)
    beta2i = 0.5/betacor*(tanh((box%y-ymgc)/w) + 1.)
    beta = 1./(beta1i+beta2i)
    do i=origin+1,iy
        den(i) = den(i-1) * ((1.+1./beta(i-1))*temp(i-1) + 0.5*box%con%gam*box%con%dy*box%con%gy)&
                          / ((1.+1./beta(i))*temp(i) - 0.5*box%con%gam*box%con%dy*box%con%gy)
    end do
    do i=origin-1,1,-1
        den(i) = den(i+1) * ((1.+1./beta(i+1))*temp(i+1) - 0.5*box%con%gam*box%con%dy*box%con%gy)&
                          / ((1.+1./beta(i))*temp(i) + 0.5*box%con%gam*box%con%dy*box%con%gy)
    end do
    
    pre = pre(origin) * (den/den(origin)) * (temp/temp(origin))
    b = sqrt(2.*pre/beta)
    
    do i=1,iy
        if (box%y(i)<yinv) then
            phi(i)=0.
        else
            phi(i)=phicor
        end if
    end do

    open(24,file="initial.dat",status="replace")
    do i=1,iy
       write (24,*) box%y(i), den(i), pre(i), temp(i), 1./beta(i), b(i), phi(i)
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
            if (y>yfsl .and. y<yfsu .and. x>0.5*wid-0.25*lp .and. x<0.5*wid+0.25*lp) then
                box%rovx(i,j) = box%ro(i,j)*amp*sin(8.*atan(1.d0)*(x-0.5*wid)/lp)
            end if
        end do
    end do

    box%bx = spread(b*cos(phi),1,ix)
    box%by = spread(b*sin(phi),1,ix)
    box%bz = 0.
    box%pr = spread(pre,1,ix)  
    box%e = 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
            + box%pr/(box%con%gam-1.) &
            + 0.5*(box%bx**2 + box%by**2 + box%bz**2)

    box%bpot(1,1)=0.
    do i=2,ix
        box%bpot(i,1) = box%bpot(i-1,1) &
                        - 0.5*box%con%dx*(box%by(i,1)+box%by(i-1,1))
    end do
    do i=1,ix
        do j=2,iy
            box%bpot(i,j) = box%bpot(i,j-1) &
                            + 0.5*box%con%dy*(box%bx(i,j)+box%bx(i,j-1))
        end do
    end do

    uboundary(1,1:marg) = den(iy-marg+1:iy)
    uboundary(2:7,1:marg) = 0
    uboundary(8,1:marg) = box%e(1,iy-marg+1:iy)
    uboundary(9,1:marg) = pre(iy-marg+1:iy)
    

end subroutine
