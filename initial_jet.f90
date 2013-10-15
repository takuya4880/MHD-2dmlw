module ic
implicit none 
contains
subroutine initial(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box
    double precision :: uboundary(9,marg)

    integer :: i,j,m,origin
    double precision :: gami               !inberse of gamma
    double precision :: wid
    double precision :: amp, tpt, tpho, tcor, x, z, a, ad, lp, phicor
    double precision :: w, ztr, zfsl, zfsu, zmgc, zinv, betafs, betacor
    double precision :: den(iz), pre(iz),temp(iz),phi(iz)
    double precision :: beta(iz),beta1i(iz),beta2i(iz),b(iz),ee(iz) 
    m = box%con%marg
    wid = box%con%wid
    gami = 1./box%con%gam

    amp = 0.05
    lp = 20. 
    phicor = 4.*atan(1.d0) 
    tpt = 25.
    tpho = 1.
    tcor = tpho * tpt
    w = 0.5
    ztr = 8.
    zfsl = -4.
    zfsu = -2.
    zmgc = 8.
    zinv = 5.
    betafs=4.
    betacor=0.2
    a = 2.
    ad = a*(box%con%gam-1.)/box%con%gam

    origin = int(5./box%con%hig*nz)+1+m
    
    forall(i=1:ix) box%x(i)=box%con%dx*(i-m)
    forall(i=1:iz) box%z(i)=box%con%dz*(i-origin)
    
    box%con%gx = 0.
    box%con%gy = 0.
    box%con%gz = -gami 

    do i=origin+1,iz
        temp(i) = tpho + 0.5 * (tcor-tpho)*(tanh((box%z(i)-ztr)/w) + 1.)
    end do
    do i=origin,1,-1
        temp(i) = tpho - ad*box%z(i)
    end do

    den(origin) = 1.
    pre(origin) = gami*temp(origin)
    
    beta1i = 0.25/betafs*(tanh((box%z-zfsl)/w) + 1.) * (-tanh((box%z-zfsu)/w) + 1.)
    beta2i = 0.5/betacor*(tanh((box%z-zmgc)/w) + 1.)
    beta = 1./(beta1i+beta2i)
    do i=origin+1,iz
        den(i) = den(i-1) * ((1.+1./beta(i-1))*temp(i-1) + 0.5*box%con%gam*box%con%dz*box%con%gz)&
                          / ((1.+1./beta(i))*temp(i) - 0.5*box%con%gam*box%con%dz*box%con%gz)
    end do
    do i=origin-1,1,-1
        den(i) = den(i+1) * ((1.+1./beta(i+1))*temp(i+1) - 0.5*box%con%gam*box%con%dz*box%con%gz)&
                          / ((1.+1./beta(i))*temp(i) + 0.5*box%con%gam*box%con%dz*box%con%gz)
    end do
    
    pre = pre(origin) * (den/den(origin)) * (temp/temp(origin))
    b = sqrt(2.*pre/beta)
    
    do i=1,iz
        if (box%z(i)<zinv) then
            phi(i)=0.
        else
            phi(i)=phicor
        end if
    end do

    open(24,file="initial.dat",status="replace")
    do i=1,iz
       write (24,*) box%z(i), den(i), pre(i), temp(i), 1./beta(i), b(i), phi(i)
    end do
    close(24)

    box%ro = spread(den,1,ix)
    box%rovx = 0.
    box%rovy = 0.
    box%rovz = 0.
    do i=1,ix
        do j=1,iz
            x = box%x(i)
            z = box%z(j)
            if (z>zfsl .and. z<zfsu .and. x>0.5*wid-0.25*lp .and. x<0.5*wid+0.25*lp) then
                box%rovz(i,j) = box%ro(i,j)*amp*cos(8.*atan(1.d0)*(x-0.5*wid)/lp)
            end if
        end do
    end do

    box%bx = spread(b*cos(phi),1,ix)
    box%by = 0.
    box%bz = spread(b*sin(phi),1,ix)
    box%pr = spread(pre,1,ix)  
    box%e = 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
            + box%pr/(box%con%gam-1.) &
            + 0.5*(box%bx**2 + box%by**2 + box%bz**2)

    box%bpot(1,1)=0.
    do i=2,ix
        box%bpot(i,1) = box%bpot(i-1,1) &
                        - 0.5*box%con%dx*(box%bz(i,1)+box%bz(i-1,1))
    end do
    do i=1,ix
        do j=2,iz
            box%bpot(i,j) = box%bpot(i,j-1) &
                            + 0.5*box%con%dz*(box%bx(i,j)+box%bx(i,j-1))
        end do
    end do

    !uboundary(1,1:marg) = den(iz-marg+1:iz)
    !uboundary(2,1:marg) = box%rovx(100,iz-marg+1:iz)
    !uboundary(3,1:marg) = box%rovy(100,iz-marg+1:iz)
    !uboundary(4,1:marg) = box%rovz(100,iz-marg+1:iz)
    !uboundary(5,1:marg) = box%bx(100,iz-marg+1:iz)
    !uboundary(6,1:marg) = box%by(100,iz-marg+1:iz)
    !uboundary(7,1:marg) = box%bz(100,iz-marg+1:iz)
    !uboundary(8,1:marg) = box%e(100,iz-marg+1:iz)
    !uboundary(9,1:marg) = pre(iz-marg+1:iz)
    
    uboundary(1,1:marg) = den(iz-marg+1:iz) - den(iz-marg:iz-1)
    uboundary(2,1:marg) = box%rovx(100,iz-marg+1:iz) - box%rovx(100,iz-marg:iz-1)
    uboundary(3,1:marg) = box%rovy(100,iz-marg+1:iz) - box%rovy(100,iz-marg:iz-1)
    uboundary(4,1:marg) = box%rovz(100,iz-marg+1:iz) - box%rovz(100,iz-marg:iz-1)
    uboundary(5,1:marg) = box%bx(100,iz-marg+1:iz) - box%bx(100,iz-marg:iz-1)
    uboundary(6,1:marg) = box%by(100,iz-marg+1:iz) - box%by(100,iz-marg:iz-1)
    uboundary(7,1:marg) = box%bz(100,iz-marg+1:iz) - box%bz(100,iz-marg:iz-1)
    uboundary(8,1:marg) = box%e(100,iz-marg+1:iz) - box%e(100,iz-marg:iz-1)
    uboundary(9,1:marg) = pre(iz-marg+1:iz) - pre(iz-marg:iz-1)

end subroutine
end module

