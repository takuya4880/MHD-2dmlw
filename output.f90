module op
implicit none
contains
subroutine outp(box,t)
    use defstruct
    use cansio
    implicit none
    type(cell) :: box
    double precision :: t

    double precision :: az(ix,iz)
    integer :: i,j
    !double precision :: v2(ix,iz), dt(ix,iz)
    !m = box%con%marg
    
    az(1,1)=0.
    do i=2,ix
        az(i,1) = az(i-1,1) &
                        - 0.5*box%con%dx*(box%bz(i,1)+box%bz(i-1,1))
    end do
    do i=1,ix
        do j=2,iz
            az(i,j) = az(i,j-1) &
                            + 0.5*box%con%dz*(box%bx(i,j)+box%bx(i,j-1))
        end do
    end do

    write(box%op%mf_t) t
    write(box%op%mf_ro) box%ro
    write(box%op%mf_pr) box%pr
    write(box%op%mf_vx) box%rovx/box%ro
    write(box%op%mf_vy) box%rovy/box%ro
    write(box%op%mf_bx) box%bx
    write(box%op%mf_by) box%by
    write(box%op%mf_az) az

    !do j=1,iz
    !    do i=1,ix
    !        write(23, *)  box%x(i), box%z(j), box%ro(i,j), box%pr(i,j), box%rovx(i,j), box%rovz(i,j), box%bpot(i,j)
    !    end do
    !    write(23,*) " "
    !end do
    !write(23,*) " "

    !i=100
    !do j=1,iz
    !    write(23, *)  box%z(j), box%ro(i,j), box%pr(i,j), box%rovx(i,j), box%rovz(i,j), box%bpot(i,j)
    !end do
    !write(23,*) " "
    !write(23,*) " "

end subroutine 

subroutine outpinit(box)
    use defstruct
    use cansio
    implicit none
    type (cell) :: box

    box%op%mf_t=10
    box%op%mf_x=11
    box%op%mf_y=12
    box%op%mf_ro=20
    box%op%mf_pr=21
    box%op%mf_vx=22
    box%op%mf_vy=23
    box%op%mf_bx=25
    box%op%mf_by=26
    box%op%mf_az=28
    
    call dacdefparam(box%op%mf_params,'params.txt')
    call dacdef0s(box%op%mf_t,'t.dac',6)
    call dacdef2s(box%op%mf_ro,'ro.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_pr,'pr.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_vx,'vx.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_vy,'vy.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_bx,'bx.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_by,'by.dac',6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_az,'az.dac',6,box%con%ix,box%con%iz)

    call dacputparami(box%op%mf_params,'ix',box%con%ix)
    call dacputparami(box%op%mf_params,'jx',box%con%iz)
    call dacputparami(box%op%mf_params,'margin',box%con%marg)


    call dacdef1d(box%op%mf_x,'x.dac',6,box%con%ix)
    write(box%op%mf_x) box%x
    call dacdef1d(box%op%mf_y,'y.dac',6,box%con%iz)
    write(box%op%mf_y) box%z
    call dacputparamd(box%op%mf_params,'gm',box%con%gam)

    close(box%op%mf_x)
    close(box%op%mf_y)
     
end subroutine

end module 
