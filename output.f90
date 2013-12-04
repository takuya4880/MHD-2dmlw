module op
implicit none
contains
subroutine outp(box,t)
    use defstruct
    use cansio
    implicit none
    type(cell) :: box
    double precision :: t

    !integer :: i,j,m
    !double precision :: v2(ix,iz), dt(ix,iz)
    !m = box%con%marg
    
    open(box%op%mf_t,file="t.dac",form="unformatted",position="append")

    write(box%op%mf_t) t
    write(box%op%mf_ro) box%ro
    write(box%op%mf_pr) box%pr
    write(box%op%mf_vx) box%rovx/box%ro
    write(box%op%mf_vy) box%rovz/box%ro
    write(box%op%mf_bx) box%bx
    write(box%op%mf_by) box%by
    write(box%op%mf_az) box%bpot

    close(box%op%mf_t)


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
    close(box%op%mf_params)
     
end subroutine

end module 
