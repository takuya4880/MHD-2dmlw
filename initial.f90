subroutine initial(box)
    use defstruct
    implicit none
    type(cell) :: box

    integer :: i
    double precision :: b      !mag field
    double precision :: gami               !inberse of gamma
    double precision :: alp = 0.5
    double precision :: ome2i = 2500/4.    !inverse of ome^2
    double precision :: gauss = 0.1     !amplitude of gaussian
    double precision :: theta = 3.141592/4.     !inclination of B
    type(constants) :: con
    con = box%con

    forall(i=1:ix) box%x(i)=con%dx*(i-con%marg-0.5)
    forall(i=1:iy) box%y(i)=con%dy*(i-con%marg-0.5)
    
    box%con%gx = 0.
    box%con%gx = 0.
    box%con%gy = -1./box%con%gam
    box%con%gz = 0.

    gami = 1./con%gam
    b=sqrt(alp*gami)

    box%ro = 1.
    box%rovx = 0.
    box%rovy = 0.
    box%rovz = 0.
    !box%bx = b*cos(theta)
    !box%by = b*sin(theta)
    box%bx = b
    box%by = 0.
    box%bz = 0.
    
    box%pr = (1. + gauss*exp(-( (spread(box%x,2,iy)-con%wid*0.5)**2 &
                    + (spread(box%y,1,ix)-con%hig*0.5)**2 )*ome2i)) * gami
    box%e = box%pr/(con%gam-1.) + 0.5*b*b

end subroutine
