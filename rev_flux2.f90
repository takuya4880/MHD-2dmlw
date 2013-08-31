subroutine flux(box, fx, fy)
    use defstruct
    implicit none
    type(cell) :: box, fx, fy
    
    !double precision, allocatable :: b2(:,:),roi(:,:), roh(:,:)
    !allocate(b2(ix,iy))
    !allocate(roi(ix,iy))
    !allocate(roh(ix,iy))

    double precision b2,roi,roh
    integer i,j

    do j=1,iy
        do i=1,ix

    b2 = box%bx(i,j)**2 + box%by(i,j)**2 + box%bz(i,j)**2
    roi = 1./box%ro(i,j)
    roh = 0.5*(box%rovx(i,j)**2+box%rovy(i,j)**2+box%rovz(i,j)**2)*roi &
            + box%pr(i,j)*box%con%gam/(box%con%gam-1.) + b2

    fx%ro(i,j) = box%rovx(i,j)
    fx%rovx(i,j) = box%rovx(i,j)*box%rovx(i,j)*roi - box%bx(i,j)*box%bx(i,j) + box%pr(i,j) + 0.5*b2
    fx%rovy(i,j) = box%rovx(i,j)*box%rovy(i,j)*roi - box%bx(i,j)*box%by(i,j)
    fx%rovz(i,j) = box%rovx(i,j)*box%rovz(i,j)*roi - box%bx(i,j)*box%bz(i,j)
    fx%by(i,j) = (box%rovx(i,j)*box%by(i,j) - box%rovy(i,j)*box%bx(i,j)) * roi
    fx%bz(i,j) = (box%rovx(i,j)*box%bz(i,j) - box%rovz(i,j)*box%bx(i,j)) * roi
    fx%e(i,j) = (roh*box%rovx(i,j) - box%bx(i,j)*(box%rovx(i,j)*box%bx(i,j) + box%rovy(i,j)*box%by(i,j) &
            + box%rovz(i,j)*box%bz(i,j)) )*roi 

    fy%ro(i,j) = box%rovy(i,j)
    fy%rovx(i,j) = box%rovy(i,j)*box%rovx(i,j)*roi - box%by(i,j)*box%bx(i,j)
    fy%rovy(i,j) = box%rovy(i,j)*box%rovy(i,j)*roi - box%by(i,j)*box%by(i,j) + box%pr(i,j) + 0.5*b2
    fy%rovz(i,j) = box%rovy(i,j)*box%rovz(i,j)*roi - box%by(i,j)*box%bz(i,j)
    fy%bx(i,j) = (box%rovy(i,j)*box%bx(i,j) - box%rovx(i,j)*box%by(i,j)) * roi
    fy%bz(i,j) = (box%rovy(i,j)*box%bz(i,j) - box%rovz(i,j)*box%by(i,j)) * roi 
    fy%e(i,j) = (roh*box%rovy(i,j) - box%by(i,j)*(box%rovx(i,j)*box%bx(i,j) + box%rovy(i,j)*box%by(i,j) &
            + box%rovz(i,j)*box%bz(i,j)) )*roi 
    
        end do
    end do
    fx%bx = 0.
    fy%by = 0.

    !deallocate(b2,roi,roh)

end subroutine
