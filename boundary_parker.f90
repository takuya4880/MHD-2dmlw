subroutine boundary(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box
    double precision :: uboundary(9,2)

    call cppmbc(box%ro, 1, uboundary)
    call cppmbc(box%rovx, 2, uboundary)
    call cppmbc2(box%rovy, 3, uboundary)
    call cppmbc(box%rovz, 4, uboundary)
    call cppmbc(box%bx, 5, uboundary)
    call cppmbc2(box%by, 6, uboundary)
    call cppmbc(box%bz, 7, uboundary)
    call cppmbc(box%e, 8, uboundary)
    call cppmbc(box%pr, 9, uboundary)

end subroutine 
!contains
subroutine allfreebc(arr)
    use defstruct
    implicit none
    double precision ::  arr(ix,iy)

    integer i

    do i=1,marg
        arr(i,:) = arr(marg+1,:)   
        arr(ix-marg+i,:) = arr(ix-marg,:)
        arr(:,i) = arr(:,marg+1) 
        arr(:,iy-marg+i) = arr(:,iy-marg)
    end do
end subroutine

subroutine fppmbc(arr)  !free, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iy)
    integer i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = arr(:,2*marg+1-i)   
        arr(:,iy-marg+i) = arr(:,iy-marg)
    end do
end subroutine

subroutine fppmbc2(arr)
    use defstruct
    implicit none
    double precision :: arr(ix,iy)
    integer i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = - arr(:,2*marg+1-i)   
        arr(:,iy-marg+i) = arr(:,iy-marg)
    end do
end subroutine

subroutine cppmbc(arr, k, ub)   !constant, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iy)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = arr(:,2*marg+1-i)   
        !arr(:,iy-marg+i) = arr(:,iy-marg)
        arr(:,iy-marg+i) = ub(k,i)
    end do
end subroutine

subroutine cppmbc2(arr, k, ub)   !constant, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iy)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = - arr(:,2*marg+1-i)   
        !arr(:,iy-marg+i) = arr(:,iy-marg)
        arr(:,iy-marg+i) = ub(k,i)
    end do
end subroutine
!end subroutine
