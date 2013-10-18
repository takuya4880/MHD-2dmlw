module bc
implicit none 
contains
subroutine boundary(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box
    double precision :: uboundary(9,marg)

    call gppmbc(box%ro, 1, uboundary)
    call gppmbc(box%rovx, 2, uboundary)
    call gppmbc(box%rovy, 3, uboundary)
    call gppmbc2(box%rovz, 4, uboundary)
    call gppmbc(box%bx, 5, uboundary)
    call gppmbc(box%by, 6, uboundary)
    call gppmbc2(box%bz, 7, uboundary)
    call gppmbc(box%e, 8, uboundary)
    call gppmbc(box%pr, 9, uboundary)

    !call cppmbc(box%ro, 1, uboundary)
    !call cppmbc(box%rovx, 2, uboundary)
    !call cppmbc(box%rovy, 3, uboundary)
    !call cppmbc2(box%rovz, 4, uboundary)
    !call cppmbc(box%bx, 5, uboundary)
    !call cppmbc(box%by, 6, uboundary)
    !call cppmbc2(box%bz, 7, uboundary)
    !call cppmbc(box%e, 8, uboundary)
    !call cppmbc(box%pr, 9, uboundary)

    !call fppmbc(box%ro)
    !call fppmbc(box%rovx)
    !call fppmbc(box%rovy)
    !call fppmbc2(box%rovz)
    !call fppmbc(box%bx)
    !call fppmbc(box%by)
    !call fppmbc2(box%bz)
    !call fppmbc(box%e)
    !call fppmbc(box%pr)
    !call fppmbc(box%bpot)


end subroutine 
!contains
subroutine allfreebc(arr)
    use defstruct
    implicit none
    double precision ::  arr(ix,iz)

    integer i

    do i=1,marg
        arr(i,:) = arr(marg+1,:)   
        arr(ix-marg+i,:) = arr(ix-marg,:)
        arr(:,i) = arr(:,marg+1) 
        arr(:,iz-marg+i) = arr(:,iz-marg)
    end do
end subroutine

subroutine fppmbc(arr)  !free, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = arr(:,2*marg+1-i)   
        arr(:,iz-marg+i) = arr(:,iz-marg)
    end do
end subroutine

subroutine fppmbc2(arr)
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = - arr(:,2*marg+1-i)   
        arr(:,iz-marg+i) = arr(:,iz-marg)
    end do
end subroutine

subroutine cppmbc(arr, k, ub)   !constant, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = arr(:,2*marg+1-i)   
        !arr(:,iz-marg+i) = arr(:,iz-marg)
        arr(:,iz-marg+i) = ub(k,i)
    end do
end subroutine

subroutine cppmbc2(arr, k, ub)   !constant, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = - arr(:,2*marg+1-i)   
        !arr(:,iz-marg+i) = arr(:,iz-marg)
        arr(:,iz-marg+i) = ub(k,i)
    end do
end subroutine

subroutine gppmbc(arr, k, ub)   !gradient, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = arr(:,2*marg+1-i)   
        !arr(:,iz-marg+i) = arr(:,iz-marg)
        arr(:,iz-marg+i) = arr(:,iz-marg) + ub(k,i)
    end do
end subroutine

subroutine gppmbc2(arr, k, ub)   !gradient, periodic, mirror
    use defstruct
    implicit none
    double precision :: arr(ix,iz)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(i,:) = arr(ix-2*marg+i,:) 
        arr(ix-marg+i,:) = arr(marg+i,:)
        arr(:,i) = - arr(:,2*marg+1-i)   
        !arr(:,iz-marg+i) = arr(:,iz-marg)
        arr(:,iz-marg+i) = arr(:,iz-marg) + ub(k,i)
    end do
end subroutine


end module
