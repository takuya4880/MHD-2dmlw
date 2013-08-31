subroutine boundary(box)
    use defstruct
    implicit none
    type(cell) :: box

    call fppmbc(box%ro)
    call fppmbc(box%rovx)
    call fppmbc2(box%rovy)
    call fppmbc(box%rovz)
    call fppmbc(box%bx)
    call fppmbc2(box%by)
    call fppmbc(box%bz)
    call fppmbc(box%e)
    call fppmbc(box%pr)

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

subroutine fppmbc(arr)
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
!end subroutine
