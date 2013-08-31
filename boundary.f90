subroutine boundary(box)
    use defstruct
    implicit none
    type(cell) :: box

    call allfreebc(box%ro)
    call allfreebc(box%rovx)
    call allfreebc(box%rovy)
    call allfreebc(box%rovz)
    call allfreebc(box%bx)
    call allfreebc(box%by)
    call allfreebc(box%bz)
    call allfreebc(box%e)
    call allfreebc(box%pr)


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

!end subroutine
