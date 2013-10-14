module pr
implicit none

contains

subroutine pressure(box)
    use defstruct
    implicit none
    type(cell) :: box

    box%pr = (box%con%gam-1.) &
             * (box%e - 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
                      - 0.5*(box%bx**2 + box%by**2 + box%bz**2) )

    box%pr = 0.5*(box%pr + abs(box%pr))

end subroutine

end module
