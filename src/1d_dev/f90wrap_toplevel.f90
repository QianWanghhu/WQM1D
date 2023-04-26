subroutine f90wrap_wq1d(istl_, is2tl_)
    implicit none
    external wq1d
    
    integer :: istl_
    integer :: is2tl_
    call wq1d(istl_, is2tl_)
end subroutine f90wrap_wq1d

