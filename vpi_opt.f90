module vpi_opt
  use vpi_defines
  use vpi_aziz
  implicit none

contains
   SUBROUTINE madr(n,p,x,nf,r,nl2sol_iteration,been_rst,stop_opt,x_h)
    USE machine_constants
    INTEGER,INTENT(in)    :: n,p,nl2sol_iteration,been_rst
    INTEGER,INTENT(inout) :: nf
    REAL(kind=kind(0.d0)),INTENT(in)          :: x(:)
    REAL(kind=kind(0.d0)),INTENT(in),OPTIONAL :: x_h(:)
    REAL(kind=kind(0.d0)),INTENT(out)         :: r(:)
    LOGICAL,INTENT(inout)                     :: stop_opt
   END SUBROUTINE madr
end module vpi_opt
