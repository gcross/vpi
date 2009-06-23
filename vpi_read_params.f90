program vpi_read_params

  implicit none

  character, dimension(10) :: pname
  integer :: pval

  call vpi_init_sim_control(pname,pval)
  print *, pname,pval

contains 

subroutine vpi_init_sim_control(pname,pval)
  character, dimension(:) :: pname
  integer :: pval
  
  open(12, file="vpi.init")
  read (12,"(a5,i)") pname, pval
  close(12)

end subroutine vpi_init_sim_control

end program vpi_read_params
