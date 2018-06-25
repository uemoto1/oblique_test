program main
  use fdtd1d1
  use inputoutput
  implicit none
  integer :: i
  
  call init_fdtd()
  
  call run_fdtd()
  stop "bye"
end program main
