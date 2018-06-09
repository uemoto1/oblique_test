program main
  use fdtd1d1
  implicit none
  integer :: i
  
  call init_fdtd()
  
  call write_ac()
  
  do i = 1, (nt / 1000)
      call update(1000)
      call write_ac()
  end do

  
  stop "bye"
end program main
