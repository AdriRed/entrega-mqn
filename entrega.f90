!0.00433 700 0.015 1000000 0.00005 0.3 70000
!0.00433 400 0.015 10000 0.0001 0.8 40000
!0.00433 600 0.015 100000 0.0001 0.4 50000
!0.00433 400 0.020 1000 0.0001 0.5 50000
!0.00433 300 0.020 100 0.0001 0.5 50000 

program practica9
   implicit none
   interface 
      function initialize_wave_function(step, r, r2, width) result(psi)
         implicit none
         double precision, dimension(:), allocatable :: psi
         double precision, intent(in):: step
         double precision, dimension(width), intent(in):: r, r2
         integer, intent(in):: width
      end function
   end interface

   common /grosspita/alpha, alpha2, a0, cvar
   common /simulation/step, time_step, atom_numbers
   common /constants/pi
   double precision :: alpha, alpha2, cvar, pi, step, a0, time_step, interaction, density_param
   double precision, dimension(:), allocatable :: r, r2, prev_psi, next_psi
   integer :: width, i, atom_numbers, iterations

   write(*,*) "a0, width, step, atom_numbers, time_step, alpha, iterations"
   read(*,*) a0, width, step, atom_numbers, time_step, alpha, iterations

   allocate(r(width), r2(width))
   
   pi = 4.0d0*datan(1.0d0)

   !calcular constantes
   alpha2 = alpha*alpha
   cvar = 2.0d0*dsqrt(alpha)**3/dsqrt(sqrt(pi))
   interaction = a0*atom_numbers
   density_param = atom_numbers*a0*a0*a0

   ! precompute spacial fields
   do i = 1, width   
      r(i) = dfloat(i-1)*step
      r2(i) = r(i)*r(i)
   end do

   !inicializar funcion de onda (input_matrix)
   prev_psi = initialize_wave_function(step, r, r2, width)

   open(1, file="wave-function.dat")
   write(1, *) "# x, phi"
   do i = 1, width
      write(1, "(E20.12, E20.12)") r(i), prev_psi(i)
   end do
   close(1)

end program practica9

function initialize_wave_function(step, r, r2, width) result(psi)
   implicit none
   double precision, intent(in):: step
   double precision, dimension(width), intent(in):: r, r2
   integer, intent(in):: width
   double precision, dimension(:), allocatable :: psi
   integer :: i
   common /grosspita/alpha, alpha2, a0, cvar
   double precision :: alpha, alpha2, a0, cvar

   allocate(psi(width))
   do i = 1, width
      psi(i) = cvar*r(i)*dexp(-0.5d0*alpha2*r2(i))
   end do

end function initialize_wave_function

subroutine poisson_equation_algorithm(method_used, length, input_matrix, steps, rho, h, output_matrix)
   implicit none

   interface
      double precision function rho(x1)
         integer, intent(in) :: x1
      end function
      subroutine method_used(width1, rho, h1, input_matrix1, output_matrix1)
         implicit none
         interface
            double precision function rho(x1)
               integer, intent(in) :: x1
            end function
         end interface
         integer, intent(in) :: width1
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(width1) :: input_matrix1
         double precision, intent(out), dimension(width1) :: output_matrix1
      end subroutine method_used
   end interface

   integer, intent(in) :: length, steps
   double precision, intent(in), dimension(length) :: input_matrix
   double precision, intent(in) :: h
   double precision, intent(out), dimension(length) :: output_matrix

   double precision, dimension(length) :: current_matrix
   integer :: i

   ! call method_used(length, rho, h, input_matrix, output_matrix)
   current_matrix = input_matrix
   output_matrix = 0

   call method_used(length, rho, h, current_matrix, output_matrix)

   do i = 1, steps
      current_matrix = output_matrix
      call method_used(length, rho, h, current_matrix, output_matrix)
   end do
end subroutine poisson_equation_algorithm

subroutine jacobi_method(width, rho, h, input_matrix, output_matrix)
   implicit none

   interface
      double precision function rho(x, width1, input_matrix1)
         integer, intent(in) :: x, width1
         double precision, intent(in), dimension(width1) :: input_matrix1
      end function
   end interface

   integer, intent(in) :: width
   double precision, intent(in) :: h
   double precision, intent(in), dimension(width) :: input_matrix
   double precision, intent(out), dimension(width) :: output_matrix

   integer :: i
   ! output_matrix = input_matrix

   do i = 2, width-1
      output_matrix(i) = (input_matrix(i+1) + input_matrix(i-1) - 2*input_matrix(i))/(h*h)
   end do

   output_matrix(width) = (input_matrix(width-1) - 2.0d0*input_matrix(width))/(h*h)

end subroutine jacobi_method


