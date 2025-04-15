program practica9
   implicit none
   interface 
      function initialize_wave_function(step1, r1, r21, width1) result(psi)
         implicit none
         double precision, dimension(:), allocatable :: psi
         double precision, intent(in):: step1
         double precision, dimension(width1), intent(in):: r1, r21
         integer, intent(in):: width1
      end function

      subroutine compute_wave_equation(length, r1, r21, psi_pre, steps, h, psi_next)
         implicit none
      
         integer, intent(in) :: length, steps
         double precision, intent(in), dimension(length) :: psi_pre, r1, r21
         double precision, intent(in) :: h
         double precision, intent(out), dimension(length) :: psi_next
      end subroutine compute_wave_equation
   end interface

   common /grosspita/alpha, alpha2, a0, cvar
   common /simulation/step, time_step, atom_numbers, interaction
   common /constants/pi
   double precision :: alpha, alpha2, cvar, pi, step, a0, time_step, interaction, density_param
   double precision, dimension(:), allocatable :: r, r2, prev_psi, next_psi
   integer :: width, i, atom_numbers, iterations

   write(*,*) "a0, width, step, atom_numbers, time_step, alpha, iterations"
   read(*,*) a0, width, step, atom_numbers, time_step, alpha, iterations

   allocate(r(width), r2(width), next_psi(width))
   
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

   open(6, file="energy-per-iteration.dat")
   call compute_wave_equation(width, r, r2, prev_psi, iterations, step, next_psi)
   close(6)
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

subroutine compute_wave_equation(length, r, r2, psi_pre, steps, h, psi_next)
   implicit none

   integer, intent(in) :: length, steps
   double precision, intent(in), dimension(length) :: psi_pre, r, r2
   double precision, intent(in) :: h
   double precision, intent(out), dimension(:), allocatable :: psi_next


   interface
      subroutine second_derivative(width, h1, input_matrix1, output_matrix1)
         implicit none
      
         integer, intent(in) :: width
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(width) :: input_matrix1
         double precision, intent(out), dimension(width) :: output_matrix1
      end subroutine second_derivative   
   end interface


   double precision, dimension(length) :: psi_curr, psi_dx2, local_mu
   integer :: i, j

   double precision :: normalization_const, energy
   common /simulation/step, time_step, atom_numbers, interaction
   double precision :: step, time_step, atom_numbers, interaction

   ! call method_used(length, rho, h, input_matrix, output_matrix)
   psi_curr = psi_pre
   psi_dx2 = 0


   do j = 1, steps
      call second_derivative(length, h, psi_curr, psi_dx2)

      normalization_const = 0.0d0
      energy = 0.0d0
      
      ! Cálculo de energía y potencial químico local
      do i = 1, length       
        if (i == 1) then
         local_mu(i) = 0.0d0  ! Condición de contorno en r=0
        else
          ! Contribuciones a la energía
          energy = energy - psi_curr(i)*psi_dx2(i)*0.5d0 + 0.5d0*r2(i)*psi_curr(i)*psi_curr(i) + &
                 0.5d0*interaction*r2(i)*(psi_curr(i)/r(i))**4
          
          ! Potencial químico local
          local_mu(i) = -0.5d0*psi_dx2(i)/psi_curr(i) + 0.5d0*r2(i) + interaction*(psi_curr(i)/r(i))**2
        end if
        
        ! Evolución temporal de la función de onda
        psi_next(i) = psi_curr(i) - time_step*local_mu(i)*psi_curr(i)
        
        ! Acumulación para normalización
        normalization_const = normalization_const + psi_next(i)*psi_next(i)
      end do
      
      ! Normalización y actualización
      normalization_const = sqrt(normalization_const*step)
      energy = energy*step
      
      ! Salida de información periódica
      if (mod(j,200) == 0) then
        write(6,*) 'Iteration:', j, ' Energy:', energy
      end if
      
      ! Normalización y actualización de la función de onda
      psi_curr = psi_next/normalization_const


   end do
end subroutine compute_wave_equation



subroutine second_derivative(width, h, input_matrix, output_matrix)
   implicit none

   integer, intent(in) :: width
   double precision, intent(in) :: h
   double precision, intent(in), dimension(width) :: input_matrix
   double precision, intent(out), dimension(width) :: output_matrix

   integer :: i
   ! output_matrix = input_matrix

   do i = 2, width-1
      output_matrix(i) = (input_matrix(i+1) + input_matrix(i-1) - 2d0*input_matrix(i))/(h*h)
   end do

   output_matrix(width) = (input_matrix(width-1) - 2.0d0*input_matrix(width))/(h*h)

end subroutine second_derivative


