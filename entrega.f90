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

      subroutine compute_wave_equation(length, r1, r21, psi_pre, steps, h, psi_next, local_mu)
         implicit none
      
         integer, intent(in) :: length, steps
         double precision, intent(in), dimension(length) :: psi_pre, r1, r21
         double precision, intent(in) :: h
         double precision, intent(out), dimension(length) :: psi_next, local_mu
      end subroutine compute_wave_equation

      subroutine calculate_final_properties(length, r1, final_psi, local_mu)
         implicit none
      
         integer, intent(in) :: length
         double precision, intent(in), dimension(length) :: final_psi, r1, local_mu
      end subroutine calculate_final_properties
   end interface
   
   common /grosspita/alpha, alpha2, a0, cvar
   common /simulation/step, time_step, atom_numbers, interaction
   common /constants/pi
   double precision :: alpha, alpha2, cvar, pi, step, a0, time_step, interaction, density_param
   double precision, dimension(:), allocatable :: r, r2, prev_psi, next_psi, local_mu
   integer :: width, i, atom_numbers, iterations

   write(*,*) "a0, width, step, atom_numbers, time_step, alpha, iterations"
   read(*,*) a0, width, step, atom_numbers, time_step, alpha, iterations

   allocate(r(width), r2(width), next_psi(width), local_mu(width))
   
   pi = 4.0d0*datan(1.0d0)

   !calcular constantes
   alpha2 = alpha*alpha
   cvar = 2.0d0*dsqrt(alpha)**3/dsqrt(sqrt(pi))
   ! interaction = 0
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
   do i = 1, width
      write(1, "(E20.12, E20.12)") r(i), prev_psi(i)
   end do
   close(1)

   open(6, file="energy-per-iteration.dat")
   call compute_wave_equation(width, r, r2, prev_psi, iterations, step, next_psi, local_mu)
   close(6)


   call calculate_final_properties(width, r, next_psi, local_mu)
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

subroutine compute_wave_equation(length, r, r2, psi_pre, steps, h, psi_next, local_mu)
   implicit none

   integer, intent(in) :: length, steps
   double precision, intent(in), dimension(length) :: psi_pre, r, r2
   double precision, intent(in) :: h
   double precision, intent(out), dimension(length) :: psi_next, local_mu

   interface
      subroutine second_derivative(width, h1, input_matrix1, output_matrix1)
         implicit none
      
         integer, intent(in) :: width
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(width) :: input_matrix1
         double precision, intent(out), dimension(width) :: output_matrix1
      end subroutine second_derivative   
   end interface


   double precision, dimension(length) :: psi_curr, psi_dx2
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
      ! if (mod(j,200) == 0) then
      write(6,*) j, energy
      ! end if
      
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


subroutine calculate_final_properties(length, r, final_psi, local_mu)
   implicit none

   integer, intent(in) :: length
   double precision, intent(in), dimension(length) :: final_psi, r, local_mu
   
   interface
      subroutine second_derivative(width, h1, input_matrix1, output_matrix1)
         implicit none
      
         integer, intent(in) :: width
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(width) :: input_matrix1
         double precision, intent(out), dimension(width) :: output_matrix1
      end subroutine second_derivative   
   end interface


   common /grosspita/alpha, alpha2, a0, cvar
   common /simulation/step, time_step, atom_numbers, interaction
   double precision :: alpha, alpha2, cvar, a0, time_step, interaction, density_param, step

   integer :: i, atom_numbers
   double precision :: r2, total_energy, accumulted_radius
   double precision :: mean_squared_radius, kinetic_energy, armonic_potential, interaction_energy
   double precision :: average_chem_potential, squared_radius, density_normalization
   double precision, dimension(length) :: psi_dxd2, density, particle_potential
   
   write(*, *) "Calculating final property"

   ! Recalcular derivada segunda para propiedades finales
   call second_derivative(length, step, final_psi, psi_dxd2)
   
   ! Inicialización de acumuladores
   accumulted_radius = 0.0d0
   kinetic_energy = 0.0d0
   armonic_potential = 0.0d0
   interaction_energy = 0.0d0
   average_chem_potential = 0.0d0
   density_normalization = 0.0d0
   
   ! Cálculo de integrales para propiedades físicas
   do i = 2, length
     r2 = r(i)*r(i)
     
     ! Contribuciones a las diferentes propiedades
     accumulted_radius = accumulted_radius + r2*final_psi(i)*final_psi(i)
     kinetic_energy = kinetic_energy + final_psi(i)*psi_dxd2(i)
     armonic_potential = armonic_potential + r2*final_psi(i)*final_psi(i)
     interaction_energy = interaction_energy + r2*(final_psi(i)/r(i))**4
     average_chem_potential = average_chem_potential + local_mu(i)*final_psi(i)*final_psi(i)
     
     ! Densidad y potencial de partícula única
     particle_potential(i) = 0.5d0*r2 + interaction*(final_psi(i)/r(i))**2
     density(i) = (final_psi(i)/r(i))**2
     
     ! Normalizaciones
     density_normalization = density_normalization + density(i)*r2
   end do
   
   ! Normalización de resultados
   squared_radius = accumulted_radius*step
   mean_squared_radius = sqrt(accumulted_radius*step)
   average_chem_potential = average_chem_potential*step
   kinetic_energy = -kinetic_energy*step*0.5d0
   armonic_potential = 0.5d0*armonic_potential*step
   interaction_energy = interaction_energy*step*interaction*0.5d0
   total_energy = kinetic_energy + armonic_potential + interaction_energy  ! Energía total
   write(*,'(/5X,"* ener",e12.5,"   average chemical=", e12.5/ &
   & 5X,"* kin-ener=",e12.5,"     total-pot= ", e12.5,/ &
   & 5X,"* potho=",E12.4,2X,"     potint =",E12.4/ &
   & 5x,"* radious =",e12.4, 5x, " radious2 =",e15.7, /)') &
   total_energy, average_chem_potential, kinetic_energy, (armonic_potential+interaction_energy), &
      armonic_potential, interaction_energy, mean_squared_radius, squared_radius

end subroutine calculate_final_properties