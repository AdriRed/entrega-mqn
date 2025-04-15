!-----------------------------------------------------------------------
! PROGRAM: BOSONIC_TRAP_SIMULATION
! Descripción: Simulación numérica de bosones interactuantes en una trampa armónica esférica
!              utilizando la aproximación de Gross-Pitaevskii.
! Métodos: Diferencias finitas, integración en rejilla radial, iteración temporal
!-----------------------------------------------------------------------

module constants
    implicit none
    ! Constantes matemáticas y físicas
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)       ! Pi
    real(8), parameter :: pi_inv = 1.0d0/(4.0d0*pi)   ! 1/(4*pi)
    real(8), parameter :: pi2_inv = sqrt(pi_inv)      ! sqrt(1/(4*pi))
  end module constants
  
  module variables
    use constants
    implicit none
    ! Parámetros de simulación
    real(8) :: a0          ! Longitud de dispersión en unidades del oscilador armónico
    integer :: n1          ! Número de puntos en la rejilla radial
    real(8) :: step        ! Paso de integración en espacio r
    real(8) :: aa          ! Número de átomos
    real(8) :: time_step   ! Paso temporal
    real(8) :: alpha       ! Parámetro inicial para la función del oscilador armónico
    integer :: iter        ! Número de iteraciones
    
    ! Arrays para almacenar valores en la rejilla
    real(8), dimension(:), allocatable :: xr       ! Coordenada radial
    real(8), dimension(:), allocatable :: frev      ! Función de onda en iteración anterior
    real(8), dimension(:), allocatable :: freo      ! Función de onda actual
    real(8), dimension(:), allocatable :: fred      ! Derivada segunda de la función de onda
    real(8), dimension(:), allocatable :: xmu       ! Potencial químico local
    real(8), dimension(:), allocatable :: fren      ! Nueva función de onda
    real(8), dimension(:), allocatable :: den       ! Densidad de partículas
    real(8), dimension(:), allocatable :: u         ! Potencial de partícula única
    
    ! Variables físicas calculadas
    real(8) :: ene0        ! Energía total
    real(8) :: radious     ! Radio cuadrático medio
    real(8) :: radious2    ! Radio al cuadrado
    real(8) :: xkin        ! Energía cinética
    real(8) :: poth0       ! Energía potencial armónica
    real(8) :: potself     ! Energía de interacción
    real(8) :: chem        ! Potencial químico promedio
    real(8) :: xaver       ! Valor promedio adicional
    real(8) :: xnormden    ! Normalización de la densidad
  end module variables
  
  !-----------------------------------------------------------------------
  ! PROGRAMA PRINCIPAL
  !-----------------------------------------------------------------------
  program bosonic_trap
    use variables
    implicit none
    integer :: i, it, itw
    real(8) :: xr2, xnorm, alpha2, cvar, cequ, as3n
    
    ! Inicialización: leer parámetros de entrada
    call read_input_parameters()
    
    ! Asignación de memoria para los arrays
    allocate(xr(n1), frev(n1), freo(n1), fred(n1), xmu(n1), fren(n1), den(n1), u(n1))
    
    ! Constantes derivadas
    alpha2 = alpha*alpha
    cvar = 2.0d0*sqrt(alpha)**3/sqrt(sqrt(pi))
    
    ! Construcción de la función de onda inicial (oscilador armónico)
    do i = 1, n1
      xr(i) = step*dble(i-1)
      xr2 = xr(i)*xr(i)
      frev(i) = cvar*xr(i)*exp(-0.5d0*alpha2*xr2)
      freo(i) = frev(i)
    end do
    
    ! Configuración de parámetros de interacción
    cequ = a0*aa           ! Término de interacción
    as3n = aa*a0*a0*a0     ! Parámetro de densidad
    
    ! Proceso iterativo principal
    itw = 0
    do it = 1, iter
      itw = itw + 1
      
      ! Cálculo de la derivada segunda (laplaciano)
      call calculate_second_derivative(freo, fred, step, n1)
      
      ! Inicialización de normas y energías
      xnorm = 0.0d0
      ene0 = 0.0d0
      
      ! Cálculo de energía y potencial químico local
      do i = 1, n1
        xr2 = xr(i)*xr(i)
        
        if (i == 1) then
          xmu(i) = 0.0d0  ! Condición de contorno en r=0
        else
          ! Contribuciones a la energía
          ene0 = ene0 - freo(i)*fred(i)*0.5d0 + 0.5d0*xr2*freo(i)*freo(i) + &
                 0.5d0*cequ*xr2*(freo(i)/xr(i))**4
          
          ! Potencial químico local
          xmu(i) = -0.5d0*fred(i)/freo(i) + 0.5d0*xr2 + cequ*(freo(i)/xr(i))**2
        end if
        
        ! Evolución temporal de la función de onda
        fren(i) = freo(i) - time_step*xmu(i)*freo(i)
        
        ! Acumulación para normalización
        xnorm = xnorm + fren(i)*fren(i)
      end do
      
      ! Normalización y actualización
      xnorm = sqrt(xnorm*step)
      ene0 = ene0*step
      
      ! Salida de información periódica
      if (mod(it,200) == 0) then
        write(6,*) 'Iteration:', it, ' Energy:', ene0
        itw = 0
      end if
      
      ! Normalización y actualización de la función de onda
      freo = fren/xnorm
      
      ! Salida final si es la última iteración
      if (it == iter) then
        open(unit=10, file='potential.dat')
        do i = 2, n1
          write(10,'(2e15.5)') xr(i), xmu(i)
        end do
        close(10)
      end if
    end do
    
    ! Cálculo de propiedades físicas finales
    call calculate_final_properties()
    
    ! Salida de resultados
    call output_results()
    
    ! Liberar memoria
    deallocate(xr, frev, freo, fred, xmu, fren, den, u)
    
  contains
  
    !---------------------------------------------------------------------
    ! Subrutina para leer parámetros de entrada
    !---------------------------------------------------------------------
    subroutine read_input_parameters()
      implicit none
      read(5,*) a0, n1, step, aa, time_step, alpha, iter
      
      ! Verificación de parámetros
      if (n1 > 1000) then
        write(6,*) 'Error: n1 too large. Maximum is 1000'
        stop
      end if
      
      ! Información de la simulación
      write(6,'(/5X,"* ",F10.0," BOSONS IN A SPHERICAL TRAP *"/ &
           & 5X,"* r-GRID IN N1=",I4," POINTS", " r-step",f8.3/ &
           & 5X,"* A0=",E12.4,2X," ALPHA=",E12.4/ &
           & 5x,"* time=",e12.4,"  number-iter =",i8/)') &
           aa, n1, step, a0, alpha, time_step, iter
    end subroutine read_input_parameters
  
    !---------------------------------------------------------------------
    ! Calcula la derivada segunda usando diferencias finitas
    !---------------------------------------------------------------------
    subroutine calculate_second_derivative(f, df, h, n)
      implicit none
      real(8), intent(in) :: f(:)
      real(8), intent(out) :: df(:)
      real(8), intent(in) :: h
      integer, intent(in) :: n
      integer :: i
      
      ! Puntos interiores
      do i = 2, n-1
        df(i) = (f(i-1) + f(i+1) - 2.0d0*f(i))/(h*h)
      end do
      
      ! Condición de contorno en el extremo
      df(n) = (f(n-1) - 2.0d0*f(n))/(h*h)
    end subroutine calculate_second_derivative
  
    !---------------------------------------------------------------------
    ! Calcula propiedades físicas finales
    !---------------------------------------------------------------------
    subroutine calculate_final_properties()
      implicit none
      integer :: i
      real(8) :: xr2
      
      ! Recalcular derivada segunda para propiedades finales
      call calculate_second_derivative(freo, fred, step, n1)
      
      ! Inicialización de acumuladores
      radious = 0.0d0
      xkin = 0.0d0
      poth0 = 0.0d0
      potself = 0.0d0
      chem = 0.0d0
      xaver = 0.0d0
      xnormden = 0.0d0
      
      ! Cálculo de integrales para propiedades físicas
      do i = 2, n1
        xr2 = xr(i)*xr(i)
        
        ! Contribuciones a las diferentes propiedades
        radious = radious + xr2*freo(i)*freo(i)
        xkin = xkin + freo(i)*fred(i)
        poth0 = poth0 + xr2*freo(i)*freo(i)
        potself = potself + xr2*(freo(i)/xr(i))**4
        chem = chem + xmu(i)*freo(i)*freo(i)
        
        ! Densidad y potencial de partícula única
        u(i) = 0.5d0*xr2 + cequ*(freo(i)/xr(i))**2
        den(i) = (freo(i)/xr(i))**2
        
        ! Normalizaciones
        xnormden = xnormden + den(i)*xr2
        xaver = xaver + freo(i)*freo(i)*as3n*den(i)
      end do
      
      ! Normalización de resultados
      radious2 = radious*step
      radious = sqrt(radious*step)
      xaver = xaver*step
      chem = chem*step
      xkin = -xkin*step*0.5d0
      poth0 = 0.5d0*poth0*step
      potself = potself*step*cequ*0.5d0
      ene0 = xkin + poth0 + potself  ! Energía total
    end subroutine calculate_final_properties
  
    !---------------------------------------------------------------------
    ! Escribe los resultados finales
    !---------------------------------------------------------------------
    subroutine output_results()
      implicit none
      integer :: i
      
      ! Archivo con la densidad
      open(unit=9, file='density.dat')
      do i = 2, n1
        write(9,'(2e15.5)') xr(i), den(i)
      end do
      close(9)
      
      ! Resumen de resultados
      write(6,'(/5X,"* ener",e12.5,"   average chemical=", e12.5/ &
           & 5X,"* kin-ener=",e12.5,"     total-pot= ", e12.5,/ &
           & 5X,"* potho=",E12.4,2X,"     potint =",E12.4/ &
           & 5x,"* radious =",e12.4, 5x, " radious2 =",e15.7, /)') &
           ene0, chem, xkin, (poth0+potself), poth0, potself, radious, radious2
      
      write(6,*) ' Density normalization = ', xnormden
    end subroutine output_results
  
  end program bosonic_trap
