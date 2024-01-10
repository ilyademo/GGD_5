Program Main
  use Functions
  implicit none
  include 'omp_lib.h'
  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt', SolverParameters='solvedata.txt' ! names of input and output files
  character MeshFile*30        ! name of file with computational mesh
  character GradType*30
  character gradT_calc*30
  integer NI,NJ,i,j, Niter_GG, mode, iter, Niter, scheme
  integer, parameter:: IO = 12 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume, T, Res, dt, T_flos! scalar arrays
  real,allocatable,dimension(:,:):: divV, divVExact, divVError, LapP, LapPExact, LapPError, rotV, rotVExact, rotVError, T_error !Errors arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,GradP, GradPExact, GradPError, V, Vel_eq, GradT ! vector arrays
  real Re, Pr, CFL, dx, T_cold, T_hot, T_0, eps, Max_Res, old_res, T1, T2
!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  read(IO,*) GradType  ! Метод расчета градиента (при расчете пространственных операторов)
  read(IO,*) mode      ! Выбор схемы аппроксимации скалярной величины на грань (при расчете дивергенции)
  NIter_GG = 5         ! Число итераций Грина-Гаусса
  CLOSE(IO)
!=== Read Solver data ===
  write(*,*) 'Read solver data file: ', SolverParameters
  OPEN(IO,FILE=SolverParameters)
  read(IO,*) CFL
  read(IO,*) T_0
  read(IO,*) Re
  read(IO,*) Pr
  read(IO,*) Niter
  read(IO,*) T_cold
  read(IO,*) T_hot
  read(IO,*) eps
  read(IO,*) scheme
  read(IO,*) gradT_calc
  CLOSE(IO)
  
!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(GradP(0:NI, 0:NJ, 2)) !Gradient of pressure
  allocate(GradPExact(0:NI, 0:NJ, 2))
  allocate(GradPError(0:NI, 0:NJ, 2))
  allocate(V(0:NI,0:NJ,2))
  allocate(Vel_eq(0:NI,0:NJ,2))
  allocate(divV(0:NI,0:NJ))
  allocate(divVExact(0:NI,0:NJ))
  allocate(divVError(0:NI,0:NJ))
  allocate(rotV(0:NI,0:NJ))
  allocate(rotVExact(0:NI,0:NJ))
  allocate(rotVError(0:NI,0:NJ))
  allocate(LapP(0:NI,0:NJ))
  allocate(LapPExact(0:NI,0:NJ))
  allocate(LapPError(0:NI,0:NJ))
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces
  allocate(T(0:NI,0:NJ))   ! T
  allocate(T_flos(0:NI,0:NJ))
  allocate(T_error(0:NI,0:NJ))
  allocate(GradT(0:NI, 0:NJ, 2))
  allocate(Res(0:NI,0:NJ))   ! Res
  allocate(dt(0:NI,0:NJ))
  
!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)

  dx = 1
!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector, dx) 
  !write (*,*) "dx = ", dx
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      P(I,J) = Pressure(CellCenter(i,j,1),CellCenter(i,j,2))
      Call CalcGradExact(CellCenter(i,j,1), CellCenter(i,j,2), GradPExact(i,j,:))
      Call Velocity(CellCenter(i,j,1), CellCenter(i,j,2), V(i,j,:))
      divVExact(i,j) = divVelocityPExact(CellCenter(i,j,1), CellCenter(i,j,2))
      LapPExact(i,j) = LapPressureExact(CellCenter(i,j,1), CellCenter(i,j,2))
      rotVExact(i,j) = RotVelExact(CellCenter(i,j,1), CellCenter(i,j,2))
    ENDDO
  ENDDO
!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'
  if (GradType .EQ. "GG") then
  GradP = 0.0
  do i = 1, NIter_GG
  Call B_CalcGradient(GradP, NI, NJ, P, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  GradPError = Abs(GradPExact-GradP)/GradPExact
  write(*,*) 'Iteration: ', i, 'maximum error ', maxval(GradPError(1:NI-1, 1:NJ-1, :))
  end do
  !Метод наименьших квадратов
  else if (GradType .EQ. "LS") then
  GradP = 0.0
  Call LeastSquareGradCalc(GradP, NI, NJ, P, CellCenter)
  GradPError = Abs(GradPExact-GradP)/GradPExact
  write(*,*) 'Maximum error LeastSquareCalc:', maxval(GradPError(1:NI-1, 1:NJ-1, :))
  end if
  write(*,*) 'Maximum Gradx-error', maxval(GradPError(1:NI-1, 1:NJ-1, 1))
  write(*,*) 'Maximum Grady-error', maxval(GradPError(1:NI-1, 1:NJ-1, 2))

  !=== CALCULATE DIVERGENCE ===
  write(*,*) "Calculate Divergence"
  divV = 0.0
  call B_CalcDivergence(mode, NI, NJ, V, divV, P, GradP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  divVError = abs(divVExact - divV)/divVExact
  write (*,*) 'Maximum error of divV', maxval(divVError(1:NI-1, 1:NJ-1))
  
  !=== CALCULATE LAPLACIAN ====
  write (*,*) "Calculate Laplacian"
  LapP = 0.0
  Call B_CalcLaplacian(NI, NJ, P, GradP, LapP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  LapPError = abs(LapPExact - LapP)/LapPExact
  write (*,*) 'Maximum error of Laplacian', maxval(LapPError(1:NI-1, 1:NJ-1))
  
  !=== CALCULATE ROTOR ====
  write (*,*) "Calculate Rotor"
  rotV = 0.0
  Call BCalcRotor(mode, NI, NJ, V, rotV, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  rotVError = abs(rotVExact - rotV)/rotVExact
  write (*,*) 'Maximum error of rotV', maxval(rotVError(1:NI-1, 1:NJ-1))
  
  
  !Чтение файлов FLOS
  if (MeshFile == 'caverna1.msh') then
  !Vel_eq = 0
  OPEN(IO,FILE='v_field.dat')
  READ(IO,*) ((Vel_eq(I,J,1),Vel_eq(I,J,2),I=0,NI),J=0,NJ)
  CLOSE(IO)
  if (scheme == 1) then
  OPEN(IO,FILE='t_field.dat')
  else if (scheme == 2) then
  OPEN(IO,FILE='t_field_SOU.dat')
  end if
  READ(IO,*) ((T_flos(I,J),I=0,NI),J=0,NJ)
  CLOSE(IO)
  P = 0.0
  OPEN(IO,FILE='p_field.dat')
  READ(IO,*) ((P(I,J),I=0,NI),J=0,NJ)
  CLOSE(IO)
  GradP = 0.0
  LapP = 0.0
  !Вычисление пространственных операторов для исследуемой задачи
  Call LeastSquareGradCalc(GradP, NI, NJ, P, CellCenter)
  Call B_CalcLaplacian(NI, NJ, P, GradP, LapP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  divV = 0.0
  rotV = 0.0
  call B_CalcDivergence(mode, NI, NJ, V, divV, P, GradP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  Call BCalcRotor(mode, NI, NJ, V, rotV, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  T = T_0
 
  !=== Solve conv-diff equation===
  T = T_0
  iter = 1
  !Вызов подпрограммы для вычисления локального шага по времени    
  call Localdt(CFL, dt, NI, NJ, CellCenter, IFaceCenter, JFaceCenter)
  !Вывод невязок
  open(IO, FILE='Res.plt')  
  T1 = omp_get_wtime()
  do while (iter <= Niter)
  Res = 0.0
  old_res = maxval(Res(1:NI-1, 1:NJ-1))/Max_Res
  T(0, :) = T_cold
  T(NI, :) = T_hot
  if (scheme .NE. 1) then
      if (gradT_calc == "GG_iter") then
          do i = 1, NIter_GG
                Call B_CalcGradient(GradT, NI, NJ, T, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
          end do
      else if (gradT_calc == "GG") then
          Call B_CalcGradient(GradT, NI, NJ, T, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
      else if (GradT_calc == "LS") then
          Call LeastSquareGradCalc(GradT, NI, NJ, T, CellCenter)
      end if
  end if
  Call B_CalcResiduals(CFL, dt, NI, NJ, T, GradT, Vel_eq, scheme, Re, Pr, Res, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  if (iter == 1) then
      Max_Res = maxval(Res(1:NI-1, 1:NJ-1))
  end if
  write (*,*) iter, ': ', maxval(Res(1:NI-1, 1:NJ-1))/Max_Res
  write (IO, *) iter, maxval(Res(1:NI-1, 1:NJ-1))/Max_Res
  if (abs(maxval(Res(1:NI-1, 1:NJ-1))/Max_Res - old_res) < eps) then !maxval(Res(1:NI-1, 1:NJ-1))/Max_Res < eps)  
      exit
  end if
  iter = iter + 1
  end do
  close(IO)
  end if
  T_error = abs(T-T_flos)/T_flos
  write (*,*) 'Maximum error of T', maxval(T_error(1:NI-1, 1:NJ-1))
  T2 = omp_get_wtime()
  !Вывод времени работы программы
  print *, 'Computation time: ', T2 - T1, ' s'
  !=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP, GradPError, V, divV, divVError, LapP, LapPError, rotV, rotVError, T, Vel_eq, T_error)
  Close(IO)
  
  
END PROGRAM Main  

