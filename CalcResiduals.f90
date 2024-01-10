subroutine B_CalcResiduals(CFL, dt, NI, NJ, T, GradT, Vel_eq, scheme, Re, Pr, Res, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
use Functions
include 'omp_lib.h'
integer:: NI,NJ, iFace, scheme
real,dimension(0:NI,0:NJ, 2):: CellCenter, Vel_eq, GradT
real,dimension(0:NI,0:NJ):: T, Res, dt
real CellVolume(NI-1, NJ-1), IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real GF(2), VOL, RF(2), SF(2), RC(2), RN(2), DC, DN, DNC, dTdn, NF(2), RNC(2), Re, Pr, CFL, VF(2), TF, TN
!Параллелизация программы
!$OMP PARALLEL PRIVATE(i, j, iFace, IN, JN, RF, SF, VOL, RC, RN, DC, DN, NF, k, VF, GFLUX, TF, TN, GC, GB, GN, DNC, dTdn)
!$OMP DO SCHEDULE(DYNAMIC)
    DO i = 1, NI-1
        DO j = 1, NJ-1
            DO iFace = 1, 4
            select case (iFace)
            case(1)
                IN=I-1
                JN=J
                RF(:)=IFaceCenter(i,j,:) !Radus Vector center 
                SF(:)=-IFaceVector(i,j,:)
            case(2)
                IN=I+1
                JN=J
                RF(:)=IFaceCenter(i+1,j,:)
                SF(:)=IFacevector(i+1,j,:)
            case(3)
                IN=I
                JN=J-1
                RF(:)=JFaceCenter(i,j,:)
                SF(:)=-JFaceVector(i,j,:)
            case(4)
                IN=I
                JN=J+1
                RF(:)=JFaceCenter(i,j+1,:)
                SF(:)=JFaceVector(i,j+1,:)
            end select
        
            VOL=CellVolume(i,j)
            RC=CellCenter(i,j,:)
            RN=CellCenter(IN,JN,:)
        
            DC=Norm2(RF(:)-RC(:))
            DN=Norm2(RF(:)-RN(:))
        
            NF = SF/Norm2(SF)
        
            do k = 1, 2
            VF(k) = RLinearInterp(DC, DN, Vel_eq(i,j,k), Vel_eq(IN, JN, k))
            end do
        
            GFLUX = dot_product(SF, VF)
            
            select case(scheme)
            case(1) !Центральная схема
                TF = RLinearInterp(DC, DN, T(i,j), T(IN, JN))
            case(2) !Противопоточная схема 2-го порядка
                if(GFLUX >= 0) then
                    TF = T(i,j) + dot_product(GradT(i,j,:), RF - CellCenter(i,j,:))
                else
                    TF = T(IN,JN) + dot_product(GradT(IN,JN,:), RF - CellCenter(IN,JN,:))
                    if (DN < 1e-6) then
                        TN = 2*T(IN,JN) - T(i,j)
                        GC = dot_product(GradT(i,j,:), CellCenter(i,j,:) - RF)
                        GB = T(i,j) - T(IN,JN)
                        GN = 4*GB - 3*GC
                        TF = TN + GN
                    end if
                end if
            end select
        
        
            DNC = Norm2(CellCenter(IN,JN,:) - CellCenter(i,j,:))
        
            dTdn = (T(IN, JN) - T(i,j))/DNC 
            !Вычисление невязки
            if ((JN .NE. 0) .and. (JN .NE. NJ)) then 
            Res(i,j) = Res(i,j) + (dTdn /(Re*Pr) - dot_product(VF, NF)*TF) * norm2(SF) 
            end if
            end do
        end do
    end do 
!$OMP END DO    
!$OMP END PARALLEL
T(1:NI-1, 1:NJ-1) = T(1:NI-1, 1:NJ-1) + dt(1:NI-1, 1:NJ-1)/CellVolume(1:NI-1, 1:NJ-1)*Res(1:NI-1, 1:NJ-1)
End Subroutine 

!Подпрограмма для вычисления локального шага по времени    
subroutine Localdt(CFL, dt, NI, NJ, CellCenter, IFaceCenter, JFaceCenter)
include 'omp_lib.h'
integer:: NI, NJ, iFace
real,dimension(0:NI,0:NJ, 2):: CellCenter
real,dimension(0:NI,0:NJ):: dt
real IFaceCenter(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), dx_av, dx_i, dx_j, RF(2)
!$OMP PARALLEL PRIVATE(i, j, iFace, RF, dx_i, dx_j, dx_av)
!$OMP DO SCHEDULE (DYNAMIC)
do i = 1, NI-1
    do j = 1, NJ-1
        do iFace = 1, 4
        select case (iFace)
        case(1)
            IN=I-1
            JN=J
            RF(:)=IFaceCenter(i,j,:) !Radus Vector center 
        case(2)
            IN=I+1
            JN=J
            RF(:)=IFaceCenter(i+1,j,:)
        case(3)
            IN=I
            JN=J-1
            RF(:)=JFaceCenter(i,j,:)
        case(4)
            IN=I
            JN=J+1
            RF(:)=JFaceCenter(i,j+1,:)
        end select
        
        if (iFace == 2) then
            dx_i = Norm2(RF - IFaceCenter(i,j,:))
        else if (iFace == 4) then
            dx_j = Norm2(RF - JFaceCenter(i,j,:))
        end if
        
        end do
        
        dx_av = (dx_i + dx_j)/2
        dt(i,j) = CFL*dx_av
    end do
end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine
