subroutine B_CalcDivergence(mode, NI,NJ,V,divV, P, GradP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
use Functions
implicit none
integer:: NI, NJ,i,j,k,iFace, IN, JN, mode
real,dimension(0:NI,0:NJ):: divV, P
real,dimension(0:NI,0:NJ,2):: CellCenter, V, GradP
real CellVolume(NI-1, NJ-1), IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real VOL, RF(2), SF(2), RC(2), RN(2), VF(2), DC, DN, PF, GFLUX, PN, GC, GB, GN
do i = 1, NI-1
    do j = 1, NJ-1
        do iFace = 1, 4
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
        
        
        do k = 1, 2
        VF(k) = RLinearInterp(DC, DN, V(i,j,k), V(IN, JN,k))
        end do
        !Скалярное произведение вектора площади грани на скорсть
        GFLUX = dot_product(SF, VF)
        !Выбор схемы
        select case(mode)
        case(1) !Центральная
            PF = RLinearInterp(DC, DN, P(i, j), P(IN, JN))
        case(2) !FOU
            
            if(GFLUX >= 0) then
                PF = P(i,j)
            else
                PF = P(IN, JN)
                if (DN < 1e-6) PF = 2*P(IN,JN) - P(i,j)
            end if
        case (3) !SOU
            if(GFLUX >= 0) then
                PF = P(i,j) + dot_product(GradP(i,j,:), RF - CellCenter(i,j,:))
            else
                PF = P(IN,JN) + dot_product(GradP(IN,JN,:), RF - CellCenter(IN,JN,:))
                if (DN < 1e-6) then
                    PN = 2*P(IN,JN) - P(i,j)
                    GC = dot_product(GradP(i,j,:), CellCenter(i,j,:) - RF)
                    GB = P(i,j) - P(IN,JN)
                    GN = 4*GB - 3*GC
                    PF = PN + GN
                end if
            end if
        end select
        divV(i,j) = divV(i,j) + dot_product(PF*VF,SF)
        
        end do
        divV(i,j) = divV(i,j)/VOL
    end do
end do 

end subroutine