Subroutine B_CalcGradient(GradP, NI, NJ, P, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
use Functions
integer:: NI,NJ, iFace
real,dimension(0:NI,0:NJ, 2)::GradP, CellCenter
real,dimension(0:NI,0:NJ):: P
real CellVolume(NI-1, NJ-1), IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real GP(2), VOL, RF(2), SF(2), RC(2), RN(2), PE, DC, DN, RE(2), GPE(2), PF
do i = 1, NI-1
    do j = 1, NJ-1
        GP = 0.0
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
        VOL=CellVolume(i,j) !Объем ячейки
        RC=CellCenter(i,j,:)!Радиус вектор центра ячейки
        RN=CellCenter(IN,JN,:)!Радиус вектор соседней ячейки
        
        DC=Norm2(RF(:)-RC(:)) !Расстояние до грани
        DN=Norm2(RF(:)-RN(:)) ! Расстояние до центра соседней ячейки
        
        
        PE = RLinearInterp(DC, DN, P(i,j), P(IN, JN)) !Линейная интерполяция
        
        RE(1) = RLinearInterp(DC, DN, CellCenter(i,j,1), CellCenter(IN,JN,1))
        RE(2) = RLinearInterp(DC, DN, CellCenter(i,j,2), CellCenter(IN,JN,2))
        GPE(1) = RLinearInterp(DC, DN, GradP(i,j,1), GradP(IN, JN, 1))
        GPE(2) = RLinearInterp(DC, DN, GradP(i,j,2), GradP(IN, JN, 2))
        
        PF = PE + DOT_PRODUCT(RF(:)-RE(:), GPE(:))
        GP(:) = GP(:) + PF*SF(:)
        end do
        GradP(i,j,:) = GP(:)/VOL
    end do
end do 
End Subroutine 
