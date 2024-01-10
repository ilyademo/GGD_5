subroutine B_CalcLaplacian(NI, NJ, P, GradP, LapP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
use Functions
integer:: NI,NJ, iFace
real,dimension(0:NI,0:NJ, 2)::GradP, CellCenter
real,dimension(0:NI,0:NJ):: P, LapP
real CellVolume(NI-1, NJ-1), IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real GF(2), VOL, RF(2), SF(2), RC(2), RN(2), PE, DC, DN, RE(2), GPE(2), PF, DNC, dpdn, dpdn_c, NF(2), RNC(2)
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
        VOL=CellVolume(i,j)
        RC=CellCenter(i,j,:)
        RN=CellCenter(IN,JN,:)
        
        DC=Norm2(RF(:)-RC(:))
        DN=Norm2(RF(:)-RN(:))
        
        DNC = Norm2(CellCenter(IN,JN,:) - CellCenter(i,j,:))
        NF = SF/Norm2(SF)
        
        RNC = (CellCenter(IN,JN,:) - CellCenter(i,j,:))/DNC
        GF(1) = RLinearInterp(DC,DN,GradP(i,j,1), GradP(IN,JN,1))
        GF(2) = RLinearInterp(DC,DN,GradP(i,j,2), GradP(IN,JN,2))
        
        dpdn = (P(IN, JN) - P(i,j))/DNC
        !¬ычисление производной по нормали на границе расчетной области
        if (DN .lt. 1e-5) then
            dpdn_c =dot_product(GradP(i,j,:), NF(:))
            dpdn = 5.0/3.0*dpdn - 2.0/3.0*dpdn_c
            GF(:) = GradP(i,j,:)
        end if
        
        dpdn = dpdn + dot_product(NF - RNC, GF) !ѕоправка на неортогональность
        
        LapP(i,j) = LapP(i,j) + dpdn*Norm2(SF)
        end do
        LapP(i,j) = LapP(i,j) / Vol
    end do
end do 
End Subroutine 
