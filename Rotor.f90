subroutine BCalcRotor(mode, NI,NJ,V,rotV, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
use Functions
implicit none
integer:: NI, NJ,i,j,k,iFace, IN, JN, mode
real,dimension(0:NI,0:NJ):: rotV
real,dimension(0:NI,0:NJ,2):: CellCenter, V
real CellVolume(NI-1, NJ-1), IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real VOL, RF(2), SF(2), RC(2), RN(2), VF(2), DC, DN, PF, GFLUX, PN, GC, GB, GN, rot
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
        
        rotV(i,j) = rotV(i,j) + (SF(1)*VF(2) - SF(2)*VF(1))
        
        end do
        rotV(i,j) = rotV(i,j)/VOL
    end do
end do 

end subroutine