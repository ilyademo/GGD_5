Subroutine LeastSquareGradCalc(GradP, NI, NJ, P, CellCenter)
use Functions
integer:: NI,NJ, iFace
real,dimension(0:NI,0:NJ, 2)::GradP, CellCenter
real,dimension(0:NI,0:NJ):: P
real RC(2), Coeff(2,3)
do i = 1, NI-1
    do j = 1, NJ-1
        Coeff = 0.0
        do iFace = 1, 4
        select case (iFace)
        case(1)
            IN=I-1
            JN=J
            RC(:)=CellCenter(IN,JN,:)-CellCenter(i,j,:)
        case(2)
            IN=I+1
            JN=J
            RC(:)=CellCenter(IN,JN,:)-CellCenter(i,j,:)
        case(3)
            IN=I
            JN=J-1
            RC(:)=CellCenter(IN,JN,:)-CellCenter(i,j,:)
        case(4)
            IN=I
            JN=J+1
            RC(:)=CellCenter(IN,JN,:)-CellCenter(i,j,:)
        end select
            Coeff(1, 1) = Coeff(1, 1) + RC(1)**2
            Coeff(1, 2) = Coeff(1, 2) + RC(1)*RC(2)
            Coeff(1, 3) = Coeff(1, 3) + RC(1)*(P(IN,JN) - P(i,j))
            Coeff(2, 1) = Coeff(2, 1) + RC(1)*RC(2)
            Coeff(2, 2) = Coeff(2, 2) + RC(2)**2
            Coeff(2, 3) = Coeff(2, 3) + RC(2)*(P(IN,JN) - P(i,j))
        end do
        call Gauss(GradP, Coeff, i, j, NI, NJ)
    end do
end do
end subroutine