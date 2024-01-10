Module Functions
    contains
Real Function Pressure(X,Y) 
    implicit none
    real x,y
    Pressure = x**4 + y**4 - x**2 + 5*y**2!exp(x*x) + exp(-y) + x**2 + 2*y**2 + 1!!exp(x*x) + exp(-y) + x + 2*y!!x + y + 1.0! !
End Function

real function RLinearInterp(d1, d2, x1, x2)
    implicit none
    real x1, x2, d1, d2
    RLinearInterp = (x1*d2 + x2*d1)/(d1+d2)
end function

real function DivVelocityPExact(x, y)
implicit none
real x, y
DivVelocityPExact = -2*x*y*exp(x**2) - y - x*exp(-y) + 2*x - 3*exp(-y) + 6!(3*x*y-10)*(y-x) + 9*y**2 + 30!5*y+x+3!-y + x + 3!
end function

real function LapPressureExact(x, y)
implicit none
real x, y
LapPressureExact = 12*x**2 + 12*y**2 + 8.0!2*exp(x*x)*(1.0 + 2*x**2) + exp(-y) + 6.0!
end function

real function RotVelExact(x, y)
implicit none
real x, y
RotVelExact = 3.0*(x**2+y**2) + 1!cos(x)-2*y*exp(y**2)+1.0!
end function

end module
    
Subroutine CalcGradExact(x, y, GPE)
real x, y, GPE(2)
GPE(1) = 2*x!2*x*exp(x*x) + 1.0!1.0!2*x + 1!!3*x**2 + 10!
GPE(2) = 2*y!-exp(-y) + 2.0!2*y + 1!!3*y**2 + 10!
end subroutine    
    
subroutine Velocity(x, y, V)
implicit none
real x, y, V(2)
V(1) = -y**3 + 2*y!exp(y**2) + 2*y!
V(2) = x**3 + 3*x!sin(x) + 3*x!
end subroutine


!Метод Гаусса
subroutine Gauss(GradP, Coeff, ind_1, ind_2, NI, NJ)
implicit none
integer NI, NJ
real,dimension(0:NI,0:NJ, 2)::GradP
real Coeff(2, 3), m
integer:: n = 2
integer:: k, i, j, l, ind_1, ind_2
do k = 2, n
    do j = k, n
        m = Coeff(j, k - 1)/ Coeff(k-1,k-1)
        do i = 1, n + 1
            Coeff(j, i) = Coeff(j, i) - m * Coeff(k-1, i)
        end do
    end do
end do
do i = n, 1, -1
    GradP(ind_1, ind_2, i) = Coeff(i, n + 1) / Coeff(i, i)
    do l = n, i + 1, -1
        GradP(ind_1, ind_2, i) = GradP(ind_1, ind_2, i) - Coeff(i, l) * GradP(ind_1, ind_2, l)/ Coeff(i, i)
    end do
end do
end subroutine
    
      
    