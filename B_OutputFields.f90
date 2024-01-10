Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP, GradPError, V, divV, divVError, LapP, LapPError, rotV, rotVError, T, Vel_eq, T_error)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P, divV, divVError, LapP, LapPError, rotV, rotVError, T, T_error
  Real,Dimension(0:NI,0:NJ,2)::GradP, GradPError, V, Vel_eq

  Write(IO, '(a190)') 'VARIABLES = "X", "Y", "P", "GradPx", "GradPy", "GradPErrorX", "GradPErrorY", "Vx", "Vy", "divV", "divVError", "LapP", "LapPError", "RotV", "RotVError", "Vel_eq_x", "Vel_eq_y", "T", "T_error"' 
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.7)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.7)') GradPError(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.7)') GradPError(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.7)') V(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.7)') divV(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') divVError(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') LapP(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') LapPError(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') rotV(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') rotVError(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') Vel_eq(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.7)') Vel_eq(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.7)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.7)') T_error(1:NI-1,1:NJ-1)
End Subroutine 
