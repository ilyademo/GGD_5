SUBROUTINE B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector, dx) 
  REAL X(NI,NJ),Y(NI,NJ),&                               !        input: nodes coordinates
       CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !        output: cell centers and volumes
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !        face centers and vectors for I-faces
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)   !        face centers and vectors for J-faces
  REAL r(2), dx
  
  !=== FACE CENTERS AND FACE VECTORS ===
  ! I-DIRECTION
  DO J = 1,NJ-1
    DO I = 1,NI
      r(1) = X(I,J+1) - X(I,J)  ! r = vector from one node to another
      r(2) = Y(I,J+1) - Y(I,J)
      IFaceVector(I,J,1) = r(2) ! IFaceVector = r rotated on 90 degree
      IFaceVector(I,J,2) =-r(1) ! IFaceVector directed to increasing I-index
      IFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i,j+1))
      IFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i,j+1))
      if (dx > r(2)) then
          dx = r(2)
      end if
    ENDDO
  ENDDO

  ! J-DIRECTION
  DO J = 1,NJ
    DO I = 1,NI-1
      r(1) = X(I+1,J) - X(I,J)  ! r = vector from one node to another
      r(2) = Y(I+1,J) - Y(I,J)
      JFaceVector(I,J,1) =-r(2) ! JFaceVector = r rotated on -90 degree
      JFaceVector(I,J,2) = r(1) ! JFaceVector directed to increasing J-index 
      JFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i+1,j))
      JFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i+1,j))
      if (dx > r(1)) then
          dx = r(1)
      end if
    ENDDO
  ENDDO


 !=== CELL VOLUMES ===
  DO J = 1,NJ-1
    DO I = 1,NI-1
      r(1)=X(I+1,J+1) - X(I,J)
      r(2)=Y(I+1,J+1) - Y(I,J)
      CellVolume(I,J) = 0.5*DOT_PRODUCT(IFaceVector(I,J,:),r)& ! sum surfaces of two triangles
                      + 0.5*DOT_PRODUCT(JFaceVector(I,J,:),r)
    ENDDO
  ENDDO


  !=== CELL CENTERS ===
  ! FOR INNER CELLS: CENTER OF CONTOUR (sum of FaceCenter*FaceLength/Perimeter)
  DO J = 1,NJ-1
    DO  I = 1,NI-1
      CellCenter(I,J,:) = ( IFaceCenter(I  ,J,:)*Norm2(IFaceVector(I  ,J,:))+&
                            IFaceCenter(I+1,J,:)*Norm2(IFaceVector(I+1,J,:))+&
                            JFaceCenter(I,J  ,:)*Norm2(JFaceVector(I,J  ,:))+&
                            JFaceCenter(I,J+1,:)*Norm2(JFaceVector(I,J+1,:)) )&
                         /( Norm2(IFaceVector(I,J,:))+Norm2(IFaceVector(I+1,J,:))+&
                            Norm2(JFaceVector(I,J,:))+Norm2(JFaceVector(I,J+1,:)) )
    ENDDO
  ENDDO

  ! FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
  ! I-BOUNDARIES -----------------------------------------------------
  DO NBOUND = 1,2
    IF (NBOUND.EQ.1) THEN
      IBOUND =  1; IOUT =  0
    ELSE 
      IBOUND = NI; IOUT =  NI
    ENDIF
    DO J = 1,NJ-1
      CellCenter(IOUT,J,:) = IFaceCenter(IBOUND,J,:)
    ENDDO
  ENDDO

  ! J-BOUNDARIES -----------------------------------------------------
  DO NBOUND = 1,2
    IF (NBOUND.EQ.1) THEN
      JBOUND = 1;  JOUT =  0
    ELSE 
      JBOUND = NJ; JOUT =  NJ
    ENDIF
    DO  I = 1,NI-1
      CellCenter(I,JOUT,:) = JFaceCenter(I,JBOUND,:) 
    ENDDO
  ENDDO

END SUBROUTINE


