
!===================================================================================================================

SUBROUTINE RS (NM, N, A, W, Z, FV, IERR)

!     BEGIN PROLOGUE  RS
!     PURPOSE   Compute the eigenvalues and, optionally, the eigenvectors of a real symmetric matrix.
!     LIBRARY   SLATEC (EISPACK)
!     CATEGORY  D4A1
!     TYPE      SINGLE PRECISION (RS-S, CH-C)
!     KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!     AUTHOR    Smith, B. T., et al.
!     DESCRIPTION
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional array parameters, A and Z, as declared in the calling
!           program dimension statement.  NM is an INTEGER variable.
!        N is the order of the matrix A.  N is an INTEGER variable. N must be less than or equal to NM.
!        A contains the real symmetric matrix.  A is a two-dimensional array, dimensioned A(NM,N).
!        FV is one-dimensional array used for temporary storage, dimensioned FV1(N).
!        MATZ is an INTEGER variable set equal to zero if only eigenvalues are desired.  Otherwise, it is set to any
!           non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        A is unaltered.
!        W contains the eigenvalues in ascending order.  W is a one-dimensional DOUBLE PRECISION array, dimensioned W(N).
!        Z contains the eigenvectors if MATZ is not zero.  The eigenvectors are orthonormal.  Z is a two-dimensional
!           array, dimensioned Z(NM,N).
!        IERR is an INTEGER flag set to
!           0          for normal return,
!           10*N       if N is greater than NM,
!           J          if the J-th eigenvalue has not been determined after 30 iterations.
!                      The eigenvalues, and eigenvectors if requested, should be correct for indices 1, 2, ..., IERR-1.
!
!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     --------------------------------------------------------------------------------------------------------------------
!
!     REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag, 1976.
!     ROUTINES CALLED  TQL2, TQLRAT,  TRED2
!     REVISION HISTORY  (YYMMDD)
!     760101  DATE WRITTEN
!     890831  Modified array declarations.  (WRB)
!     890831  REVISION DATE from Version 3.2
!     891214  Prologue converted to Version 4.0 format.  (BAB)
!     920501  Reformatted the REFERENCES section.  (WRB)
!     END PROLOGUE  RS

      use MyPrec
      IMPLICIT NONE
      INTEGER N, NM, IERR
      REAL(RP) A(NM, N), Z(NM, N), W(N), FV(N)

!     FIRST EXECUTABLE STATEMENT  RS
      IERR = 10 * N
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      CALL  TRED2(NM, N, A, W,  FV, Z)
      CALL  TQL2 (NM, N, W, FV, Z,  IERR)

END SUBROUTINE RS

!===================================================================================================================

SUBROUTINE TRED2 (NM, N, A, D, E, Z)

!     BEGIN PROLOGUE  TRED2
!     PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations.
!     LIBRARY   SLATEC (EISPACK)
!     CATEGORY  D4C1B1
!     TYPE      SINGLE PRECISION (TRED2-S)
!     KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!     AUTHOR  Smith, B. T., et al.
!     DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED2, NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional array parameters, A and Z, as declared in the calling
!           program dimension statement.  NM is an INTEGER variable.
!        N is the order of the matrix A.  N is an INTEGER variable. N must be less than or equal to NM.
!        A contains the real symmetric input matrix.  Only the lower triangle of the matrix need be supplied.  A is a two-
!           dimensional real array, dimensioned A(NM,N).
!
!     On Output
!
!        D contains the diagonal elements of the symmetric tridiagonal matrix.  D is a one-dimensional real array, dimensioned D(N).
!        E contains the subdiagonal elements of the symmetric tridiagonal matrix in its last N-1 positions.  E(1) is set
!           to zero.  E is a one-dimensional real array, dimensioned E(N).
!        Z contains the orthogonal transformation matrix produced in the reduction.  Z is a two-dimensional real array, dimensioned Z(NM,N).
!        A and Z may coincide.  If distinct, A is unaltered.
!
!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     --------------------------------------------------------------------------------------------------------------------
!
!     REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag, 1976.
!     ROUTINES CALLED  (NONE)
!     REVISION HISTORY  (YYMMDD)
!     760101  DATE WRITTEN
!     890831  Modified array declarations.  (WRB)
!     890831  REVISION DATE from Version 3.2
!     891214  Prologue converted to Version 4.0 format.  (BAB)
!     920501  Reformatted the REFERENCES section.  (WRB)
!     END PROLOGUE  TRED2

      use MyPrec
      IMPLICIT NONE
      INTEGER NM, N
      INTEGER I, J, K, L, II, JP1
      REAL(RP) A(NM, N), Z(NM, N), D(N), E(N)
      REAL(RP) F, G, H, HH, SCALE

!     FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N

         DO 100 J = 1, I
            Z(I, J) = A(I, J)
  100 CONTINUE

      IF (N .EQ. 1) GO TO 320
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0._RP
         SCALE = 0._RP
         IF (L .LT. 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I, K))

         IF (SCALE .NE. 0._RP) GO TO 140
  130    E(I) = Z(I, L)
         GO TO 290

  140    DO 150 K = 1, L
            Z(I, K) = Z(I, K) / SCALE
            H = H + Z(I, K) * Z(I, K)
  150    CONTINUE

         F = Z(I, L)
         G = -SIGN(SQRT(H), F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I, L) = F - G
         F = 0._RP

         DO 240 J = 1, L
            Z(J, I) = Z(I, J) / H
            G = 0._RP
!     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J, K) * Z(I, K)

            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220

            DO 200 K = JP1, L
  200       G = G + Z(K, J) * Z(I, K)
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I, J)
  240    CONTINUE

         HH = F / (H + H)
!     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I, J)
            G = E(J) - HH * F
            E(J) = G

            DO 260 K = 1, J
               Z(J, K) = Z(J, K) - F * E(K) - G * Z(I, K)
  260    CONTINUE

  290    D(I) = H
  300 CONTINUE

  320 D(1) = 0._RP
      E(1) = 0._RP
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0._RP) GO TO 380

         DO 360 J = 1, L
            G = 0._RP

            DO 340 K = 1, L
  340       G = G + Z(I, K) * Z(K, J)

            DO 360 K = 1, L
               Z(K, J) = Z(K, J) - G * Z(K, I)
  360    CONTINUE

  380    D(I) = Z(I, I)
         Z(I, I) = 1._RP
         IF (L .LT. 1) GO TO 500

         DO 400 J = 1, L
            Z(I, J) = 0._RP
            Z(J, I) = 0._RP
  400    CONTINUE

  500 CONTINUE

END SUBROUTINE TRED2

!===================================================================================================================

subroutine TQL2(NM,N,D,E,Z,IERR)

!**************************************************************************************************    
! PROGRAM:  TQL2(NM,N,D,E,Z,IERR)
! TYPE:     subroutine
! PURPOSE:  Compute the eigenvalues and eigenvectors of symmetric tridiagonal matrix.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4A5, D4C2A
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQL2, NUM. MATH. 
!     11, 293-306(1968) by Bowdler, Martin, Reinsch, and Wilkinson. HANDBOOK FOR 
!     AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors of a SYMMETRIC TRIDIAGONAL 
!     matrix by the QL method. The eigenvectors of a FULL SYMMETRIC matrix can also be 
!     found if  TRED2  has been used to reduce this full matrix to tridiagonal form.
!
!     On Input
!
!        NM: must be set to the row dimension of the two-dimensional array parameter, 
!            Z, as declared in the calling program dimension statement.  NM is an INTEGER variable.
!        N:  the order of the matrix.  N is an INTEGER variable. N must be less than or equal to NM.
!        D:  contains the diagonal elements of the symmetric tridiagonal matrix.  
!            D is a one-dimensional REAL array, dimensioned D(N).
!        E:  contains the subdiagonal elements of the symmetric tridiagonal matrix in its last N-1 positions. 
!            E(1) is arbitrary.  E is a one-dimensional REAL array, dimensioned E(N).
!        Z:  contains the transformation matrix produced in the reduction by  TRED2, if performed. 
!            If the eigenvectors of the tridiagonal matrix are desired, Z must contain the identity matrix. 
!            Z is a two-dimensional REAL array, dimensioned Z(NM,N).
!
!      On Output
!
!        D: contains the eigenvalues in ascending order.  If an error exit is made, 
!           the eigenvalues are correct but unordered for indices 1, 2, ..., IERR-1.
!        E: has been destroyed.
!        Z: contains orthonormal eigenvectors of the symmetric tridiagonal (or full) matrix. 
!           If an error exit is made, Z contains the eigenvectors associated with the stored eigenvalues.
!        IERR is an INTEGER flag set to
!           Zero       for normal return,
!           J          if the J-th eigenvalue has not been determined after 30 iterations.
!
!        Calls PyThagC(A,B) for sqrt(A**2+B**2).
!
!     Questions and comments should be directed to B. S. Garbow, 
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,Y. Ikebe,
!              V. C. Klema and C. B. Moler, Matrix Eigen-system Routines - EISPACK Guide, Springer-Verlag,1976.
! SUBROUTINES CALLED:  PyThagC
! REVISION HISTORY  (YYMMDD)
!     760101  DATE WRITTEN
!     890831  Modified array declarations.  (WRB)
!     890831  REVISION DATE from Version 3.2
!     891214  Prologue converted to Version 4.0 format.  (BAB)
!     920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  TQL2
!**************************************************************************************************  

        use MyPrec
        implicit none

        integer NM,N,IERR
        integer I,J,K,L,M,II,L1,L2,MML
        real(rp) B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
        real(rp) PyThag
        real(rp) Z(NM,N)
        real(rp) D(N),E(N)
		
!***FIRST EXECUTABLE STATEMENT  TQL2
		IERR=0
		IF (N .EQ. 1) GO TO 1001
!
		DO 100 I = 2, N
100   E(I-1) = E(I)
!
		F = 0.0_rp
		B = 0.0_rp
		E(N) = 0.0_rp
!
		DO 240 L = 1, N
			J = 0
			H = ABS(D(L)) + ABS(E(L))
			IF (B .LT. H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
			DO 110 M = L, N
				IF (B + ABS(E(M)) .EQ. B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
110      CONTINUE
!
120      IF (M .EQ. L) GO TO 220
130      IF (J .EQ. 30) GO TO 1000
			J = J + 1
!     .......... FORM SHIFT ..........
			L1 = L + 1
			L2 = L1 + 1
			G = D(L)
			P = (D(L1) - G) / (2.0_rp * E(L))
			R = PyThag(P,1.0_rp)
			D(L) = E(L) / (P + SIGN(R,P))
			D(L1) = E(L) * (P + SIGN(R,P))
			DL1 = D(L1)
			H = G - D(L)
			IF (L2 .GT. N) GO TO 145
!
			DO 140 I = L2, N
140      D(I) = D(I) - H
!
145      F = F + H
!     .......... QL TRANSFORMATION ..........
			P = D(M)
			C = 1.0_rp
			C2 = C
			EL1 = E(L1)
			S = 0.0_rp
			MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
			DO 200 II = 1, MML
				C3 = C2
				C2 = C
				S2 = S
				I = M - II
				G = C * E(I)
				H = C * P
				IF (ABS(P) .LT. ABS(E(I))) GO TO 150
				C = E(I) / P
				R = SQRT(C*C+1.0_rp)
				E(I+1) = S * P * R
				S = C / R
				C = 1.0_rp / R
				GO TO 160
150         C = P / E(I)
				R = SQRT(C*C+1.0_rp)
				E(I+1) = S * E(I) * R
				S = 1.0_rp / R
				C = C * S
160         P = C * D(I) - S * G
				D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
				DO 180 K = 1, N
					H = Z(K,I+1)
					Z(K,I+1) = S * Z(K,I) + C * H
					Z(K,I) = C * Z(K,I) - S * H
180         CONTINUE
!
200      CONTINUE
!
			P = -S * S2 * C3 * EL1 * E(L) / DL1
			E(L) = S * P
			D(L) = C * P
			IF (B + ABS(E(L)) .GT. B) GO TO 130
220      D(L) = D(L) + F
240   CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
		DO 300 II = 2, N
			I = II - 1
			K = I
			P = D(I)
!
			DO 260 J = II, N
				IF (D(J) .GE. P) GO TO 260
				K = J
				P = D(J)
260      CONTINUE
!		
			IF (K .EQ. I) GO TO 300
			D(K) = D(I)
			D(I) = P
!
			DO 280 J = 1, N
				P = Z(J,I)
				Z(J,I) = Z(J,K)
				Z(J,K) = P
280      CONTINUE
!
300   CONTINUE
!
		GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
1000  IERR=L
1001  RETURN

end subroutine

!==================================================================================================
    
function pythag(a, b)
!	find sqrt(a ** 2 + b ** 2) without overflow or destructive underflow

	use MyPrec
	implicit none
	real(rp) pythag, a, b, p, r, s, t

	p = max(abs(a), abs(b))
	
	if (p == 0._rp) then
		pythag = p
		return
	end if

	r = (min(abs(a), abs(b)) / p) ** 2
	do while (4 + r /= 4)
		s = r / (4._rp + r)
		t = 1 + 2 * s
		p = t * p
		r = (s / t) ** 2 * r
	end do

	pythag = p

end function

!==================================================================================================