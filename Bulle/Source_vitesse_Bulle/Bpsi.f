      function Bpsi(c2,n,zetap)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      double precision, intent(in) :: c2,zetap
      integer, intent(in) :: n

c     Declaration des variables en sortie
c     -----------------------------------
      double precision :: Bpsi

c     Declaration des variables locales
c     ---------------------------------
      double precision :: kn,xn
c     ==================================================================
c2345678912345678912345678912345678912345678912345678912345678912
      xn=dble(n)
      kn=c2*(xn*(xn+1.d0))/(dsqrt(2.D0)*(2.d0*xn+1.d0)
     &   *(2.d0*xn-1.d0))

c     Coefficient pour la particule solide
c      Bpsi=-kn*(2.D0*dexp((2*n+1)*zetap)+(2*n+1)*dexp(-2.D0*zetap)
c     &     -dble(2*n-1))/(2.D0*dsinh(dble(2*n+1)*(zetap))
c     &     -dble(2*n+1)*dsinh(2.D0*zetap))
c      write(*,*) 'n= ',n,'Bn=',Bpsi
c       write(*,*) 'n= ',n,'kn= ',kn    
c     Coefficient pour la bulle
      Bpsi=kn*(exp(-2.D0*zeta0)-exp(dble(2*n+1)*zeta0))/
     &     (cosh(dble(2*n+1)*zeta0)-cosh(2.D0*zeta0))

c     =======================
c     Fin de la fonction Bpsi
c     =======================
      end function Bpsi
