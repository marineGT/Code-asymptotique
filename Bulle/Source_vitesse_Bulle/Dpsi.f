      function Dpsi(c2,n,zetap)
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
      double precision :: Dpsi

c     Declaration des variables locales
c     ---------------------------------
      double precision :: kn
c     ==================================================================

      kn=c2*dble(n*(n+1))/(sqrt(2.D0)*dble(2*n+1)*dble(2*n+3))

c     Coefficient par la particule solide
c     -----------------------------------
c      Dpsi=kn*(2.D0*exp(dble(2*n+1)*(zetap))-dble(2*n+1)
c     &     *dexp(2.D0*zetap)+2*n+3)/(2.D0*sinh(dble(2*n+1)*(zetap))
c     &     -dble(2*n+1)*sinh(2.D0*zetap))
c      write(*,*) 'n= ',n,'Dn=',Dpsi
c       write(*,*) 'n= ',n,'kn= ',kn
c     Coefficient pour la bulle
c     -------------------------
      Dpsi=kn*(exp(dble(2*n+1)*zeta0)-exp(2.D0*zeta0))/
     &     (cosh(dble(2*n+1)*zeta0)-cosh(2.D0*zeta0))
c     =======================
c     Fin de la fonction Dpsi
c     =======================
      end function Dpsi
