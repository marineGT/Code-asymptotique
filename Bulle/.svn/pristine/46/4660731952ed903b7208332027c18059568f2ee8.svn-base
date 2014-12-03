      function flegendre(n,x)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      integer, intent(in) :: n
      double precision, intent(in) :: x

c     Declaration des variables en sortie
c     -----------------------------------
      double precision :: flegendre

c     Declaration des variables locales
c     ---------------------------------
      integer :: m
      double precision :: Pmm1,Pmm2

c     ==================================================================
     
      if (n.eq.0) then
         flegendre=1.D0
      elseif (n.eq.1) then
         flegendre=x
      else         
         Pmm2=1.D0
         Pmm1=x
         do m=2,n
c           Calcul de Pm
            flegendre=(dble(2*m-1)*x*Pmm1-dble(m-1)*
     &                Pmm2)/dble(m)
c            write(*,*) 'flegendre(',m,')',flegendre
c           Sauvegarde des Pmn1 et Pmn2
            Pmm2=Pmm1
            Pmm1=flegendre                    
         enddo
         
      endif
      
c     ===========================
c     Fin de la fonction legendre
c     ===========================
      end function flegendre
