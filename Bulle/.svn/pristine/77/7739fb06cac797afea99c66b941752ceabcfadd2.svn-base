      function secondmembre(Nt,c,c2,eps,tabs,tabs0,x)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      integer,intent(in) :: Nt 
      double precision, intent(in) :: eps, x, tabs0, c, c2
      double precision, dimension(Nt+1),intent(in) :: tabs

c     Declaration des variables en sortie
c     -----------------------------------
      double precision :: secondmembre

c     Declaration des variables locales
c     ---------------------------------
      integer :: n
      double precision :: cn

c     Declaration des interfaces
c     --------------------------
      interface
         function flegendre(n,x)
            integer :: n
            double precision :: x
            double precision :: flegendre
         end function flegendre
     
      end interface
c     ==================================================================

c     Calcul du coefficient pour n=0
c     ------------------------------
      secondmembre=tabs0*flegendre(0,x)
      
c      write(*,*)'tabs0=', tabs0
c      write(*,*)'flegendre0=',flegendre(0,x)
c      write(*,*) 'secondmembre0=', secondmembre

c     Calcul de la somme jusqu'à Nt
c     -----------------------------
      do 1 n=1,Nt

c        Calcul du coefficient de la série pour n
c        ----------------------------------------
         cn=tabs(n)*flegendre(n,x)
c         write(*,*)'cn(',n,')=', cn
c         write(*,*)'tabs=', tabs(n)
c         write(*,*)'flegendre(',n,')=',flegendre(n,x)

c        Incrémentation de la somme
c        --------------------------
         secondmembre=secondmembre+cn
 1	 continue      
  
c     Multiplication de la somme par le préfacteur
c     --------------------------------------------
      secondmembre=secondmembre/c2
c      write(*,*) 'secondmembre=', secondmembre
c     ===============================
c     Fin de la fonction secondmembre
c     ===============================
      end function secondmembre
