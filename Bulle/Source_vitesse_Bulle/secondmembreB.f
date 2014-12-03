      function secondmembreB(Nt,zetap,c,c2,eps,lambda0,gammahat,
     &                        tabsB,tabsB0,xi,Rn,R0)
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
      double precision, intent(in) :: c, c2, eps, lambda0,zetap
      double precision, intent(in) :: gammahat, tabsB0, xi 
      double precision, intent(in) :: R0
      double precision, dimension(Nt+1),intent(in) :: tabsB
      double precision, dimension(Nt),intent(in) :: Rn

c     Declaration des variables en sortie
c     -----------------------------------
      double precision :: secondmembreB

c     Declaration des variables locales
c     ---------------------------------
      integer :: n
      double precision :: cn,ch,sh
      double precision :: pn,pi,aux,aux3,aux2

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
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)

c     calcul de ch et sh
c     ------------------
c      ch=dcosh(zetap)
c      sh=dsinh(zetap)
      sh=0.d0
      ch=1.d0

c     Calcul du coefficient pour n=0
c     ------------------------------
c      secondmembreB=(-(1.d0/(c*c2))*tabs0+dsqrt(2.d0)*3.d0*lambda0
c     &              *c2*R0)*flegendre(0,xi) !autre expression du second membre
c       secondmembreB=+(1.d0/(c*c2))*tabs0*flegendre(0,xi)
       secondmembreB=tabsB0*flegendre(0,xi)

       aux3=R0*flegendre(0,xi)

c     Calcul du coefficient pour n=1
c     ------------------------------      
c      cn=(1.d0/(c*c2))*tabs(1)*flegendre(1,xi)
      cn=tabsB(1)*flegendre(1,xi)
      aux3=aux3+Rn(1)*flegendre(1,xi)
      secondmembreB=secondmembreB+cn
      
c     Calcul de la somme jusqu'à Nt
c     -----------------------------
      do 1 n=2,Nt

c        Calcul du coefficient de la série pour n
c        ----------------------------------------
c         cn=+(1.d0/(c*c2))*tabs(n)*flegendre(n,xi)
         cn=tabsB(n)*flegendre(n,xi)
         
         aux3=aux3+Rn(n)*flegendre(n,xi)
c        Incrémentation de la somme
c        --------------------------
         secondmembreB=secondmembreB+cn
 1	 continue      
          
        
c     Ajout de la contribution des force de gravites
c     ----------------------------------------------
c      aux=(3.D0)*lambda0*sh/(ch-xi)
      aux2=(3.D0)*lambda0*c*dsqrt(2.d0)*aux3
      aux2=((ch-xi)**(1.d0/2.d0))*aux2
      
      aux3=(3.D0)*lambda0*c2*dsqrt(2.d0)*aux3

      secondmembreB=(((ch-xi)**(1.d0/2.d0))/(c2*c))
     &              *secondmembreB  

      write(16,*) xi,secondmembreB,aux3            
c     ===============================
c     Fin de la fonction secondmembre
c     ===============================
      end function secondmembreB
