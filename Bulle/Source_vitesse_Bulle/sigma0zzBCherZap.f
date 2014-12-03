c****************************************************************
c                 fonction qui calcule la composante de
c                     la contrainte normale sigma0zz
c                        sur la surface de la bulle
c
c     Appelle par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force d'ordre 1
c
c     Programmeur: M.Guemas; Octobre 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function sigma0zzBCherZap(Nt,c,zetap,xi,alpha,alpha0,
     &                          Pn,Un,dUn)          
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: c, zetap 
      double precision, intent(in) :: xi, alpha0
      double precision, dimension(Nt), intent(in) :: alpha           
      double precision, dimension(Nt+3), intent(in) :: Pn
      double precision, dimension(Nt+1), intent(in) :: Un  
      double precision, dimension(Nt+1), intent(in) :: dUn         

c     Declaration des variables en sorties
c     ------------------------------------      
      double precision :: sigma0zzBCherZap


c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2, aux3
      double precision :: xn, xNt, sinEta, invsinEta
      double precision :: xi2, sh, ch, aux4, xn2
      double precision :: termesA, A1, A2
      double precision :: SigmaZetaB

c     ==================================================================       
c     on calcule la composante de la contrainte tangentielle
c     en utilisant la forme de la contrainte definie par
c     Chervenivanova et Zapryanov (1988).

c     Donc la contrainte totale est 
c
c      sigmazzB = (ch-x)^(-1/2)/c2   
c                 sum { - alphan*ch[(n+1/2)zeta] 
c                       + 2(2n+1)[sh/2*Un + ch*U'n] 
c                       - ((2n-1)/2)*U'n-1
c                       - ((2n+3)/2)*U'n+1
c                      }   
c     ==================================================================    
c      write(*,*) 'Appel de la fonction sigma0zzB'
c
c     Calcul de cosh(zetap) et sinh(zetap)
c     ----------------------------------- 
      ch=dcosh(zetap)
      sh=dsinh(zetap)
c      ch=1.d0
c      sh=0.d0
 
c     Calcul de la somme sigma0zzBCherZap
c     -------------------------------
       aux=0.d0
       aux0=0.d0
       aux1=0.d0
       aux2=0.d0
       aux3=0.d0
       aux4=(3.d0*c*c)/(4.d0*dsqrt(2.d0))

c     on sait que Pn0=1.D0 
c     quand n= 0 -> U'0(-zetap) = 0 et U0(-zetap) = 0

c     Alors si n = 0
c     sigmazzB = -alpha0*ch(zetap/2) -(3/2)*U'1


c     si n = 1 
c     sigmazzB = [-alpha1*ch((3/2)*zetap) + (6*ch*U'1 + 3*sh*U1) -(5/2)*U'2]*P1
c                
c     si n = 1 
c     sigmazzB = [-alpha1*ch((3/2)*zetap) + (6*ch*U'1 + 3*sh*U1) -(5/2)*U'2]*P1
c         


c     Ligne pour n=0
c     --------------
      aux0=-(3.d0)*dUn(1) !tau0zz
c       aux0=+alpha0*dcosh(zetap/2.d0)-3.d0*dUn(1)
      aux3=-alpha0*dcosh((1.d0/2.d0)*zetap) !alphan

      A2=(2.d0/3.d0)*(dexp(-(5.d0/2.d0)*zetap)
     &        -dexp(-0.5d0*zetap))
     
c     Ligne pour n=1
c     --------------
      aux1=+(6.d0*ch*dUn(1)+3.d0*sh*Un(1) !tau0zz
     &     -(5.d0)*dUn(2))*Pn(1)

   

      aux3=aux3-alpha(1)*Pn(1)*dcosh((3.d0/2.d0)*zetap)!alphan
c$$$      aux1=+alpha(1)*Pn(1)*dcosh((3.d0/2.d0)*zetap) !vrai forme
c$$$     &     +(6.d0*ch*dUn(1)+3.d0*sh*Un(1)
c$$$     &     -5.d0*dUn(2))*Pn(1)
       
       A2=A2+aux4*(6.d0/5.d0)*(dexp(-(7.d0/2.d0)*zetap)
     &        -dexp(-(5.d0/2.d0)*zetap))*Pn(1)
      
c     Ligne pour n=2,Nt-1
c     -------------------
      do 51 n=2,Nt
       xn=dble(n)
       xn2=xn*xn
       A1=A1+aux4*(xn*(xn-1.d0)/(2.d0*xn-1.d0))*
     &       dexp(-(xn-3.d0/2.d0)*zetap)
       A2=A2+aux4*((xn+2.d0)*(xn+1.d0)/(2.d0*xn+3.d0)*
     &       dexp(-(xn+5.d0/2.d0)*zetap)
     &       -(2.d0*(2.d0*xn+1.d0)*(xn2-xn-1.d0)
     &       /((2.d0*xn+3.d0)*(2.d0*xn-1.d0)))
     &       *dexp(-(xn+1.d0/2.d0)*zetap))

       aux=aux+((2.d0*xn+1.d0)*(sh*Un(n)+2.d0*ch*dUn(n))
     &     -((2.d0*xn-1.d0))*dUn(n-1)
     &     -((2.d0*xn+3.d0))*dUn(n+1))*Pn(n)   
       aux3=aux3-Pn(n)*alpha(n)*dcosh((xn+(1.d0/2.d0))*zetap)
 51   continue 

c$$$      do 51 n=2,Nt
c$$$       xn=dble(n)
c$$$       aux=aux+Pn(n)*alpha(n)*dcosh((xn+(1.d0/2.d0))*zetap)
c$$$     &     +((2.d0*xn+1.d0)*(sh*Un(n)+2.d0*ch*dUn(n))
c$$$     &     -(2.d0*xn-1.d0)*dUn(n-1)
c$$$     &     -(2.d0*xn+3.d0)*dUn(n+1))*Pn(n)
c$$$ 51   continue 

      SigmaZetaB=aux+aux1+aux0
      
c      sigma0zzB=aux+aux1+aux0
c      write(*,*) 'avant scaling sigma0zz=', sigma0zzB

c     scaling par (ch-x)^{1/2}/c^2 les alphan
c     ---------------------------------------
      aux3=aux3*(1.d0/(c*c*c))
     &     *(ch-xi)**(1.d0/2.d0)

c     scaling par  (ch-x)^{1/2}/c^2 tau0zzB
c     -------------------------------------
      SigmaZetaB=(1.d0/(c*c*c))
     &           *(ch-xi)**(1.d0/2.d0)
     &           *SigmaZetaB

c     scaling par  (ch-x)^{1/2}/c2 
c     ------------------------------
      sigma0zzBCherZap=aux3+SigmaZetaB

c      sigma0zzBCherZap=(1.d0/(c*c))*(ch-xi)
c     &           *SigmaZetaB     
c      write(*,*) 'dans fonction sigma0zz=', sigma0zzB
      write(15,*) xi, sigma0zzBCherZap, SigmaZetaB, aux3,
     &            A1,A2
c      write(8,*) xi, sigma0zzB, zetap,c,ch

c     =============================
c     Fin de la fonction sigma0zzB
c     ==============================
      end function sigma0zzBCherZap

