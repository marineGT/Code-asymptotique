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

      function sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                   Un,dUn)          
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: c, zetap 
      double precision, intent(in) :: xi
      double precision, dimension(Nt), intent(in) :: alpha  
      double precision, intent(in) :: alpha0 
      double precision, dimension(Nt+3), intent(in) :: Pn
      double precision, dimension(Nt+1), intent(in) :: Un  
      double precision, dimension(Nt+1), intent(in) :: dUn         

c     Declaration des variables en sorties
c     ------------------------------------      
      double precision :: sigma0zzB


c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2, aux3
      double precision :: xn, xNt, sinEta, inv2sinEta
      double precision :: xi2, sh, ch
      double precision, dimension(Nt+3) :: Vn
      double precision, dimension(Nt+3) :: dVn
      double precision :: SigmaZetaB

c     ==================================================================       
c     on calcule la composante de la contrainte tangentielle

c     soit sigmazzB = -P + Tau0zzB
c
c       Tau0zzB(zeta,x)= -(ch-x)^(-1/2)/(2c**3)
c                       { (sh*Un [3*Vn + 2(ch-x)V'n ]
c                         +2*(ch-x)*U'n [Vn + 2(ch-x)V'n ]
c                       } 
c     Avec 

c     Vn = [Pn-1 - Pn+1]
c     V'n = (1/(1-xi**2))*[nxPn-1 - nPn -(n+2)xPn+1 +(n+2)Pn+2]
 
c     De plus, en coordonneees bispherique la pression est

c      P = (ch-x)^(1/2)/(c³) sum{alphan*ch[(n+1/2)zeta]*Pn} 

c     avec les betan = 0 car les An=Cn=0 en raison des conditions 
c     limites sur la particule
c     ==================================================================    
c      write(*,*) 'Appel de la fonction sigma0zzB'
c
c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=-dsinh(zetap)

c     Calcul de sin(eta) et 1/sin(eta)
c     --------------------------------   
      sinEta = dsqrt((1.d0-xi)*(1.d0+xi)) 
      inv2sinEta=1.d0/(sinEta**2.d0)
      xi2=xi*xi
      
c$$$       if (xi.eq.0.d0) then
c$$$          write(*,*) aux2, c, ch
c$$$       endif      

c     On definit les tableaux Vn et dVn
c     --------------------------------  
c     on sait que Pn0=1.D0 
      Vn(1)=1.d0-Pn(2)
      dVn(1)=(xi-Pn(1)-3.d0*xi*Pn(2)+3.d0*Pn(3))
      dVn(1)=dVn(1)/(1.d0-xi2)

c$$$       if (xi.eq.(-0.999d0)) then
c$$$          write(*,*) 1, Vn(1), dVn(1)
c$$$       endif         

      do n=2,Nt
         xn=dble(n)
         Vn(n)=Pn(n-1)-Pn(n+1)
         dVn(n)=(xn*xi*Pn(n-1)-xn*Pn(n)
     &   -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
         dVn(n)=dVn(n)/(1.d0-xi2)
      enddo

c     Calcul de la somme dSigma_barre
c     -------------------------------
       aux=0.d0
       aux0=0.d0
       aux1=0.d0
       aux2=0.d0
       aux3=0.d0
       aux2=inv2sinEta/2.d0
       

      !Test (verif avec Mathematica OK)
c      n=1
c$$$       aux1=sh*Un(1)*(3.d0*Vn(1)+2.d0*(ch-xi)*dVn(1))
c$$$     &      +2.d0*(ch-xi)*dUn(1)*(Vn(1)
c$$$     &      +2.d0*(ch-xi)*dVn(1))
c       aux1=sh*Un(1)*2.d0*(ch-xi)*dVn(1)

c$$$      do n=2,Nt 
c$$$       xn=dble(n)
c$$$       aux=aux+(sh*Un(n)*(3.d0*Vn(n)+2.d0*(ch-xi)*dVn(n)) 
c$$$     &          +2.d0*(ch-xi)*dUn(n)*(Vn(n)
c$$$     &          +2.d0*(ch-xi)*dVn(n)))
c$$$c       aux=aux+sh*Un(n)*2.d0*(ch-xi)*dVn(n)  
c$$$      enddo

c     Ligne pour n=0
c     --------------
      aux3=-alpha0*dcosh((1.d0/2.d0)*zetap)

c     Ligne pour n=1
c     --------------
       aux3=aux3-alpha(1)*Pn(1)*dcosh((3.d0/2.d0)*zetap)

       aux1=sh*Un(1)*(3.d0*Vn(1)+2.d0*(ch-xi)*dVn(1))
     &      +2.d0*(ch-xi)*dUn(1)*(Vn(1)
     &      +2.d0*(ch-xi)*dVn(1))

c     Ligne pour n=2,Nt-1
c     -------------------
      do 51 n=2,Nt
       xn=dble(n)      
       aux=aux+(sh*Un(n)*(3.d0*Vn(n)+2.d0*(ch-xi)*dVn(n)) 
     &          +2.d0*(ch-xi)*dUn(n)*(Vn(n)
     &          +2.d0*(ch-xi)*dVn(n)))
       aux3=aux3-Pn(n)*alpha(n)*dcosh((xn+(1.d0/2.d0))*zetap)
 51   continue   

      SigmaZetaB=aux+aux1

c     scaling par  -(ch-x)^{-1/2}/c^3 tau0zzB
c     -------------------------------------
c$$$      SigmaZetaB=-(1.d0/(2.d0*c*c*c))
c$$$     &           *(ch-xi)**(-1.d0/2.d0)
c$$$     &           *SigmaZetaB
      SigmaZetaB=-(1.d0/(2.d0*c*c*c))
     &           *(ch-xi)**(-1.d0/2.d0)
     &           *SigmaZetaB
     
c     scaling par (ch-x)^{1/2}/c^3 les alphan
c     ---------------------------------------
      aux3=aux3*(1.d0/(c*c*c))
     &     *(ch-xi)**(1.d0/2.d0)

       sigma0zzB=aux3+SigmaZetaB
c      sigma0zzB=aux3+SigmaZetaB+3.3493734285132062
c      write(*,*) 'dans fonction sigma0zz=', sigma0zzB
c      write(15,*) xi, sigma0zzB,aux3,SigmaZetaB
c     =============================
c     Fin de la fonction sigma0zzB
c     ==============================
      end function sigma0zzB

