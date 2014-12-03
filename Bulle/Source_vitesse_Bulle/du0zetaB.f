c****************************************************************
c               fonction qui calcule la composante de
c                    la vitesse duzeta sur la 
c                        surface libre
c     
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Septembre 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function du0zetaB(Nt,c,zetap,xi,Pn,Un,dUn)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt
      double precision, intent(in) :: xi, c, zetap
      double precision, dimension(Nt+3), intent(in) :: Pn 
      double precision, dimension(Nt+1), intent(in) :: Un
      double precision, dimension(Nt+1), intent(in) :: dUn  

c     declaration des variables locales
c     ---------------------------------
      integer :: n   
      double precision :: aux, aux0, aux1, aux2
      double precision :: aux3, aux4, aux5 
      double precision :: xn, xNt, xi2, ch, sh
      double precision :: duZetaB
     
c     Declaration des variables en sorties
c     ------------------------------------      
      double precision :: du0zetaB

c     ==================================================================       
c      on calcule la composante de la vitesse du0zetaB 
c
c      du0zetaB = {(-(ch-x)^(-3/2))/4c²}
c                   {3Vn[-sh*Un +2(ch-x)U'n] 
c                    +2(ch-x)V'n[2(ch-x)*U'n+sh*Un]}

c      ou bien suivant les Pn

c      du0zetaB = {(-(ch-x)^(-3/2))/4c²}
c                   {3[Pn-1 - Pn+1] (-sh*Un +2(ch-x)U'n) 
c                    +(2(ch-x)/(1-x²))*[nxPn-1 -nPn -(n+2)xPn+1  
c                    +(n+2)Pn+2] (2(ch-x)*U'n+sh*Un)}

c     Calcul de cosh(zetap) et sinh(zetap)
c     -----------------------------------      
      ch=dcosh(zetap)
      sh=dsinh(zetap)

c     On initialise duZeta_barre à O
c     --------------------------------      
      duZetaB=0.d0

c     On note que la somme de duZeta est bien de  n=1 à n=Nt+1   
c     Calcul de la somme dSigma_barre
c     -------------------------------
       aux=0.d0       
       aux1=0.d0
       aux2=0.d0
       aux3=0.d0
       aux4=0.d0
       aux5=0.d0

       xi2=xi*xi
       aux4=(2.d0*(ch-xi))/(1.d0-xi2)

c      on sait que Pn0=1.D0 
c      quand n= 0 U'0(-zetap) = 0 et U0(-zetap) = 0
c       du0zetaB = 0
c      quand n= 1,  
c      du0zetaB ={(-(ch-x)^(-3/2))/4c²}
c                {3[1 - P(2)] (-sh*U(1)+2(ch-x)U'(1)) 
c                 +(2(ch-x)/(1-x²))*[+x -P(1) -3xP(2)  
c                 +3P(3)] (2(ch-x)*U'(1)+sh*U(1))}

      !Test
c$$$       n=1
c$$$       aux1=3.d0*(1.d0-Pn(2))*(2.d0*(ch-xi)*dUn(1)
c$$$     &     -sh*Un(1))+aux4*(xi-Pn(1)-3.d0*xi*Pn(2)
c$$$     &     +3.d0*Pn(3))*(2.d0*(ch-xi)*dUn(1)+sh*Un(1))
c$$$
c$$$c$$$       n=3
c$$$       do 7 n=2,Nt+1 
c$$$       xn=dble(n)
c$$$       aux=aux+3.d0*(Pn(n-1)-Pn(n+1))*(2.d0*(ch-xi)*dUn(n)
c$$$     &     -sh*Un(n))+aux4*(xn*xi*Pn(n-1)-xn*Pn(n)
c$$$     &     -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
c$$$     &     *(2.d0*(ch-xi)*dUn(n)+sh*Un(n))
c$$$ 7      continue

c$$$      n=Nt+1
c$$$      xn=dble(n)
c$$$       aux2=3.d0*(Pn(n-1)-Pn(n+1))*(2.d0*(ch-xi)*dUn(n)
c$$$     &     -sh*Un(n))+aux4*(xn*xi*Pn(n-1)-xn*Pn(n)
c$$$     &     -(xn+2.d0)*xi*Pn(n+1))
c$$$     &     *(2.d0*(ch-xi)*dUn(n)+sh*Un(n))

     

c     Ligne pour n=1
c     --------------
      aux1=3.d0*(1.d0-Pn(2))*(2.d0*(ch-xi)*dUn(1)
     &     -sh*Un(1))+aux4*(xi-Pn(1)-3.d0*xi*Pn(2)
     &     +3.d0*Pn(3))*(2.d0*(ch-xi)*dUn(1)+sh*Un(1))

c     Ligne pour n=2,Nt+1
c     -------------------
      do 7 n=2,Nt+1 
       xn=dble(n)    
       aux=aux+3.d0*(Pn(n-1)-Pn(n+1))*(2.d0*(ch-xi)*dUn(n)
     &     -sh*Un(n))+aux4*(xn*xi*Pn(n-1)-xn*Pn(n)
     &     -(xn+2.d0)*xi*Pn(n+1)+(xn+2.d0)*Pn(n+2))
     &     *(2.d0*(ch-xi)*dUn(n)+sh*Un(n))
 7    continue

      duZetaB=aux+aux1

c     scaling par  -1/(4c²(ch-x)^(3/2))
c     ---------------------------------
      du0zetaB=-1.d0/(4*c*c
     &         *(ch-xi)**(3.d0/2.d0))
     &         *duZetaB

c      write(8,*) xi,du0zetaB,duZetaB

c     ===========================
c     Fin de la fonction du0zetaB
c     ===========================
      end function du0zetaB

