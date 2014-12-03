      subroutine diagonalessysteme1(Nt,Nx,precision,lsa,c,c2,dx,Bo,
     &                   zetap,tabs,tabs0,tabxi,ainf,aprin,asup,d)
c     ==================================================================
c     Description :
c     ------------
c     Ce programme cree la matrice AX=B en utilisant une methode de 
c     difference finie centree.
c     Soit f'(x) = [f(x+h) - f(x-h)]/2h ou bien f'n = [fn+1 - fn-1]/2dx
c     Et de même pour la derivee seconde
c          f"n = [fn+1 - 2fn + fn-1]/h²

c     Si l'on a l'eq differentielle telle que :

c                a(xn)f"(xn)+b(xn)f'(xn)+c(xn)f=d(xn)

c     On peut donc remplacer les derivees premiere et seconde par leur 
c     expressions, d'ou
c      an{[fn+1 - 2fn + fn-1]/h²} + bn{[fn+1 - fn-1]/2h} + cn.fn = dn

c     Mais encore
c                fn+1.[(an/h²) + (bn/2h)] + fn.[(-2an/h²) + cn] 
c                    + fn-1.[(an/h²) - (bn/2h)] = dn

c     avec les parties superieur, asup, inferieur, ainf et la diagonale
c     aprinc, de la matrice A telle que 
      
c                 fn+1.asup + fn.aprinc + fn-1.ainf = dn
c
c     Soit        |asup   = (an/h²) + (bn/2h)
c                 |aprinc = (-2an/h²) + cn
c                 |ainf   = (an/h²) - (bn/2h)
c
c
c     Condition aux limites :
c     ----------------------
c     Sachant qu'aux bords de la matrice les termes f(x=-1) et f(x=N+1)
c     n'existe pas ( taille de la matrice A [0,N])
c     On utilise la methode au difference finies decentrees a droite d'ordre
c     deux.
c     Ainsi, on exprime la derivee premiere en fonction des trois point 
c     suivant, f0,f1,f2.

c     f(x+h)  = f(x) + hf'(x) + (h²/2)f"(x) + 0(h**3)
c     f(x+2h) = f(x) + 2hf'(x) + (2*h²)f"(x) + 0(h**3) 
c     On elimine le terme de derivee seconde pour n'avoir que la derivee 
c     premiere en fonction des termes de f
c     4*f(x+h) - f(x+2h) = 3f(x) + 2hf'(x) + 0(h**3)

c     Donc    f'(x) = [4*f(x+h) - f(x+2h) - 3f(x)]/2h + 0(h**3)

c     Notre condition au limite etant, pour le cas de la fonction h(xi)
c         h'0-h0 = 0 en Xi = -1

c     donc [4*h1 - h2 - 3h0]/2h - h0 = 0
c     Soit -(1/2h)*h2 + (2/h)*h1 - [1 + (3/2h)]*h0 = 0  
c
c     D'ou les coefficients a0=1+(3/2h), b0=(2/h) et c0=-(1/2h)

c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      integer, intent(in) :: Nt, Nx
      double precision, intent(in) :: precision, lsa
      double precision, intent(in) :: c, c2, dx, zetap 
      double precision, intent(in) :: Bo, tabs0
      double precision, dimension(Nt+1) :: tabs

c     Declaration des variables en sortie
c     -----------------------------------
      double precision, dimension(0:Nx) :: tabxi      
      double precision, dimension(0:Nx) :: aprin, asup, ainf, d

c     Declaration des variables locales
c     ---------------------------------
      integer :: i
      double precision :: invdx,invdx2,xi,a0,b0,c0,secondmembre0

c     Declaration des interfaces
c     --------------------------
      interface

         function secondmembre(Nt,c,c2,eps,tabs,tabs0,x)
            integer,intent(in) :: Nt 
            double precision, intent(in) :: eps, x, tabs0, c, c2
            double precision, dimension(Nt+1), intent(in):: tabs
            double precision :: secondmembre            
         end function secondmembre
     
      end interface
c     
c     ==================================================================

c     Calcul des inverses de dx et dx2
c     --------------------------------
      invdx=1.D0/dx
      invdx2=invdx*invdx
c      write(*,*) 'dx=', dx
c      open(unit=13,file='secondmembre_res.dat')
      tabxi(0)=-1.d0
      
c      write(*,*) 'tabxi0=', tabxi(0) 
c     Determination des coefficients du systeme tridiagonal pour i=1 à Nx-1
c     ---------------------------------------------------------------------
      do i=1,Nx-1
c        Determination de xi
c        -------------------
         xi=-1.D0+dx*dble(i)
          tabxi(i)=xi
cccc            write(*,*) 'tabxi=', tabxi(i)

c        Calcul du coefficient de la diagonale inferieure
c        ------------------------------------------------
c         ainf(i)=((1.d0-xi)**3)*(1.D0+xi)*invdx2-((1.d0-xi)**2)
c     &            *invdx
         ainf(i)=((1.d0-xi)**3)*(1.d0+xi)*invdx2-(1.d0/2.d0)
     &           *((1.d0-xi)**2)*(2.d0-3.d0*(1.d0+xi))*invdx

c        Calcul du coefficient de la diagonale inferieure
c        ------------------------------------------------
c         aprin(i)=2.D0*(1.D0-xi)-Bo*c2-2.D0*((1.d0-xi)**3)
c     &           *(1.D0+xi)*invdx2
         aprin(i)=(1.D0-xi)*((3.d0/4.d0)*(1.d0+xi)-1.d0)-Bo*c2
     &           -2.d0*((1.d0-xi)**3)*(1.d0+xi)*invdx2
c        Calcul du coefficient de la diagonale inferieure
c        ------------------------------------------------
c         asup(i)=((1.d0-xi)**3)*(1.D0+xi)*invdx2+((1.d0-xi)**2)
c     &           *invdx
         asup(i)=((1.d0-xi)**3)*(1.d0+xi)*invdx2+(1.d0/2.d0)
     &           *((1.d0-xi)**2)*(2.d0-3.d0*(1.d0+xi))*invdx

c        Calcul du coefficient du second membre
c        --------------------------------------
         
c         d(i)=c*((1.d0-xi)**(3.d0/2.d0))*secondmembre(Nt,c,c2,
c     &                  precision,tabs,tabs0,xi)
         d(i)=secondmembre(Nt,c,c2,precision,tabs,tabs0,xi)


c$$$         if (i.eq.10) then
c$$$            stop
c$$$         endif 
c         write(13,*) xi, d(i) 
      enddo
       tabxi(Nx)=1.d0
cccc       write(*,*) 'tabxi=', tabxi(Nx) 
      close(13)
c     Calcul des coefficients pour i=0
c     --------------------------------
c     On utilise une approximation d'ordre deux pour approcher la condition
c     de Robin en xi=0. Afin de rendre le système tridiagonale, on compose
c     cette équation avec celle en i=1.

c     Calcul des coefficients tenant compte uniquement de la condition
c     ----------------------------------------------------------------
c     aux limites
c     -----------
c      a0=1.D0-3.D0*invdx !cas pour zeta(x)
c      b0=4.D0*invdx
c      c0=-invdx
c      
      a0=-1.D0-6.D0*invdx !condition en zero pour h(x)
      b0=8.D0*invdx
      c0=-2.d0*invdx

c       a0=-(12.d0*invdx+2.d0+Bo*c2)!pas de condition en xi=-1
c       c0=-4.d0*invdx
c       b0=16.d0*invdx

c     Calcul des coefficients en i=0 par combination des relations en i=0 et 1
c     ------------------------------------------------------------------------
c     si CL en xi=-1
      ainf(0)=0.D0
      aprin(0)=a0-c0*ainf(1)/asup(1)
      asup(0)=b0-c0*aprin(1)/asup(1)
      d(0)=-c0*d(1)/asup(1) 
c$$$c     Si pas CL en xi=-1
c$$$      ainf(0)=0.D0
c$$$      aprin(0)=a0-c0*ainf(1)/asup(1)
c$$$      asup(0)=b0-c0*aprin(1)/asup(1)
c$$$      secondmembre0=secondmembre(Nt,c,c2,precision,tabs,tabs0,-1.d0)
c$$$      d(0)=secondmembre0-c0*d(1)/asup(1) 

      
c     Calcul des coefficients pour i=Nx
c     ---------------------------------
      xi=1.d0
      ainf(Nx)=0.D0
      aprin(Nx)=1.D0
      asup(Nx)=0.D0
c      d(Nx)=0.D0
c      d(Nx)=1.d0
      d(Nx)=-secondmembre(Nt,c,c2,precision,
     &       tabs,tabs0,xi)/Bo*c2
      write(*,*)'secondmembre(xi=1)=',d(Nx)
c     =======================================
c     Fin du sous-programme diagonalessysteme1
c     =======================================
      end subroutine diagonalessysteme1
