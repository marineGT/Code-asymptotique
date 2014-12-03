c****************************************************************
c                 Programme de calcul de la Force de 
c             trainee pour une deformation de l'ordre 1
c             lorsqu'une Bulle remonte vers la surface
c            plane en utilisant les coordonnees bipolaires
c               dans le cas Ca petit et Bond = 1
c                     
c   
c
c     spline, gaussint1dim : appel de sous programme
c
c
c     Programmeur: M.Guemas; Octobre 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      subroutine Force1(lsa,Nt,Nx,hxi,dhxi,tabxi,Bn,Dn,
     &                  alpha,alpha0,F1)

c     ===========================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none

c     ===========================================================
c     Declaration des variables en entree
c     -----------------------------------      
      integer :: Nt, Nx
      double precision :: lsa, alpha0      
      double precision, dimension(Nt) :: alpha
      double precision, dimension(Nt+1) :: Bn, Dn
      double precision, dimension(Nx) :: hxi, dhxi      
      double precision, dimension(0:Nx) :: tabxi 
   
c     Declaration des variables locales
c     ---------------------------------
      integer :: n, k, i
      double precision :: prec, re, a, b, xi, r
      double precision :: xn, xNt, pi, c, resultat, aux 
      double precision :: zetap
      double precision, dimension(Nx) :: coefb1, coefb2 
      double precision, dimension(Nx) :: coefc1, coefc2
      double precision, dimension(Nx) :: coefd1, coefd2


c     Variables pour le test de spline.f90 et ispline.f90
c     ---------------------------------------------------
      integer :: test
      integer, parameter :: np=11
      integer, parameter :: nint=21
      double precision :: xmax, xmin, step, x, y
      double precision :: dy1, dyN, ys1
      double precision :: ys, error
      double precision, dimension(np) :: xa, ya
      double precision, dimension(np) :: dy, dyExa 
      double precision, dimension(np) :: coefb, coefe, coefd
      double precision, dimension(nint) :: xderiv 



c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: F1
c      double precision :: F1, zeta0
c      double precision, dimension(Nx) :: zeta

c     Declaration des interfaces
c     ---------------------------
      interface

          subroutine spline(x, y, b, c, d, n)
            integer  n
            double precision  x(n),y(n),b(n),c(n),d(n)
          end subroutine spline
       
         function fonctionb(eta,Nt,Nx,c,zetap,hxi,dhxi,tabxi,
     &            Bn,Dn,alpha,alpha0,b1,e1,d1,b2,e2,d2)
          integer, intent(in)  ::  Nt,Nx      
          double precision, intent(in)  :: eta, c, zetap 
          double precision, dimension(Nx) :: hxi, dhxi
          double precision, dimension(0:Nx) :: tabxi
          double precision, dimension(Nt+1) :: Bn,Dn
          double precision, dimension(Nt) :: alpha
          double precision, intent(in) :: alpha0        
          double precision, dimension(Nx) :: b1, e1, d1  
          double precision, dimension(Nx) :: b2, e2, d2     
          double precision :: fonctionb 
         end function fonctionb
c$$$
c$$$         function gaussint1dim(fonctionb,a,b,Nt,Nx,hxi,dhxi,
c$$$     &                      Bn,Dn,alpha,alpha0,b1,e1,d1,
c$$$     &                      b2,e2,d2,precision)
c$$$           integer Nt           
c$$$           character name*(*)
c$$$           parameter (name = 'gaussint1dim')
c$$$           external fonctionb           
c$$$         end function gaussint1dim
               
      end interface 

      double precision :: gaussint1dim
      external gaussint1dim

c     ===========================================================
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)
c      write (*,*) 'pi=',pi      

c     calcul du reel c 
c     -----------------
      c=dsqrt(lsa**2-1.d0)

c     calcul de la constante zetap
c     ----------------------------
      zetap=dlog(lsa+dsqrt(lsa**2-1.d0))

c     Initialisation de la Force
c     --------------------------
      F1=0.d0

c     Tests :
c     -------
      test=0   
c     ==========================================================
c     Test spline
c     On teste avec la fonction f(x)=sin(x)

      if (test.eq.1) then
      xmin=0.d0
      xmax=2.d0
c     On genere xa, ya 
      write(*,*) 'calcul des derivees de sin(x) '
      write(*,*) 
      step=(xmax-xmin)/(np-1)

c     Calcul de y=sin(x) et dyExa=cos(x)
c     ---------------------------------
      do 76 i=1,np
         xa(i)=xmin+step*dble(i-1)
         ya(i)=dsin(xa(i))
         dyExa(i)=dcos(xa(i))
c         write(*,*) i,xa(i),dyExa(i)
 76   continue

c     Calcul de la derivee de sin(x) suivant la methode des 
c     difference finies centrees
c     Calcul pour n=2,np-1
c     --------------------
      do 77 i=2,np-1
         xa(i)=xmin+step*dble(i-1)
         dy(i)=(ya(i+1)-ya(i-1))/(2*step) 
        write(*,*) i,xa(i),dyExa(i),dy(i)
 77   continue

c     Calcul au bornes selon la methode difference finie decentrees
c     df/dx=[4*f(x+h)-f(x+2*h)-3*f(x)]/(2*f(x)) + O(h²)
c     ----  ----------------------------------------------
      dy(1)=(4*ya(2)-ya(3)-3*ya(1))/(2*step) 
      dy(np)=(ya(np-2)-4*ya(np-1)+3*ya(np))/(2*step) 
c      write(*,*) 'dy(1) =',dy(1),'dy(np) =',dy(np)
c      write(*,*) dyExa(1),dyExa(np)

c     Interpolation de la derivee methode 1 :
c     ---------------------------------------
      write(*,*) 'appel de la fonction spline pour cos(x) '
      write(*,*) 
       call spline(xa,dyExa,coefb,coefe,coefd,np)
c       do i=1,np
c          write(*,*) i,b(i),e(i),d(i)
c          write(*,*) i,e(i)
c          write(*,*) i,d(i)
c       end do  
       
c     Interpolation de la derivee methode 2 :
c     ---------------------------------------  
      write(*,*) 'appel de la fonction spline pour methode diff finie '
      write(*,*) 
       call spline(xa,dy,coefb,coefe,coefd,np)
c       do i=1,np
c       write(*,*) i,b(i),e(i),d(i)
c          write(*,*) i,e(i)
c          write(*,*) i,d(i)
c       end do  

      endif
c     Fin du Test Spline
c     ==========================================================

c     Calcul des coefficients d'interpolation 
c     ---------------------------------------
c     write(*,*) 'appel de la fonction spline pour la fonction h'
c      write(*,*) 

       call spline(tabxi,hxi,coefb1,coefc1,coefd1,Nx)
       call spline(tabxi,dhxi,coefb2,coefc2,coefd2,Nx)


c$$$        open(unit=7,file='coef_1.dat')
c$$$        write(7,*)'Coefficients 1 :'
c$$$        write(7,*)
c$$$         do 1 i=1,Nx           
c$$$         write(7,*) i, coefb1(i), coefc1(i), coefd1(i)
c$$$1       continue
c$$$        write(7,*)
c$$$        write(7,*)'Coefficients 2 :'
c$$$        write(7,*)
c$$$         do 10 i=1,Nx           
c$$$         write(7,*) i, coefb2(i), coefc2(i), coefd2(i)
c$$$ 10      continue
c$$$        close(7)

c       open(unit=7,file='tabxi_F1.dat')
c      On range les valeurs de zeta à chaque pas de temps 
c      dans le fichier de sortie position
c       do 11 i=0,Nx 
c          write(7,*) tabxi(i)
c 11    continue
            
c$$$        r=c*sqrt((1.d0+tabxi(0))/(1.d0-tabxi(0)))
c$$$        aux=(1.d0-tabxi(0))**(3.d0/2.d0) 
c$$$        zeta0=aux*hxi(1)
c$$$        write(numfich1,*) tabxi(0),r, (c*c)*zeta0/(1.d0-tabxi(0))
c$$$
c$$$        do 11 i=2,Nx          
c$$$         aux=(1.d0-tabxi(i))**(3.d0/2.d0) 
c$$$         r=c*sqrt((1.d0+tabxi(i))/(1.d0-tabxi(i)))
c$$$         zeta(i)=aux*hxi(i-1)
c$$$         write(numfich1,*) tabxi(i), r, (c*c)*zeta(i)/(1.d0-tabxi(i))
c$$$ 11     continue
c$$$        write(numfich1,*) 
c       close(7)

c     Calcul de la Force sur la sphere solide
c     ---------------------------------------
      prec=1.d-5 ! precision de l'integrale
      a=0.d0
      b=pi

c     On calcule l'integrale de fonctionb
c     -----------------------------------

c      open(unit=8,file='values_fonctionb.dat')
       F1=gaussint1dim(fonctionb,a,b,Nt,Nx,c,zetap,
     &    hxi,dhxi,tabxi,Bn,Dn,alpha,alpha0,coefb1,
     &    coefc1,coefd1,coefb2,coefc2,coefd2,prec)
c       write (*,*) 'Force1 dans sous programme ',F1 
c      F1=1.d-3       
c     On scale F1 avec 2pi
c     -----------------------------------         
      F1=c*c*2.d0*pi*F1
c      write(*,*)'2pi dasn Force1.f', 2.d0*pi
 
        end subroutine Force1   
c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

