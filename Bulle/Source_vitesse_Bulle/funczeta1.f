c****************************************************************
c                 Sous-programme calculant les 
c		   fonctions zeta1(xi) et zeta1'(xi) 
c                      pour la surface libre
c
c	diagonalessystem,                     
c	thomas1d : appel de sous programme
c
c     Programmeur: M.Guemas; Janvier 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      subroutine funczeta1(Nt,Nx,precision,lsa,c,c2,Bo,Ca,zetap,
     &                 tabs,tabs0,alpha0,R0,Rn,hxi,dhxi,tabxi)
c     ===========================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none

c     ===========================================================
c     Declaration des variables en entrees
c     ------------------------------------
      integer, intent(in) :: Nt, Nx
      double precision, intent(in) :: lsa, c, c2, Bo, Ca,zetap
      double precision :: precision, tabs0, alpha0,R0
      double precision, dimension(Nt),intent(in) :: Rn
      double precision,dimension(Nt+1) :: tabs
      
c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, k, erreur 
      double precision :: aux, xk, dx, y 
      double precision, dimension(0:Nx) :: ainf, aprin
      double precision, dimension(0:Nx) :: asup, d
      
c     Declaration des variables en sorties
c     ------------------------------------ 
      double precision, dimension(Nx) :: hxi, dhxi     
      double precision, dimension(0:Nx) :: tabxi

c     Declaration des interfaces
c     ---------------------------
      interface

         subroutine diagonalessysteme1(Nt,Nx,precision,lsa,
     &                  c,c2,dx,Bo,zetap,tabs,tabs0,
     &                  tabxi,ainf,aprin,asup,d)
	   integer, intent(in) :: Nt, Nx
           double precision, intent(in) :: precision, lsa
	   double precision, intent(in) :: c, c2, dx, zetap
	   double precision, intent(in) :: Bo, tabs0           
           double precision, dimension(Nt+1) ::	tabs
	   double precision, dimension(0:Nx) :: tabxi
           double precision, dimension(0:Nx) :: ainf, aprin
	   double precision, dimension(0:Nx) :: asup, d
         end subroutine diagonalessysteme1	

         subroutine thomas1d(imin,imax,a,b,c,d,x,erreur)
            integer, intent(in) :: imin, imax
            double precision, dimension(imin:imax) :: a,b,c,d
c            dimension(imin:imax) :: a,b,c,d
            double precision, dimension(imin:imax) :: x
            integer :: erreur
         end subroutine thomas1d	      
	
      end interface

c     ===========================================================

c     Lecture de la tronquature Nt
c     ----------------------------
      write(*,*) 'la troncature Nt =', Nt

c     Determination de dx
c     -------------------
      dx=2.D0/dble(Nx)
      
c      write(*,*)'dx',dx
c	 ______________________________________________________
c	 ------ On construit le tableau des h(xi) -------------


c     Calcul des diagonales du syst�me
c     --------------------------------*
c$$$       do 11 k=1,Nt         
c$$$         write(*,*) 'tabs(',k,') =',tabs(k) 
c$$$ 11    continue

      call diagonalessysteme1(Nt,Nx,precision,lsa,c,c2,dx,Bo,
     &              zetap,tabs,tabs0,tabxi,ainf,aprin,asup,d)
c      write(*,*) 'apres diagonalessysteme'
c      open(unit=12,file='secondmembre_CLX1.dat')
c      do 22 n=0,Nx
c       write(*,*) tabxi(n)
c       write(*,*) d(n)
c 22    continue
c     Resolution du systeme tridiagonal
c     ---------------------------------
c      call thomas1d(0,Nx,ainf,aprin,asup,d,zeta1,erreur)
      call thomas1d(0,Nx,ainf,aprin,asup,d,hxi,erreur)
c      write(*,*) 'apres thomas1d'

c$$$c     Ecriture du tableau hxi
c$$$c     -----------------------
c$$$c      write(*,*) 'hxi0',hxi(1)
c$$$      open(unit=12,file='hxi_res1.dat')
c$$$
c$$$      do 11 i=0,Nx-1
c$$$
c$$$c       Calcul de la valeur de x en i
c$$$c       -----------------------------
c$$$        y=-1.D0+dx*dble(i)  
c$$$c        write(*,*) 'i=',i, 'xi=',xi
c$$$
c$$$c       Impression dans le fichier de sortie
c$$$c       ------------------------------------
c$$$        write(12,*) y,c*dsqrt((1.d0+y)/(1.d0-y)),hxi(i+1)
c$$$
c$$$ 11   continue
c$$$      close(12)



c$$$c     Calcul de la deformee zeta1 = (1-xi)**(3/2)
c$$$c     -------------------------------------------       
c$$$      do 1 k=0,Nx-1
c$$$         aux=(1.d0-tabxi(k))**(3.d0/2.d0)
c$$$         zeta1(k+1)=aux*hxi(k+1)     
c$$$1     continue   

      
c	 ______________________________________________________
c	 ---------- Calcul de la deriv�e dzeta1 --------------

    
       dx=2.d0/dble(Nx-1)  !Calcul du nouveaux pas
c       write(*,*) 'le pas dx = ', dx 
c       write(*,*) 'Nx = ', Nx 
c$$$      dzeta1(1)=(4*zeta1(2)-zeta1(3)+3*zeta1(1))/(2.d0*dx)
c$$$c      write(*,*) 'hxi(2) = ', hxi(2)
c$$$c      write(*,*) 'hxi(3) = ', hxi(3)
c$$$c      write(*,*) 'hxi(1) = ', hxi(1)
c$$$      write(*,*) 'dzeta1(1) = ', dzeta1(1)
c$$$      do 2 k=2,Nx-1        
c$$$       dzeta1(k)=(zeta1(k+1)-zeta1(k-1))/(2.d0*dx)
c$$$ 2    continue     
c$$$      dzeta1(Nx)=(3*zeta1(Nx)-4*zeta1(Nx-1)
c$$$     &           +zeta1(Nx-2))/(2.d0*dx)

       dhxi(1)=(4*hxi(2)-hxi(3)-3*hxi(1))/(2.d0*dx)
c$$$       write(*,*) 'hxi(2) = ', hxi(2)
c$$$       write(*,*) 'hxi(3) = ', hxi(3)
c$$$       write(*,*) 'hxi(1) = ', hxi(1)
c$$$       write(*,*) 'dhxi(1) = ', dhxi(1)
       do 2 k=2,Nx-1        
       dhxi(k)=(hxi(k+1)-hxi(k-1))/(2.d0*dx)
 2     continue     
       dhxi(Nx)=(4*hxi(Nx-1)-hxi(Nx-2)-3*hxi(Nx))
     &          /(2.d0*dx)
      
      
c$$$      do 3 k=0,Nx-1
c$$$         aux=(1.d0-tabxi(k))
c$$$         dzeta1(k+1)=((aux)**(3.d0/2.d0))*dhxi(k+1)
c$$$     &        -(3.d0/2.d0)*(aux)**(1.d0/2.d0)*hxi(k+1)
c$$$3     continue   

 
      end subroutine funczeta1

c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

