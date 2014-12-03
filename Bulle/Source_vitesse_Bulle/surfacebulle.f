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

      subroutine surfacebulle(Nt,Nx,Ng,precision,lsa,c,c2,
     &	           zetap,lambda0,gammahat,tabsB,tabsB0,alpha,
     &             alpha0,Bn,Dn,zeta1B,dzeta1B,tabx,
     &             Rn,R0,phi,tabteta)
c     ===========================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none

c     ===========================================================
c     Declaration des variables en entrees
c     ------------------------------------
      integer, intent(in) :: Nt, Nx, Ng
      double precision, intent(in) :: lsa, c, c2, zetap
      double precision, intent(in) :: lambda0, gammahat
      double precision :: precision, tabsB0, alpha0      
      double precision,dimension(Nt) :: alpha,Rn
      double precision :: R0
      double precision, dimension(Nt+1) :: Bn,Dn
      double precision,dimension(Nt+1) :: tabsB
      
      
c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, k, erreur 
      double precision :: aux, xk, dx, y  
      double precision, dimension(Nx+1) :: tabsigma0zz

c     Declaration des variables en sorties
c     ------------------------------------      
      double precision, dimension(0:Nx) :: zeta1B,dzeta1B
      double precision, dimension(0:Nx) :: tabx
      double precision, dimension(0:Nx) :: phi
      double precision, dimension(0:Nx) :: tabteta

c     Declaration des interfaces
c     ---------------------------
      interface

         subroutine resolutionsystLU(Nt,Nx,Ng,precision,lsa,c,
     &	                   c2,dx,zetap,lambda0,gammahat,tabsB,
     &                     tabsB0,alpha,alpha0,Bn,Dn,tabx,
     &                     zeta1B,Rn,R0,tabsigma0zz)
	   integer, intent(in) :: Nt, Nx, Ng
           double precision, intent(in) :: precision, lsa
	   double precision, intent(in) :: c, c2, dx, zetap
	   double precision, intent(in) :: lambda0, gammahat
           double precision, intent(in) :: tabsB0, alpha0,R0           
           double precision,dimension(Nt) :: alpha,Rn
           double precision, dimension(Nt+1) :: Bn,Dn
	   double precision, dimension(Nt+1) ::	 tabsB
           double precision, dimension(0:Nx) :: tabx,zeta1B
           double precision, dimension(Nx+1) :: tabsigma0zz
         end subroutine resolutionsystLU
	
        subroutine  TestCL(Nt,Nx,c,zetap,lambda0,alpha,
     &              alpha0,Bn,Dn,tabsigma0zz,phi,tabteta)
          integer, intent(in) :: Nt, Nx
          double precision, intent(in) :: c, zetap
          double precision, intent(in) :: lambda0, alpha0
          double precision,dimension(Nt) :: alpha
          double precision, dimension(Nt+1) :: Bn,Dn
          double precision, dimension(Nx+1) :: tabsigma0zz
          double precision, dimension(0:Nx) :: phi
          double precision, dimension(0:Nx) :: tabteta
        end subroutine TestCL

      end interface

c     ===========================================================

c     Lecture de la tronquature Nt
c     ----------------------------
      write(*,*) 'la troncature Nt =', Nt

c     Determination de dx
c     -------------------
      dx=2.D0/dble(Nx)      
      write(*,*)'dx',dx
c	 ______________________________________________________
c	 ------ On construit le tableau des h(xi) -------------
   
c     Resolution du systeme lineaire
c     ------------------------------
c       write(*,*)'alpha0',alpha0
c       open(unit=15,file='S0_part2.dat')
c       open(unit=16,file='S0_sigma0zz.dat')
      call resolutionsystLU(Nt,Nx,Ng,precision,lsa,c,c2,
     &	                dx,zetap,lambda0,gammahat,tabsB,
     &                  tabsB0,alpha,alpha0,Bn,Dn,tabx,
     &                  zeta1B,Rn,R0,tabsigma0zz)
      write(*,*) 'apres resolutionsystLU '

     
      
      write(*,*)'On teste les conditions limites'

       call TestCL(Nt,Nx,c,zetap,lambda0,alpha,alpha0,
     &             Bn,Dn,tabsigma0zz,phi,tabteta)
      write(*,*) 'apres TestCL'

c      close(15)
c      close(16)
c	 ______________________________________________________
c	 ---------- Calcul de la derivée dzeta1B --------------

    
      dx=2.d0/dble(Nx-1)  !Calcul du nouveaux pas
      write(*,*) 'le pas dx = ', dx 
      write(*,*) 'Nx = ', Nx 
      dzeta1B(0)=(4*zeta1B(1)-zeta1B(2)+3*zeta1B(0))/(2.d0*dx)

      do 2 k=1,Nx-1        
       dzeta1B(k)=(zeta1B(k+1)-zeta1B(k-1))/(2.d0*dx)
 2    continue
  
      dzeta1B(Nx)=(3*zeta1B(Nx)-4*zeta1B(Nx-1)
     &           +zeta1B(Nx-2))/(2.d0*dx)
      
      
      end subroutine surfacebulle

c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

