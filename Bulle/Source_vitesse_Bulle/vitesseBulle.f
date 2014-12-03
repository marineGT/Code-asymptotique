c****************************************************************
c                 Programme de calcul de la vitesse                
c             d'une Bulle remontant vers une surface
c            plane en utilisant les coordonnees bipolaires
c               dans la cas Ca petit et Bond = 1
c                     
c
c     Conventions:
c     Bond  : nombre de Bond
c     lsa   : altitude centre de la sphere/rayon a de la sphÃ¨re
c     Nt    : indice de troncature 
c     coefA : coefficiant A appartenant Ã  la solution generale de 
c            l'equation differentielle voir fonction zeta 
c     zeta  : fonction liee a la deformation de l'interface
c     eta   : coordonnee bipsherique (bipolaire), est un angle 
c             compris entre 0 et pi      
c     xi    : correspond a cos(eta) 
c     c     : parametre liee a la distance l entre l'interface et la 
c            sphere, ie sqrt(lÂ²-1)
c     Neta: 
c     xNeta, xk : les parametres Neat et k sont en simple precision, on veut 
c                 les utiliser aussi en double pecision, on en cree des 
c                 nouveaux de mÃªme valeur en double precision
c
c     spline, ispline, Force1, Force0 : appel de sous programme
c
c
c     Programmeur: M.Guemas; Janvier 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      program vitesseBulle
c     ===========================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ===========================================================
c     Declaration des parametres
c     --------------------------
c      integer, parameter :: Nx=200
      character (len=50), parameter :: nomconfig='config' 
      integer, parameter :: numfichdeformee=1 
      integer, parameter :: numfichbulle=2 
      integer, parameter :: numfichhmin=11 

c     Declaration des variables
c     -------------------------
      character (len=100) :: nomfichentree
      character (len=100) :: fichdeformee, fichbulle
      integer :: Nx, Nt, n, k, i, NxInt, Ntinitial, Ng
      double precision :: lsa, Bo, Ca, pi, c, c2, eta, r      
      double precision :: eps, g, dt, t, h, nu, aux
      double precision :: precision, zetap, ch, rb, zb
      double precision :: zdeformee
      double precision :: xb1, zb1, testzetap
      double precision :: alpha0, tabs0, tabsB0, R0
      double precision, dimension(5000) :: Bn,Dn,tabs,tabsB
      double precision, dimension(5000) :: alpha, Rn
      double precision, dimension(5000) :: zeta1
      double precision, dimension(5000) :: dzeta1
      double precision, dimension(:), allocatable :: hxi,dhxi 
      double precision, dimension(:), allocatable :: zeta1B 
      double precision, dimension(:), allocatable :: dzeta1B 
      double precision, dimension(:), allocatable :: tabxi,tabx
      double precision, dimension(:), allocatable :: phi
      double precision, dimension(:), allocatable :: tabteta
c      double precision, dimension(2*Nx-1) :: hInt, dhInt 
      double precision :: F0, F1, force, lambda, lambda0 
      double precision :: u, u1, u0, poids, gammahat


c     Declaration des interfaces
c     ---------------------------
      interface
         subroutine lecturedata(nomfichentree,Nx,Ntinitial,Ng,
     &                         eps,precision,lsa,Bo,gammahat,
     &                         fichdeformee,fichbulle)
            character (len=50), intent(in) :: nomfichentree
            integer :: Nx,Ntinitial,Ng
            double precision :: eps,precision,lsa,Bo,gammahat
            character (len=100) :: fichdeformee
            character (len=100) :: fichbulle
         end subroutine lecturedata

         subroutine Ntoptimal(lsa,Ntinitial,precision,Bn,
     &     Dn,alpha,alpha0,tabs,tabs0,tabsB,tabsB0,R0,Rn,Nt)
           integer, intent(in) :: Ntinitial
           double precision, intent(in) :: lsa, precision
           integer :: Nt
           double precision, intent(out) :: alpha0,tabs0
           double precision, intent(out) :: tabsB0,R0
           double precision,dimension(Nt+1) :: Bn,Dn
           double precision,dimension(Nt+1) :: tabs,tabsB           
           double precision,dimension(Nt) :: alpha,Rn
         end subroutine Ntoptimal

	 subroutine funczeta1(Nt,Nx,precision,lsa,c,c2,Bo,Ca,
     &	  		zetap,tabs,tabs0,alpha0,R0,Rn,hxi,
     &                  dhxi,tabxi)
	   integer :: Nt, Nx
           double precision :: lsa, c, c2, Bo, Ca, zetap
	   double precision :: precision, tabs0, alpha0, R0 
	   double precision,dimension(Nt+1) :: tabs    
           double precision,dimension(Nt) :: Rn
	   double precision, dimension(Nx) :: hxi, dhxi
           double precision, dimension(0:Nx) :: tabxi 
	 end subroutine funczeta1

        subroutine surfacebulle(Nt,Nx,Ng,precision,lsa,c,c2,
     &                zetap,lambda0,gammahat,tabsB,tabsB0,
     &                alpha,alpha0,Bn,Dn,zeta1B,dzeta1B,
     &                tabx,Rn,R0,phi,tabteta)
           integer, intent(in) :: Nt, Nx, Ng
           double precision, intent(in) :: lsa, c, c2, zetap
           double precision, intent(in) :: lambda0, gammahat
           double precision :: precision, tabsB0, alpha0,R0           
           double precision,dimension(Nt) :: alpha,Rn
           double precision, dimension(Nt+1) :: Bn,Dn
           double precision,dimension(Nt+1) :: tabsB           
           double precision, dimension(0:Nx) :: zeta1B,dzeta1B          
           double precision, dimension(0:Nx) :: tabx
           double precision, dimension(0:Nx) :: phi
           double precision, dimension(0:Nx) :: tabteta
        end subroutine surfacebulle
	       
         subroutine Force0(Nt,c,Bn,Dn,lambda0,F0)
           integer,intent(in) :: Nt
           double precision,intent(in) :: c
           double precision, dimension(Nt+1):: Bn, Dn
           double precision :: F0, lambda0
         end subroutine Force0
 
      end interface

c     ===========================================================
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)

c     Recuperation du nom du fichier d'entree
c     ---------------------------------------
      open(unit=1,file=nomconfig)
      read(1,*) nomfichentree
      close(1)

c     Lecture des donnees
c     -------------------
      call lecturedata(nomfichentree,Nx,Ntinitial,Ng,eps,
     &                 precision,lsa,Bo,gammahat,
     &                 fichdeformee,fichbulle)

c     Allocation de la memoire
c     ------------------------      
c      allocate(zeta1(Nx),dzeta1(Nx))
      allocate(hxi(Nx),dhxi(Nx))
      allocate(zeta1B(0:Nx),dzeta1B(0:Nx))      
      allocate(tabxi(0:Nx),tabx(0:Nx))
      allocate(phi(0:Nx),tabteta(0:Nx))
   
c     calcul du reel c
c     ----------------
      c=dsqrt((lsa**2)-1.d0)
      c2=c*c 
      zetap=dlog(lsa+dsqrt((lsa**2)-1.d0))      
      ch=dcosh(zetap)
      write (*,*) 'c=',c,'zetap=',zetap
      write (*,*) 'ch=',ch
c     Initialisation de h et du temps
c     -------------------------------
      h=2.d0
      t=0.d0 ! addimensionnne td=(a/U)*t
      write (*,*) 'h=',h
c      write (*,*) 2.d0*0.7506104411350174
c     Boucle sur lsa
c     --------------
      open(unit=numfichdeformee,file=fichdeformee)
      open(unit=numfichbulle,file=fichbulle)
      open(unit=numfichhmin,file='hmin_bulle.dat')
c      open(unit=12,file='h(xi)_avec_Cl_xi=1.dat')
c$$$      open(unit=12,file='bulle_coor-sphe.dat')
c      open(unit=16,file='Rn.dat')
c$$$      open(unit=16,file='UndUn_resolutionLU.dat')
c      open(unit=15,file='sol-part_Bo01l12.dat')
c      open(unit=17,file='sol-spherique_Bo01l12.dat')
c      open(unit=18,file='sol_totale_Bo01l12.dat')
c$$$      open(unit=19,file='bulle-non-deformee.dat')
c$$$      open(unit=20,file='tesf.dat')

      write (*,*) 'debut de la boucle sur le temps'
c      do while (lsa.ne.(1.d0+eps)) 
         
        ! Calcul de la troncature optimale
        !---------------------------------
         call Ntoptimal(lsa,Ntinitial,precision,Bn,Dn,alpha,
     &             alpha0,tabs,tabs0,tabsB,tabsB0,R0,Rn,Nt)
         
        !Calcul de la Force a l'ordre 0 et du lambda0 
        !--------------------------------------------
         call Force0(Nt,c,Bn,Dn,lambda0,F0)

         write(*,*)'Force0=', F0
         write(*,*)'lambda0=', lambda0
         if(gammahat.eq.0) then
            write(*,*) 'pb division par zero pour rapport R
     &                  tension surface'
         endif         
         Bo=Bo/(gammahat)
         Ca=Bo/(3.d0*lambda0)!On definit la valeur de Ca 

c         Ca=1.d0/5.d0 !test avec Berdan et Leal
c         Ca=5.d-2 !test gravity-only
         write(*,*)'Ca=',Ca
        
        !Construction de zeta1 pour la Surface libre
        !--------------------------------------------
c        Bo=Bo/(gammahat)
         call funczeta1(Nt,Nx,precision,lsa,c,c2,Bo,Ca,zetap,
     &        tabs,tabs0,alpha0,R0,Rn,hxi,dhxi,tabxi)
          write (*,*) 'FIN calcul de funczeta1'
          

        !Calcul de la deformee et sauvegarde dans un fichier de sortie
        !-------------------------------------------------------------
        !ici on trace z/Ca qui est le profil de la deformee 
        !ce profil nous montre de combien la surface se deforme et
        !nous permet de calculer la force
        
        tabxi(0)=-1.d0
        do 11 i=0,Nx-1   
           aux=(1.d0-tabxi(i))**(3.d0/2.d0) 
           r=c*sqrt((1.d0+tabxi(i))/(1.d0-tabxi(i)))
           zeta1(i+1)=aux*hxi(i+1)
           zdeformee=c*dsinh(Ca*zeta1(i+1))/(
     &               (dcosh(Ca*zeta1(i+1))-tabxi(i)))
c           write(*,*) tabxi(i),zeta1(i+1)
           write(numfichdeformee,*)r,zdeformee,zdeformee/Ca
 11     continue
 	            
        !Calcul de la derivée
        !--------------------
         do 111 k=0,Nx-1 
           aux=(1.d0-tabxi(k))
           dzeta1(k+1)=((aux)**(3.d0/2.d0))*dhxi(k+1)
     &                -(3.d0/2.d0)*(aux)**(1.d0/2.d0)*hxi(k+1)
 111      continue  
          
           stop
        !Construction de zeta1 pour la surface de la bulle
        !-------------------------------------------------
         call  surfacebulle(Nt,Nx,Ng,precision,lsa,c,c2,
     &	           zetap,lambda0,gammahat,tabsB,tabsB0,
     &             alpha,alpha0,Bn,Dn,zeta1B,dzeta1B,
     &             tabx,Rn,R0,phi,tabteta)
       
         do 12 i=0,Nx-1 
          !coordonnees bipolaires 
c$$$          rb=c*dsin(tabteta(i))/
c$$$     &       (dcosh(-zetap+Ca*zeta1B(i))-dcos(tabteta(i)))
c$$$          zb=c*dsinh(-zetap+Ca*zeta1B(i))/
c$$$     &       (dcosh(-zetap+Ca*zeta1B(i))-dcos(tabteta(i)))
           !coordonnees spheriques
           xb1=(1.d0+Ca*phi(i))*dsin(tabteta(i))
           zb1=(1.d0+Ca*phi(i))*dcos(tabteta(i))
c           write(12,*) zb1, xb1 
           write(numfichbulle,*) zb1, xb1 
           xb1=(1.d0)*dsin(tabteta(i))
           zb1=(1.d0)*dcos(tabteta(i))
c            write(12,*) tabx(i),zeta1B(i)  
c            write(19,*) zb1, xb1 
            
c            write(numfichbulle,*) rb,zb
 12      continue
        
	!Calcul et Sauvegarde du drainage au cours du temps
        !--------------------------------------------------
c         h=lsa-1.d0+Ca*zeta(1)
         h=lsa-1.d0
          !on inscrit l'evolution de l'epaisseur du film h 
          !en fonction du temps
         write(numfichhmin,*) t,h     
c         write(*,*) t,h
c         write(*,*) 'Ca*zeta1',Ca*zeta(1)
         
        !Calcul de la Force a l'ordre 1 :
        !--------------------------------
c         write(*,*)'c avant Force1', c
         call Force1(lsa,Nt,Nx,hxi,dhxi,tabxi,Bn,Dn,
     &        alpha,alpha0,F1)
c          write(*,*)'2pi*Force1=',2.d0*pi*F1 
c           write(*,*)'2pi=',2.d0*pi0
          write(*,*)'Force1=', F1    
         
c         write(*,*)'c apres Force1', c
	!Calcul des vitesses u1 et u0 :
        !------------------------------	         
c         F1=-F1
c         u0=1.d0/lambda0        !calcul de la vitesse a l'ordre 0

c        utilisation du coefficiant de Franck 
c        dans la relation:  h=h-4.d0*dt/(3.d0*lambda)
         u0=4.d0/(3.d0*lambda0)
         write(*,*) 'u0=',u0
c         write(15,*) lsa,zetap,h,lambda0
c          write(*,*) 'F0=',F0
          
         u1=-(F1*u0)/F0     !calcul de la vitesse a l'ordre 1
c         write(15,*) u1
         write(*,*) 'u1=',u1

         u=u0+Ca*u1 !vitesse totale apres petite deformation
c          u=u0 !vitesse totale apres petite deformation
         write(*,*) 'u=', u

        !Calcul u coeff lambda et de h:
        !------------------------------
        !on rescale les forces par les vitesses
         F0=u0*F0
         write(*,*)'Force0 rescalée', F0
         F1=u0*F1
         write(*,*)'Force1 rescalée', F1

         force=(F0+Ca*F1)  !force totale exercee sur la bulle
         write(*,*) 'force=', force
c         lambda=-force/(6*pi)

        !Calcul du nouveau lsa, c, zetap :
        !---------------------------------
         dt=1.d-3               !le pas de temps
         lsa=lsa-u*dt 
         c=dsqrt(lsa**2-1.d0)
         c2=c*c 
         zetap=-dlog(lsa+dsqrt(lsa**2-1.d0)) 
c         write(12,*) zetap,lsa
c         write(*,*) zetap,lsa
         t=t+dt !incrementation du temps
         
c      enddo
        
      close(numfichdeformee)
      close(numfichbulle)
      close(numfichhmin)
      close(12)
      close(16)
      close(15) 
      close(17) 
      close(18)  
      close(19)  
      close(20)  

c     Liberation de la memoire
c     ------------------------  
      write(*,*)' LIBERATION DE LA MEMOIRE'
c      deallocate(zeta1,dzeta1)
      deallocate(hxi,dhxi)
      deallocate(zeta1B,dzeta1B)      
      deallocate(tabxi,tabx)
      deallocate(phi,tabteta)

      write(*,*)'fin LIBERATION DE LA MEMOIRE'
      end program vitesseBulle 
c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912



