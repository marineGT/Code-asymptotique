c****************************************************************
c                 Programme de calcul de la Force de 
c             trainee a l'ordre 0, i.e sans deformation
c             lorsqu'une bulle remonte vers la surface
c            plane en utilisant les coordonnees bipolaires
c               dans la cas Ca petit et Bond = 1
c                     
c
c     Connaissant les tableaux Bn et Dn, on calcule directement 
c     le coefficiant de trainee Lambda relie à la force, 
c     sachant que les tableaux An et Cn sont nuls dans notre cas.
c
c     external : Bpsi.f, Dpsi.f
c
c     sachant que  F = -4*pi*nu*U*lambda pour une bulle
c
c     Programmeur: M.Guemas; Janvier 2013
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      subroutine Force0(Nt,c,Bn,Dn,lambda0,F0)
c     ===========================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ===========================================================
c     declaration des variables en entrees
c     ------------------------------------
      integer, intent(in) :: Nt    
      double precision,intent(in) :: c
      double precision, dimension(Nt):: Bn, Dn    

c     Declaration des variables locales
c     --------------------------------- 
      integer :: n
      double precision :: xn, pi   

c     Declaration des variables en sorties
c     ------------------------------------    
      double precision :: F0, lambda0

c      external Bpsi,Dpsi
c     ==================================================================
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)

c     Boucle calculant la Force 
c     -------------------------- 
      F0=0.d0

      do 6 n=1,Nt+1
        xn=dble(n)
        F0=F0-(Bn(n)+Dn(n))*(2.d0*xn+1.d0)
 6    continue
      F0=(2*pi*sqrt(2.d0)/c)*F0
c      write(*,*) 'Force0=',F0 



c     Boucle calculant le coefficient de trainee
c     ------------------------------------------ 
c     initialisation du tableau Lambda     
      lambda0=0.d0     
      lambda0=-F0/(4*pi)
c      write(*,*) 'lambda0=',lambda0
 
        

       end subroutine  Force0  
c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
