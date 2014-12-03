      subroutine thomas1d(imin,imax,a,b,c,d,x,erreur)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Nom : thomas1d.
c     -----
c     Description : Resolution d'un systeme tridiagonal par un algorith-
c     -------------
c                   de Thomas. Il s'agit de faire un double balayage.
c     Langage : FORTRAN 77.
c     ---------
c     Auteur : S.G.R., EV/F.P.
c     --------
c     ==================================================================
c     Sous-programmes appeles : Aucun.
c     ==================================================================
c     Variables en entree :
c     ---------------------
c     -------------------------------------------------------------------
c     |   Nom   |            Description               | unite |  Type  |
c     |---------|--------------------------------------|-------|--------|
c     |  imin   | indice inf. des matrices a,b,c,d,x   |   -   | entier |
c     |  imax   | indice sup. des matrices a,b,c,d,x   |   -   | entier |
c     |  a      | matrice de la diag. inf.             |   -   |  reel  |
c     |  b      | matrice de la diag.                  |   -   |  reel  |
c     |  c      | matrice de la diag. sup.             |   -   |  reel  |
c     |  d      | matrice du second membre             |   -   |  reel  |
c     |  x      | matrice du vecteur des inconnues     |   -   |  reel  |
c     -------------------------------------------------------------------

c     Variables en entree/sortie :
c     ----------------------------
c     Variables en sortie :
c     ---------------------
c     -------------------------------------------------------------------
c     |   Nom   |            Description               | unite |  Type  |
c     |---------|--------------------------------------|-------|--------|
c     |  x      | matrice du vecteur des inconnues     |   -   |  reel  |
c     -------------------------------------------------------------------

c     Variables internes :
c     --------------------
c     ==================================================================
c     Code retour :
c     -------------
c     erreur = 0, execution normale.
c     erreur = 1, division par zero.
c     ==================================================================
c     Version : V1.0,15/02/01.
c     ---------
c     ==================================================================
c     Historique : V1.0,15/02/01, creation.
c     ------------
c     ==================================================================

c     Declaration des variables en entree
c     -----------------------------------
      integer, intent(in) :: imin,imax
      double precision, dimension(imin:imax), intent(in) :: a,b,c,d

c     Declaration des variables en sortie
c     -----------------------------------
      double precision, dimension(imin:imax) :: x
      integer :: erreur

c     Declaration des variables locales
c     ---------------------------------
      integer i,ip,im
      double precision :: inv,denom
      double precision, dimension(:), allocatable :: cp
      double precision, dimension(:), allocatable :: dp
c     ==================================================================

c     Allocation de la memoire
c     ------------------------
      allocate(cp(imin:imax),dp(imin:imax))

c         Premier balayage, suppression de la diagonale inferieure :
c         ----------------------------------------------------------

c     Initialisation de la diagonale superieure et du second membre :
c     ---------------------------------------------------------------
      if (b(imin).eq.0.D0) then
c        Erreur execution dans le deroulement du sous-programme
         goto 1000
      endif
      inv=1./b(imin)
      cp(imin)=c(imin)*inv
      dp(imin)=d(imin)*inv
c      write(*,*)'dp(0)=',dp(imin)
c     Determination de la diagonale superieure et du second membre
c     ------------------------------------------------------------
c     pour i=imin+1,imax :
c     --------------------
      do 10 i=imin+1,imax
         im=i-1
         denom=b(i)-a(i)*cp(i-1)
         if (denom.eq.0) then
c           Erreur execution dans le deroulement du sous-programme
            goto 1000
         endif          
         inv=1./denom
         cp(i)=c(i)*inv
         dp(i)=(d(i)-a(i)*dp(i-1))*inv
10    continue
c      Deuxieme balayage, suppression de la diagonale superieure :
c      -----------------------------------------------------------      
      x(imax)=dp(imax)
c      write(*,*)'x(Nt)=',x(imax)
      do 20 i=imax-1,imin,-1
         ip=i+1
         x(i)=dp(i)-x(ip)*cp(i)          
20    continue      
c     ==================================================================
c     Fin normale d'execution du sous-programme thomas1d
      erreur=0
c     Liberation de la memoire
      deallocate(cp,dp)
      return
c     ==================================================================
c     Traitement des erreurs :
c     ------------------------
 1000 write(*,*)'Resolution impossible dans le sous-programme thomas1d:'
      write(*,*)'Division par zero.'
      erreur=1
c     Liberation de la memoire
      deallocate(cp,dp)
      return
c     ==================================================================
      end subroutine thomas1d
