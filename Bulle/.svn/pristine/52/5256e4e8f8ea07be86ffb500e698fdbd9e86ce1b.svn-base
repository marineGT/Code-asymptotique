      subroutine lecturedata(nomfichentree,Nx,Ntinitial,Ng,eps,
     &                       precision,lsa,Bo,gammahat,
     &                       fichdeformee,fichbulle)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------
      integer, parameter :: numfich=1

c     Declaration des variables en entree
c     -----------------------------------
      character (len=50), intent(in) :: nomfichentree

c     Declaration des variables en sortie
c     -----------------------------------
      integer :: Nx,Nbitmax,Ntinitial,Ng
      double precision :: eps,precision,lsa,Bo,gammahat
      character (len=100) :: fichdeformee
      character (len=100) :: fichbulle

c     Declaration des variables locales
c     ---------------------------------
      character (len=20) :: nomvar

c     ==================================================================

      open(unit=numfich,file=nomfichentree)

c     Lecture du nombre d'espace sur l'intervalle [-1;1]
c     --------------------------------------------------
      read(numfich,*) nomvar,Nx
      print *,'Nx=',Nx

c     Lecture de la troncature intiale imposee ( au minimum = 20)
c     ----------------------------------------
      read(numfich,*) nomvar,Ntinitial
      print *,'Ntinitial=',Ntinitial

c     Lecture du nombre de point de Gauss sur l'intervalle [-1;1]
c     ----------------------------------------------------------
      read(numfich,*) nomvar,Ng
      print *,'Ng=',Ng

c     Lecture de la valeur minimale pour lsa
c     --------------------------------------
      read(numfich,*) nomvar,eps
      print *,'eps=',eps

c     Lecture de la precision des erreurs de calculs
c     ----------------------------------------------
      read(numfich,*) nomvar,precision
      print *,'precision=',precision

c     Lecture de la position de la particule
c     --------------------------------------
      read(numfich,*) nomvar,lsa
      print *,'lsa=',lsa

c     Lecture du nombre de Bond
c     -------------------------
      read(numfich,*) nomvar,Bo
      print *,'Bo=',Bo

c     Lecture du rapport des tension de surface
c     -----------------------------------------
      read(numfich,*) nomvar,gammahat
      print *,'gammahat=',gammahat

c     Lecture du nom de fichier forme de la surface libre 
c     ---------------------------------------------------
      read(numfich,*) nomvar,fichdeformee
      print *,'fichdeformee=',fichdeformee

c     Lecture du nom de fichier forme de la bulle
c     -------------------------------------------
      read(numfich,*) nomvar,fichbulle
      print *,'fichbulle=',fichbulle

c     Fermeture du fichier de donnees
c     -------------------------------
      close(numfich)

c     ====================================
c     Fin du sous-programme lecturedata
c     ====================================
      end subroutine lecturedata
