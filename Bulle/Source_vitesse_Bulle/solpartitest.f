c****************************************************************
c                Solution particuliere de l'equation 
c                   differentielle pour la surface
c                           de la bulle
c     
c
c     Appelée par la fonction gaussint1dim3 pour le calcul
c     de l'integrale cette solution pour differente valeur
c     de teta
c
c     Programmeur: M.Guemas; Fevrier 2014
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function solpartitest(teta)
c     ============================================================
c     declaration des variables en entree
c     -----------------------------------
      double precision, intent(in)  :: teta

c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: solpartitest

c     ================================================================== 

c     Calcul de la solution totale pour une valeur de teta
c     ----------------------------------------------------
      solpartitest=(dcos(3.d0*teta)-dcos(teta))*dsin(teta)

c      write(16,*) teta, solpartitest

      end function solpartitest
c     -------- Fin du programme----------------------------------
c****************************************************************

