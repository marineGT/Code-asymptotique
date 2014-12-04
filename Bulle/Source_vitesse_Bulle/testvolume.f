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

      function testvolume(teta,Nx,g,tabg,b,e,d)
c     ============================================================
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nx
      double precision, intent(in)  :: teta
      double precision, dimension(Nx) :: b, e, d
      double precision, dimension(0:Nx) :: g
      double precision, dimension(0:Nx) :: tabg

c     declaration des variables locales
c     ---------------------------------
      double precision :: gInt

c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: testvolume

c     Declaration des interfaces
c     ---------------------------
      interface

         function ispline(u, x, y, b, e, d, n)
           integer n
           double precision  u,x(n),y(n),b(n),e(n),d(n)
           double precision ispline
         end function ispline

      end interface
c     ================================================================== 

c     On interpole la valeur de phi correspondant au teta
c     --------------------------------------------------- 
      gInt=ispline(teta,tabg,g,b,e,d,Nx) 
c      write (*,*) 'gInt=', gInt

c     Calcul de la constante du CM pour une valeur de teta
c     ----------------------------------------------------
     
      testvolume=fInt*dsin(teta)

c      write(16,*) teta, testvolume

      end function testvolume
c     -------- Fin du programme----------------------------------
c****************************************************************

