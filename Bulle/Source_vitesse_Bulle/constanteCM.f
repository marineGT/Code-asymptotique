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

      function constanteCM(teta,Nx,phi,tabteta,b1,e1,d1)
c     ============================================================
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nx
      double precision, intent(in)  :: teta
      double precision, dimension(Nx) :: b1, e1, d1
      double precision, dimension(0:Nx) :: phi
      double precision, dimension(0:Nx) :: tabteta

c     declaration des variables locales
c     ---------------------------------
      double precision :: phiInt

c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: constanteCM

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
      phiInt=ispline(teta,tabteta,phi,b1,e1,d1,Nx) 

c     Calcul de la constante du CM pour une valeur de teta
c     ----------------------------------------------------
     
      constanteCM=phiInt*dsin(teta)*dcos(teta)   

c      write(16,*) teta, constanteCM

      end function constanteCM
c     -------- Fin du programme----------------------------------
c****************************************************************

