
c****************************************************************
c                 Sous-programme optimisant la 
c               valeur du Nt en fonction du choix
c                       de lsa de départ
c
c
c
c     CalculTabs : appel de sous programme
c
c
c     Programmeur: M.Guemas; Janvier 2013.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      subroutine Ntoptimal(lsa,Ntinitial,precision,Bn,Dn,alpha,
     &             alpha0,tabs,tabs0,tabsB,tabsB0,R0,Rn,Nt)

c     ===========================================================
      implicit none 

c     declaration des variables en entrees
c     ------------------------------------
      integer, intent(in) :: Ntinitial
      double precision, intent(in) :: lsa, precision

c     Declaration des variables locales
c     ---------------------------------
      integer :: i, ind, n, k, Neta, pi, c, Nt2
      double precision :: eta, xn, xi, error
      double precision :: tabs01, tabs02, alpha01, alpha02 
      double precision :: tabsB01, tabsB02, R01, R02
      double precision,dimension(5000) :: tabD1
      double precision,dimension(5000) :: tabD2 
      double precision,dimension(5000) :: diffD2D1
      double precision,dimension(5000) :: tabs1, tabs2
      double precision,dimension(5000) :: diffs2s1
      double precision,dimension(5000) :: tabsB1, tabsB2
      double precision,dimension(5000) :: diffsB2sB1
      double precision,dimension(5000) :: tabB1, tabB2
      double precision,dimension(5000) :: diffB2B1
      double precision,dimension(5000) :: alpha1, alpha2
      double precision,dimension(5000) :: diffal2al1
      double precision,dimension(5000) :: Rn1, Rn2
      double precision,dimension(5000) :: diffRn1Rn2

c     declaration des variables en sorties
c     ------------------------------------
      integer :: Nt
      double precision, intent(out) :: alpha0,tabs0
      double precision, intent(out) :: tabsB0,R0
      double precision,dimension(5000) :: Bn,Dn
      double precision,dimension(5000) :: tabs,tabsB
      double precision,dimension(5000) :: alpha,Rn      
      
c     Declaration des interfaces
c     --------------------------     
      interface

         subroutine CalculTabsB(lsa,Nt,tabs,tabs0,tabsB,tabsB0,
     &                  Bn,Dn,alpha,alpha0,R0,Rn)
           integer, intent(in) :: Nt
           double precision, intent(in) :: lsa
           double precision, intent(out) :: alpha0,tabs0
           double precision, intent(out) :: tabsB0,R0
           double precision,dimension(Nt+1) :: Bn,Dn
           double precision,dimension(Nt+1) :: tabs,tabsB
           double precision,dimension(Nt) :: alpha, Rn
         end subroutine CalculTabsB

         function fmax(t,ind,n)
           integer::ind,n
           double precision,dimension(n)::t
           double precision ::fmax
         end function fmax

      end interface
c     ===========================================================      
      
 
c     Donnée de l'entier de tronquature Nt
c     -----------------------------------
      write(*,*)'Ntinitial=',Ntinitial

c     variable d'erreur minimale
c     --------------------------
      write(*,*)'precision=',precision

      open(unit=7,file='alphan_Ntoptimal.dat')
 
c     Boucle sur la valeur de Nt
c     -------------------------
      do 1 k=1,4
         Nt=(2**(k-1))*Ntinitial
c         write(*,*) 'Nt=', Nt
c         write(*,*) 'k=', k
         if (k.eq.1) then
c            write (*,*) '1er appel de Calcultabs'
            call CalculTabsB(lsa,Nt,tabs1,tabs01,tabsB1,
     &         tabsB01,tabB1,tabD1,alpha1,alpha01,R01,Rn1) 
c             write (*,*) '2ème appel de Calcultabs'
           Nt2=Nt+2 
           write(*,*) 'Nt2=', Nt2 
           call CalculTabsB(lsa,Nt2,tabs2,tabs02,tabsB2,
     &       tabsB02,tabB2,tabD2,alpha2,alpha02,R02,Rn2) 
c           write(*,*)'alpha01=',alpha01         
c$$$           do 101 n=1,Nt
c$$$             write(*,*)'(tabB1',n,')=',tabB1(n)
c$$$ 101       continue 
           write(*,*)
c           write(*,*)'alpha02=',alpha02
c          write(*,*)'tabs02=',tabs02
c$$$           do 102 n=1,Nt2 
c$$$            write(*,*)'(tabs2',n,')=',tabs2(n)
c$$$ 102       continue    
         else              
c            write (*,*) '3ème appel de Calcultabs'   
c            Nt=20
            
           call CalculTabsB(lsa,Nt,tabs2,tabs02,tabsB2,
     &          tabsB02,tabB2,tabD2,alpha2,alpha02,R02,Rn2)  
         write(*,*)
         write(*,*) 'Nt=', Nt, 'Nt2=',Nt2
          
c$$$           write(*,*)'tabsB02=',tabsB02
c$$$           do 122 n=1,Nt 
c$$$              write(*,*)'(tabB2',n,')=',tabB2(n)
c$$$ 122       continue    

          
c$$$          do 122 n=1,Nt 
c$$$              write(*,*)'(tabB2',n,')=',tabB2(n)
c$$$ 122       continue    
c$$$           
c$$$           do 123 n=1,Nt 
c$$$              write(*,*)'(tabD2',n,')=',tabD2(n)
c$$$ 123       continue    



         endif  

         do 11 n=1,Nt/2
              diffB2B1(n)=tabB2(n)-tabB1(n)
              diffD2D1(n)=tabD2(n)-tabD1(n)
              diffs2s1(n)=tabs2(n)-tabs1(n)
              diffsB2sB1(n)=tabsB2(n)-tabsB1(n)
              diffal2al1(n)=alpha2(n)-alpha1(n)
              diffRn1Rn2(n)=Rn2(n)-Rn1(n)
c             write(*,*) 'diffB2B1(',n,')=',diffB2B1(n)
c             write(*,*) 'diffD2D1(',n,')=',diffD2D1(n)
c              write(*,*) 'diffsB2sB1(',n,')=',diffsB2sB1(n)
c           write(*,*) 'alpha2 et alpha1',alpha2(n), alpha1(n)
 11      continue 

c$$$         write(*,*) 'fmax(diffRn1Rn2)=',fmax(diffRn1Rn2,Nt,Nt)
c$$$         write(*,*) 'fmax(diffal2al1)=',fmax(diffal2al1,Nt,Nt)     
c$$$         write(*,*) 'fmax(diffs2s1)=',fmax(diffs2s1,Nt,Nt)
c$$$         write(*,*) 'fmax(diffB2B1)=',fmax(diffB2B1,Nt,Nt)
c$$$         write(*,*) 'fmax(diffD2D1)=',fmax(diffD2D1,Nt,Nt)
c$$$         write(*,*) 'tabs02-tabs01',tabs02-tabs01
c$$$         write(*,*) 'alpha02-alpha01',alpha02-alpha01
c$$$         write(*,*) 'R02-R01',R02-R01
c$$$         write(*,*) 'fmax(diffsB2sB1)=',fmax(diffsB2sB1,Nt,Nt)

      error= max(abs(fmax(diffB2B1,Nt,Nt)),
     &       abs(fmax(diffD2D1,Nt,Nt)),
     &       abs(fmax(diffal2al1,Nt,Nt)),
     &       abs(fmax(diffs2s1,Nt,Nt)),tabs02-tabs01,
     &       alpha02-alpha01,abs(fmax(diffRn1Rn2,Nt,Nt)),
     &       R02-R01,abs(fmax(diffsB2sB1,Nt,Nt)),
     &       tabsB02-tabsB01)            

c         write(*,*) 'error =',error

         if (error < precision) then
            write(*,*) 'error < precision'
             
            ! initialisation des tableaux de sortie
             do 123 n=1,Nt+1
                 Bn(n)=0.d0
                 Dn(n)=0.d0
                 tabs(n)=0.d0
                 tabsB(n)=0.d0                
 123          continue

          do 12 n=1,Nt
              Bn(n) = tabB2(n)              
c             write(*,*) 'Bn(',n,')=',Bn(n)
              Dn(n) = tabD2(n)
            tabs(n) = tabs2(n)
	   tabsB(n) = tabsB2(n)            
 12      continue
         
          do 13 n=1,Nt
            alpha(n) = alpha2(n)
            write(7,*) n, alpha(n) 
               Rn(n) = Rn2(n)              
 13       continue  
           tabsB0 = tabsB02
           alpha0 = alpha02
            tabs0 = tabs02
               R0 = R02
            write(*,*) 'error =',error            
            write(*,*)'erreur acceptable, on sort de Ntoptimal' 
          go to 999
         else 
            write(*,*) 'error > precision'
            tabB1  = tabB2            
            tabD1  = tabD2
            tabs1  = tabs2
            alpha1 = alpha2
           alpha01 = alpha02
           tabs01  = tabs02
               R01 = R02
               Rn1 = Rn2
            tabsB1 = tabsB2
           tabsB01 = tabsB02
        endif

 1    continue

 999  continue
       close(7)
        end subroutine Ntoptimal   
c      -------- Fin du programme--------------------------------
****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

c     fmax : Fonction calculant le maximum d'un tableau

      function fmax(t,ind,n)

      implicit none 
c     ===========================================================
        integer::ind,n,i
        double precision, dimension(n)::t
        double precision ::fmax
c     ===========================================================

        fmax=t(1)
        do i=1,ind
         if (abs(t(i))>abs(fmax)) then
c            write(*,*) 'abs(t(i))=',abs(t(i))
            fmax=t(i)
c            write(*,*) 'fmax=',t(i)
         end if   
        end do
      end function fmax
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
