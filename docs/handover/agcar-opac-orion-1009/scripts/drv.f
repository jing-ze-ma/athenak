      program drv
      integer nz,ne,is,k,nk
      real e,s,st,emin,emax
      nk=3000
      emin=0.5
      emax=300.
      open(10,file='verner_ground.txt')
      do nz=1,30
      do ne=1,nz
      do k=1,nk
        e=emin*(emax/emin)**(real(k-1)/real(nk-1))
        st=0.
        do is=1,7
          call phfit2(nz,ne,is,e,s)
          st=st+s
        enddo
        write(10,'(2i3,1p2e13.5)') nz,ne,e,st
      enddo
      enddo
      enddo
      close(10)
      end
