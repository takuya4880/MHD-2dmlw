      module cansio
      contains

c======================================================================|
      subroutine dacdef0s(idf,file,mtype)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      write(idf) 1

c   version
      write(idf) 0

c   for data type: integer=4, real=5, double=6
      write(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      write(idf) 1

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      write(idf) -1

      return
      end subroutine


c======================================================================|
      subroutine dacdef1d(idf,file,mtype,ix)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      write(idf) 1

c   version
      write(idf) 0

c   for data type: integer=4, real=5, double=6
      write(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      write(idf) 1

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      write(idf) ix

      return
      end subroutine 



c======================================================================|
      subroutine dacdef2s(idf,file,mtype,ix,jx)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      write(idf) 1

c   version
      write(idf) 0

c   for data type: integer=4, real=5, double=6
      write(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      write(idf) 3

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      write(idf) ix,jx,-1

      return
      end subroutine 



c======================================================================|
      subroutine dacdefparam(idf,file)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open (idf,file=file,form='formatted')

      return
      end subroutine 



c======================================================================|
      subroutine dacputparamd(mf_params,name,value)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
      double precision value
c----------------------------------------------------------------------|

      write(mf_params,9020) 5,name,value
9020  format('#dac ',i4,a32,' ',e20.13)

      return
      end subroutine 



c======================================================================|
      subroutine dacputparami(mf_params,name,mvalue)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
c----------------------------------------------------------------------|

      write(mf_params,9010) 4,name,mvalue
9010  format('#dac ',i4,a32,' ',i20)

      return
      end subroutine 


c======================================================================|
      subroutine dacopnr0s(idf,file,mtype,nx)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      read(idf) munity
      if (munity.ne.1) then
        write(6,*) 'Error:: endian different'
        stop
      endif

c   version
      read(idf) mversion
      if (mversion.ne.0) then
        write(6,*) 'Error:: version different. version=', mversion
        stop
      endif

c   for data type: integer=4, real=5, double=6
      read(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      read(idf) mdim
      if (mdim.ne.1) then
        write(6,*) 'Error:: dimension different. mdim=', mdim
        stop
      endif

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      read(idf) nx

      return
      end subroutine

c======================================================================|
      subroutine dacopnr2s(idf,file,mtype,ix,jx,nx)
c======================================================================|
c
c NAME  ncdefss
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      read(idf) munity
      if (munity.ne.1) then
        write(6,*) 'Error:: endian different'
        stop
      endif

c   version
      read(idf) mversion
      if (mversion.ne.0) then
        write(6,*) 'Error:: version different. version=', mversion
        stop
      endif

c   for data type: integer=4, real=5, double=6
      read(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      read(idf) mdim
      if (mdim.ne.3) then
        write(6,*) 'Error:: dimension different. mdim=', mdim
        stop
      endif

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      read(idf) ix,jx,nx

      return
      end subroutine


      end module
