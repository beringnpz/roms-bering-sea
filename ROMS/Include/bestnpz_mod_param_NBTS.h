#  ifdef FEAST
        integer, parameter :: NBTS = 155 ! May not be able to do FEAST w/o CARBON?
#  elif defined CARBON
        integer, parameter :: NBTS = 150
#  else
        integer, parameter :: NBTS = 148
#  endif
