c
c	Parameter file for RFI subroutines
c
        integer maxlayer, maxsublayer, max_nd
        integer maxmoddim, maxdata, maxwave
        parameter       (maxlayer = 15,
     &                   maxsublayer = 100,
     &                   max_nd = maxlayer * 4,
     &                   maxmoddim = max_nd,
     &                   maxdata = 1024,
     &                   maxwave = 10 )

