OBJS     = sdpmain.o sdpallo.o sdpbasi.o sdpchec.o sdpdata.o\
           sdpdire.o sdpinit.o sdpinvs.o sdpmatx.o sdpmmdg.o\
           sdpnfac.o sdpperm.o sdppque.o sdpprin.o sdpproc.o\
           sdpshut.o sdpsymb.o sdpvalu.o sdpvect.o sdpupda.o\
           sdp2vec.o sdpvhat.o sdpmult.o

CC       = gcc
CFLAGS   = -O

sdp: sdpmain.o $(OBJS)
	${CC} $(CFLAGS) $(OBJS) -o sdp -lm
