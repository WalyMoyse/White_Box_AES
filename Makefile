CFLAGS = -Wall -Wextra -Werror
CPPFLAGS = -I../include/
LDFLAGS = -lm

all :
	(cd src; make all; mv AES ..)

clean :
	(rm AES; cd src; make clean)

help : 
	(cd src; make help)

.PHONY : all clean help