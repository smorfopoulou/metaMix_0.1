SRC = concatReads.c svn_version.c
HEADERS = concatReads.h svn_version.h
OBJECTS = concatReads.o svn_version.o

all:concatReads

concatReads : $(OBJECTS)
	gcc -o concatReads $(OBJECTS)

## taken from: http://subversion.apache.org/faq.html#version-value-in-source
##
## on every build, record the working copy revision string
##
#svn_version.c: FORCE
	#printf 'const char* svn_version(void) { const char* SVN_Version = "' > svn_version.c
	#svnversion -n .                   >> svn_version.c
	#echo '"; return SVN_Version; }'   >> svn_version.c

%.o : %.c %.h
	gcc -std=gnu99 -g -Wall -O2 -c $<

clean:
	rm -f concatReads *.o
.PHONY: all clean FORCE
