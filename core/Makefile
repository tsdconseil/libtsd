# Makefile pour libdsp/core
#   Copyright 2017 - 2020 J.A. / http://www.tsdconseil.fr

# Nom du projet
NOM_PROJET = libtsd

# Fichiers sources


 
SOURCES := hadamard ct hb frat tableau tests-gen moniteur-cpu  doa stats tod cqt limits 
SOURCES += clock-rec modulateur demod-dec demod-ndec carrier-rec itrp wav 
SOURCES += fourier  prbs telecom ecp bitstream 
SOURCES += estimation-delais detection emetteur
SOURCES += filtre-plot filtre-analyse ra fenetres divers
SOURCES += figure tsd axes filtrage filtre-rt stdo image freetype
SOURCES += axes  canva test-figure
SOURCES += goertzel polyphase hilbert egalisation freqestim
SOURCES += rif-eq rif-cs rif-freq rif-fen cic rii analogique  
SOURCES += modulations
SOURCES += etalement-spectre canalisation transpo-bb simulation
SOURCES += image-bmp cmaps unites recepteur temps geometrie
SOURCES += kalman modele-imu  ssm-plot misc-plot

# test-kalman 
SOURCES_TEST += test-tab test-dsp test-geometrie test-temps test-telecom test-stats test-poly test-image 
SOURCES_TEST += test-cqt test-tod test-ra test-fenetres test-fourier test-filtres test-crec test-tsd 
SOURCES_TEST += test-filtre-fft test-detecteur test-wav test-itrp test-filtrage-analyse


# Dossiers où trouver les fichiers sources (autres que ./src and ./include)
AUTRES_DOSSIERS_SOURCES += src/apps src/vue src/misc src/fourier src/filtrage tests src/reechan
AUTRES_DOSSIERS_SOURCES += src/telecom src/telecom/codes
AUTRES_DOSSIERS_SOURCES += src/telecom/canalisation src/kalman

# Optional dependencies:
include build-options



#ifeq ($(LIBTSD_USE_OMP),1)
#	OPTIONS_COMPIL += -DLIBTSD_USE_OMP=1 -fopenmp 
#	OPTIONS_LINK	 += -fopenmp 
#endif


OPTIONS_COMPIL += -I/mingw64/include/freetype2

# Makefile standard
include std-makefile-defs


# PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif


install: all-nocls
	install -d /usr/share/tsd/
	install -m 644 ./data/fonts/OpenSans-Regular.ttf /usr/share/tsd
	install -d $(DESTDIR)$(PREFIX)/lib/
	cp $(BUILD_DIR)/lib.a $(BUILD_DIR)/libtsd.a
	install -m 644 -D $(BUILD_DIR)/libtsd.a $(DESTDIR)$(PREFIX)/lib
	(cd ./include && find tsd -type f -exec install -Dm 644 "{}" "$(DESTDIR)$(PREFIX)/include/{}" \;)
	(cd ./include && find dsp -type f -exec install -Dm 644 "{}" "$(DESTDIR)$(PREFIX)/include/{}" \;)



docs:
	sh -c "cd ../doxy; make"

docs-en:
	sh -c "cd ../doxy-en; make"

