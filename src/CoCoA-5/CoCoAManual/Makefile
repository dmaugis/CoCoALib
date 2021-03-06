# Makefile for CoCoALib/src/CoCoA-5/CoCoAManual/ directory.

# COCOA_ROOT=../../..
# include $(COCOA_ROOT)/configuration/autoconf.mk
include ../release-files/version
COCOA5_VERSION_DEFINES=-DCOCOA5_VER_MAJ=$(COCOA5_VER_MAJ)  -DCOCOA5_VER_MIN=$(COCOA5_VER_MIN)  -DCOCOA5_VER_MINMIN=$(COCOA5_VER_MINMIN)
CWD=src/CoCoA-5/CoCoAManual/

COCOA5_DOC_DIR=CoCoAManual
COCOA5_DOC_XML=CoCoAHelp.xml
COCOA5_DOC_AUX_DIR=aux-files
COCOA5_DOC_HTML_DIR=html
COCOA5_DOC_TEX_DIR=tex

# Default target
.PHONY: pdf-and-html-doc
pdf-and-html-doc: $(COCOA5_DOC_HTML_DIR)/TimeStamp  $(COCOA5_DOC_TEX_DIR)/TimeStamp

# Target for HTML files: uses a TimeStamp to tell whether needs to be rebuilt
$(COCOA5_DOC_HTML_DIR)/TimeStamp: CoCoAHelp.xml
	/bin/mkdir -p $(COCOA5_DOC_HTML_DIR)
	cd $(COCOA5_DOC_HTML_DIR); \
	touch WillBecomeTimeStamp; \
	java -jar ../aux-files/saxon.jar  ../$(COCOA5_DOC_XML) ../aux-files/GUI_help.xsl  CoCoAVersion="$(COCOA5_VER_MAJ).$(COCOA5_VER_MIN).$(COCOA5_VER_MINMIN)" ; \
	/bin/mv WillBecomeTimeStamp TimeStamp



$(COCOA5_DOC_TEX_DIR)/TimeStamp: CoCoAHelp.xml
	/bin/mkdir -p $(COCOA5_DOC_TEX_DIR)
	cd $(COCOA5_DOC_TEX_DIR); \
	touch WillBecomeTimeStamp; \
	java -jar ../aux-files/saxon.jar  ../$(COCOA5_DOC_XML) ../aux-files/TeX.xsl CoCoAVersion="$(COCOA5_VER_MAJ).$(COCOA5_VER_MIN).$(COCOA5_VER_MINMIN)"  > CoCoAManual.tex
	/bin/cp -r $(COCOA5_DOC_AUX_DIR)/TeX-extra-files/*  $(COCOA5_DOC_TEX_DIR)/.
	@cd $(COCOA5_DOC_TEX_DIR); \
	echo "LaTeX --> PDF: pass 1 of 3"; pdflatex -halt-on-error CoCoAManual > /dev/null; if [ $$? -ne 0 ]; then exit 1; fi; \
	echo "LaTeX --> PDF: pass 2 of 3"; pdflatex -halt-on-error CoCoAManual > /dev/null; if [ $$? -ne 0 ]; then exit 1; fi; \
	echo "LaTeX --> PDF: final pass";  pdflatex -halt-on-error CoCoAManual > /dev/null; if [ $$? -ne 0 ]; then exit 1; fi; \
	fgrep "LaTeX Warning:" CoCoAManual.log; \
	/bin/mv WillBecomeTimeStamp TimeStamp

.PHONY: clean
clean:
	@echo "Cleaned CoCoALib/$(CWD)"

.PHONY: veryclean
veryclean:
	@/bin/rm -rf  $(COCOA5_DOC_HTML_DIR)  $(COCOA5_DOC_TEX_DIR)
	@echo "Verycleaned CoCoALib/$(CWD)"
