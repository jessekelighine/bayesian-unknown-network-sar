# Latex Makefile

MAIN := main
NAME := bayesian-unknown-network-sar
TIKZ_DIR := figures-tikz
FIGURES_TIKZ := $(shell find $(TIKZ_DIR) -type f -name '*.tex' | sed 's/\.tex$$/.pdf/')
FIGURES := figures
FILES := $(MAIN).tex \
	 *.tex \
	 *.sty \
	 references.bib

.PHONY: all clean

all: $(NAME).pdf

clean:
	latexmk -C
	cd $(TIKZ_DIR) ; latexmk -C ; cd ..
	rm -rf $(NAME).pdf

$(NAME).pdf: $(MAIN).pdf
	cp $(MAIN).pdf $(NAME).pdf

$(MAIN).pdf: $(FILES) $(FIGURES)/* $(FIGURES_TIKZ)
	latexmk -xelatex -synctex=1 $(MAIN).tex

$(TIKZ_DIR)/%.pdf: $(TIKZ_DIR)/%.tex
	latexmk -outdir=$(<D) -xelatex $<
