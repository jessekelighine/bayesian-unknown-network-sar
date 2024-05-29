FORMATTER := format.vim
TRIALS := $(shell ls . | grep -E 'txt$$' | sed 's/.txt$$//')

.PHONY: all
all: $(TRIALS)

.PHONY: clean
clean:
	rm -rf *.tex

%: %.txt
	nvim --headless +'source $(FORMATTER) | write! $*.tex | q!' $< >/dev/null
