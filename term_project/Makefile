all: build/presentation.pdf build/report.pdf

texoptions = \
	     --lualatex \
	     --interaction=nonstopmode \
	     --halt-on-error \
	     --output-directory=build

build/presentation.pdf: FORCE | build
	latexmk $(texoptions) presentation.tex

build/report.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	report.tex

FORCE:

build:
	mkdir -p build

clean:
	rm -r build
