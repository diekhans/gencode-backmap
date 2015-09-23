ROOT = .
.SECONDARY:

docs += gencode-backmap

docsPdf = ${docs:%=output/%.pdf}
docsFontchk = ${docs:%=output/%.fontchk}
drawings = output/transMap2.pdf
images = $(wildcard images/*.pdf) $(wildcard images/*.jpg) $(wildcard images/*.png)
depends = ${drawings} ${images}

all: ${docsPdf} ${docsFontchk} ${depends}

ltxopts = -file-line-error-style -output-directory=output
ltxmode =\\nonstopmode\\input

output/%.pdf: %.ltx ${depends}
	@mkdir -p $(dir $@)
	pdflatex ${ltxopts} ${ltxmode} $< </dev/null || (rm -f $@; false)
	pdflatex ${ltxopts} ${ltxmode} $< </dev/null || (rm -f $@; false)

output/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

output/%.fontchk: output/%.pdf
	pdfFontChk $<
	touch $@

output/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

output/%.pdf: output/%.eps
	epspdf $< $@.tmp.pdf
	mv -f $@.tmp.pdf $@

output/%.eps: drawings/%.gnuplot
	@mkdir -p $(dir $@)
	gnuplot $< >$@.tmp
	mv -f $@.tmp $@

clean:
	rm -rf output
