#
gnuplot homogeneousoverT.plt
pdflatex tex2pdf.tex 
pdfcrop tex2pdf.pdf homogeneousoverT.pdf
rm plot.tex plot.eps plot-eps-converted-to.pdf
rm tex2pdf.log tex2pdf.aux tex2pdf.pdf

gnuplot homogeneousoverE.plt
pdflatex tex2pdf.tex 
pdfcrop tex2pdf.pdf homogeneousoverE.pdf
rm plot.tex plot.eps plot-eps-converted-to.pdf
rm tex2pdf.log tex2pdf.aux tex2pdf.pdf

gnuplot kappaoverTlfixedsmall.plt
pdflatex tex2pdf.tex 
pdfcrop tex2pdf.pdf kappaoverTlfixedsmall.pdf
rm plot.tex plot.eps plot-eps-converted-to.pdf
rm tex2pdf.log tex2pdf.aux tex2pdf.pdf

gnuplot kappaoverETfixedlsmall.plt
pdflatex tex2pdf.tex 
pdfcrop tex2pdf.pdf kappaoverETfixedlsmall.pdf
rm plot.tex plot.eps plot-eps-converted-to.pdf
rm tex2pdf.log tex2pdf.aux tex2pdf.pdf

gnuplot kappaoverlEfixed.plt
pdflatex tex2pdf.tex 
pdfcrop tex2pdf.pdf kappaoverlEfixed.pdf
rm plot.tex plot.eps plot-eps-converted-to.pdf
rm tex2pdf.log tex2pdf.aux tex2pdf.pdf
