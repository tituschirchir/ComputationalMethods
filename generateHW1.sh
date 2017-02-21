file=$1
echo "$file"
Rscript -e "rmarkdown::render('$file".Rmd"')"
evince $file".pdf"
