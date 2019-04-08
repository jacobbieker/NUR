#!/bin/bash

echo "Seed is 5227"

echo "Run handin template"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist create it!"
  mkdir plots
fi

# Script that returns a plot
echo "Run main script ..."
python3 main.py

echo "Generating the pdf"

pdflatex template.tex
bibtex template.aux
pdflatex template.tex
pdflatex template.tex


