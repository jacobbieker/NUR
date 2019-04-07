#!/bin/bash

echo "Run handin template"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist create it!"
  mkdir plots
fi

echo "Check if the sine movie exist"
if [ -e sinemovie.mp4 ]; then
  echo "Remove mp4 file"
  rm sinemovie.mp4
fi

# Script that returns a plot
echo "Run main script ..."
python3 main.py

echo "Generating the pdf"

pdflatex template.tex
bibtex template.aux
pdflatex template.tex
pdflatex template.tex


