#!/bin/sh

jupyter nbconvert --to webpdf --allow-chromium-download --no-input showcase/README.ipynb
