#!/bin/sh

jupyter nbconvert --ClearMetadataPreprocessor.enabled=True --ClearOutput.enabled=True --to markdown showcase/README.ipynb --stdout > README.md
