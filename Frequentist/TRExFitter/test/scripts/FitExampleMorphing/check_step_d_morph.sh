#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing" LOG_MORPH_d test/logs/FitExampleMorphing/LOG_MORPH_d && for file in `ls FitExampleMorphing/Tables/*.txt`; do diff -w $file test/reference/$file; done
