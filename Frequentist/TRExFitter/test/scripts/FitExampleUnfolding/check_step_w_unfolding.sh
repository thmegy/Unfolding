 #!/bin/bash
diff -w -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleUnfolding" LOG_UNFOLDING_w test/logs/FitExampleUnfolding/LOG_UNFOLDING_w && diff FitExampleUnfolding/PruningText.txt test/reference/FitExampleUnfolding/PruningText.txt
