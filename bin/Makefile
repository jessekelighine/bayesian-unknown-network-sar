###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Description: make inference files                                           #
#                                                                             #
# Last Modified: 2024-05-26                                                   #
###############################################################################

INFERENCE_FILE := inference.R
TRIALS_P := P10 P11 P12 P13 P14 P20 P21 P22 P23 P24
TRIALS_S := S10 S11 S12 S20 S21 S22
TRIALS_K := K10 K11 K12 K20 K21 K22
TRIALS := ${TRIALS_P} ${TRIALS_S} ${TRIALS_Z}

.PHONY: inference ${TRIALS}

inference: ${TRIALS}

P%:
	Rscript ${INFERENCE_FILE} $@ 10

S%:
	Rscript ${INFERENCE_FILE} $@ 12

K%:
	Rscript ${INFERENCE_FILE} $@ 15
